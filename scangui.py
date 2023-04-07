"""
scangui.py
Author: Thomas Sutter
Description: 
    GUI for the UED data acquisition program
"""


"""
The goal of this script is to provide a user interface for data collection.
It should be able to:

1.) Live display of images from phosphor camera
2.) Live display of intensity within ROI
3.) Live display of scan results
4.) Scan settings configuration with tkinter gui
5.) "smart" scans, a scan that tries to fill in high error data points

"""

import tkinter as tk
from tkinter import filedialog as fd
from tkinter import ttk
import random as rand
import matplotlib
matplotlib.use('TkAgg')
# ~ from matplotlib.ticker import FormatStrFormatter
from matplotlib.backends.backend_tkagg import FigureCanvasTk, NavigationToolbar2Tk
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import matplotlib.widgets
# ~ import matplotlib.animation as animation
# ~ from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk)
import numpy as np
import time
import datetime
import threading
import tifffile
import cv2
import json


def update_errorbar(errobj, x, y, xerr=None, yerr=None):
    ln, caps, bars = errobj


    if len(bars) == 2:
        assert xerr is not None and yerr is not None, "Your errorbar object has 2 dimension of error bars defined. You must provide xerr and yerr."
        barsx, barsy = bars  # bars always exist (?)
        try:  # caps are optional
            errx_top, errx_bot, erry_top, erry_bot = caps
        except ValueError:  # in case there is no caps
            pass

    elif len(bars) == 1:
        assert (xerr is     None and yerr is not None) or\
               (xerr is not None and yerr is     None),  \
               "Your errorbar object has 1 dimension of error bars defined. You must provide xerr or yerr."

        if xerr is not None:
            barsx, = bars  # bars always exist (?)
            try:
                errx_top, errx_bot = caps
            except ValueError:  # in case there is no caps
                pass
        else:
            barsy, = bars  # bars always exist (?)
            try:
                erry_top, erry_bot = caps
            except ValueError:  # in case there is no caps
                pass

    ln.set_data(x,y)

    try:
        errx_top.set_xdata(x + xerr)
        errx_bot.set_xdata(x - xerr)
        errx_top.set_ydata(y)
        errx_bot.set_ydata(y)
    except NameError:
        pass
    try:
        barsx.set_segments([np.array([[xt, y], [xb, y]]) for xt, xb, y in zip(x + xerr, x - xerr, y)])
    except NameError:
        pass

    try:
        erry_top.set_xdata(x)
        erry_bot.set_xdata(x)
        erry_top.set_ydata(y + yerr)
        erry_bot.set_ydata(y - yerr)
    except NameError:
        pass
    try:
        barsy.set_segments([np.array([[x, yt], [x, yb]]) for x, yt, yb in zip(x, y + yerr, y - yerr)])
    except NameError:
        pass

class RoiRectangle():
    def __init__(self,cx,cy,w,h,ax):
        self.cx = cx
        self.cy = cy
        self.w = w
        self.h = h
        self.patch = patches.Rectangle((self.cx, self.cy), self.w, self.h, linewidth=1, edgecolor='r', facecolor='none', alpha=0.8)
        self.live = False
        
        self.cRadius = 1
        # ~ self.cpatch = patches.Circle((self.cx,self.cy),self.cRadius, edgecolor=(1, 0, 0, 0.5), facecolor=(1, 0, 0, 1),fill=False)
        
        ax.add_patch(self.patch) # Add the patch to the axis
        # ~ ax.add_patch(self.cpatch) # Add the center patch to the axis
        
    def goLive(self):
        if self.w < 1:
            self.w = 1 # Minimum width is a single pixel
        if self.h < 1:
            self.h = 1 # Minimum height is a single pixel
        tl_x = self.cx - self.w/2
        tl_y = self.cy - self.h/2
        self.patch.set_xy((tl_x,tl_y)) # Update patch top left corner
        self.patch.set_width(self.w) # Update the patch with new width
        self.patch.set_height(self.h) # Update the patch with new height
        self.live = True
        
class Camera():
    #SCAF
    def __init__(self,exposure=1,gain=30):
        self.exposure = exposure
        self.gain = gain
    def setExposure(self, exposure):
        self.exposure = exposure
    def setGain(self, gain):
        self.gain = gain
    def collectImage(self):
        time.sleep(self.exposure)
        
        x, y = np.linspace(-3,3,100), np.linspace(-3,3,100)
        X, Y = np.meshgrid(x,y)
        
        camImage = np.exp(-(X**2/4+Y**2)/.2)*np.cos((X**2+Y**2)*10)**2
        camImage = np.array(100*(np.random.random(np.shape(camImage))),dtype=np.uint16)
        return camImage

class DelayStage():
    def __init__(self, dsPos=0):
        self.dsPos = dsPos
    def setPos(self, dsPos):
        self.dsPos = dsPos
        print(self.dsPos)

class SetupApp():
    def __init__(self, wait = .033):        
        # Setup the tkinter GUI window
        self.root = tk.Tk()
        self.root.geometry("1400x1000")
        self.root.configure(bg='white')
        self.root.wm_title("Smart Scan")
        self.root.iconbitmap("kogar_0.ico")
        
        # Connect to the camera, collect the initial image and plot the image
        # SCAF
        self.cam = Camera()
        self.camImage = self.cam.collectImage()
        self.h, self.w = self.camImage.shape
        gs = gridspec.GridSpec(2, 2,width_ratios=[self.w,self.w*.1], height_ratios=[self.h,self.h*.1])
        self.camFig = plt.figure()
        self.camAx = [plt.subplot(gs[0]),]
        self.camAx.append(plt.subplot(gs[1],sharey=self.camAx[0]))
        self.camAx.append(plt.subplot(gs[2],sharex=self.camAx[0]))
        self.camAx[0].xaxis.set_tick_params(labelbottom=False)
        self.camAx[0].yaxis.set_tick_params(labelleft=False)
        self.camAx[0].set_xticks([])
        self.camAx[0].set_yticks([])
        
        # Connect to the delay stage
        self.ds = DelayStage()
        
        # Data display settings and variables
        self.wait = int(1000*wait) # Convert to ms
        self.ppw = 100
        self.rm = [] # list of rois
        self.roiActive = False
        
        # Some general formatting vars
        padyControls = 2
        padySep = 5
        entryWidth = 30
        sepWidth = 190
        sepHeight = 3
        cntrlBg = '#2F8D35'
        
        # Create intensity history arrays for data plots
        self.resetHistory = False
        self.resetLongHistory = False
        self.intensityHistory = [0]
        self.avgIntensity = np.array([0],dtype='float')
        self.stdIntensity = np.array([0],dtype='float')
        self.timeHistory = [0]
        self.t0 = time.time()
        self.t1 = self.t0
        
        # Make a Tkinter frame for the data plots
        self.dataFrame = tk.Frame(master=self.root, width=100, height=100, padx=5, pady=5)
        self.dataFrame.pack(side=tk.RIGHT)
        
        # Make a Tkinter frame for the ROI selection buttons
        self.controls = tk.Frame(master=self.root, width=100, height=100, padx=5, pady=5, bg=cntrlBg)
        self.controls.pack(side=tk.RIGHT)

        # This is the GUI controls section
        
        # Adding the ROI controls
        self.roiSelection = tk.StringVar(self.controls)
        self.roiSelection.set("line profile") # default value
        add_roi_button = tk.OptionMenu(self.controls, self.roiSelection, "line profile", "roi rectangle")
        def saveRoisButtonFunc():
            print('Saved ROIS')
            roiDict = {}
            roiDict['#'] = ['cx','cy','w','h']
            for i in range(len(self.rm)):
                roiDict[i] = [self.rm[i].cx,self.rm[i].cy,self.rm[i].w,self.rm[i].h]
            with open('rois.txt','w') as rf:
                for key, value in roiDict.items():
                    rf.write('%s:%s\n' % (key, value))
        self.saveRoisButton = tk.Button(
            master=self.controls,
            text="Save ROIS",
            width=18,
            height=1,
            fg='black', bg='white',
            command=saveRoisButtonFunc)
        def loadRoisButtonFunc():
            # Clear all existing ROIS
            for i in range(len(self.rm)):
                self.rm[i].patch.remove() # Remove the patch from the plot
                del self.rm[i].patch # Delete the patch
            self.rm = [] # make rm an empty list
            # Load the new ROIS from the file
            with open('rois.txt','r') as rf:
                for line in rf:
                    if line[:line.find(":")] != '#':
                        cx,cy,w,h = list(map(int,json.loads(line[line.find(":")+1:].replace('\n',''))))
                        self.rm.append(RoiRectangle(cx,cy,w,h,self.camAx[0]))
                        self.rm[-1].live = True
        self.loadRoisButton = tk.Button(
            master=self.controls,
            text="Load ROIS",
            width=18,
            height=1,
            fg='black', bg='white',
            command=loadRoisButtonFunc)
        def clearRoisButtonFunc():
            # Clear all existing ROIS
            for i in range(len(self.rm)):
                self.rm[i].patch.remove() # Remove the patch from the plot
                del self.rm[i].patch # Delete the patch
            self.rm = [] # make rm an empty list
        self.clearRoisButton = tk.Button(
            master=self.controls,
            text="Clear ROIS",
            width=18,
            height=1,
            fg='black', bg='white',
            command=clearRoisButtonFunc)
        
        # Adding the camera controls
        self.saveImageButton = tk.Button(
            master=self.controls,
            text="Save Cam Image",
            width=18,
            height=1,
            fg='black', bg='white',
            command=self.saveImage)
        self.gain_label_string = tk.StringVar()
        self.gain_label_string.set(f'Gain = {self.cam.gain} dB')
        self.gain_label = tk.Label(master=self.controls,textvariable=self.gain_label_string,bg=cntrlBg)
        self.gainEntry = tk.Entry(master=self.controls,width=entryWidth)
        self.gainEntry.insert(-1,self.cam.gain)
        self.exposure_label_string = tk.StringVar()
        self.exposure_label_string.set(f'Exposure = {self.cam.exposure} s')
        self.exposure_label = tk.Label(master=self.controls,textvariable=self.exposure_label_string,bg=cntrlBg)
        self.exposureEntry = tk.Entry(master=self.controls,width=entryWidth)
        self.exposureEntry.insert(-1,self.cam.exposure)
        def camSettingsButtonFunc():
            if self.gainEntry.get() != '':
                gain = float(self.gainEntry.get())
                self.cam.setGain(gain)
                self.gain_label_string.set(f'Gain = {gain} dB')
            if self.exposureEntry.get() != '':
                exposure = float(self.exposureEntry.get())
                self.cam.setExposure(exposure)
                self.exposure_label_string.set(f'Exposure = {exposure} s')
        self.camSettingsButton = tk.Button(
            master=self.controls,
            text="Update Cam Settings",
            width=18,
            height=1,
            fg='black', bg='white',
            command=camSettingsButtonFunc)
        
        # Adding the ds controls
        self.currentDsPos_label_string = tk.StringVar()
        self.currentDsPos_label_string.set(f'Current dsPos = {self.ds.dsPos}')
        self.currentDsPos_label = tk.Label(master=self.controls,textvariable=self.currentDsPos_label_string,bg=cntrlBg)
        dsPos1_label = tk.Label(master=self.controls,text="dsPos1 (mm)",bg=cntrlBg)
        def dsPos1ButtonFunc():
            dsPos = float(self.dsPos1Entry.get())
            self.ds.setPos(dsPos)
            self.currentDsPos_label_string.set(f'Current dsPos = {dsPos}')
        self.dsPos1Entry = tk.Entry(master=self.controls,width=entryWidth)
        self.dsPos1Button = tk.Button(
            master=self.controls,
            text="Move dsPos1",
            width=18,
            height=1,
            fg='black', bg='white',
            command=dsPos1ButtonFunc)
        dsPos2_label = tk.Label(master=self.controls,text="dsPos2 (mm)",bg=cntrlBg)
        def dsPos2ButtonFunc():
            dsPos = float(self.dsPos2Entry.get())
            self.ds.setPos(dsPos)
            self.currentDsPos_label_string.set(f'Current dsPos = {dsPos}')
        self.dsPos2Entry = tk.Entry(master=self.controls,width=entryWidth)
        self.dsPos2Button = tk.Button(
            master=self.controls,
            text="Move dsPos2",
            width=18,
            height=1,
            fg='black', bg='white',
            command=dsPos2ButtonFunc)
            
        # Adding the data plot controls
        def newPointButtonFunc():
            self.resetHistory = True
        def clearDataButtonFunc():
            self.resetHistory = True
            self.resetLongHistory = True
            
        self.newPointButton = tk.Button(
            master=self.controls,
            text="New Point",
            width=18,
            height=1,
            fg='black', bg='white',
            command=newPointButtonFunc)
        self.clearDataButton = tk.Button(
            master=self.controls,
            text="Clear Data",
            width=18,
            height=1,
            fg='black', bg='white',
            command=clearDataButtonFunc)
        
        # Adding the scan controls
        scanDsPos_label = tk.Label(master=self.controls,text="dsPos Scan List: (start,end,step)",bg=cntrlBg)
        self.scanDsPosEntry = tk.Entry(master=self.controls,width=entryWidth)
        scanBatchSize_label = tk.Label(master=self.controls,text="Batch Size",bg=cntrlBg)
        self.batchSizeEntry = tk.Entry(master=self.controls,width=entryWidth)
        scanDirLoc_label = tk.Label(master=self.controls,text="Scan Directory",bg=cntrlBg)
        self.scanDirLocEntry = tk.Entry(master=self.controls,width=entryWidth)
        varTransCorr = tk.IntVar()
        self.transCorrBox = tk.Checkbutton(master=self.controls, text='Translation Correction',variable=varTransCorr, onvalue=1, offvalue=0,bg=cntrlBg)
        varSmartScan = tk.IntVar()
        self.smartScanBox = tk.Checkbutton(master=self.controls, text='Smart Scan ®',variable=varSmartScan, onvalue=1, offvalue=0,bg=cntrlBg)
        self.startScanButton = tk.Button(
            master=self.controls,
            text="Start Scan",
            width=18,
            height=1,
            fg='black', bg='white',
            command=lambda:print('test'))
        
        # Pack controls
        add_roi_button.pack(pady=padyControls)
        self.saveRoisButton.pack(pady=padyControls)
        self.loadRoisButton.pack(pady=padyControls)
        self.clearRoisButton.pack(pady=padyControls)
        tk.Frame(master=self.controls, bd=100, relief='flat',height=sepHeight,width=sepWidth,bg='black').pack(side='top', pady=padySep)
        self.saveImageButton.pack(pady=padyControls)
        self.gain_label.pack(pady=padyControls)
        self.gainEntry.pack(pady=padyControls)
        self.exposure_label.pack(pady=padyControls)
        self.exposureEntry.pack(pady=padyControls)
        self.camSettingsButton.pack(pady=padyControls)
        tk.Frame(master=self.controls, bd=100, relief='flat',height=sepHeight,width=sepWidth,bg='black').pack(side='top', pady=padySep)
        self.currentDsPos_label.pack()
        dsPos1_label.pack(pady=padyControls)
        self.dsPos1Entry.pack(pady=padyControls)
        self.dsPos1Button.pack(pady=padyControls)
        dsPos2_label.pack(pady=padyControls)
        self.dsPos2Entry.pack(pady=padyControls)
        self.dsPos2Button.pack(pady=padyControls)
        tk.Frame(master=self.controls, bd=100, relief='flat',height=sepHeight,width=sepWidth,bg='black').pack(side='top', pady=padySep)
        self.newPointButton.pack(pady=padyControls)
        self.clearDataButton.pack(pady=padyControls)
        tk.Frame(master=self.controls, bd=100, relief='flat',height=sepHeight,width=sepWidth,bg='black').pack(side='top', pady=padySep)
        scanDsPos_label.pack(pady=padyControls)
        self.scanDsPosEntry.pack(pady=padyControls)
        scanBatchSize_label.pack(pady=padyControls)
        self.batchSizeEntry.pack(pady=padyControls)
        scanDirLoc_label.pack(pady=padyControls)
        self.scanDirLocEntry.pack(pady=padyControls)
        self.transCorrBox.pack(pady=padyControls)
        self.smartScanBox.pack(pady=padyControls)
        self.startScanButton.pack(pady=padyControls)
        
        # Build sliders for adjusting the brightness and contrast of the image
        ax_setVals = [plt.axes([0.15, 0.06, 0.5, 0.02]), plt.axes([0.15, 0.02, 0.5, 0.02])]
        imageMin, imageMax = np.min(self.camImage), np.max(self.camImage)
        slider_vmax = matplotlib.widgets.Slider(ax_setVals[0], r'$v_{max}$', imageMin, imageMax, valinit=imageMax)
        slider_vmin = matplotlib.widgets.Slider(ax_setVals[1], r'$v_{min}$', imageMin, imageMax, valinit=imageMin)
        def sliderUpdateVmax(val):
            if val > self.camImageObj.get_clim()[0]:
                self.camImageObj.set_clim(self.camImageObj.get_clim()[0], val)
            else:
                slider_vmax.set_val(self.camImageObj.get_clim()[0]+1)
        def sliderUpdateVmin(val):
            if val < self.camImageObj.get_clim()[1]:
                self.camImageObj.set_clim(val, self.camImageObj.get_clim()[1])
            else:
                slider_vmin.set_val(self.camImageObj.get_clim()[1]-1)
        slider_vmax.on_changed(sliderUpdateVmax)
        slider_vmin.on_changed(sliderUpdateVmin)
        
        # Make the line profile crosshairs
        self.lineProfile_x = self.w//2
        self.lineProfile_y = self.h//2
        self.lineProfile_v = self.camAx[0].axvline(x=self.lineProfile_x,c='black')
        self.lineProfile_h = self.camAx[0].axhline(y=self.lineProfile_y,c='black')
        
        # Plot the cam image and the line profile crosshairs
        # SCAF
        self.camImageObj = self.camAx[0].imshow(self.camImage)
        self.lineProfilePlot_v, = self.camAx[1].plot(self.camImage[:,int(self.lineProfile_x)],range(self.h))
        self.lineProfilePlot_h, = self.camAx[2].plot(range(self.w),self.camImage[int(self.lineProfile_y),:])
        
        pos0 = self.camAx[0].get_position()
        pos1 = self.camAx[1].get_position()
        pos2 = self.camAx[2].get_position()
        self.camAx[1].set_position([pos1.x0,pos0.y0,pos1.width,pos0.height])
        self.camAx[2].set_position([pos0.x0,pos2.y0,pos0.width,pos2.height])

        # Connect the cam figure to key press events
        self.camFig.canvas.mpl_connect('button_press_event', self.camonclick)
        self.camFig.canvas.mpl_connect('scroll_event', self.camonscroll)
        self.camFig.canvas.mpl_connect('motion_notify_event', self.camonmove)

        # Attach the matplotlib image figures to a tkinter gui window
        self.canvas = FigureCanvasTkAgg(self.camFig, master=self.root)
        # ~ self.toolbar = NavigationToolbar2Tk(self.canvas, self.root)
        # ~ self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.RIGHT, fill=tk.BOTH, expand=1)
        self.label = tk.Label(text="")
        self.label.pack()
        
        # Create the data plot
        self.dataFig, self.dataAx = plt.subplots(2,1,figsize=(5, 7))
        self.liveIntensityPlot, = self.dataAx[0].plot([],[],'go-') # live intensity plot
        self.dataAx[0].set_xlabel('Time (s)')
        self.avgIntensityPlot = self.dataAx[1].errorbar(range(len(self.avgIntensity)),self.avgIntensity,self.stdIntensity,\
                                                        marker='o',color='green',elinewidth=1,capsize=5) # average intensiity plot with error bars at 1 sigma
        self.dataAx[1].set_xlabel('Data Point #')
        
        # Attach the matplotlib data figures to a tkinter gui window
        self.datacanvas = FigureCanvasTkAgg(self.dataFig, master=self.dataFrame)
        # ~ self.datatoolbar = NavigationToolbar2Tk(self.datacanvas, self.dataFrame)
        # ~ self.datatoolbar.update()
        self.datacanvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.datalabel = tk.Label(text="")
        self.datalabel.pack()
        
        def on_closing():
            # ~ if tk.messagebox.askokcancel("Quit", "Do you want to quit?"):
            self.root.quit()
            
        # Activate the camera thread - this will take images on repeat asynchronously
        thr = threading.Thread(target=self.collectCamImageLive, args=(), kwargs={})
        thr.start() # Will run takeImageTestRepeat
        
        self.update()
        self.root.protocol("WM_DELETE_WINDOW", on_closing)
        self.root.mainloop()
    
    def camonclick(self,event):
        ax = self.camAx[0]
        # Find the location of the mouse click and update the crosshair variable to this value
        selection = self.roiSelection.get()
        if event.inaxes == ax:
            if event.button == 1:
                if selection == 'line profile':
                    self.lineProfile_y = event.ydata
                    self.lineProfile_x = event.xdata
                if selection == 'roi rectangle':
                    if self.roiActive:
                        self.roiActive = False
                        self.rm[-1].goLive() # Set the roi that was just placed to live status
                    else:
                        self.roiActive = True
                        self.rm.append(RoiRectangle(event.xdata,event.ydata,0,0,ax))
            if event.button == 2:
                cur_xlim = ax.get_xlim()
                cur_ylim = ax.get_ylim()
                # compute new limits
                xlims = [event.xdata-(cur_xlim[1]-cur_xlim[0])/2, event.xdata+(cur_xlim[1]-cur_xlim[0])/2]
                ylims = [event.ydata-(cur_ylim[1]-cur_ylim[0])/2, event.ydata+(cur_ylim[1]-cur_ylim[0])/2]
                # make sure the limits don't fall outside the image size
                if xlims[0] < 0:
                    xlims[1] += -xlims[0]
                    xlims[0] += -xlims[0]
                if xlims[1] > self.w:
                    xlims[0] += self.w-xlims[1]
                    xlims[1] += self.w-xlims[1]
                if ylims[1] < 0:
                    ylims[0] += -ylims[1]
                    ylims[1] += -ylims[1]
                if ylims[0] > self.h:
                    ylims[1] += self.h-ylims[0]
                    ylims[0] += self.h-ylims[0]
                # set the new limits
                ax.set_xlim(*xlims)
                ax.set_ylim(*ylims)
            if event.button == 3:
                for i in range(len(self.rm)):
                    if abs(event.xdata-self.rm[i].cx) < self.rm[i].w/2 and abs(event.ydata-self.rm[i].cy) < self.rm[i].h/2:
                        self.rm[i].patch.remove() # Remove the patch from the plot
                        del self.rm[i].patch # Delete the patch
                        self.rm.pop(i) # Pop this item from the rm list
                        break # Exits the loop on the first occurance of a ROI to remove.
                
    def camonscroll(self, event, base_scale = 1.15):
        # ~ print(event.button)
        ax = self.camAx[0]
        # get the current x and y limits
        if event.inaxes == ax:
            cur_xlim = list(ax.get_xlim())
            cur_ylim = list(ax.get_ylim())
            xdata = event.xdata # get event x location
            ydata = event.ydata # get event y location
            if event.button == 'up':
                # deal with zoom in
                scale_factor = base_scale
            elif event.button == 'down':
                # deal with zoom out
                scale_factor = 1/base_scale
            # compute new limits
            xlims = [xdata - (xdata-cur_xlim[0]) / scale_factor, xdata + (cur_xlim[1]-xdata) / scale_factor]
            ylims = [ydata - (ydata-cur_ylim[0]) / scale_factor, ydata + (cur_ylim[1]-ydata) / scale_factor]
            # make sure the limits aren't too big
            if (xlims[1]-xlims[0] > self.w) or (ylims[0]-ylims[1] > self.h):
                xlims = cur_xlim
                ylims = cur_ylim
            # make sure the limits don't fall outside the image size
            if xlims[0] < 0:
                xlims[1] += -xlims[0]
                xlims[0] += -xlims[0]
            if xlims[1] > self.w:
                xlims[0] += self.w-xlims[1]
                xlims[1] += self.w-xlims[1]
            
            if ylims[1] < 0:
                ylims[0] += -ylims[1]
                ylims[1] += -ylims[1]
            if ylims[0] > self.h:
                ylims[1] += self.h-ylims[0]
                ylims[0] += self.h-ylims[0]
            # set the new limits
            ax.set_xlim(*xlims)
            ax.set_ylim(*ylims)
    
    def camonmove(self, event):
        ax = self.camAx[0]
        if event.inaxes == ax:
            if self.roiActive:
                self.rm[-1].w = abs(event.xdata - self.rm[-1].cx)*2
                self.rm[-1].h = abs(event.ydata - self.rm[-1].cy)*2
                tl_x = self.rm[-1].cx - self.rm[-1].w/2
                tl_y = self.rm[-1].cy - self.rm[-1].h/2
                self.rm[-1].patch.set_xy((tl_x,tl_y)) # Update patch top left corner
                self.rm[-1].patch.set_width(self.rm[-1].w) # Update the patch with new width
                self.rm[-1].patch.set_height(self.rm[-1].h) # Update the patch with new height
    
    def saveImage(self):
        now = datetime.datetime.now()
        nowString = now.strftime("%Y_%m_%d-%H_%M_%S")
        defaultFileName = nowString+'expo'+ str(self.exposureEntry.get()) + 'gain' + str(self.exposureEntry.get()) + '.tiff'
        self.saveImageFileName = tk.filedialog.asksaveasfile(mode='w',initialfile=defaultFileName,defaultextension=".*",filetypes = [("image files", ".tiff")])
        if self.saveImageFileName != None:
            tifffile.imwrite(self.saveImageFileName.name,self.camImage)
    
    def collectCamImageLive(self):
        while True:
            deltaTString = '%.3f' % round(time.time()-self.t1,5)
            self.camAx[0].set_title(f"Image Rate =  {deltaTString} s, ROI # = {len(self.rm)}")
            self.t1 = time.time()
            
            self.camImage = self.cam.collectImage() # Take a new image from the camera
            
            self.camImageObj.set_data(self.camImage)
            
            if len(self.rm) > 0:
                roiPixelSum = 0
                roiTotalArea = 0
                # Loop through all the ROIS
                for r in self.rm:
                    # Only count an ROI if it is active
                    if r.live:
                        roiImage = self.camImage[int(r.cy-r.h/2):int(r.cy+r.h/2),int(r.cx-r.w/2):int(r.cx+r.w/2)]
                        roiPixelSum += np.sum(roiImage)
                        roiTotalArea += roiImage.shape[0]*roiImage.shape[1]
                        
                if roiTotalArea > 0:
                    self.timeHistory.append(self.t1 - self.t0)
                    self.intensityHistory.append(roiPixelSum/roiTotalArea)
                    self.avgIntensity[-1] = np.mean(self.intensityHistory)
                    self.stdIntensity[-1] = np.std(self.intensityHistory)/np.sqrt(len(self.intensityHistory))
                    
                    self.dataAx[0].set_xlim([min(self.timeHistory)-0.1,max(self.timeHistory)+0.1])
                    self.dataAx[0].set_ylim([min(self.intensityHistory)-0.1,max(self.intensityHistory)+0.1])
                    self.dataAx[1].set_xlim([-1,len(self.avgIntensity)])
                    self.dataAx[1].set_ylim([min(self.avgIntensity-self.stdIntensity)-0.1,max(self.avgIntensity+self.stdIntensity)+0.1])
    
            self.liveIntensityPlot.set_xdata(self.timeHistory)
            self.liveIntensityPlot.set_ydata(self.intensityHistory)
                
            update_errorbar(self.avgIntensityPlot,range(len(self.avgIntensity)),self.avgIntensity,yerr=self.stdIntensity)
            
            if self.resetHistory:
                self.resetHistory = False
                self.intensityHistory = []
                self.timeHistory = []
                self.avgIntensity = np.append(self.avgIntensity,[0])
                self.stdIntensity = np.append(self.stdIntensity,[0])
                self.t0 = time.time()
            if self.resetLongHistory:
                self.resetLongHistory = False
                self.avgIntensity = np.array([0],dtype='float')
                self.stdIntensity = np.array([0],dtype='float')

    def update(self):
        self.lineProfile_v.set_xdata(self.lineProfile_x)
        self.lineProfile_h.set_ydata(self.lineProfile_y)
        self.lineProfilePlot_v.set_xdata(self.camImage[:,int(self.lineProfile_x)])
        self.camAx[1].set_xlim(min(self.camImage[:,int(self.lineProfile_x)]),max(self.camImage[:,int(self.lineProfile_x)]))
        self.lineProfilePlot_h.set_ydata(self.camImage[int(self.lineProfile_y),:])
        self.camAx[2].set_ylim(min(self.camImage[int(self.lineProfile_y),:]),max(self.camImage[int(self.lineProfile_y),:]))
        self.canvas.draw()
        self.datacanvas.draw()
        self.root.after(self.wait, self.update)
        
        

app = SetupApp()
