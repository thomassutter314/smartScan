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
import random as rand
import matplotlib
matplotlib.use('TkAgg')
# ~ from matplotlib.ticker import FormatStrFormatter
from matplotlib.backends.backend_tkagg import FigureCanvasTk, NavigationToolbar2Tk
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
# ~ import matplotlib.animation as animation
# ~ from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk)
import numpy as np
import time
import datetime
import threading

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
        
        self.cRadius = 1
        # ~ self.cpatch = patches.Circle((self.cx,self.cy),self.cRadius, edgecolor=(1, 0, 0, 0.5), facecolor=(1, 0, 0, 1),fill=False)
        
        ax.add_patch(self.patch) # Add the patch to the axis
        # ~ ax.add_patch(self.cpatch) # Add the center patch to the axis

class SetupApp():
    def __init__(self, wait = .033):
        self.rm = [] # list of rois
        self.roiActive = False
        self.wait = int(1000*wait) # Convert to ms
        self.root = tk.Tk()
        self.root.geometry("1400x800")
        self.root.configure(bg='white')
        self.root.wm_title("Smart Scan")
        
        # Create intensity history arrays for data plots
        self.resetHistory = False
        self.resetLongHistory = False
        self.intensityHistory = [0]
        self.avgIntensity = np.array([0])
        self.stdIntensity = np.array([0])
        self.timeHistory = [0]
        self.t0 = time.time()
        
        # Make a Tkinter frame for the data plots
        self.dataFrame = tk.Frame(master=self.root, width=100, height=100, padx=15, pady=15)
        self.dataFrame.pack(side=tk.RIGHT)
        
        # Make a Tkinter frame for the ROI selection buttons
        self.controls = tk.Frame(master=self.root, width=100, height=100, padx=15, pady=15)
        self.controls.pack(side=tk.RIGHT)

        # Adding the ROI controls
        self.roiSelection = tk.StringVar(self.controls)
        self.roiSelection.set("line profile") # default value
        add_roi_button = tk.OptionMenu(self.controls, self.roiSelection, "line profile", "roi rectangle")
        
        # Adding the camera controls
        self.saveImageButton = tk.Button(
            master=self.controls,
            text="Save Cam Image",
            width=18,
            height=2,
            fg='black', bg='white',
            command=lambda:print('test'))
        gain_label = tk.Label(master=self.controls,text="Gain")
        self.gainEntry = tk.Entry(master=self.controls)
        exposure_label = tk.Label(master=self.controls,text="Exposure (s)")
        self.exposureEntry = tk.Entry(master=self.controls)
        self.camSettingsButton = tk.Button(
            master=self.controls,
            text="Update Cam Settings",
            width=18,
            height=2,
            fg='black', bg='white',
            command=lambda:print('test'))
        
        # Adding the ds controls
        dsPos1_label = tk.Label(master=self.controls,text="dsPos1 (mm)")
        self.dsPos1Entry = tk.Entry(master=self.controls)
        self.dsPos1Button = tk.Button(
            master=self.controls,
            text="Move dsPos1",
            width=18,
            height=2,
            fg='black', bg='white',
            command=lambda:print('test'))
        dsPos2_label = tk.Label(master=self.controls,text="dsPos2 (mm)")
        self.dsPos2Entry = tk.Entry(master=self.controls)
        self.dsPos2Button = tk.Button(
            master=self.controls,
            text="Move dsPos2",
            width=18,
            height=2,
            fg='black', bg='white',
            command=lambda:print('test'))
            
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
            height=2,
            fg='black', bg='white',
            command=newPointButtonFunc)
        self.clearDataButton = tk.Button(
            master=self.controls,
            text="Clear Data",
            width=18,
            height=2,
            fg='black', bg='white',
            command=clearDataButtonFunc)
            
        # Pack controls
        add_roi_button.pack()
        self.saveImageButton.pack()
        gain_label.pack()
        self.gainEntry.pack()
        exposure_label.pack()
        self.exposureEntry.pack()
        self.camSettingsButton.pack()
        dsPos1_label.pack()
        self.dsPos1Entry.pack()
        self.dsPos1Button.pack()
        dsPos2_label.pack()
        self.dsPos2Entry.pack()
        self.dsPos2Button.pack()
        self.newPointButton.pack()
        self.clearDataButton.pack()

        # Make the cam plot window for the diffraction image
        # SCAF
        self.takeImageTest()
        self.h, self.w = self.image_0.shape
        gs = gridspec.GridSpec(2, 2,width_ratios=[self.w,self.w*.1], height_ratios=[self.h,self.h*.1])
        self.camFig = plt.figure()
        self.camAx = [plt.subplot(gs[0]),]
        self.camAx.append(plt.subplot(gs[1],sharey=self.camAx[0]))
        self.camAx.append(plt.subplot(gs[2],sharex=self.camAx[0]))
        self.camAx[0].xaxis.set_tick_params(labelbottom=False)
        self.camAx[0].yaxis.set_tick_params(labelleft=False)
        self.camAx[0].set_xticks([])
        self.camAx[0].set_yticks([])
        
        # Make the line profile crosshairs
        self.lineProfile_x = self.w//2
        self.lineProfile_y = self.h//2
        self.lineProfile_v = self.camAx[0].axvline(x=self.lineProfile_x,c='black')
        self.lineProfile_h = self.camAx[0].axhline(y=self.lineProfile_y,c='black')
        
        # Plot the cam image and the line profile crosshairs
        # SCAF
        self.camImage = self.camAx[0].imshow(self.image_0)
        self.lineProfilePlot_v, = self.camAx[1].plot(self.image_0[:,int(self.lineProfile_x)],range(self.h))
        self.lineProfilePlot_h, = self.camAx[2].plot(range(self.w),self.image_0[int(self.lineProfile_y),:])
        
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
        self.liveIntensityPlot, = self.dataAx[0].plot([],[],'ro-') # live intensity plot
        self.avgIntensityPlot = self.dataAx[1].errorbar(range(len(self.avgIntensity)),self.avgIntensity,self.stdIntensity,\
                                                        marker='o',elinewidth=1,capsize=5) # average intensiity plot with error bars at 1 sigma
        print("self.avgIntensityPlot",self.avgIntensityPlot)
        
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
        thr = threading.Thread(target=self.takeImageTestRepeat, args=(), kwargs={})
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
    
    def takeImageTest(self):
        x, y = np.linspace(-3,3,1000), np.linspace(-3,3,800)
        X, Y = np.meshgrid(x,y)
        time.sleep(1)
        self.image_0 = np.exp(-(X**2/4+Y**2)/.2)*np.cos((X**2+Y**2)*10)**2
        print('ok')
        
    def takeImageTestRepeat(self):
        while True:
            time.sleep(.1)
                        
            self.timeHistory.append(time.time() - self.t0)
            self.image_0 = (np.random.random(np.shape(self.image_0)))
            self.camImage.set_data(self.image_0)
            
            roiPixelSum = 0
            for r in self.rm:
                roiPixelSum += np.sum(self.image_0[int(r.cy-r.h/2):int(r.cy+r.h/2),int(r.cx-r.w/2):int(r.cx+r.w/2)])
                
            self.intensityHistory.append(roiPixelSum)
            self.avgIntensity[-1] = np.mean(self.intensityHistory)
            self.stdIntensity[-1] = np.std(self.intensityHistory)/np.sqrt(len(self.intensityHistory))
            
            self.dataAx[0].set_xlim([min(self.timeHistory),max(self.timeHistory)])
            self.dataAx[0].set_ylim([min(self.intensityHistory)-0.1,max(self.intensityHistory)+0.1])
            
            self.dataAx[1].set_xlim([-1,len(self.avgIntensity)])
            self.dataAx[1].set_ylim([min(self.avgIntensity-self.stdIntensity)-0.1,max(self.avgIntensity+self.stdIntensity)+0.1])

            self.liveIntensityPlot.set_xdata(self.timeHistory)
            self.liveIntensityPlot.set_ydata(self.intensityHistory)
            
            update_errorbar(self.avgIntensityPlot,range(len(self.avgIntensity)),self.avgIntensity,yerr=self.stdIntensity)
            print('New image')
            
            if self.resetHistory:
                self.resetHistory = False
                self.intensityHistory = []
                self.timeHistory = []
                self.avgIntensity = np.append(self.avgIntensity,[0])
                self.stdIntensity = np.append(self.stdIntensity,[0])
                self.t0 = time.time()
            if self.resetLongHistory:
                self.resetLongHistory = False
                self.avgIntensity = np.array([0])
                self.stdIntensity = np.array([0])

                
            
    def update(self):
        self.lineProfile_v.set_xdata(self.lineProfile_x)
        self.lineProfile_h.set_ydata(self.lineProfile_y)
        self.lineProfilePlot_v.set_xdata(self.image_0[:,int(self.lineProfile_x)])
        self.camAx[1].set_xlim(min(self.image_0[:,int(self.lineProfile_x)]),max(self.image_0[:,int(self.lineProfile_x)]))
        self.lineProfilePlot_h.set_ydata(self.image_0[int(self.lineProfile_y),:])
        self.camAx[2].set_ylim(min(self.image_0[int(self.lineProfile_y),:]),max(self.image_0[int(self.lineProfile_y),:]))
        self.canvas.draw()
        self.datacanvas.draw()
        self.root.after(self.wait, self.update)
        
        

app = SetupApp()
