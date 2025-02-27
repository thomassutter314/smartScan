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
# ~ import random as rand
import matplotlib
# ~ matplotlib.use('TkAgg')
# ~ from matplotlib.ticker import FormatStrFormatter
from matplotlib.backends.backend_tkagg import FigureCanvasTk, NavigationToolbar2Tk
import matplotlib.pyplot as plt
import matplotlib.patches as patches
# ~ import matplotlib.gridspec as gridspec
import matplotlib.widgets
# ~ import matplotlib.animation as animation
# ~ from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk)
    
# Use a matplotlib backend that doesn't have an associated gui
plt.switch_backend('agg')

import numpy as np
import scipy
from scipy import fft
from scipy.ndimage import median_filter
from scipy import optimize as opt
import time
import datetime
import threading
import tifffile
import cv2

import json
import csv
from itertools import zip_longest
import os
import io
import sys
import psutil
osprocess = psutil.Process(os.getpid())

import translationcorr as tc
import lab_instruments
import data_analysis


def parsePosStr(s):
    dsPositions = np.array([])
    
    while len(s) > 0:
        i1 = 1 + s.find('(')
        i2 = s[i1:].find(')')
        j1 = 1 + s.find('[')
        j2 = s[j1:].find(']')
        
        if i1 != 0 and (i1 < j1 or j1 == 0):
            # ~ print(s)
            # ~ print(s[i1:i1+i2].split(','))
            start, stop, step = list(map(float,s[i1:i1+i2].split(',')))
            newPosArr = np.arange(start,stop+step,step)
            dsPositions = np.append(dsPositions,newPosArr)
            s = s[i1+i2+1:]
            # ~ print(s)
        if j1 != 0 and (j1 < i1 or i1 == 0):
            newPosArr = list(map(float,s[j1:j1+j2].split(',')))
            dsPositions = np.append(dsPositions,newPosArr)
            s = s[j1+j2+1:]
        
    dsPositions.sort()
    dsPositions_unique = np.unique(dsPositions.round(4)) # rounds the values to 4 decimal places and removes duplicate positions
    
    return dsPositions_unique
        
def popupmsg(msg):
    popup = tk.Tk()
    popup.wm_title("!")
    label = ttk.Label(popup, text=msg)
    label.pack(side="top", fill="x", pady=10)
    B1 = ttk.Button(popup, text="Okay", command = popup.destroy)
    B1.pack()
    
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

def incrementBatchMetaDataWeight(fileDir, pi, increment = 1):
    with open(fileDir, 'r') as f:
        metaData = f.readlines()
    
    idString_start = f'pi = '
    idString_finish = f', dsPos'
    weightString = f'weight = '
    for l in range(len(metaData)):
        if idString_start in metaData[l]:
            index_start = metaData[l].find(idString_start) + len(idString_start)
            index_finish = metaData[l].find(idString_finish)
            if pi == int(metaData[l][index_start:index_finish]):
                indexWeight = metaData[l].find(weightString) + len(weightString)
                oldWeight = int(metaData[l][indexWeight:])
                newLine = metaData[l][:indexWeight] + str(oldWeight + increment) + '\n'
                metaData[l] = newLine # Update the list storing the lines
                break # exit the for loop after the correct position index has been found
                
    # Write the lines back to the text file
    with open(fileDir, 'w') as f:
        f.writelines(metaData)

def readBatchMetaDataWeight(fileDir, pi):
    with open(fileDir, 'r') as f:
        metaData = f.readlines()

    idString_start = f'pi = '
    idString_finish = f', dsPos'
    weightString = f'weight = '
    for l in range(len(metaData)):
        if idString_start in metaData[l]:
            index_start = metaData[l].find(idString_start) + len(idString_start)
            index_finish = metaData[l].find(idString_finish)
            if pi == int(metaData[l][index_start:index_finish]):
                indexWeight = metaData[l].find(weightString) + len(weightString)
                weight = int(metaData[l][indexWeight:])
                return weight

def incrementAverageMetaDataWeight(fileDir, pi, fi, increment = 1):
    # function runs on the average metadata too
    with open(fileDir, 'r') as f:
        metaData = f.readlines()

    idString_1_start = f'pi = '
    idString_1_finish = f', dsPos'
    idString_2_start = f'fi = '
    idString_2_finish = f', fluence'
    weightString = f'weight = '
    for l in range(len(metaData)):
        if weightString in metaData[l]:
            index_1_start = metaData[l].find(idString_1_start) + len(idString_1_start)
            index_1_finish = metaData[l].find(idString_1_finish)
            index_2_start = metaData[l].find(idString_2_start) + len(idString_2_start)
            index_2_finish = metaData[l].find(idString_2_finish)
            if pi == int(metaData[l][index_1_start:index_1_finish]) and fi == int(metaData[l][index_2_start:index_2_finish]):
                indexWeight = metaData[l].find(weightString) + len(weightString)
                oldWeight = int(metaData[l][indexWeight:])
                newLine = metaData[l][:indexWeight] + str(oldWeight + increment) + '\n'
                metaData[l] = newLine # Update the list storing the lines
                break # exit the for loop after the correct position index has been found
                
    # Write the lines back to the text file
    with open(fileDir, 'w') as fileObj:
        fileObj.writelines(metaData)

def readAverageMetaDataWeight(fileDir, pi, fi):
    # function runs on the average metadata too
    with open(fileDir, 'r') as f:
        metaData = f.readlines()

    idString_1_start = f'pi = '
    idString_1_finish = f', dsPos'
    idString_2_start = f'fi = '
    idString_2_finish = f', fluence'
    weightString = f'weight = '
    for l in range(len(metaData)):
        if weightString in metaData[l]:
            index_1_start = metaData[l].find(idString_1_start) + len(idString_1_start)
            index_1_finish = metaData[l].find(idString_1_finish)
            index_2_start = metaData[l].find(idString_2_start) + len(idString_2_start)
            index_2_finish = metaData[l].find(idString_2_finish)
            if pi == int(metaData[l][index_1_start:index_1_finish]) and fi == int(metaData[l][index_2_start:index_2_finish]):
                indexWeight = metaData[l].find(weightString) + len(weightString)
                weight = int(metaData[l][indexWeight:])
                return weight

 
class RoiRectangle():
    def __init__(self,cx,cy,w,h,ax):
        self.cx = cx
        self.cy = cy
        self.w = w
        self.h = h
        tl_x = self.cx - w/2
        tl_y = self.cy - h/2
        self.patch = patches.Rectangle((tl_x, tl_y), self.w, self.h, linewidth=1, edgecolor='r', facecolor='none', alpha=0.8)
        self.live = False
        
        self.cRadius = 1
        # ~ self.cpatch = patches.Circle((self.cx,self.cy),self.cRadius, edgecolor=(1, 0, 0, 0.5), facecolor=(1, 0, 0, 1),fill=False)
        
        ax.add_patch(self.patch) # Add the patch to the axis
        # ~ ax.add_patch(self.cpatch) # Add the center patch to the axis
    
    def updateShape(self, w, h,):
        self.w = w
        self.h = h
        tl_x = self.cx - w/2
        tl_y = self.cy - h/2
        self.patch.set_xy((tl_x,tl_y)) # Update patch top left corner
        self.patch.set_width(w) # Update the patch with new width
        self.patch.set_height(h) # Update the patch with new height
    
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
        
    def updateAxis(self, ax):
        del self.patch
        tl_x = self.cx - self.w/2
        tl_y = self.cy - self.h/2
        self.patch = patches.Rectangle((tl_x, tl_y), self.w, self.h, linewidth=1, edgecolor='r', facecolor='none', alpha=0.8)
        ax.add_patch(self.patch) # Add the patch to the axis

class SetupApp():
    def __init__(self, wait = .033):
        # Load the persistent_settings
        persistent_settings_dict = {}
        with open('persistent_settings.txt','r') as f:
            for line in f:
                if line[:line.find(":")] != '#':
                    persistent_settings_dict[line[:line.find(":")]] = line[line.find(":")+1:].replace('\n','')
        
        # Define the pump malus curve
        self.imalus = lambda x: data_analysis.malusFunc(x, float(persistent_settings_dict['pump_malus_A']), float(persistent_settings_dict['pump_malus_phi']), float(persistent_settings_dict['pump_malus_offset']),invert=True)
        self.malus = lambda x: data_analysis.malusFunc(x, float(persistent_settings_dict['pump_malus_A']), float(persistent_settings_dict['pump_malus_phi']), float(persistent_settings_dict['pump_malus_offset']),invert=False)
        
        # Setup the tkinter GUI window
        self.root = tk.Tk()
        self.root.geometry("1000x700")
        self.root.configure(bg='white')
        self.root.wm_title("Smart Scan")
        self.root.iconbitmap("kogar_0.ico")
        #self.root.iconify()
        

        self.camera = lab_instruments.BlackFlyCamera()
        self.gain = 30
        self.camera.setGain(self.gain)
            
            
        self.exposure = 1
        self.camera.setExposure(self.exposure)
        
        time.sleep(0.5)
        self.camSettingsUpdateInProgress = False
        
        # take a sample image
        self.camera.trigger()
        self.camImage = self.camera.grabImage()
        self.displayImage = self.camImage # This can be different than the direct cam image if filters or averaging is applied.
        self.h, self.w = self.camImage.shape
        
        # ~ self.camFig = plt.figure()
        self.camFig = plt.figure(figsize=(12,10))
        # ~ self.camFig = plt.figure(figsize=(8,10))
        self.camAx = [self.camFig.add_gridspec(top=0.75, right=0.9, hspace = 0).subplots()] 
        # ~ self.camAx = [self.camFig.add_gridspec.subplots()]

        self.camAx.append(self.camAx[0].inset_axes([0, 1.05, 1, 0.25], sharex=self.camAx[0]))
        
        # Connect to the delay stage
        self.ds = lab_instruments.DelayStage(com = 'COM3')
        
        # Connect to the half wave plate
        self.hwp = lab_instruments.HalfWavePlate(com = 'COM4') # COM4 is pump, COM5 is heaterline
        
        # Data display settings and variables
        self.wait = wait
        self.ppw = 100
        self.rm = [] # list of rois
        self.roiActive = False
        
        # Some general formatting vars
        padyControls = 1
        padySep = 1
        entryWidth = 10
        largeEntryWidth = 18
        sepWidth = 150
        sepHeight = 1
        buttonWidth = 16
        cntrlBg = '#2F8D35'
        
        # Create intensity history arrays for data plots
        self.resetHistory = False
        self.resetLongHistory = False
        self.intensityHistory = []
        self.avgIntensity = np.array([0],dtype='float')
        self.stdIntensity = np.array([0],dtype='float')
        self.timeHistory = []
        self.t0 = time.time()
        self.t1 = self.t0
        
        # global variables for mouse location
        self.mouse_x = 0
        self.mouse_y = 0

        # Make a Tkinter frame for the detector
        self.detectorFrame = tk.Frame(master=self.root, width=0, height=0, padx=0, pady=0, bg=cntrlBg)
        self.detectorFrame.pack(side=tk.LEFT)
        
        # Make a Tkinter frame for the ROI selection buttons
        self.controls = tk.Frame(master=self.root, width=0, height=0, padx=0, pady=0, bg=cntrlBg)
        self.controls.pack(side=tk.LEFT)
        
        # Make a Tkinter frame for the scan controls
        self.scanFrame = tk.Frame(master=self.root, width=0, height=0, padx=0, pady=0, bg=cntrlBg)
        self.scanFrame.pack(side=tk.LEFT)
        
        # Make a Tkinter frame for the data plots
        self.dataFrame = tk.Frame(master=self.root, width=0, height=0, padx=0, pady=0)
        self.dataFrame.pack(side=tk.LEFT)
                
        # This is the GUI controls section
        # Adding the ROI controls
        self.roiSelection = tk.StringVar(self.controls)
        self.roiSelection.set("roi rectangle") # default value
        add_roi_button = tk.OptionMenu(self.controls, self.roiSelection, "line profile", "roi rectangle")
        
        # Adding a label for the total number of ROI's
        self.roi_number_label_string = tk.StringVar()
        self.roi_number_label_string.set('ROI # = 0')
        self.roi_number_label = tk.Label(master=self.controls,textvariable=self.roi_number_label_string, bg=cntrlBg)
        
        # Adding labels for the x, y, and intensity at the mouse positions
        self.cursor_xycoord_label_string = tk.StringVar()
        self.cursor_xycoord_label_string.set('(x, y) = (0, 0)')
        self.cursor_xycoord_label = tk.Label(master=self.controls,textvariable=self.cursor_xycoord_label_string, bg=cntrlBg)
        
        self.cursor_intensity_label_string = tk.StringVar()
        self.cursor_intensity_label_string.set('pixel value = 0')
        self.cursor_intensity_label = tk.Label(master=self.controls,textvariable=self.cursor_intensity_label_string, bg=cntrlBg)

        # Adding the camera controls
        self.saveImageButton = tk.Button(
            master=self.controls,
            text="Save Cam Image",
            width=buttonWidth,
            height=1,
            fg='black', bg='white',
            command=self.saveImage)
        self.gain_label_string = tk.StringVar()
        cgs = '%.2f' % self.camera.getGain()
        self.gain_label_string.set(f'Gain = {cgs} dB')
        self.gain_label = tk.Label(master=self.controls,textvariable=self.gain_label_string,bg=cntrlBg)
        self.gainEntry = tk.Entry(master=self.controls,width=entryWidth)
        self.gainEntry.insert(-1,self.camera.getGain())
        self.exposure_label_string = tk.StringVar()
        ces = '%.2f' % self.camera.getExposure()
        self.exposure_label_string.set(f'Exposure = {ces} s')
        self.exposure_label = tk.Label(master=self.controls,textvariable=self.exposure_label_string,bg=cntrlBg)
        self.exposureEntry = tk.Entry(master=self.controls,width=entryWidth)
        self.exposureEntry.insert(-1,self.camera.getExposure())
        
        def camSettingsButtonFunc():
            if self.gainEntry.get() != '':
                self.gain = float(self.gainEntry.get())
            else:
                self.gain = self.camera.getGain()
                cgs = '%.2f' % self.gain
                self.gainEntry.insert(-1,cgs)
                
            if self.exposureEntry.get() != '':
                self.exposure = float(self.exposureEntry.get())
            else:
                self.exposure = self.camera.getExposure()
                ces = '%.2f' % self.exposure
                self.exposureEntry.insert(-1,ces)
            
            # If the update is substantive, then write the new values to the GUI and set a ticket for the physical camera settings to be updated at next opportunity.
            if (abs(self.gain - self.camera.getGain()) > .001) or (abs(self.exposure - self.camera.getExposure()) > .001):
                print('Updating Camera Settings')
                cgs = '%.2f' % self.gain
                self.gain_label_string.set(f'Gain = {cgs} dB')
                ces = '%.2f' % self.exposure
                self.exposure_label_string.set(f'Exposure = {ces} s')
                self.camSettingsUpdateInProgress = True
                    
        # Adding the ds controls
        self.currentDsPos_label_string = tk.StringVar()
        dps = '%.2f' % self.ds.getPos()
        self.currentDsPos_label_string.set(f'dsPos = {dps} mm')
        self.currentDsPos_label = tk.Label(master=self.controls,textvariable=self.currentDsPos_label_string,bg=cntrlBg)
        
        dsPos1_label = tk.Label(master=self.controls,text="dsPos1 (mm)",bg=cntrlBg)
        self.dsPos1Entry = tk.Entry(master=self.controls,width=entryWidth)
        
        dsPos2_label = tk.Label(master=self.controls,text="dsPos2 (mm)",bg=cntrlBg)
        self.dsPos2Entry = tk.Entry(master=self.controls,width=entryWidth)
        
        def dsPos1ButtonFunc():
            dsPos = float(self.dsPos1Entry.get())
            self.ds.setPos(dsPos)
            dps = '%.2f' % self.ds.getPos()
            self.currentDsPos_label_string.set(f'dsPos = {dps}')
            
        def dsPos2ButtonFunc():
            dsPos = float(self.dsPos2Entry.get())
            self.ds.setPos(dsPos)
            dps = '%.2f' % self.ds.getPos()
            self.currentDsPos_label_string.set(f'dsPos = {dps}')
        

        # Adding the hwp pump controls
        self.currentFluence_label_string = tk.StringVar()
        currentHwpAngle = self.hwp.getPos()
        cfs = '%.2f' % self.malus(currentHwpAngle)
        self.currentFluence_label_string.set(f'Fluence = {cfs} Ko')
        self.currentFluence_label = tk.Label(master=self.controls,textvariable=self.currentFluence_label_string,bg=cntrlBg)
        
        self.currentHwpPos_label_string = tk.StringVar()
        chs = '%.2f' % currentHwpAngle
        self.currentHwpPos_label_string.set(f'HWP Angle = {chs} deg')
        self.currentHwpPos_label = tk.Label(master=self.controls,textvariable=self.currentHwpPos_label_string,bg=cntrlBg)
        
        fluence1_label = tk.Label(master=self.controls,text="Fluence (Ko)",bg=cntrlBg)
        self.fluence1Entry = tk.Entry(master=self.controls,width=entryWidth)
        
        def fluence1ButtonFunc():
            fluence = float(self.fluence1Entry.get())
            hwp_angle = self.imalus(fluence)
            
            # Check to make sure this command makes sense
            if hwp_angle == None:
                print('HWP motion rejected: input fluence outside possible range')
                return
                     
            if hwp_angle < 0:
                print('HWP motion rejected: angle less than zero')
                return
                
            if fluence < 0:
                print('HWP motion rejected: fluence less than zero')
                return
                
            # Move the hwp to the new fluence
            self.hwp.moveAbsolute(hwp_angle)
            currentHwpAngle = self.hwp.getPos()
            cfs = '%.2f' % self.malus(currentHwpAngle)
            chs = '%.2f' % currentHwpAngle
            self.currentFluence_label_string.set(f'Fluence = {cfs} Ko')
            self.currentHwpPos_label_string.set(f'HWP Angle = {chs} deg')

        hwp1_label = tk.Label(master=self.controls,text="HWP Angle (deg)",bg=cntrlBg)
        self.hwp1Entry = tk.Entry(master=self.controls,width=entryWidth)
        
        def hwp1ButtonFunc():
            hwp_angle = float(self.hwp1Entry.get())
            
            # Check to make sure this command makes sense
            if hwp_angle < 0:
                print('HWP motion rejected: angle less than zero')
                return
                
            if hwp_angle >= 45:
                print('HWP motion rejected: angle greater than or equal to 45')
                return
                
            # Move the hwp to the new value
            self.hwp.moveAbsolute(hwp_angle)
            currentHwpAngle = self.hwp.getPos()
            cfs = '%.2f' % self.malus(currentHwpAngle)
            chs = '%.2f' % currentHwpAngle
            self.currentFluence_label_string.set(f'Fluence = {cfs} Ko')
            self.currentHwpPos_label_string.set(f'HWP Angle = {chs} deg')

        # Adding the data plot controls
        def newPointButtonFunc():
            if self.resetHistory == False:
                self.resetHistory = True
                self.newPointButton.config(text='Continue',bg='yellow')
            else:
                self.resetHistory = False
                self.newPointButton.config(text='New Point',bg='white')
            
        def clearDataButtonFunc():
            self.resetLongHistory = True
        
        def saveDataButtonFunc():
            defaultFileName = 'data.csv'
            fileLoc = tk.filedialog.asksaveasfile(mode='w',initialfile=defaultFileName,defaultextension=".*",filetypes = [("csv files", ".csv")])
            if fileLoc == None:
                return None
            # Your list of lists
            listOlists = [['timeHistory'] + list(self.timeHistory),
                          ['intensityHistory'] + list(self.intensityHistory),
                          ['avgIntensity'] + list(self.avgIntensity),
                          ['stdIntensity'] + list(self.stdIntensity)]
            
            # Determine the maximum length of the sublists
            max_length = max(len(sublist) for sublist in listOlists)
            
            # Fill shorter sublists with placeholders (empty strings)
            listOlists = [sublist + [''] * (max_length - len(sublist)) for sublist in listOlists]
            
            # Transpose the list of lists to get separate columns
            transposed_list = list(map(list, zip(*listOlists)))
            
            # Open the CSV file in write mode
            with open(fileLoc.name, mode='w', newline='') as file1:
                # Create the CSV writer object
                writer = csv.writer(file1, delimiter = ",", quotechar = '"')
                
                # Write each row (column) from the transposed list
                for row in transposed_list:
                    writer.writerow(row)
                    
        self.newPointButton = tk.Button(
            master=self.controls,
            text="New Point",
            width=buttonWidth,
            height=1,
            fg='black', bg='white',
            command=newPointButtonFunc)
        self.clearDataButton = tk.Button(
            master=self.controls,
            text="Clear Data",
            width=buttonWidth,
            height=1,
            fg='black', bg='white',
            command=clearDataButtonFunc)
        self.saveDataButton = tk.Button(
            master=self.controls,
            text="Save Data",
            width=buttonWidth,
            height=1,
            fg='black', bg='white',
            command=saveDataButtonFunc)
        
        # Adding the scan controls
        scanSweepSelection_label = tk.Label(master=self.scanFrame,text='Scan Sweep Selection',bg=cntrlBg)
        self.scanSweepSelection = tk.StringVar(self.controls)
        self.scanSweepSelection.set(persistent_settings_dict['scanSweepSelection']) # default value
        scanSweepSelection_button = tk.OptionMenu(self.scanFrame, self.scanSweepSelection, "one direction", "back and forth", "quasi random", 'full random')
        scanDsPos_label = tk.Label(master=self.scanFrame,text="dsPos Scan List (mm):\n(start,end,step)",bg=cntrlBg)
        self.scanDsPosEntry = tk.Entry(master=self.scanFrame,width=largeEntryWidth)
        self.scanDsPosEntry.insert(-1,persistent_settings_dict['dsPositions'])
        scanDsWait_label = tk.Label(master=self.scanFrame,text="ds Wait (s)",bg=cntrlBg)
        self.scanDsWaitEntry = tk.Entry(master=self.scanFrame,width=largeEntryWidth)
        self.scanDsWaitEntry.insert(-1,persistent_settings_dict['dsWait'])

        scanFluence_label = tk.Label(master=self.scanFrame,text="Fluence Scan List (Ko):\n(start,end,step)",bg=cntrlBg)
        self.scanFluenceEntry = tk.Entry(master=self.scanFrame,width=largeEntryWidth)
        self.scanFluenceEntry.insert(-1,persistent_settings_dict['fluences'])
        scanHwpWait_label = tk.Label(master=self.scanFrame,text="HWP Wait (s)",bg=cntrlBg)
        self.scanHwpWaitEntry = tk.Entry(master=self.scanFrame,width=largeEntryWidth)
        self.scanHwpWaitEntry.insert(-1,persistent_settings_dict['hwpWait'])
        
        scanBatchSize_label = tk.Label(master=self.scanFrame,text="Batch Size",bg=cntrlBg)
        self.batchSizeEntry = tk.Entry(master=self.scanFrame,width=largeEntryWidth)
        self.batchSizeEntry.insert(-1,persistent_settings_dict['batchSize'])
        scanDirLoc_label = tk.Label(master=self.scanFrame,text="Scan Directory",bg=cntrlBg)
        self.scanDirLocEntry = tk.Entry(master=self.scanFrame,width=largeEntryWidth)
        self.varTransCorr = tk.IntVar()
        self.varTransCorr.set(int(persistent_settings_dict['tcq'] == 'True'))
        self.transCorrBox = tk.Checkbutton(master=self.scanFrame, text='Translation\nCorrection',variable=self.varTransCorr, onvalue=1, offvalue=0,bg=cntrlBg)
        self.varScanMedianFilter = tk.IntVar()
        self.varScanMedianFilter.set(int(persistent_settings_dict['smfq'] == 'True'))
        self.scanMedianFilterBox = tk.Checkbutton(master=self.scanFrame, text='Median Filter\nDuring Scan',variable=self.varScanMedianFilter, onvalue=1, offvalue=0,bg=cntrlBg)        
        self.varSmartScan = tk.IntVar()
        self.varSmartScan.set(int(persistent_settings_dict['ssq'] == 'True'))
        self.smartScanBox = tk.Checkbutton(master=self.scanFrame, text='Smart Scan ®',variable=self.varSmartScan, onvalue=1, offvalue=0,bg=cntrlBg)
        
        def startScan():
            scanDir = str(self.scanDirLocEntry.get())
            if not os.path.isdir(scanDir):
                popupmsg('Please select a valid scan directory for storing data')
                return
            if len(self.rm) == 0:
                print('No ROI selections, setting full image as ROI')
                self.rm.append(RoiRectangle(self.w//2,self.h//2,self.w,self.h,self.camAx[0]))
                self.rm[-1].goLive() # Set the roi that was just generated to live status
            self.camLive = 2
            
        self.startScanButton = tk.Button(
            master=self.scanFrame,
            text="Start Scan",
            width=buttonWidth,
            height=1,
            fg='black', bg='white',
            command=startScan)
        
        # Add a image progress bar
        self.imProgBar = ttk.Progressbar(master=self.controls,length=140)
        
        # Image set and averaging controls
        self.numAvgDis = 0
        self.varAvgDis = tk.IntVar()
        self.varAvgDis_label_string = tk.StringVar()
        self.varAvgDis_label_string.set(f'Average Images (N = {self.numAvgDis})')
        self.avgDisBox = tk.Checkbutton(master=self.controls, textvariable=self.varAvgDis_label_string,variable=self.varAvgDis, onvalue=1, offvalue=0,bg=cntrlBg)
        self.varSaveDir = tk.IntVar()
        
        self.avgDirLoc_label = tk.Label(master=self.controls,text='Image Set Directory',bg=cntrlBg)
        self.avgDirLocEntry = tk.Entry(master=self.controls,width=largeEntryWidth)
        
        # Despeckle controls
        self.varMedianFilter = tk.IntVar()
        self.medianFilterBox = tk.Checkbutton(master=self.controls, text='Median Filter (Despeckle)',variable=self.varMedianFilter, onvalue=1, offvalue=0,bg=cntrlBg)
        
        self.saveImageSet = False
        def imageSetButtonFunc():
            if self.saveImageSet == False:
                self.saveImageSet = True
                self.imageSetButton.config(bg='yellow')
            else:
                self.saveImageSet = False
                self.imageSetButton.config(bg='white')
                
        self.imageSetButton = tk.Button(
            master=self.controls,
            text="Collect Image Set",
            width=buttonWidth,
            height=1,
            fg='black', bg='white',
            command=imageSetButtonFunc)
        
        # Pack controls
        add_roi_button.pack(pady=padyControls)
        self.roi_number_label.pack(pady=padyControls)
        self.cursor_xycoord_label.pack(pady=padyControls)
        self.cursor_intensity_label.pack(pady=padyControls)
        tk.Frame(master=self.controls, bd=100, relief='flat',height=sepHeight,width=sepWidth,bg='black').pack(side='top', pady=padySep)
        self.gain_label.pack(pady=padyControls)
        self.gainEntry.pack(pady=padyControls)
        self.exposure_label.pack(pady=padyControls)
        self.exposureEntry.pack(pady=padyControls)
        self.saveImageButton.pack(pady=padyControls)
        self.imProgBar.pack(pady=padyControls)
        tk.Frame(master=self.controls, bd=100, relief='flat',height=sepHeight,width=sepWidth,bg='black').pack(side='top', pady=padySep)
        self.currentDsPos_label.pack()
        dsPos1_label.pack(pady=padyControls)
        self.dsPos1Entry.pack(pady=padyControls)
        dsPos2_label.pack(pady=padyControls)
        self.dsPos2Entry.pack(pady=padyControls)
        tk.Frame(master=self.controls, bd=100, relief='flat',height=sepHeight,width=sepWidth,bg='black').pack(side='top', pady=padySep)
        self.currentFluence_label.pack(pady=padyControls)
        self.currentHwpPos_label.pack(pady=padyControls)
        fluence1_label.pack(pady=padyControls)
        self.fluence1Entry.pack(pady=padyControls)
        
        hwp1_label.pack(pady=padyControls)
        self.hwp1Entry.pack(pady=padyControls)
        
        tk.Frame(master=self.controls, bd=100, relief='flat',height=sepHeight,width=sepWidth,bg='black').pack(side='top', pady=padySep)
        self.newPointButton.pack(pady=padyControls)
        self.clearDataButton.pack(pady=padyControls)
        self.saveDataButton.pack(pady=padyControls)
        self.avgDisBox.pack(pady=padyControls)
        self.medianFilterBox.pack(pady=padyControls)
        self.avgDirLoc_label.pack(pady=padyControls)
        self.avgDirLocEntry.pack(pady=padyControls)
        self.imageSetButton.pack(pady=padyControls)
        # ~ tk.Frame(master=self.controls, bd=100, relief='flat',height=sepHeight,width=sepWidth,bg='black').pack(side='top', pady=padySep)
        
        scanSweepSelection_label.pack()
        scanSweepSelection_button.pack(pady=padyControls)
        scanDsPos_label.pack(pady=padyControls)
        self.scanDsPosEntry.pack(pady=padyControls)
        scanDsWait_label.pack(pady=padyControls)
        self.scanDsWaitEntry.pack(pady=padyControls)
        scanFluence_label.pack(pady=padyControls)
        self.scanFluenceEntry.pack(pady=padyControls)
        scanHwpWait_label.pack(pady=padyControls)
        self.scanHwpWaitEntry.pack(pady=padyControls)
        scanBatchSize_label.pack(pady=padyControls)
        self.batchSizeEntry.pack(pady=padyControls)
        scanDirLoc_label.pack(pady=padyControls)
        self.scanDirLocEntry.pack(pady=padyControls)
        self.transCorrBox.pack(pady=padyControls)
        self.scanMedianFilterBox.pack(pady=padyControls)
        self.smartScanBox.pack(pady=padyControls)
        self.startScanButton.pack(pady=padyControls)
        
        def userPressReturn(event):
            entry = self.root.focus_get() # Get the entry that the user has selected
            if entry == self.gainEntry or entry == self.exposureEntry:
                camSettingsButtonFunc()
            if entry == self.dsPos1Entry:
                dsPos1ButtonFunc()
            if entry == self.dsPos2Entry:
                dsPos2ButtonFunc()
            if entry == self.fluence1Entry:
                fluence1ButtonFunc()
            if entry == self.hwp1Entry:
                hwp1ButtonFunc()
                
        # Bind the enter key to activate selected entry
        self.root.bind('<Return>', userPressReturn)
        
        # Build sliders for adjusting the brightness and contrast of the image
        ax_setVals = [plt.axes([0.15, 0.06, 0.5, 0.02]), plt.axes([0.15, 0.02, 0.5, 0.02])]
        
        # Check whether the image is 8bit or 16bit
        if type(self.camImage[0,0]) == np.dtype('uint8'):
            imageTypeMin, imageTypeMax = 0, 255
            print('Image Type is 8 bit')
        elif type(self.camImage[0,0]) == np.dtype('uint16'):
            imageTypeMin, imageTypeMax = 0, 2**16-1
            print('Image Type is 16 bit')
        else:
            print(f'Invalid Image Type: {type(self.camImage[0,0])}')
            
        imageCurrentMin, imageCurrentMax = np.min(self.camImage), np.max(self.camImage)
 
        slider_vmax = matplotlib.widgets.Slider(ax_setVals[0], r'$v_{max}$', imageTypeMin, imageTypeMax, valinit=imageCurrentMax)
        slider_vmin = matplotlib.widgets.Slider(ax_setVals[1], r'$v_{min}$', imageTypeMin, imageTypeMax, valinit=imageCurrentMin)
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
        self.lineProfile_x1 = self.w//2
        self.lineProfile_y1 = self.h//2
        
        # 2nd point te line profile crosses through
        self.lineProfile_x2 = self.lineProfile_x1 + 1
        self.lineProfile_y2 = self.lineProfile_y1
        
        # ~ self.lineProfile_h = self.camAx[0].axhline(y=self.lineProfile_y,c='black')
        self.lineProfile_v = self.camAx[0].axline((self.lineProfile_x1, self.lineProfile_y1), (self.lineProfile_x1-(self.lineProfile_y2 - self.lineProfile_y1), self.lineProfile_y1 + (self.lineProfile_x2 - self.lineProfile_x1)), c='black', linestyle = 'dashed')
        self.lineProfile_upper_v = self.camAx[1].axvline(x = self.lineProfile_x1, c='black', linestyle = 'dashed')
        self.lineProfile_h = self.camAx[0].axline((self.lineProfile_x1, self.lineProfile_y1), (self.lineProfile_x2, self.lineProfile_y2), c='black')
        
        # Plot the cam image and the line profile crosshairs
        self.camImageObj = self.camAx[0].imshow(self.camImage)
        self.lineProfilePlot_h, = self.camAx[1].plot(self.displayImage[self.displayImage.shape[0]//2,:])
        self.updateLineProfilePlot()
        
        # Connect the cam figure to key press events
        self.camFig.canvas.mpl_connect('button_press_event', self.camonclick)
        self.camFig.canvas.mpl_connect('scroll_event', self.camonscroll)
        self.camFig.canvas.mpl_connect('motion_notify_event', self.camonmove)

        # Attach the matplotlib image figures to a tkinter gui window
        self.canvas = FigureCanvasTkAgg(self.camFig, master=self.detectorFrame)
        # ~ self.toolbar = NavigationToolbar2Tk(self.canvas, self.root)
        # ~ self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.label = tk.Label(text="")
        self.label.pack()
        
        # Create the data plot
        self.dataFig, self.dataAx = plt.subplots(2,1,figsize=(4, 6))
        self.dataFig.subplots_adjust(left=0.2, hspace=0.4)
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
            if tk.messagebox.askokcancel("Quit", "Do you want to quit?"):
                self.camLive = 0
                # ~ self.root.quit()
            
        # Activate the camera thread - this will take images on repeat asynchronously
        self.camLive = 1
        self.camFree = False # Checks whether we have successfully disconnected from the camera at the end of the program
        self.thrDataAcquisition = threading.Thread(target=self.dataAcquisitionLoop, args=(), kwargs={})
        self.thrDataAcquisition.start() # Will run takeImageTestRepeat
        
        self.root.protocol("WM_DELETE_WINDOW", on_closing)
        self.root.mainloop()
        
        # If 2 we start a scan
        if self.camLive == 2:
            print('Starting Scan')
            
            scanDir = str(self.scanDirLocEntry.get())
            
            # Make sure this directory isn't a previous scan's directory
            while os.path.isdir(f'{scanDir}//average'):
                print('Warning! scanDir has files in it, making a new scan dir')
                scanDir = scanDir + '_new'
                original_umask = os.umask(0)
                os.makedirs(f'{scanDir}')
                os.umask(original_umask)
                print(f'New dir is {scanDir}')
            
            scansweepselection =  str(self.scanSweepSelection.get())
            print(f'scanDir: {scanDir} \n')
            
            # Read the user entries
            dsPositions = parsePosStr(str(self.scanDsPosEntry.get()))
            dsWait = float(self.scanDsWaitEntry.get())
            fluences = parsePosStr(str(self.scanFluenceEntry.get()))
            hwpWait = float(self.scanHwpWaitEntry.get())
            batchSize = int(self.batchSizeEntry.get())
            tcq = bool(self.varTransCorr.get())
            smfq = bool(self.varScanMedianFilter.get())
            ssq = bool(self.varSmartScan.get())
            
            # Compute the array of hwp values from the fluences
            hwp_arr = np.empty(len(fluences))
            for i in range(len(hwp_arr)):
                hwp_val = self.imalus(fluences[i]) 
                if hwp_val == None:
                    print(f'WARNING! Fluence outside acceptable range: {fluences[i]}')
                hwp_arr[i] = hwp_val
            
            # Write the metadata file
            nowString = datetime.datetime.now().strftime("%Y_%m_%d-%H_%M_%S")
            scanParamsDict = { \
                            'scanStartTime':nowString,\
                            'gain':str(self.gainEntry.get()),\
                            'exposure':str(self.exposureEntry.get()),\
                            'dsPositions entry':str(self.scanDsPosEntry.get()),\
                            'dsPositions_array':str(dsPositions),\
                            'ds wait':str(self.scanDsWaitEntry.get()),\
                            'fluences entry':str(self.scanFluenceEntry.get()),\
                            'hwp_array':str(hwp_arr),\
                            'hwp wait':str(self.scanHwpWaitEntry.get()),\
                            'batchSize':str(self.batchSizeEntry.get()),\
                            'translationCorrection':str(self.varTransCorr.get()),\
                            'smartScan':str(self.varSmartScan.get()),\
                            'scansweepselection':scansweepselection,\
                            'directory':scanDir}
            with open(f'{scanDir}//scanMetaData_{nowString}.txt','w') as smdf:
                for key, value in scanParamsDict.items():
                    smdf.write('%s:%s\n' % (key, value))
                    
            # Update the persistent settings file
            psParamsDict = { \
                'dsPositions':str(self.scanDsPosEntry.get()),\
                'fluences':str(self.scanFluenceEntry.get()),\
                'dsWait':str(dsWait),\
                'hwpWait':str(hwpWait),\
                'batchSize':str(batchSize),\
                'scanSweepSelection':str(self.scanSweepSelection.get()),\
                'tcq':str(tcq),\
                'smfq':str(smfq),\
                'ssq':str(ssq),\
                'pump_malus_A':persistent_settings_dict['pump_malus_A'],\
                'pump_malus_phi':persistent_settings_dict['pump_malus_phi'],\
                'pump_malus_offset':persistent_settings_dict['pump_malus_offset']}
                                
            with open(f'persistent_settings.txt','w') as psf:
                for key, value in psParamsDict.items():
                    psf.write('%s:%s\n' % (key, value))

            
            self.root.destroy() # Totally destroy the setup gui so we can move into the scan gui
            
            scanapp = ScanApp(camera=self.camera,ds=self.ds,hwp=self.hwp,dsPositions=dsPositions,fluences=fluences,hwpWait=hwpWait,dsWait=dsWait,imalus=self.imalus,batchSize=batchSize,\
                              translationCorrectionQ=tcq, scanMedianFilterQ=smfq, smartScanQ=ssq, scanDir=scanDir, rm=self.rm, scansweepselection = scansweepselection)
        else:
            print('Disconnecting camera and delay stage')
            self.ds.disconnect() # disconnect from delay stage
            #  Ending acquisition appropriately helps ensure that devices clean up
            #  properly and do not need to be power-cycled to maintain integrity.
            self.camera.setExposure(1) # Reset exposure to 1.0 s and disconnect from the camera.
            self.camera.disconnect() # disconnect from camera
      
    def camonclick(self,event):
        ax = self.camAx[0]
        # Find the location of the mouse click and update the crosshair variable to this value
        selection = self.roiSelection.get()
        if event.inaxes == ax:
            if selection == 'line profile':
                if event.button == 1:
                    # compute the vector from the 1 point to the 2 point of the line
                    delta_x, delta_y = self.lineProfile_x2 - self.lineProfile_x1, self.lineProfile_y2 - self.lineProfile_y1
                    
                    # update the 1st point coordinates
                    self.lineProfile_x1, self.lineProfile_y1 = event.xdata, event.ydata
                    
                    # update the 2nd point coordinates
                    self.lineProfile_x2, self.lineProfile_y2 = delta_x + self.lineProfile_x1, delta_y + self.lineProfile_y1
                    
                    # Set the center position for the horizontal line
                    # ~ self.lineProfile_h.set_xy1(x=self.lineProfile_x1, y=self.lineProfile_y1)
                    self.lineProfile_h._xy1 = (self.lineProfile_x1, self.lineProfile_y1)
                    
                    # Set the 2nd position of the horizontal line 
                    # ~ self.lineProfile_h.set_xy2(x=self.lineProfile_x2, y=self.lineProfile_y2)
                    self.lineProfile_h._xy2 = (self.lineProfile_x2, self.lineProfile_y2)
                    
                    # Set the center position for the vertical line in both upper and lower plots
                    # ~ self.lineProfile_v.set_xy1(x=self.lineProfile_x1, y=self.lineProfile_y1)
                    self.lineProfile_v._xy1 = (self.lineProfile_x1, self.lineProfile_y1)
                    self.lineProfile_upper_v.set_xdata([self.lineProfile_x1])

                if event.button == 3:
                    if (self.lineProfile_x1 != event.xdata or self.lineProfile_y1 != event.ydata):
                        self.lineProfile_x2, self.lineProfile_y2 = event.xdata, event.ydata
                        # ~ self.lineProfile_h.set_xy2(x=self.lineProfile_x2, y=self.lineProfile_y2)
                        self.lineProfile_h._xy2 = (self.lineProfile_x2, self.lineProfile_y2)
                                     
                # Keep the vertical line vertical
                # ~ self.lineProfile_v.set_xy2(x=self.lineProfile_x1, y=self.lineProfile_y1 - 1)
                self.lineProfile_v._xy2 = (self.lineProfile_x1, self.lineProfile_y1 - 1)               

                self.updateLineProfilePlot()

            if selection == 'roi rectangle':
                if event.button == 1:
                    if self.roiActive:
                        self.roiActive = False
                        self.rm[-1].goLive() # Set the roi that was just placed to live status
                    else:
                        self.roiActive = True
                        self.rm.append(RoiRectangle(event.xdata,event.ydata,0,0,ax))
                        
                if event.button == 3:
                    for i in range(len(self.rm)):
                        if abs(event.xdata-self.rm[i].cx) < self.rm[i].w/2 and abs(event.ydata-self.rm[i].cy) < self.rm[i].h/2:
                            self.rm[i].patch.remove() # Remove the patch from the plot
                            del self.rm[i].patch # Delete the patch
                            self.rm.pop(i) # Pop this item from the rm list
                            break # Exits the loop on the first occurance of a ROI to remove.
                
                self.roi_number_label_string.set(f'ROI # = {len(self.rm)}')        
                
            
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
            
    
            self.canvas.draw()
            
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
                xlims = [0, self.w]
                ylims = [self.h, 0]
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
            
            self.canvas.draw()
    
    def camonmove(self, event):
        ax = self.camAx[0]
        if event.inaxes == ax:
            if (int(event.xdata) >= 0) and (int(event.ydata) >= 0) and (int(event.xdata) < self.displayImage.shape[1]) and (int(event.ydata) < self.displayImage.shape[0]):
                self.mouse_x = int(event.xdata)
                self.mouse_y = int(event.ydata)
            if self.roiActive:
                self.rm[-1].updateShape(abs(self.mouse_x - self.rm[-1].cx)*2,abs(self.mouse_y - self.rm[-1].cy)*2) # width, height
                self.canvas.draw()
            
            self.cursor_xycoord_label_string.set(f'(x, y) = ({self.mouse_x}, {self.mouse_y})')
            self.cursor_intensity_label_string.set(f'pixel value = {int(self.displayImage[self.mouse_y, self.mouse_x])}')
            
    def saveImage(self):
        now = datetime.datetime.now()
        nowString = now.strftime("%Y_%m_%d-%H_%M_%S")
        defaultFileName = nowString+'expo'+ str(self.camera.getExposure()) + 'gain' + str(self.camera.getGain()) + '.tiff'
        self.saveImageFileName = tk.filedialog.asksaveasfile(mode='w',initialfile=defaultFileName,defaultextension=".*",filetypes = [("image files", ".tiff")])
        
        # If image averaging is on, the type will be a 64 bit float. We convert to 32 bit float here so that imageJ can open the saved image file
        if type(self.displayImage[0,0]) == np.dtype('float64'):
            self.displayImage = np.array(self.displayImage, dtype = 'float32')
            
        print(type(self.displayImage[0,0]))
        
        if self.saveImageFileName != None:
            tifffile.imwrite(self.saveImageFileName.name,self.displayImage)
            
        del self.saveImageFileName
    
    def updateLineProfilePlot(self):
        # Put the data from the display image to the line profile
        line_profile_data = data_analysis.getLineProfile(self.displayImage, (self.lineProfile_x1, self.lineProfile_y1), (self.lineProfile_x2, self.lineProfile_y2))
        self.lineProfilePlot_h.set_ydata(line_profile_data)
        if len(line_profile_data[line_profile_data!=0]) > 0:
            self.camAx[1].set_ylim(min(line_profile_data[line_profile_data!=0])-1,max(line_profile_data)+1) # set the scales for us to see the data
            self.camAx[1].set_yticks([min(line_profile_data[line_profile_data!=0]),(0.5*max(line_profile_data)+0.5*min(line_profile_data[line_profile_data!=0])),max(line_profile_data)])
        
    def dataAcquisitionLoop(self):
        processPrev = False
        newPointSwitch = False
        
        while self.camLive == 1:
            # ~ print(time.time(),' ',osprocess.memory_info().rss / 2 ** 20) # Print this to check for memory leaks
            if not self.resetHistory:
                # Sends the software trigger, camera starts collecting async from this loop
                self.camera.trigger()
                t_open = time.time() # Moment exposure started
            else:
                newPointSwitch = True
            
            # Process the previous image (not the one that was just triggered)
            if processPrev:
                # Check whether we are saving image sets, saves the raw image (camImage)
                if self.saveImageSet and self.avgDirLocEntry.get() != '':
                    tifffile.imwrite(f'{self.avgDirLocEntry.get()}//{time.time()}.tiff', self.camImage)
                
                # Check whether despeckle is applied
                if self.varMedianFilter.get() == 1:
                    self.camImage = median_filter(self.camImage, size = 3)
                
                # Check whether image averaging is selected in which the display will be updated to an average of all the images taken    
                if self.varAvgDis.get() == 0:
                    self.displayImage = self.camImage
                    self.numAvgDis = 0
                else:
                    self.displayImage = (self.camImage + self.numAvgDis*self.displayImage)/(self.numAvgDis+1)
                    self.numAvgDis += 1
                
                # Update the label to indicate the number of images that have been taken to create the displayed image
                self.varAvgDis_label_string.set(f'Average Images (N = {self.numAvgDis})')
                
                # Compute the image rate
                deltaTString = '%.3f' % round(time.time()-self.t1,5)
                self.camAx[0].set_title(f"Image Rate =  {deltaTString} s")
                
                # Update the intensity value for the cursor given the new image
                self.cursor_intensity_label_string.set(f'pixel value = {int(self.displayImage[self.mouse_y, self.mouse_x])}')
                
                self.t1 = time.time()
                
                self.camImageObj.set_data(self.displayImage) # Set image data to the plot
                self.updateLineProfilePlot()
            
                if len(self.rm) > 0:
                    roiPixelSum = 0
                    roiTotalArea = 0
                    # Loop through all the ROIS
                    for r in self.rm:
                        # Only count an ROI if it is active
                        if r.live:
                            roiImage = np.array(self.camImage[int(r.cy-r.h/2):int(r.cy+r.h/2),int(r.cx-r.w/2):int(r.cx+r.w/2)], dtype = 'float64') # The roi sum is taken from single images even when averaging is on
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
                
                self.dataAx[0].set_title(f"N = {len(self.intensityHistory)}")
                self.liveIntensityPlot.set_xdata(self.timeHistory)
                self.liveIntensityPlot.set_ydata(self.intensityHistory)
                
                update_errorbar(self.avgIntensityPlot,range(len(self.avgIntensity)),self.avgIntensity,yerr=self.stdIntensity)
                
                # Update the figures
                self.canvas.draw()
                self.datacanvas.draw()
            
            if newPointSwitch:
                processPrev = False
                self.intensityHistory = []
                self.timeHistory = []
                self.avgIntensity = np.append(self.avgIntensity,[0]) # Don't worry about the zeros, they are rewritten and don't through off the average or stdev
                self.stdIntensity = np.append(self.stdIntensity,[0])
                
                # Zero the display image (necessary for averaging but no reason not to just always do it)
                self.displayImage = 0*self.displayImage
                self.numAvgDis = 0
            
                # Freeze here until user is ready to move on
                while newPointSwitch and self.resetHistory:
                    time.sleep(self.wait)
                    self.t0 = time.time()
                newPointSwitch = False
            else:
                t_rem = 0.7*self.exposure - (time.time() - t_open) # amount of time remaining in exposure
                N = int(30*t_rem) # 30 frames per second
                if N >= 4:
                    # Waits here while updating the progress bar
                    spf = t_rem/N
                    self.imProgBar['value'] = 100/N
                    for n in range(N):
                        self.imProgBar.step(100/N)
                        time.sleep(spf)
                else:
                    self.imProgBar['value'] = 100
                
                # Grabs the image
                self.camImage = self.camera.grabImage()
                processPrev = True # Tell the data acquisition loop that there is an image that needs processing
            
            # Here is where we send exposure and gain updates to the camera
            if self.camSettingsUpdateInProgress == True:
                # Send the set values to the camera
                self.camera.setGain(self.gain)
                self.camera.setExposure(self.exposure)
                
                # Take an image here to ensure update has taken place, camera doesn't update exposure until after this image has been taken.
                self.camera.trigger()
                self.camera.grabImage()
                
                # Turn off the update in progress variable
                self.camSettingsUpdateInProgress = False
            
            # Check for a total reset of roi data plots
            if self.resetLongHistory:
                self.resetLongHistory = False
                                
                self.intensityHistory = []
                self.timeHistory = []
                
                # Zero the display image (necessary for averaging but no reason not to just always do it)
                self.displayImage = 0*self.displayImage
                self.numAvgDis = 0
                
                self.avgIntensity = np.array([0],dtype='float')
                self.stdIntensity = np.array([0],dtype='float')
    
    
        self.root.quit() # Close the prep GUI window when the scan is started
    
        
class ScanApp():
    def __init__(self, camera, ds, hwp, dsPositions, fluences, hwpWait, dsWait, imalus, batchSize, translationCorrectionQ, scanMedianFilterQ, smartScanQ, scanDir, rm, scansweepselection, wait = .033):
        self.camera=camera
        self.ds=ds
        self.hwp=hwp
        self.exposure=self.camera.getExposure()
        self.gain=self.camera.getGain()
        self.dsPositions=dsPositions
        self.fluences=fluences
        self.hwpWait=hwpWait
        self.imalus=imalus
        self.batchSize=batchSize
        self.translationCorrectionQ=translationCorrectionQ
        self.scanMedianFilterQ = scanMedianFilterQ
        self.smartScanQ=smartScanQ
        self.scanDir=scanDir
        self.dsWait = dsWait
        self.rm = rm
        self.scansweepselection = scansweepselection
        
        now = datetime.datetime.now()
        self.runLog = f'scan app started at {now.strftime("%Y_%m_%d-%H_%M_%S")} \n'
        self.runLogFileName = f'runLog_{now.strftime("%Y_%m_%d-%H_%M_%S")}.txt'
        
        print('imalus test',self.imalus(3))
        
        # Some general formatting vars
        padyControls = 2
        padySep = 2
        entryWidth = 30
        sepWidth = 190
        sepHeight = 3
        cntrlBg = '#2F8D35'
        self.dirPermission = 777
        
        # Pause and stop variables
        self.paused = False
        self.gentleStop = False
        self.hardStop = False
        
        # Setup the tkinter GUI window
        self.root = tk.Tk()
        #self.root.attributes("-fullscreen", True)
        self.root.geometry("1000x500")
        self.root.configure(bg='white')
        self.root.wm_title("Smart Scan")
        self.root.iconbitmap("kogar_0.ico")
        #self.root.iconify()
        
        # Make a Tkinter frame for the data plots
        self.dataFrame = tk.Frame(master=self.root, width=200, height=100, padx=5, pady=5)
        self.dataFrame.pack(side=tk.RIGHT)
        
        # Make a Tkinter frame for controls and indicators
        self.controls = tk.Frame(master=self.root, width=100, height=100, padx=5, pady=5, bg=cntrlBg)
        self.controls.pack(side=tk.RIGHT)
        
        # Setup the camera window
        self.camFig, self.camAx = plt.subplots()
        self.camAx.xaxis.set_tick_params(labelbottom=False)
        self.camAx.yaxis.set_tick_params(labelleft=False)
        self.camAx.set_xticks([])
        self.camAx.set_yticks([])
        
        # Put the rois on the new cam fig
        for r in self.rm:
            r.updateAxis(self.camAx)
        
        # Data type for saving images
        self.image_dtype = np.float32
        
        # Sends the software trigger to the camera
        self.camera.trigger()
        
        initImage = self.camera.grabImage()
            
        self.zeroImage = np.zeros(np.shape(initImage), dtype=self.image_dtype) # A zero image that can be used later to generate empty initial files
        
        # If translation correction is active, find the initial reference location for later comparison
        self.tcorr_0 = np.zeros([len(self.rm),2]) # array that stores the reference positions for translation correction
        self.tcorr = np.zeros([len(self.rm),2])
        self.tcorr_log = np.zeros([1,2])
        self.roi_total_area = 0
        
        for ri in range(len(self.rm)):
            r = self.rm[ri]
            # Get pixel data in the roi
            roi_image = initImage[int(r.cy-r.h/2):int(r.cy+r.h/2),int(r.cx-r.w/2):int(r.cx+r.w/2)]
            self.roi_total_area += roi_image.shape[0]*roi_image.shape[1]
            if self.translationCorrectionQ:
                popt, rms, self.guess_prms = tc.fitImageToGaussian(roi_image)
                self.tcorr_0[ri, :] = popt[:2]
        
        # Data display settings and variables
        self.wait = wait
        self.ppw = 100
        
        # Create intensity history arrays for data plots
        self.intensityLabTime = []
        self.intensityScanTime = np.zeros(len(self.dsPositions))
        self.stdIntensityScanTime = np.zeros(len(self.dsPositions))
        self.roiScanData = np.zeros([1,len(self.dsPositions)]) # full record of intensity in the roi
        self.timeHistory = np.zeros([1,len(self.dsPositions)]) # coincident record of when data was taken
        
        self.scanProgress_label_string = tk.StringVar()
        self.scanProgress_label_string.set(f'Scan Progress = 0/{len(self.dsPositions)}')
        self.scanProgress_label = tk.Label(master=self.controls,textvariable= self.scanProgress_label_string,bg=cntrlBg)
        
        self.batchSize_label_string = tk.StringVar()
        self.batchSize_label_string.set(f'Batch Size = {self.batchSize}')
        self.batchSize_label = tk.Label(master=self.controls,textvariable=self.batchSize_label_string,bg=cntrlBg)
        
        self.scanNumber_label_string = tk.StringVar()
        self.scanNumber_label_string.set(f"Scan # = ...")
        self.scanNumber_label = tk.Label(master=self.controls,textvariable=self.scanNumber_label_string,bg=cntrlBg)
        
        self.batchNumber_label_string = tk.StringVar()
        self.batchNumber_label_string.set(f"Batch # = ...")
        self.batchNumber_label = tk.Label(master=self.controls,textvariable=self.batchNumber_label_string,bg=cntrlBg)
        
        self.dsPosition_label_string = tk.StringVar()
        self.dsPosition_label_string.set(f"dsPos = ... ")
        self.dsPosition_label = tk.Label(master=self.controls,textvariable=self.dsPosition_label_string,bg=cntrlBg)
        
        self.exposure_label_string = tk.StringVar()
        ces = '%.2f' % self.camera.getExposure()
        self.exposure_label_string.set(f'Exposure = {ces} s')
        self.exposure_label = tk.Label(master=self.controls,textvariable=self.exposure_label_string,bg=cntrlBg)
        
        self.gain_label_string = tk.StringVar()
        cgs = '%.2f' % self.camera.getGain()
        self.gain_label_string.set(f'Gain = {cgs} dB')
        self.gain_label = tk.Label(master=self.controls,textvariable=self.gain_label_string,bg=cntrlBg)
        
        def pauseScanButtonFunc():
            if self.paused == False:
                self.paused = True
                self.pauseScanButton.config(text='Resume Scan',bg='yellow')
            else:
                self.paused = False
                self.pauseScanButton.config(text='Pause Scan',bg='white')
                
        self.pauseScanButton = tk.Button(
            master=self.controls,
            text="Pause Scan",
            width=18,
            height=1,
            fg='black', bg='white',
            command=pauseScanButtonFunc)
            
        def gentleStopButtonFunc():
            if self.gentleStop == False:
                self.gentleStop = True
                self.gentleStopButton.config(bg='yellow')
            else:
                self.gentleStop = False
                self.gentleStopButton.config(bg='white')
                
        self.gentleStopButton = tk.Button(
            master=self.controls,
            text="Gentle Stop",
            width=18,
            height=1,
            fg='black', bg='white',
            command=gentleStopButtonFunc)
        
        def hardStopButtonFunc():
            if self.hardStop == False:
                self.hardStop = True
                self.hardStopButton.config(bg='yellow')
            else:
                self.hardStop = False
                self.hardStopButton.config(bg='white')
                
        self.hardStopButton = tk.Button(
            master=self.controls,
            text="Hard Stop",
            width=18,
            height=1,
            fg='black', bg='white',
            command=hardStopButtonFunc)
        
        self.errorLog_label_string = tk.StringVar()
        self.errorLog_label_string.set(f'Error Log: \n')
        self.errorLog_label = tk.Label(master=self.controls,textvariable=self.errorLog_label_string,bg=cntrlBg)
        
        # Pack everything on controls frame
        self.batchSize_label.pack(pady=padySep)
        self.scanNumber_label.pack(pady=padySep)
        self.batchNumber_label.pack(pady=padySep)
        self.scanProgress_label.pack(pady=padySep)
        self.dsPosition_label.pack(pady=padySep)
        self.exposure_label.pack(pady=padySep)
        self.gain_label.pack(pady=padySep)
        self.pauseScanButton.pack(pady=padySep)
        self.gentleStopButton.pack(pady=padySep)
        self.hardStopButton.pack(pady=padySep)
        self.errorLog_label.pack(pady=padySep)
        
        # Build sliders for adjusting the brightness and contrast of the image
        ax_setVals = [plt.axes([0.15, 0.06, 0.5, 0.02]), plt.axes([0.15, 0.02, 0.5, 0.02])]
        
        # Check whether the image is 8bit or 16bit
        if type(initImage[0,0]) == np.dtype('uint8'):
            imageTypeMin, imageTypeMax = 0, 255
            print('Image Type is 8 bit')
        elif type(initImage[0,0]) == np.dtype('uint16'):
            imageTypeMin, imageTypeMax = 0, 2**16-1
            print('Image Type is 16 bit')
        else:
            print(f'Invalid Image Type: {type(initImage[0,0])}')
            
        imageCurrentMin, imageCurrentMax = np.min(initImage), np.max(initImage)
 
        slider_vmax = matplotlib.widgets.Slider(ax_setVals[0], r'$v_{max}$', imageTypeMin, imageTypeMax, valinit=imageCurrentMax)
        slider_vmin = matplotlib.widgets.Slider(ax_setVals[1], r'$v_{min}$', imageTypeMin, imageTypeMax, valinit=imageCurrentMin)
        
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
        
        # Plot the cam image
        self.camImageObj = self.camAx.imshow(initImage)
        pos0 = self.camAx.get_position()

        # Attach the matplotlib image figures to a tkinter gui window
        self.canvas = FigureCanvasTkAgg(self.camFig, master=self.root)
        self.canvas._tkcanvas.pack(side=tk.RIGHT, fill=tk.BOTH, expand=1)
        self.label = tk.Label(text="")
        self.label.pack()
        
        # Create the data plot
        self.dataFig, self.dataAx = plt.subplots(2,1)
        self.intensityLabTimePlot, = self.dataAx[0].plot([],[],'go-') # live intensity plot
        self.dataAx[0].set_xlabel('Scan #')
        self.intensityScanTimePlot = self.dataAx[1].errorbar(self.dsPositions,self.intensityScanTime,self.stdIntensityScanTime,\
                                    marker='o',color='green',elinewidth=1,capsize=5,linestyle='') # scan intensity plot with error bars at 1 sigma
        self.dataAx[1].set_xlabel('dsPositions (mm)')
        
        # Attach the matplotlib data figures to a tkinter gui window
        self.datacanvas = FigureCanvasTkAgg(self.dataFig, master=self.dataFrame)
        self.datacanvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        # ~ self.datalabel = tk.Label(text="")
        # ~ self.datalabel.pack()
        
        def on_closing():
            # ~ if tk.messagebox.askokcancel("Quit", "Do you want to quit?"):
            self.scanLive = False
            self.root.quit()
            
        # Activate the data acquisition loop
        self.scanLive = True
        self.camFree = False # Checks whether we have successfully disconnected from the camera at the end of the program
        self.dsFree = False # Checks whether we have successfully disconnected from the delay stage at the end of the program
        thr = threading.Thread(target=self.dataAcquisitionLoop, args=(), kwargs={})
        thr.start() # Will run takeImageTestRepeat
        
        self.root.protocol("WM_DELETE_WINDOW", on_closing)
        self.root.mainloop()
        
    def dataAcquisitionLoop(self):
        # Log the current time
        self.runLog += f'dataAcquisitionLoop started at {datetime.datetime.now().strftime("%Y_%m_%d-%H_%M_%S")} \n'
        
        original_umask = os.umask(0)
        self.errorLog_label_string.set(self.errorLog_label_string.get() + f'all clear')
             
        os.makedirs(f'{self.scanDir}//average',self.dirPermission) # Create a directory for the batch average files
        os.umask(original_umask)
        
        # Create the files for the total average
        for p in range(len(self.dsPositions)):
            for f in range(len(self.fluences)):
                dsPos = self.dsPositions[p] # direction doesn't matter here; we just need the positions for the files names
                flu = self.fluences[f]
                dsPosString = '%.4f' % dsPos
                fluString = '%.4f' % flu
                tifffile.imwrite(f'{self.scanDir}//average//pos={dsPosString}_flu={fluString}.tiff',self.zeroImage)

        thrBatchAverage = None # Variable that contains the moving batch average thread
        
        # Create the index variables
        b = -1
        s = -1
        
        # Make a string that records the formattted start time
        self.t_startString = datetime.datetime.now().strftime("%Y_%m_%d-%H_%M_%S")
        
        # Record the exact starting time
        self.t_start = time.time()
        t_lastImage = self.t_start # This variable will later be used to record when images finish
        
        # Create a meta data file for the total average folder, write a row for each delay stage position
        with open(f'{self.scanDir}//average//averageMetaData.txt', 'w') as fileObj:
            for p in range(len(self.dsPositions)):
                for f in range(len(self.fluences)):
                    fileObj.write(f'pi = {p}, dsPos = {self.dsPositions[p]} mm, fi = {f}, fluence = {self.fluences[f]} Ko, weight = 0 \n')
        
        # Double array makes back and forth motion automatic
        self.fluence_double_array = np.append(self.fluences,self.fluences[::-1])
        
        ### Start the experiment loop ###
        while self.scanLive:
            b += 1 # Record that a new batch is being taken, note that b is initiated at a value of -1
            self.batchNumber_label_string.set(f"Batch # = {b}") # Write the current batch number to the gui
            
            # Update the fluence
            flu = self.fluence_double_array[b%len(self.fluence_double_array)] # get the new fluence from the doubled over fluence array, it's doubled so that motion is back and forth
            hwpPos = self.imalus(flu)
            print(f'moving to hwpPow = {hwpPos}, fluence = {flu}, batch # = {b}')
            self.hwp.moveAbsolute(hwpPos) # Sends command to the HWP
            readHwpPos = self.hwp.getPos() # Read the half wave plate value back to compare to the send value
            
             # If there is more than one fluence, then wait for hwpWait many seconds so that temp can stabilize
            if len(self.fluences) > 1:
                time.sleep(self.hwpWait)
            
            # Log batch start time
            self.runLog += f'batch {b} started at {datetime.datetime.now().strftime("%Y_%m_%d-%H_%M_%S")} \n'
            # Log current set fluence and the read hwp pos
            self.runLog += f'set fluence = {flu} Ko, set hwpPos = {hwpPos} deg, read hwpPos = {readHwpPos} deg\n'
            
            # Create a directory for the new batch files
            os.makedirs(f'{self.scanDir}//batch{b}')
            
            # Create a batch meta data file, write a row for each delay stage position
            with open(f'{self.scanDir}//batch{b}//batchMetaData.txt', 'w') as f:
                f.write(f'Experiment start time  = {self.t_startString}\n')
                f.write(f'batch {b} started at {datetime.datetime.now().strftime("%Y_%m_%d-%H_%M_%S")} \n')
                f.write(f'batchSize = {self.batchSize}\n')
                f.write(f'fluence = {flu} Ko\n')
                f.write(f'set HWP Pos = {hwpPos} deg\n')
                f.write(f'read HWP Pos = {readHwpPos} deg\n')
                for p in range(len(self.dsPositions)):
                    f.write(f'pi = {p}, dsPos = {self.dsPositions[p]} mm, weight = 0 \n')

            # Create the files of a new batch
            for p in range(len(self.dsPositions)):
                dsPos = self.dsPositions[p] # direction doesn't matter here; we just need the positions for the files names
                dsPosString = '%.4f' % dsPos
                tifffile.imwrite(f'{self.scanDir}//batch{b}//pos={dsPosString}.tiff',self.zeroImage)
        
            ### Start the new batch ###
            for bi in range(self.batchSize):
                s += 1 # Record that a new scan is being taken, note that s is inititated at a value of -1
                self.scanNumber_label_string.set(f"Scan # = {s}") # Write the current scan number to the gui
                
                # Log scan start time
                self.runLog += f'\t scan {s} (bi = {bi}) started at {datetime.datetime.now().strftime("%Y_%m_%d-%H_%M_%S")} \n'
                
                # Get the position index list depending on scan sweep option
                piList = self.computePiList(s)
                
                ### Start a new scan ###
                for p in range(len(self.dsPositions)):
                    # ~ print(time.time(),' ',osprocess.memory_info().rss / 2 ** 20) # Print this to check for memory leaks
                    pi = piList[p]
                    dsPos = self.dsPositions[pi]
                    
                    # Move to the next ds pos
                    self.ds.setPos(dsPos)
                    
                    # Give a wait for the motor time
                    time.sleep(self.dsWait)
                    
                    # read back the delay stage position
                    dsReadBack = self.ds.getPos()
                    
                    # Write the current ds pos to the gui
                    dps = '%.2f' % dsReadBack
                    self.dsPosition_label_string.set(f"dsPos = {dps} mm")
                    self.scanProgress_label_string.set(f'Scan Progress = {p+1}/{len(self.dsPositions)}')
                    
                    # Record the current ds pos and the read back pos in the run log
                    self.runLog += f'\t \t Delay Stage Set: pi = {pi}, dsSet = {dsPos}, dsReadBack = {dsReadBack}| (t-t_start) =  {round(time.time()-self.t_start,2)} \n'
                    
                    # Sends the software trigger to the camera
                    self.camera.trigger()
                    
                    # Process previous image (that is the one taken before the currently running trigger), skip initial
                    if p > 0:
                        self.processImage(imageData, pi_prev, s, b, t_lastImage)
                    
                    # Get the new image from the camera and record the parameters at the time of this image
                    imageData, pi_prev = self.camera.grabImage(), pi
                    
                    # Compute the image rate
                    deltaTString = '%.3f' % round(time.time()-t_lastImage,5)
                    deltaETString = '%.3f' % round(time.time()- t_lastImage - self.exposure - self.dsWait, 5)
                    self.camTitle = self.camAx.set_title(f"Image Rate =  {deltaTString} s (et = {deltaETString} s)")
                    t_lastImage = time.time()
                    
                    # Before moving to next loop check hard stop and pause
                    if self.hardStop:
                        # Process the last image of the scan before hard stop
                        self.processImage(imageData, pi_prev, s, b, t_lastImage) 
                        
                        # Make sure all extra threads are closed out
                        while len(threading.enumerate()) > 2:
                            time.sleep(self.wait)
                            
                        # Do the final batch average
                        self.movingAverageBatch(b)
                        
                        #exits the dataAcquisitionLoop
                        return self.endScan(finalImage = imageData, finalScanNumber = s)
                    
                    while self.paused:
                        time.sleep(self.wait) # block here if paused

                # This point in the loop occurs when a scan has concluded
                self.processImage(imageData, pi_prev, s, b, t_lastImage) # Process the last image of the scan
                
                # Before moving on to a new scan, check whether a gentle stop call is active
                if self.gentleStop and not self.hardStop:
                    # Make sure all extra threads are closed out
                    while len(threading.enumerate()) > 2:
                        time.sleep(self.wait)
                        
                    # Do the final batch average
                    self.movingAverageBatch(b)
                    
                    #exits the dataAcquisitionLoop
                    return self.endScan(finalImage = imageData, finalScanNumber = s)
                    
                # Update the lab time intensity monitor
                self.intensityLabTime.append(0)
                self.intensityLabTime[s] = np.mean(self.roiScanData[s,:]) # Update the lab time plot since the scan has ended
            
            # This point in the loop occurs when a batch has concluded
            # Make sure the previous batch moving average thread has concluded
            if (thrBatchAverage != None) and (thrBatchAverage.is_alive()):
                print('Joining Batch Average Thread')
                thrBatchAverage.join()
                
            # Start the batch moving average thread
            thrBatchAverage = threading.Thread(target=self.movingAverageBatch, args=(b,), kwargs={})
            thrBatchAverage.start()
            
            # We do a temporary save of the roiScanData and timeHistory
            np.savetxt(self.scanDir + '//' + 'roiScanData.csv',self.roiScanData,delimiter=',')
            np.savetxt(self.scanDir + '//' + 'timeHistory.csv',self.timeHistory,delimiter=',')
            
            # Save the run log to a txt file
            with open(self.scanDir + '//' + self.runLogFileName, "w") as text_file:
                text_file.write(self.runLog)
            if len(self.runLog) > 1e6:
                now = datetime.datetime.now()
                self.runLog = f'log started at {now.strftime("%Y_%m_%d-%H_%M_%S")} \n'
                self.runLogFileName = f'runLog_{now.strftime("%Y_%m_%d-%H_%M_%S")}.txt'
        
        self.endScan()
    
    def processImage(self, image, pi, s, b, t_image):
        bi_calc = s%self.batchSize
        # Record image process info in the batch meta data file                            
        incrementBatchMetaDataWeight(fileDir=f'{self.scanDir}//batch{b}//batchMetaData.txt', pi = pi)        
        
        if self.scanMedianFilterQ:
            image = median_filter(image, size=3)
        
        # Update the ROI data plots
        if len(self.rm) > 0:
            if self.translationCorrectionQ:
                # Loop through all the ROIS
                for ri in range(len(self.rm)):
                    r = self.rm[ri]
                    roi_image = image[int(r.cy-r.h/2):int(r.cy+r.h/2),int(r.cx-r.w/2):int(r.cx+r.w/2)]
                    popt, rms, self.guess_prms = tc.fitImageToGaussian(roi_image, prevFit = self.guess_prms)
                    self.tcorr[ri, :] = popt[:2]
                    
                correction = np.mean(self.tcorr_0 - self.tcorr, axis = 0)
                self.tcorr_log = np.append(self.tcorr_log, correction[None, :], axis=0)
                image = scipy.ndimage.shift(image, correction[::-1], mode = 'grid-wrap') # Apply translation correction to the image with cylic boundaries.
                
            roi_pixel_sum = 0
            for r in self.rm:
                # Get pixel data in the roi
                roi_image = np.array(image[int(r.cy-r.h/2):int(r.cy+r.h/2),int(r.cx-r.w/2):int(r.cx+r.w/2)],dtype='float64')                              
                # get the sum inside all the rois and the roi total area, we will average with a normalization to area
                roi_pixel_sum += np.sum(roi_image)
                
            if (s + 1) > self.roiScanData.shape[0]:
                self.roiScanData = np.append(self.roiScanData,np.zeros([1,len(self.dsPositions)]),axis=0) # append a new scan row to the array
                self.timeHistory = np.append(self.timeHistory,np.zeros([1,len(self.dsPositions)]),axis=0)
            
            self.roiScanData[s,pi] = roi_pixel_sum/self.roi_total_area
            self.timeHistory[s,pi] = t_image - self.t_start #SCAF
                
            self.intensityScanTime[pi] = np.mean(self.roiScanData[:,pi])
            self.stdIntensityScanTime[pi] = np.std(self.roiScanData[:,pi])/np.sqrt(len(self.roiScanData[:,pi]))
            
            
            if len(self.intensityLabTime) != 0:
                self.dataAx[0].set_xlim([-0.1,s-1+0.1])
                self.dataAx[0].set_ylim([min(self.intensityLabTime)-0.01,max(self.intensityLabTime)+0.01])
            
            self.dataAx[1].set_ylim([min(self.intensityScanTime-self.stdIntensityScanTime)-0.01,max(self.intensityScanTime+self.stdIntensityScanTime)+0.01])
                
            # draw the data canvas, we don't blit with this because we actually need the plot labels to change with each new data point
            # set the data to the data plots
            self.intensityLabTimePlot.set_xdata(range(len(self.intensityLabTime)))
            self.intensityLabTimePlot.set_ydata(self.intensityLabTime)
            update_errorbar(self.intensityScanTimePlot,self.dsPositions,self.intensityScanTime,yerr=self.stdIntensityScanTime)
            self.datacanvas.draw()

        ## Draw the figures on the canvas and the data canvas
        
        # cam canvas
        self.camImageObj.set_data(image)
        self.canvas.draw()
        
        #SCAF
        # Load the image corresponding to this position from the current batch, update the moving average, save
        dsPosString = '%.4f' % self.dsPositions[pi]
        imAvg = tifffile.imread(f'{self.scanDir}//batch{b}//pos={dsPosString}.tiff')
        imAvg = np.array((image + bi_calc*imAvg)/(bi_calc+1),dtype=self.image_dtype) # moving average
        tifffile.imwrite(f'{self.scanDir}//batch{b}//pos={dsPosString}.tiff',imAvg) # save the updated image
    
    def movingAverageBatch(self, b):
        # Get the fluence of this batch using the back and forth formula (maybe make this more robust later SCAF)
        flu = self.fluence_double_array[b%len(self.fluence_double_array)]
        fluString = '%.4f' % flu
        fi = np.where(np.array(self.fluences) == flu)[0][0] # Get the fluence index
        print(f'Computing moving average on batch {b}, fluence {fluString} Ko, fi = {fi}')
        
        for p in range(len(self.dsPositions)):
            # ~ print(f'moving batch average on position index {p}')
            dsPos = self.dsPositions[p]
            dsPosString = '%.4f' % dsPos
            
            # Load the image from the b batch folder
            imNew = tifffile.imread(f'{self.scanDir}//batch{b}//pos={dsPosString}.tiff') # Load the image corresponding to this position and the latest batch
            # Load the image from the average folder
            imAvg = tifffile.imread(f'{self.scanDir}//average//pos={dsPosString}_flu={fluString}.tiff') # Load the moving average image corresponding to this position and fluence
            # Get the weight of the new image
            wNew = readBatchMetaDataWeight(fileDir=f'{self.scanDir}//batch{b}//batchMetaData.txt', pi = p)
            # get the weight of the average image
            wAvg = readAverageMetaDataWeight(fileDir=f'{self.scanDir}//average//averageMetaData.txt', pi = p, fi = fi)
            
            # ~ print('wNew, wAvg',wNew, wAvg) #SCAF
            # Compute the new average image
            if wNew != 0 or wAvg != 0:
                imFinal = np.array((wNew*imNew + wAvg*imAvg)/(wNew + wAvg), dtype=self.image_dtype) # compute moving average
            else:
                imFinal = imNew # If this is triggered, it means the scan was stopped before any photo was taken for this image
            
            # Upate the average image weight in the average metadata file
            incrementAverageMetaDataWeight(fileDir=f'{self.scanDir}//average//averageMetaData.txt', pi = p, fi = fi, increment = wNew)
            
            # Save the updated image
            tifffile.imwrite(f'{self.scanDir}//average//pos={dsPosString}_flu={fluString}.tiff', imFinal)
            
        now = datetime.datetime.now()
        self.runLog += f'movingAverageBatch thread completed for b = {b} at {now.strftime("%Y_%m_%d-%H_%M_%S")} \n'
    
    def computePiList(self, s, probability = 0.25):
        pi_array_init = np.array(range(len(self.dsPositions)))
        
        if self.scansweepselection == 'one direction':
            return list(pi_array_init)
        
        if self.scansweepselection == 'back and forth':
            if s%2 == 0:
                pi_array_final = pi_array_init
            else:
                pi_array_final = pi_array_init[::-1]
            return list(pi_array_final)
            
        if self.scansweepselection == 'full random':
            np.random.shuffle(pi_array_init)
            return list(pi_array_init)
            
        if self.scansweepselection == 'quasi random':
            if s%2 == 0:
                pi_list = list(pi_array_init)
            else:
                pi_list = list(pi_array_init[::-1])
            
            pi_list_final = []
            
            while len(pi_list) > 1:
                for i in range(len(pi_list)):
                    coinflip = np.random.choice([False, True], p=[1 - probability, probability])
                    if coinflip:
                        pi_list_final.append(pi_list.pop(i))
                        break
                        
            # No need to run the random choice on the last element of pi_list
            pi_list_final.append(pi_list.pop(0))
                    
            return pi_list_final
                
    def endScan(self, finalImage, finalScanNumber, makeDataPlots = True):
        # Disconnect from devices    
        self.ds.disconnect() # disconnect from delay stage
        #  Ending acquisition appropriately helps ensure that devices clean up
        #  properly and do not need to be power-cycled to maintain integrity.
        self.camera.setExposure(1) # Reset exposure to 1.0 s and disconnect from the camera.
        self.camera.disconnect() # disconnect from camera
        self.root.quit() # Close the gui window
                
        # Save some things
        np.savetxt(self.scanDir + '//' + 'roiScanData.csv',self.roiScanData,delimiter=',')
        np.savetxt(self.scanDir + '//' + 'timeHistory.csv',self.timeHistory,delimiter=',')
        
        if self.translationCorrectionQ:
            # save the corrections log
            np.savetxt(self.scanDir + '//' + f'tcorrLog_batch.csv',self.tcorr_log,delimiter=',')
        
        # Save the run log to a txt file
        with open(self.scanDir + '//' + self.runLogFileName, "w") as text_file:
            text_file.write(self.runLog)
                
        # Disconnect from camera and clean up
        # Disconnect from delay stage and clean up
        
        if makeDataPlots:
            # Save data plots
            fig, axs = plt.subplots(2,3)
            
            y00 = self.roiScanData.flatten()
            t00 = self.timeHistory.flatten()
            
            if self.hardStop == True:
                t00 = np.delete(t00, y00 == 0)
                y00 = np.delete(y00, y00 == 0)
            
            sorted_indices = np.argsort(t00)
            
            t00_sorted = t00[sorted_indices]
            y00_sorted = y00[sorted_indices]
            
            mean_y00 = np.mean(y00_sorted)
            std_y00 = np.std(y00_sorted)
            
            fig.suptitle(f't0={self.t_startString}')
            
            axs[0,0].plot(t00_sorted,y00_sorted)
            axs[0,0].set_xlabel('Lab Time (s)')
            # ~ axs[0,0].legend()
            # ~ axs[0,0].set_ylim([mean_y00-6*std_y00],mean_y00+6*std_y00)
            
            deltaT = t00_sorted-np.roll(t00_sorted,1)
            deltaT = deltaT[1:]
            T = np.mean(deltaT)
            deltaS = np.max(self.timeHistory)-np.min(self.timeHistory)
            S = deltaS/(finalScanNumber+1)
            N = len(y00_sorted)
            yf = fft.fft(y00_sorted)
            xf = fft.fftfreq(N, T)[:N//2]
            axs[0,1].plot(xf, 2.0/N * np.abs(yf[0:N//2]),alpha=1,c='blue')
            
            axs[0,1].axvline(x=1/S,c='red',zorder=1)
            axs[0,1].axvline(x=1/(2*S),c='green',zorder=1)
            
            axs[0,1].set_xscale('log')
            axs[0,1].set_yscale('log')
            
            axs[0,1].set_xlabel('Freq (Hz)')
            #axs[0,1].set_ylabel('FFT')
            #axs[0,1].legend()
            
            axs[0,2].plot(range(len(self.intensityLabTime)),self.intensityLabTime,'ro')
            axs[0,2].set_xlabel('Scan #')
            
            axs[1,0].errorbar(self.dsPositions,self.intensityScanTime,self.stdIntensityScanTime,\
                               marker='o',color='green',elinewidth=1,capsize=5,linestyle='') # scan intensity plot with error bars at 1 sigma
            axs[1,0].set_xlabel('dsPos (mm)')
            
            scanTimes = (55-self.dsPositions)*6.671281904 # 55 is chosen as a nominal "time zero" to keep the values reasonable
            
            axs[1,1].errorbar(scanTimes,self.intensityScanTime,self.stdIntensityScanTime,\
                       marker='o',color='green',elinewidth=1,capsize=5,linestyle='') # scan intensity plot with error bars at 1 sigma
            axs[1,1].set_xlabel('t-t0 (ps)')
            
            axs[1,2].imshow(finalImage)
            for r in self.rm:
                tl_x = r.cx - r.w/2
                tl_y = r.cy - r.h/2
                axs[1,2].add_patch(patches.Rectangle((tl_x, tl_y), r.w, r.h, linewidth=1, edgecolor='r', facecolor='none', alpha=0.8))
            
            axs[1,2].xaxis.set_tick_params(labelbottom=False)
            axs[1,2].yaxis.set_tick_params(labelleft=False)
            axs[1,2].set_xticks([])
            axs[1,2].set_yticks([])
            
            plt.tight_layout()
            fig.savefig(self.scanDir + '//' + 'roiFigure.pdf')
            
            
            # Save a plot showing the translation corrections
            fig, axs = plt.subplots(2)

                
            if self.translationCorrectionQ:
                if self.hardStop == True:
                    t00_sorted = t00_sorted[0:len(self.tcorr_log[1:,0])]
                # ~ print(len(t00_sorted),len(self.tcorr_log[1:,0]))
                axs[0].plot(t00_sorted,self.tcorr_log[1:,0])
                axs[0].set_ylabel('X shift')
                axs[0].set_xlabel('Lab Time (s)')
                axs[1].plot(t00_sorted,self.tcorr_log[1:,1])
                axs[1].set_ylabel('Y shift')
                axs[1].set_xlabel('Lab Time (s)')
                fig.suptitle(f't0={self.t_startString}')
                
                plt.tight_layout()
                fig.savefig(self.scanDir + '//' + 'translationCorrectionFigure.pdf')
        
        # Close the scan GUI window
        self.root.quit()
                  
                    
if __name__ == '__main__':
    app = SetupApp()
               
