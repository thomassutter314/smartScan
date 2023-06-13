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
    
# Use a matplotlib backend that doesn't have an associated gui
plt.switch_backend('agg')


import numpy as np
import scipy
from scipy import fft
from scipy.ndimage import median_filter
import time
import datetime
import threading
import tifffile
# ~ import cv2
import json
import os
import io
import PySpin
import sys
import clr
from clr import System
full_filename = r'C:\Windows\Microsoft.NET\assembly\GAC_64\Newport.DLS.CommandInterface\v4.0_1.0.1.0__90ac4f829985d2bf\Newport.DLS.CommandInterface.dll'
clr.AddReference(full_filename)
from CommandInterfaceDLS import *

import translationcorr as tc

def parsePosStr(s):
    print('checkpoint')
    dsPositions = np.array([])
    
    while len(s) > 0:
        i1 = 1 + s.find('(')
        i2 = s[i1:].find(')')
        print(s)
        print(s[i1:i1+i2].split(','))
        start, stop, step = list(map(float,s[i1:i1+i2].split(',')))
        newPosArr = np.arange(start,stop+step,step)
        dsPositions = np.append(dsPositions,newPosArr)
        s = s[i1+i2+1:]
        print(s)
        
    dsPositions.sort()
    dsPositions_unique = np.unique(dsPositions.round(4)) # rounds the values to 4 decimal places and removes duplicate positions
    
    
    print(dsPositions)
    return dsPositions
        
def popupmsg(msg):
    popup = tk.Tk()
    popup.wm_title("!")
    label = ttk.Label(popup, text=msg)
    label.pack(side="top", fill="x", pady=10)
    B1 = ttk.Button(popup, text="Okay", command = popup.destroy)
    B1.pack()
    popup.mainloop()

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

class DelayStage():    
    def __init__(self):
        #import System
        instrument = "COM11"
        self.myDLS = DLS()
        result = self.myDLS.OpenInstrument(instrument)
        print('OpenInstrument error status (0 is success): ', result)
        # Get current position and print it
        self.pos, self.errorStatus = 0, ''
        self.pos = self.myDLS.PA_Get(self.pos, self.errorStatus)[1]
        print('pos = %s' % str(self.pos))
   
    def disconnect(self):
        # Close the instrument
        result = self.myDLS.CloseInstrument()
        if result == 0:
            print('Closing instrument')
        else:
            print('Failed to close instrument')
    def getPos(self):
        self.pos = self.myDLS.PA_Get(self.pos, self.errorStatus)[1]
        #print('pos = %s' % str(self.pos))
        return self.pos
    def setPos(self,newPos):
        if (newPos > 0) and (newPos < 225):
            self.myDLS.PA_Set(newPos,self.errorStatus)
            #self.pos = newPos
        else:
            print('Invalid Input Position; DS motion denied')
        
class Camera():
    def __init__(self):
        # Retrieve singleton reference to system object
        self.system = PySpin.System.GetInstance()
        
        # Get current library version
        version = self.system.GetLibraryVersion()
        print('Library version: %d.%d.%d.%d' % (version.major, version.minor, version.type, version.build))
        
        # Retrieve list of cameras from the system
        self.cam_list = self.system.GetCameras()
        num_cams = self.cam_list.GetSize()
        print('Number of cameras detected: %i' % num_cams)
        
        # Finish if there are no cameras
        if num_cams == 0:
            # Clear camera list before releasing system
            self.cam_list.Clear()
            # Release system instance
            self.system.ReleaseInstance()
            print('Not enough cameras!')
            self.cam = None
        else:
            self.cam = self.cam_list[0]
            

        try:
            # Initialize camera
            self.cam.Init()
            # Retrieve GenICam nodemap
            self.nodemap = self.cam.GetNodeMap()
            self.nodemap_tldevice = self.cam.GetTLDeviceNodeMap()
        except PySpin.SpinnakerException as ex:
            print('Error: %s' % ex)
            
        self.continue_recording = False
        self.exposure = self.cam.ExposureTime.GetValue()/1e6 # Division by 1e6 to convert from us to s
        self.gain = self.cam.Gain.GetValue()
        #print('Exposure = %s s' % self.exposure)
        #print('Gain = %s dB' % self.gain)
            
    def disconnect(self):
        
        # Deinitialize camera
        self.cam.DeInit()
        
        # Release reference to camera
        del self.cam
    
        # Clear camera list before releasing system
        self.cam_list.Clear()
    
        # Release system instance
        self.system.ReleaseInstance()
        
    def handle_close(evt):
        """
        This function will close the GUI when close event happens.
    
        :param evt: Event that occurs when the figure closes.
        :type evt: Event
        """
        self.continue_recording = False
        
    def liveAcquire(self):
        self.continue_recording = True
        # Acquire images continuously and display

        ### result &= acquire_and_display_images(cam, nodemap, nodemap_tldevice)
        
        sNodemap = self.cam.GetTLStreamNodeMap()
        
        # Change bufferhandling mode to NewestOnly
        node_bufferhandling_mode = PySpin.CEnumerationPtr(sNodemap.GetNode('StreamBufferHandlingMode'))
        if not PySpin.IsAvailable(node_bufferhandling_mode) or not PySpin.IsWritable(node_bufferhandling_mode):
            print('Unable to set stream buffer handling mode.. Aborting...')
            return False
            
        # Retrieve entry node from enumeration node
        node_newestonly = node_bufferhandling_mode.GetEntryByName('NewestOnly')
        if not PySpin.IsAvailable(node_newestonly) or not PySpin.IsReadable(node_newestonly):
            print('Unable to set stream buffer handling mode.. Aborting...')
            return False
            
        # Retrieve integer value from entry node
        node_newestonly_mode = node_newestonly.GetValue()
        
        # Set integer value from entry node as new value of enumeration node
        node_bufferhandling_mode.SetIntValue(node_newestonly_mode)
        
        try:
            node_acquisition_mode = PySpin.CEnumerationPtr(self.nodemap.GetNode('AcquisitionMode'))
            if not PySpin.IsAvailable(node_acquisition_mode) or not PySpin.IsWritable(node_acquisition_mode):
                print('Unable to set acquisition mode to continuous (enum retrieval). Aborting...')
                return False
    
            # Retrieve entry node from enumeration node
            node_acquisition_mode_continuous = node_acquisition_mode.GetEntryByName('Continuous')
            if not PySpin.IsAvailable(node_acquisition_mode_continuous) or not PySpin.IsReadable(
                    node_acquisition_mode_continuous):
                print('Unable to set acquisition mode to continuous (entry retrieval). Aborting...')
                return False
                
            # Retrieve integer value from entry node
            acquisition_mode_continuous = node_acquisition_mode_continuous.GetValue()
    
            # Set integer value from entry node as new value of enumeration node
            node_acquisition_mode.SetIntValue(acquisition_mode_continuous)
    
            print('Acquisition mode set to continuous...')
    
            #  Begin acquiring images
            #
            #  *** NOTES ***
            #  What happens when the camera begins acquiring images depends on the
            #  acquisition mode. Single frame captures only a single image, multi
            #  frame catures a set number of images, and continuous captures a
            #  continuous stream of images.
            #
            #  *** LATER ***
            #  Image acquisition must be ended when no more images are needed.
            self.cam.BeginAcquisition()
                        
            print('Acquiring images...')
    
            #  Retrieve device serial number for filename
            #
            #  *** NOTES ***
            #  The device serial number is retrieved in order to keep cameras from
            #  overwriting one another. Grabbing image IDs could also accomplish
            #  this.
            
            device_serial_number = ''
            node_device_serial_number = PySpin.CStringPtr(self.nodemap_tldevice.GetNode('DeviceSerialNumber'))
            if PySpin.IsAvailable(node_device_serial_number) and PySpin.IsReadable(node_device_serial_number):
                device_serial_number = node_device_serial_number.GetValue()
                print('Device serial number retrieved as %s...' % device_serial_number)
                                
            # Close program
            print('Press enter to close the program..')
    
            # Figure(1) is default so you can omit this line. Figure(0) will create a new window every time program hits this line
            fig = plt.figure(1)
    
            # Close the GUI when close event happens
            fig.canvas.mpl_connect('close_event', self.handle_close)
            
            # Retrieve and display images
            while(self.continue_recording):
                try:
                    #  Retrieve next received image
                    #  *** NOTES ***
                    #  Capturing an image houses images on the camera buffer. Trying
                    #  to capture an image that does not exist will hang the camera.
                    #  *** LATER ***
                    #  Once an image from the buffer is saved and/or no longer
                    #  needed, the image must be released in order to keep the
                    #  buffer from filling up.
                    image_result = self.cam.GetNextImage()
                    
                    #  Ensure image completion
                    if image_result.IsIncomplete():
                        print('Image incomplete with image status %d ...' % image_result.GetImageStatus())
                    else:
                        # Getting the image data as a numpy array
                        image_data = image_result.GetNDArray()
    
                        # Draws an image on the current figure
                        plt.imshow(image_data, cmap='gray')
    
                        # Interval in plt.pause(interval) determines how fast the images are displayed in a GUI
                        # Interval is in seconds.
                        plt.pause(0.001)
    
                        # Clear current reference of a figure. This will improve display speed significantly
                        plt.clf()

                        # If user presses enter, close the program
                        if keyboard.is_pressed('ENTER'):
                            print('Program is closing...')
    
                            # Close figure
                            plt.close('all')
                            input('Done! Press Enter to exit...')
                            self.continue_recording=False
                            
                    #  Release image
                    #  *** NOTES ***
                    #  Images retrieved directly from the camera (i.e. non-converted
                    #  images) need to be released in order to keep from filling the
                    #  buffer.
                    image_result.Release()
                    
                except PySpin.SpinnakerException as ex:
                    print('Error: %s' % ex)
                    return False
                    
            #  End acquisition
            #
            #  *** NOTES ***
            #  Ending acquisition appropriately helps ensure that devices clean up
            #  properly and do not need to be power-cycled to maintain integrity.
            self.cam.EndAcquisition()
            
        except PySpin.SpinnakerException as ex:
            print('Error: %s' % ex)
            return False
            
        return True
   
    def snapShot(self, display = True, save = True, fileName = 'default.tiff'):
        # set acquisition to single-frame; take an image; save image in directory
        node_acquisition_mode = PySpin.CEnumerationPtr(self.nodemap.GetNode('AcquisitionMode'))
        if not PySpin.IsAvailable(node_acquisition_mode) or not PySpin.IsWritable(node_acquisition_mode):
            print('Unable to set acquisition mode to single-frame (enum retrieval). Aborting...')
            return False
            
        # Retrieve entry node from enumeration node
        node_acquisition_mode_singleFrame = node_acquisition_mode.GetEntryByName('SingleFrame')
        if not PySpin.IsAvailable(node_acquisition_mode_singleFrame) or not PySpin.IsReadable(
                node_acquisition_mode_singleFrame):
            print('Unable to set acquisition mode to single-frame (entry retrieval). Aborting...')
            return False
            
        # Retrieve integer value from entry node
        acquisition_mode_singleFrame = node_acquisition_mode_singleFrame.GetValue()

        # Set integer value from entry node as new value of enumeration node
        node_acquisition_mode.SetIntValue(acquisition_mode_singleFrame)

        print('Acquisition mode set to single-frame...')
        
        print('Acquiring image...')
        self.cam.BeginAcquisition()
        image_result = self.cam.GetNextImage()
        #  Ensure image completion
        if image_result.IsIncomplete():
            print('Image incomplete with image status %d ...' % image_result.GetImageStatus())
        else:
            # Getting the image data as a numpy array
            image_data = image_result.GetNDArray()
        #  Images retrieved directly from the camera (i.e. non-converted
        #  images) need to be released in order to keep from filling the
        #  buffer.
        image_result.Release()
        #  Ending acquisition appropriately helps ensure that devices clean up
        #  properly and do not need to be power-cycled to maintain integrity.
        self.cam.EndAcquisition()
        # Show the image that you collected
        print('Image Complete')
        
        if display == True:
            plt.imshow(image_data, cmap='viridis')
            plt.show()
            
        if save == True:
            cv2.imwrite(fileName,image_data)
            
        return True
        
    def getExposure(self):
        self.exposure = self.cam.ExposureTime.GetValue()/1e6 # Division by 1e6 to convert from us to s
        return self.exposure
        
    def getGain(self):
        self.gain = self.cam.Gain.GetValue()
        return self.gain
        #print('Gain = %s dB' % self.gain)
        
    def setExposure(self,setValueSecond):
        setValue = setValueSecond*1e6 # converts to us
        # Make sure that auto exposure has been disabled
        if self.cam.ExposureAuto.GetValue() == True:
            print('Cannot set exposure because auto-exposure is currently activated. Deactivate auto-exposure to use this method.')
            return False
        # Ensure desired exposure time falls between the max and min
        if setValue > self.cam.ExposureTime.GetMax():
            print('Cannot set exposure because set value %s micro-second is too large.' % setValue)
            print('Max exposure time is %s micro-second' % self.cam.ExposureTime.GetMax())
            return False
        if setValue < self.cam.ExposureTime.GetMin():
            print('Cannot set exposure because set value %s micro-second is too small.' % setValue)
            print('Min exposure time is %s micro-second' % self.cam.ExposureTime.GetMin())
            return False
                      
        # Send command to update the exposure to the camera
        try:
            self.cam.ExposureTime.SetValue(setValue)
            return True
        except:
            return False        
        
    def setGain(self,setValue):
        
        # Make sure that auto gain has been disabled
        if self.cam.GainAuto.GetValue() == True:
            print('Cannot set gain because auto-gain is currently activated. Deactivate auto-gain to use this method.')
            return False
        # Ensure desired gain falls between the max and min
        if setValue > self.cam.Gain.GetMax():
            print('Cannot set gain because set value %s dB is too large.' % setValue)
            print('Max gain is %s dB' % self.cam.Gain.GetMax())
            return False
        if setValue < self.cam.Gain.GetMin():
            print('Cannot set gain because set value %s dB is too small.' % setValue)
            print('Min gain is %s dB' % self.cam.Gain.GetMin())
            return False

        # Send command to update the gain to the camera
        try:
            self.cam.Gain.SetValue(setValue)
            return True
        except:
            return False

    def configureForContinuousCapture(self):
        sNodemap = self.cam.GetTLStreamNodeMap()
        
        # Change bufferhandling mode to NewestOnly
        node_bufferhandling_mode = PySpin.CEnumerationPtr(sNodemap.GetNode('StreamBufferHandlingMode'))
        if not PySpin.IsAvailable(node_bufferhandling_mode) or not PySpin.IsWritable(node_bufferhandling_mode):
            print('Unable to set stream buffer handling mode.. Aborting...')
            return False
            
        # Retrieve entry node from enumeration node
        node_newestonly = node_bufferhandling_mode.GetEntryByName('NewestOnly')
        if not PySpin.IsAvailable(node_newestonly) or not PySpin.IsReadable(node_newestonly):
            print('Unable to set stream buffer handling mode.. Aborting...')
            return False
            
        # Retrieve integer value from entry node
        node_newestonly_mode = node_newestonly.GetValue()
        
        # Set integer value from entry node as new value of enumeration node
        node_bufferhandling_mode.SetIntValue(node_newestonly_mode)
        
        try:
            node_acquisition_mode = PySpin.CEnumerationPtr(self.nodemap.GetNode('AcquisitionMode'))
            if not PySpin.IsAvailable(node_acquisition_mode) or not PySpin.IsWritable(node_acquisition_mode):
                print('Unable to set acquisition mode to continuous (enum retrieval). Aborting...')
                return False
    
            # Retrieve entry node from enumeration node
            node_acquisition_mode_continuous = node_acquisition_mode.GetEntryByName('Continuous')
            if not PySpin.IsAvailable(node_acquisition_mode_continuous) or not PySpin.IsReadable(
                    node_acquisition_mode_continuous):
                print('Unable to set acquisition mode to continuous (entry retrieval). Aborting...')
                return False
                
            # Retrieve integer value from entry node
            acquisition_mode_continuous = node_acquisition_mode_continuous.GetValue()
    
            # Set integer value from entry node as new value of enumeration node
            node_acquisition_mode.SetIntValue(acquisition_mode_continuous)
    
            print('Acquisition mode set to continuous...')
            
        except PySpin.SpinnakerException as ex:
            print('Error: %s' % ex)
            return False
      

class SetupApp():
    def __init__(self, wait = .033, dsWait = 0.1):
        self.exposure = 1
        self.gain = 30
        
        # Setup the tkinter GUI window
        self.root = tk.Tk()
        self.root.geometry("1000x700")
        self.root.configure(bg='white')
        self.root.wm_title("Smart Scan")
        self.root.iconbitmap("kogar_0.ico")
        #self.root.iconify()
        
        # Connect to the camera, collect the initial image and plot the image
        # initialize camera
        self.camera = Camera()
        self.camera.setExposure(self.exposure)
        self.camera.setGain(self.gain)
        self.camera.configureForContinuousCapture()
        self.camera.cam.BeginAcquisition() # start image acquisition
        print('Acquiring Images...')
        
        # take image
        sampleImage = self.camera.cam.GetNextImage()
        #  Ensure image completion
        if sampleImage.IsIncomplete():
            print('Image incomplete with image status %d ...' % sampleImage.GetImageStatus())
        else:
            # Getting the image data as a numpy array
            self.camImage = sampleImage.GetNDArray()
            self.h, self.w = self.camImage.shape
        sampleImage.Release()
        
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
        self.wait = wait
        self.ppw = 100
        self.rm = [] # list of rois
        self.roiActive = False
        
        # Wait time for delay stage motion
        self.dsWait = dsWait
        
        # Some general formatting vars
        padyControls = 1
        padySep = 1
        entryWidth = 10
        largeEntryWidth = 20
        sepWidth = 190
        sepHeight = 1
        buttonWidth = 18
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
        
        # Make a Tkinter frame for the data plots
        self.dataFrame = tk.Frame(master=self.root, width=50, height=100, padx=5, pady=5)
        self.dataFrame.pack(side=tk.RIGHT)
        
        # Make a Tkinter frame for the scan plots
        self.scanFrame = tk.Frame(master=self.root, width=50, height=100, padx=5, pady=5, bg=cntrlBg)
        self.scanFrame.pack(side=tk.RIGHT)
        
        # Make a Tkinter frame for the ROI selection buttons
        self.controls = tk.Frame(master=self.root, width=50, height=100, padx=5, pady=5, bg=cntrlBg)
        self.controls.pack(side=tk.RIGHT)

        # This is the GUI controls section
        
        # Adding the ROI controls
        self.roiSelection = tk.StringVar(self.controls)
        self.roiSelection.set("roi rectangle") # default value
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
            width=buttonWidth,
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
            width=buttonWidth,
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
            width=buttonWidth,
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
        self.gain_label_string.set(f'Gain = {self.camera.getGain()} dB')
        self.gain_label = tk.Label(master=self.controls,textvariable=self.gain_label_string,bg=cntrlBg)
        self.gainEntry = tk.Entry(master=self.controls,width=entryWidth)
        self.gainEntry.insert(-1,self.camera.getGain())
        self.exposure_label_string = tk.StringVar()
        self.exposure_label_string.set(f'Exposure = {self.camera.getExposure()} s')
        self.exposure_label = tk.Label(master=self.controls,textvariable=self.exposure_label_string,bg=cntrlBg)
        self.exposureEntry = tk.Entry(master=self.controls,width=entryWidth)
        self.exposureEntry.insert(-1,self.camera.getExposure())
        def camSettingsButtonFunc():
            if self.gainEntry.get() != '':
                self.gain = float(self.gainEntry.get())
                self.camera.setGain(self.gain)
            else:
                self.gainEntry.insert(-1,self.camera.getGain())
                
            self.gain_label_string.set(f'Gain = {self.camera.getGain()} dB')
            
            if self.exposureEntry.get() != '':
                self.exposure = float(self.exposureEntry.get())
                self.camera.setExposure(self.exposure)
            else:
                self.exposureEntry.insert(-1,self.camera.getExposure())
                
            self.exposure_label_string.set(f'Exposure = {self.camera.getExposure()} s')
            
            
        self.camSettingsButton = tk.Button(
            master=self.controls,
            text="Update Cam Settings",
            width=buttonWidth,
            height=1,
            fg='black', bg='white',
            command=camSettingsButtonFunc)
        
        # Adding the ds controls
        self.currentDsPos_label_string = tk.StringVar()
        self.currentDsPos_label_string.set(f'Current dsPos = {self.ds.getPos()}')
        self.currentDsPos_label = tk.Label(master=self.controls,textvariable=self.currentDsPos_label_string,bg=cntrlBg)
        dsPos1_label = tk.Label(master=self.controls,text="dsPos1 (mm)",bg=cntrlBg)
        
        def dsPos1ButtonFunc():
            dsPos = float(self.dsPos1Entry.get())
            self.ds.setPos(dsPos)
            time.sleep(self.dsWait)
            self.currentDsPos_label_string.set(f'Current dsPos = {self.ds.getPos()}')
            
        self.dsPos1Entry = tk.Entry(master=self.controls,width=entryWidth)
        self.dsPos1Button = tk.Button(
            master=self.controls,
            text="Move dsPos1",
            width=buttonWidth,
            height=1,
            fg='black', bg='white',
            command=dsPos1ButtonFunc)
        dsPos2_label = tk.Label(master=self.controls,text="dsPos2 (mm)",bg=cntrlBg)
        
        def dsPos2ButtonFunc():
            dsPos = float(self.dsPos2Entry.get())
            self.ds.setPos(dsPos)
            time.sleep(self.dsWait)
            self.currentDsPos_label_string.set(f'Current dsPos = {self.ds.getPos()}')
            
        self.dsPos2Entry = tk.Entry(master=self.controls,width=entryWidth)
        self.dsPos2Button = tk.Button(
            master=self.controls,
            text="Move dsPos2",
            width=buttonWidth,
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
        
        # Adding the scan controls
        scanDsPos_label = tk.Label(master=self.scanFrame,text="dsPos Scan List: (start,end,step)",bg=cntrlBg)
        self.scanDsPosEntry = tk.Entry(master=self.scanFrame,width=largeEntryWidth)
        tzpos_label = tk.Label(master=self.scanFrame,text="time zero position (mm)",bg=cntrlBg)
        self.tzposEntry = tk.Entry(master=self.scanFrame,width=largeEntryWidth)
        self.tzposEntry.insert(-1,0)
        scanBatchSize_label = tk.Label(master=self.scanFrame,text="Batch Size",bg=cntrlBg)
        self.batchSizeEntry = tk.Entry(master=self.scanFrame,width=largeEntryWidth)
        self.batchSizeEntry.insert(-1,1)
        scanDirLoc_label = tk.Label(master=self.scanFrame,text="Scan Directory",bg=cntrlBg)
        self.scanDirLocEntry = tk.Entry(master=self.scanFrame,width=largeEntryWidth)
        self.varTransCorr = tk.IntVar()
        self.transCorrBox = tk.Checkbutton(master=self.scanFrame, text='Translation Correction',variable=self.varTransCorr, onvalue=1, offvalue=0,bg=cntrlBg)
        self.varSmartScan = tk.IntVar()
        self.smartScanBox = tk.Checkbutton(master=self.scanFrame, text='Smart Scan ®',variable=self.varSmartScan, onvalue=1, offvalue=0,bg=cntrlBg)
        
        def startScan():
            scanDir = str(self.scanDirLocEntry.get())
            if not os.path.isdir(scanDir):
                popupmsg('Please select a valid scan directory for storing data')
                return
            if len(self.rm) == 0:
                popupmsg('Please select an roi before starting')
                return
            self.camLive = 2
            
        self.startScanButton = tk.Button(
            master=self.scanFrame,
            text="Start Scan",
            width=buttonWidth,
            height=1,
            fg='black', bg='white',
            command=startScan)
        
        # Add a image progress bar
        self.imProgBar = ttk.Progressbar(master=self.controls,length=160)
        
        
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
        self.imProgBar.pack(pady=padyControls)
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
        self.avgDisBox.pack(pady=padyControls)
        self.medianFilterBox.pack(pady=padyControls)
        self.avgDirLoc_label.pack(pady=padyControls)
        self.avgDirLocEntry.pack(pady=padyControls)
        self.imageSetButton.pack(pady=padyControls)
        # ~ tk.Frame(master=self.controls, bd=100, relief='flat',height=sepHeight,width=sepWidth,bg='black').pack(side='top', pady=padySep)
        
        scanDsPos_label.pack(pady=padyControls)
        self.scanDsPosEntry.pack(pady=padyControls)
        tzpos_label.pack(pady=padyControls)
        self.tzposEntry.pack(pady=padyControls)
        scanBatchSize_label.pack(pady=padyControls)
        self.batchSizeEntry.pack(pady=padyControls)
        scanDirLoc_label.pack(pady=padyControls)
        self.scanDirLocEntry.pack(pady=padyControls)
        self.transCorrBox.pack(pady=padyControls)
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
                
        
        # Bind the enter key to activate selected entry
        self.root.bind('<Return>', userPressReturn)
        
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
            if tk.messagebox.askokcancel("Quit", "Do you want to quit?"):
                self.camLive = 0
                # ~ self.root.quit()
            
        # Activate the camera thread - this will take images on repeat asynchronously
        self.camLive = 1
        self.camFree = False # Checks whether we have successfully disconnected from the camera at the end of the program
        self.thrDataAcquisition = threading.Thread(target=self.dataAcquisitionLoop, args=(), kwargs={})
        self.thrDataAcquisition.start() # Will run takeImageTestRepeat
        
        print('check point 1')
        
        self.root.protocol("WM_DELETE_WINDOW", on_closing)
        self.root.mainloop()
        
        # If 2 we start a scan
        if self.camLive == 2:
            print('Starting Scan')
            scanDir = str(self.scanDirLocEntry.get())
            print(scanDir)
            now = datetime.datetime.now()
            nowString = now.strftime("%Y_%m_%d-%H_%M_%S")
            scanParamsDict = { \
                            'scanStartTime':nowString,\
                            'gain':str(self.gainEntry.get()),\
                            'exposure':str(self.exposureEntry.get()),\
                            'dsPositions':str(self.scanDsPosEntry.get()),\
                            'batchSize':str(self.batchSizeEntry.get()),\
                            'translationCorrection':str(self.varTransCorr.get()),\
                            'smartScan':str(self.varSmartScan.get()),\
                            'directory':scanDir}
            with open(f'{scanDir}//scanMetaData_{nowString}.txt','w') as smdf:
                for key, value in scanParamsDict.items():
                    smdf.write('%s:%s\n' % (key, value))
                    
            dsPositions = parsePosStr(str(self.scanDsPosEntry.get()))
            batchSize = int(self.batchSizeEntry.get())
            tcq = bool(self.varTransCorr.get())
            ssq = bool(self.varSmartScan.get())
            tzpos = float(self.tzposEntry.get())
            
            self.root.destroy() # Totally destroy the setup gui so we can move into the scan gui
            
            scanapp = ScanApp(camera=self.camera,ds=self.ds,dsPositions=dsPositions,batchSize=batchSize, \
                              translationCorrectionQ=tcq, smartScanQ=ssq, scanDir=scanDir, tzpos=tzpos, rm=self.rm)
        else:
            print('Disconnecting camera and delay stage')
            self.ds.disconnect() # disconnect from delay stage
            #  Ending acquisition appropriately helps ensure that devices clean up
            #  properly and do not need to be power-cycled to maintain integrity.
            self.camera.setExposure(1) # Reset exposure to 1.0 s and disconnect from the camera.
            self.camera.cam.EndAcquisition()
            self.camera.disconnect() # disconnect from camera
      
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
                        
                        
            self.lineProfile_v.set_xdata(self.lineProfile_x)
            self.lineProfile_h.set_ydata(self.lineProfile_y)
            self.lineProfilePlot_v.set_xdata(self.camImage[:,int(self.lineProfile_x)])
            self.camAx[1].set_xlim(min(self.camImage[:,int(self.lineProfile_x)])-0.01,max(self.camImage[:,int(self.lineProfile_x)])+0.01)
            self.lineProfilePlot_h.set_ydata(self.camImage[int(self.lineProfile_y),:])
            self.camAx[2].set_ylim(min(self.camImage[int(self.lineProfile_y),:])-0.01,max(self.camImage[int(self.lineProfile_y),:])+0.01)   
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
            if self.roiActive:
                self.rm[-1].updateShape(abs(event.xdata - self.rm[-1].cx)*2,abs(event.ydata - self.rm[-1].cy)*2) # width, height
                
                self.canvas.draw()
    
    def saveImage(self):
        now = datetime.datetime.now()
        nowString = now.strftime("%Y_%m_%d-%H_%M_%S")
        defaultFileName = nowString+'expo'+ str(self.camera.getExposure()) + 'gain' + str(self.camera.getGain()) + '.tiff'
        self.saveImageFileName = tk.filedialog.asksaveasfile(mode='w',initialfile=defaultFileName,defaultextension=".*",filetypes = [("image files", ".tiff")])
        if self.saveImageFileName != None:
            tifffile.imwrite(self.saveImageFileName.name,self.camImage)
        
    def dataAcquisitionLoop(self):
        thrProgBar = None # Used as variable for the progress bar thread
        while self.camLive == 1:            
            image_result = self.camera.cam.GetNextImage()
            #  Ensure image completion
            if image_result.IsIncomplete():
                print('Image incomplete with image status %d ...' % image_result.GetImageStatus())
            else:
                # Getting the image data as a numpy array
                newImage = image_result.GetNDArray()
                
                
            #  Images retrieved directly from the camera (i.e. non-converted
            #  images) need to be released in order to keep from filling the
            #  buffer.
            image_result.Release()
            
            # Check whether despeckle is applied
            if self.varMedianFilter.get() == 1:
                newImage = median_filter(newImage, size=2)
            
            # Check whether we are saving image sets
            if self.saveImageSet and self.avgDirLocEntry.get() != '':
                tifffile.imwrite(f'{self.avgDirLocEntry.get()}//{time.time()}.tiff',newImage)
            
            # Check whether image averaging is selected in which the display will be updated to an average of all the images taken    
            if self.varAvgDis.get() == 0:
                self.camImage = newImage
                self.numAvgDis = 0
            else:
                self.camImage = (newImage + self.numAvgDis*self.camImage)/(self.numAvgDis+1)
                self.numAvgDis += 1
            
            # Update the label to indicate the number of images that have been taken to create the displayed image
            self.varAvgDis_label_string.set(f'Average Images (N = {self.numAvgDis})')
            
            # Compute the image rate
            deltaTString = '%.3f' % round(time.time()-self.t1,5)
            self.camAx[0].set_title(f"Image Rate =  {deltaTString} s, ROI # = {len(self.rm)}")
            
            self.t1 = time.time()
            
            # Make sure the previous progress bar thread has concluded
            if (thrProgBar != None) and (thrProgBar.is_alive()):
                eventProgBar.set()
                thrProgBar.join()
            
            eventProgBar = threading.Event()    
            thrProgBar = threading.Thread(target=self.stepProgBar, args=(eventProgBar,self.wait + self.exposure/50), kwargs={})
            thrProgBar.start()
            
            self.camImageObj.set_data(self.camImage) # Set image data to the plot
            self.lineProfile_v.set_xdata(self.lineProfile_x)
            self.lineProfile_h.set_ydata(self.lineProfile_y)
            self.lineProfilePlot_v.set_xdata(self.camImage[:,int(self.lineProfile_x)])
            self.camAx[1].set_xlim(min(self.camImage[:,int(self.lineProfile_x)]),max(self.camImage[:,int(self.lineProfile_x)]))
            self.lineProfilePlot_h.set_ydata(self.camImage[int(self.lineProfile_y),:])
            self.camAx[2].set_ylim(min(self.camImage[int(self.lineProfile_y),:]),max(self.camImage[int(self.lineProfile_y),:]))
            
            if len(self.rm) > 0:
                roiPixelSum = 0
                roiTotalArea = 0
                # Loop through all the ROIS
                for r in self.rm:
                    # Only count an ROI if it is active
                    if r.live:
                        roiImage = newImage[int(r.cy-r.h/2):int(r.cy+r.h/2),int(r.cx-r.w/2):int(r.cx+r.w/2)]
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
            
            # Update the figures
            self.canvas.draw()
            self.datacanvas.draw()
            
            if self.resetHistory:
                self.resetHistory = False
                self.intensityHistory = []
                self.timeHistory = []
                self.avgIntensity = np.append(self.avgIntensity,[0])
                self.stdIntensity = np.append(self.stdIntensity,[0])
                
                # Zero the cam image if there is averaging
                if self.varAvgDis.get() == 1:
                    self.camImage = 0*self.camImage
                    self.numAvgDis = 0
                
                self.t0 = time.time()
                
            if self.resetLongHistory:
                self.resetLongHistory = False
                self.avgIntensity = np.array([0],dtype='float')
                self.stdIntensity = np.array([0],dtype='float')
            
        # Make sure the previous progress bar thread has concluded
        if (thrProgBar != None) and (thrProgBar.is_alive()):
            eventProgBar.set()
            thrProgBar.join()
        self.root.quit() # Close the prep GUI window when the scan is started
    
    def stepProgBar(self, event, spf):
        N = int(self.exposure/spf)
        self.imProgBar['value'] = 100/N
        while not event.is_set():
            self.imProgBar.step(100/N)
            time.sleep(spf)
            
        
        
        
class ScanApp():
    def __init__(self, camera, ds, dsPositions, batchSize, translationCorrectionQ, smartScanQ, scanDir, tzpos, rm, dsWait = 0.1, wait = .033):
        self.camera=camera
        self.ds=ds
        self.exposure=self.camera.getExposure()
        self.gain=self.camera.getGain()
        self.dsPositions=dsPositions
        self.batchSize=batchSize
        self.translationCorrectionQ=translationCorrectionQ
        self.smartScanQ=smartScanQ
        self.scanDir=scanDir
        self.dsWait = dsWait
        self.tzpos = tzpos
        self.rm = rm
        
        now = datetime.datetime.now()
        self.runLog = f'scan app started at {now.strftime("%Y_%m_%d-%H_%M_%S")} \n'
        self.runLogFileName = f'runLog_{now.strftime("%Y_%m_%d-%H_%M_%S")}.txt'
        
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
        for r in rm:
            r.updateAxis(self.camAx)
        
        # SCAF
        self.s = -1 # Counts the number of scans taken, +1 each run so starts on zero
        self.b = -1 # Counts the number of batches taken, +1 each run so starts on zero
        self.pi = 0 # Position index
        
        # Data type for saving images
        self.image_dtype = np.float32
        
        # collect the initial image
        image_result = self.camera.cam.GetNextImage()
        #  Ensure image completion
        if image_result.IsIncomplete():
            print('Image incomplete with image status %d ...' % image_result.GetImageStatus())
        else:
            # Getting the image data as a numpy array
            self.camImage = image_result.GetNDArray()
            self.h, self.w = self.camImage.shape
            
        self.zeroImage = np.zeros(np.shape(self.camImage),dtype=self.image_dtype)
        
        
        # If translation correction is active, find the initial reference location for later comparison
        self.tcorr_0 = np.zeros(2) # array that stores the reference position for translation correction
        self.tcorr = np.zeros(2)
        self.tcorr_log = np.zeros([1,2])
        self.roiTotalArea = 0
        
        for ri in range(len(self.rm)):
            r = self.rm[ri]
            # Get pixel data in the roi
            roiImage = self.camImage[int(r.cy-r.h/2):int(r.cy+r.h/2),int(r.cx-r.w/2):int(r.cx+r.w/2)]
            self.roiTotalArea += roiImage.shape[0]*roiImage.shape[1]
            # ~ tifffile.imwrite(self.scanDir + '//' + f'initialRoi_{ri}.tiff',roiImage)
            if self.translationCorrectionQ:
                popt, rms, self.guess_prms = tc.fitImageToGaussian(roiImage)
                self.tcorr_0 += popt[:2]*roiImage.shape[0]*roiImage.shape[1]
            
        self.tcorr_0 = self.tcorr_0/self.roiTotalArea
        
        # Data display settings and variables
        self.wait = wait
        self.ppw = 100
        
        # Create intensity history arrays for data plots
        self.intensityLabTime = []
        self.intensityScanTime = np.zeros(len(self.dsPositions))
        self.stdIntensityScanTime = np.zeros(len(self.dsPositions))
        self.roiScanData = np.zeros([1,len(self.dsPositions)]) # full record of intensity in the roi
        self.timeHistory = np.zeros([1,len(self.dsPositions)]) # coincident record of when data was taken
        
        self.t0 = time.time()
        self.t0start = datetime.datetime.now()
        self.t0startString = self.t0start.strftime("%Y_%m_%d-%H_%M_%S")
        self.t1 = self.t0
        
        self.scanNumber_label_string = tk.StringVar()
        self.scanNumber_label_string.set(f"scan # = {self.s}")
        self.scanNumber_label = tk.Label(master=self.controls,textvariable=self.scanNumber_label_string,bg=cntrlBg)
        
        self.batchNumber_label_string = tk.StringVar()
        self.batchNumber_label_string.set(f"batch # = {self.b}")
        self.batchNumber_label = tk.Label(master=self.controls,textvariable=self.batchNumber_label_string,bg=cntrlBg)
        
        self.dsPosition_label_string = tk.StringVar()
        self.dsPosition_label_string.set(f"dsPos = {self.ds.getPos()} mm")
        self.dsPosition_label = tk.Label(master=self.controls,textvariable=self.dsPosition_label_string,bg=cntrlBg)
        
        self.exposure_label_string = tk.StringVar()
        self.exposure_label_string.set(f'Exposure = {self.camera.getExposure()} s')
        self.exposure_label = tk.Label(master=self.controls,textvariable=self.exposure_label_string,bg=cntrlBg)
        
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
        self.scanNumber_label.pack(pady=padySep)
        self.batchNumber_label.pack(pady=padySep)
        self.dsPosition_label.pack(pady=padySep)
        self.exposure_label.pack(pady=padySep)
        self.pauseScanButton.pack(pady=padySep)
        self.gentleStopButton.pack(pady=padySep)
        self.hardStopButton.pack(pady=padySep)
        self.errorLog_label.pack(pady=padySep)
        
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
        
        # Plot the cam image and the line profile crosshairs
        # SCAF
        self.camImageObj = self.camAx.imshow(self.camImage)
        pos0 = self.camAx.get_position()

        # Attach the matplotlib image figures to a tkinter gui window
        self.canvas = FigureCanvasTkAgg(self.camFig, master=self.root)
        self.canvas._tkcanvas.pack(side=tk.RIGHT, fill=tk.BOTH, expand=1)
        self.label = tk.Label(text="")
        self.label.pack()
        
        # Create the data plot
        self.dataFig, self.dataAx = plt.subplots(2,1,figsize=(5, 7))
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
            
        # Activate the camera thread - this will take images on repeat asynchronously
        self.scanLive = True
        self.camFree = False # Checks whether we have successfully disconnected from the camera at the end of the program
        self.dsFree = False # Checks whether we have successfully disconnected from the delay stage at the end of the program
        thr = threading.Thread(target=self.dataAcquisitionLoop, args=(), kwargs={})
        thr.start() # Will run takeImageTestRepeat
        
        self.root.protocol("WM_DELETE_WINDOW", on_closing)
        self.root.mainloop()
       
    def dataAcquisitionLoop(self):
        now = datetime.datetime.now()
        self.runLog += f'dataAcquisitionLoop started at {now.strftime("%Y_%m_%d-%H_%M_%S")} \n'
        
        original_umask = os.umask(0)
        while os.path.isdir(f'{self.scanDir}//average'):
            print('Warning! scanDir is already filled, making a new scan dir')
            self.errorLog_label_string.set(self.errorLog_label_string.get() + f'W: scanDir has data in it \n changing to scanDir_new \n')
            self.runLog += f'WARNING: scanDir({self.scanDir}) has data in it. Changing scanDir to {self.scanDir}_new \n'
            self.scanDir = self.scanDir + '_new'
        else:
            self.errorLog_label_string.set(self.errorLog_label_string.get() + f'all clear')
             
        os.makedirs(f'{self.scanDir}//average',self.dirPermission) # Create a directory for the batch average files
        os.umask(original_umask)
        
        # Create the files of this batch average
        for p in range(len(self.dsPositions)):
            dsPos = self.dsPositions[p] # direction doesn't matter here; we just need the positions for the files names
            dsPosString = '%.4f' % dsPos
            tifffile.imwrite(f'{self.scanDir}//average//pos={dsPosString}.tiff',self.zeroImage)
        
        thrImage = None # Variable that contains the image analysis async thread
        thrBatchAverage = None # Variable that contains the moving batch average thread
        # Start the scan loop
        while self.scanLive:
            self.b += 1 # Record that a new batch is being taken
            self.batchNumber_label_string.set(f"batch # = {self.b}") # Write the current batch number to the gui
            # SCAF
            os.makedirs(f'{self.scanDir}//batch{self.b}') # Create a directory for the new batch files
            # Create the files of a new batch
            for p in range(len(self.dsPositions)):
                dsPos = self.dsPositions[p] # direction doesn't matter here; we just need the positions for the files names
                dsPosString = '%.4f' % dsPos
                tifffile.imwrite(f'{self.scanDir}//batch{self.b}//pos={dsPosString}.tiff',self.zeroImage)
                
            now = datetime.datetime.now()
            self.runLog += f'batch {self.b} started at {now.strftime("%Y_%m_%d-%H_%M_%S")} \n'
            # Start the new batch
            for bi in range(self.batchSize):
                # Before starting new scan check whether a gentle stop call is active
                if self.gentleStop:
                    return self.endScan() #exits the dataAcquisitionLoop
                    
                self.s += 1 # Record that a new scan is being taken
                self.scanNumber_label_string.set(f"scan # = {self.s}") # Write the current scan number to the gui
                
                now = datetime.datetime.now()
                self.runLog += f'\t scan {self.s} (bi = {bi}) started at {now.strftime("%Y_%m_%d-%H_%M_%S")} \n'
                # Start a new scan
                for p in range(len(self.dsPositions)):
                    # Before collecting a new position, check whether a hard stop call is active
                    if self.hardStop:
                        return self.endScan() #exits the dataAcquisitionLoop
                    
                    while self.paused:
                        time.sleep(self.wait) # block here if paused
                    
                    # Compute the image rate
                    deltaTString = '%.3f' % round(time.time()-self.t1,5)
                    deltaETString = '%.3f' % round(time.time()-self.t1 - self.camera.exposure - self.dsWait,5)
                    self.camTitle = self.camAx.set_title(f"Image Rate =  {deltaTString} s (et = {deltaETString} s)")
                    self.t1 = time.time()

                    # Compute the next ds pos
                    if self.s%2==0:
                        self.pi = p
                        dsPos = self.dsPositions[self.pi] # if s is even, go through the list the normal way
                    else:
                        self.pi = len(self.dsPositions)-(1+p)
                        dsPos = self.dsPositions[self.pi] # if s is odd, go through the list backwards
                    
                    # Move to the next ds pos
                    self.ds.setPos(dsPos)
                    
                    # Give a wait for the motor time
                    time.sleep(self.dsWait)
                    
                    # Write the current ds pos to the gui
                    self.dsPosition_label_string.set(f"dsPos = {self.ds.getPos()} mm")
                    
                    self.runLog += f'\t \t dsPos = {dsPos} mm, pi = {self.pi}, deltaT = {deltaTString}, ET = {deltaETString} \n'
                    
                    # Make sure the previous image thread has concluded
                    if (thrImage != None) and (thrImage.is_alive()):
                            thrImage.join()
                    # Take an image and save it
                    thrImage = threading.Thread(target=self.processImage, args=(bi,self.pi,self.s,self.b), kwargs={})
                    thrImage.start()
                    time.sleep(self.exposure)
                    
                # This point in the loop occurs when a scan has concluded
                # Make sure the image analysis thread is concluded so that it is safe to go forward in the data acquisition loop
                if (thrImage != None) and (thrImage.is_alive()):
                    print('Waiting on thrImage at the end of a scan!')
                    thrImage.join()
                    
                # Update the lab time intensity monitor
                self.intensityLabTime.append(0)
                self.intensityLabTime[self.s] = np.mean(self.roiScanData[self.s,:]) # Update the lab time plot since the scan has ended
            
            # This point in the loop occurs when a batch has concluded
            # Make sure the previous batch moving average thread has concluded
            if (thrBatchAverage != None) and (thrBatchAverage.is_alive()):
                print('Joining Batch Average Thread')
                thrBatchAverage.join()
            # Start the batch moving average thread
            thrBatchAverage = threading.Thread(target=self.movingAverageBatch, args=(self.b,), kwargs={})
            thrBatchAverage.start()
            
            now = datetime.datetime.now()
            self.runLog += f'batch {self.b} concluded at {now.strftime("%Y_%m_%d-%H_%M_%S")} \n'
            self.runLog += f'Active Threads: {threading.enumerate()} \n'
            
            # Save some things
            
            # We do a temporary save of the roiScanData and timeHistory
            np.savetxt(self.scanDir + '//' + 'roiScanData.csv',self.roiScanData,delimiter=',')
            np.savetxt(self.scanDir + '//' + 'timeHistory.csv',self.timeHistory,delimiter=',')
            
            # SCAF, save the run log to a txt file
            with open(self.scanDir + '//' + self.runLogFileName, "w") as text_file:
                text_file.write(self.runLog)
            if len(self.runLog) > 1e5:
                now = datetime.datetime.now()
                self.runLog = f'log started at {now.strftime("%Y_%m_%d-%H_%M_%S")} \n'
                self.runLogFileName = f'runLog_{now.strftime("%Y_%m_%d-%H_%M_%S")}.txt'
        
        self.endScan()        
        
    def processImage(self, bi, pi, s, b):
        t0 = time.time()
        dsPos = self.dsPositions[pi]
        t_threadStart = time.time()
        
        # Take a new image from the camera
        image_result = self.camera.cam.GetNextImage()
        #  Ensure image completion
        if image_result.IsIncomplete():
            print('Image incomplete with image status %d ...' % image_result.GetImageStatus())
        else:
            # Getting the image data as a numpy array
            self.camImage = image_result.GetNDArray()
            
        #  Images retrieved directly from the camera (i.e. non-converted
        #  images) need to be released in order to keep from filling the
        #  buffer.
        image_result.Release()
            
        t_newImage = time.time()
                
        self.runLog += f'\t \t processImage thread image acquired for pi = {pi}, s = {s} | duration =  {round(time.time()-t0,2)} \n'
        
        # Update the ROI data plots
        if len(self.rm) > 0:
            if self.translationCorrectionQ:
                # Loop through all the ROIS
                for r in self.rm:
                    roiImage = self.camImage[int(r.cy-r.h/2):int(r.cy+r.h/2),int(r.cx-r.w/2):int(r.cx+r.w/2)]
                    popt, rms, self.guess_prms = tc.fitImageToGaussian(roiImage,prevFit=self.guess_prms)
                    self.tcorr += popt[:2]*roiImage.shape[0]*roiImage.shape[1]
                    
                self.tcorr = self.tcorr/self.roiTotalArea
                correction = np.array(self.tcorr_0 - self.tcorr,dtype=np.int32)
                self.tcorr_log = np.append(self.tcorr_log,np.array([correction]),axis=0)
                self.camImage = np.roll(self.camImage,correction,axis=[1,0]) # Apply translation correction to the image with cylic boundaries.
            
            roiPixelSum = 0
            for r in self.rm:
                # Get pixel data in the roi
                roiImage = self.camImage[int(r.cy-r.h/2):int(r.cy+r.h/2),int(r.cx-r.w/2):int(r.cx+r.w/2)]                                
                # get the sum inside all the rois and the roi total area, we will average with a normalization to area
                roiPixelSum += np.sum(roiImage)
                
            if (s + 1) > self.roiScanData.shape[0]:
                self.roiScanData = np.append(self.roiScanData,np.zeros([1,len(self.dsPositions)]),axis=0) # append a new scan row to the array
                self.timeHistory = np.append(self.timeHistory,np.zeros([1,len(self.dsPositions)]),axis=0)
            
            self.roiScanData[s,pi] = roiPixelSum/self.roiTotalArea
            self.timeHistory[s,pi] = t_newImage - self.t0
            
            if s > 1:
                print(abs(self.roiScanData[s,pi]-self.intensityScanTime[pi])/(self.stdIntensityScanTime[pi]*np.sqrt(len(self.roiScanData[:,pi]))))
                
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
        self.camImageObj.set_data(self.camImage)
        self.canvas.draw()
        
        #SCAF
        # Load the image corresponding to this position from the current batch, update the moving average, save
        dsPosString = '%.4f' % dsPos
        imAvg = tifffile.imread(f'{self.scanDir}//batch{b}//pos={dsPosString}.tiff')
        imAvg = np.array((self.camImage + bi*imAvg)/(bi+1),dtype=self.image_dtype) # moving average
        tifffile.imwrite(f'{self.scanDir}//batch{b}//pos={dsPosString}.tiff',imAvg) # save the updated image
        
        now = datetime.datetime.now()
        self.runLog += f'\t \t processImage thread completed for pi = {pi}, s = {s} | duration =  {round(time.time()-t0,2)} \n'
    
    def movingAverageBatch(self, b):
        print(f'Computing moving average on batch {b}')
        for p in range(len(self.dsPositions)):
            print(f'moving batch average on position index {p}')
            dsPos = self.dsPositions[p]
            dsPosString = '%.4f' % dsPos
            imNew = tifffile.imread(f'{self.scanDir}//batch{b}//pos={dsPosString}.tiff') # Load the image corresponding to this position and the latest batch
            imAvg = tifffile.imread(f'{self.scanDir}//average//pos={dsPosString}.tiff') # Load the moving average image corresponding to this position
            imAvg = np.array((imNew + b*imAvg)/(b+1),dtype=self.image_dtype) # compute moving average
            tifffile.imwrite(f'{self.scanDir}//average//pos={dsPosString}.tiff',imAvg) # save the updated image
            
        now = datetime.datetime.now()
        self.runLog += f'movingAverageBatch thread completed for b = {b} at {now.strftime("%Y_%m_%d-%H_%M_%S")} \n'
    
    def endScan(self):        
        # Make sure all extra threads are closed out
        while len(threading.enumerate()) > 2:
            time.sleep(self.wait)
        
        # Disconnect from devices    
        self.ds.disconnect() # disconnect from delay stage
        #  Ending acquisition appropriately helps ensure that devices clean up
        #  properly and do not need to be power-cycled to maintain integrity.
        self.camera.setExposure(1) # Reset exposure to 1.0 s and disconnect from the camera.
        self.camera.cam.EndAcquisition()
        self.camera.disconnect() # disconnect from camera
        self.root.quit() # Close the prep GUI window when the scan is started
        
        # Do the final batch average
        self.movingAverageBatch(self.b)
        
        # Save some things
        np.savetxt(self.scanDir + '//' + 'roiScanData.csv',self.roiScanData,delimiter=',')
        np.savetxt(self.scanDir + '//' + 'timeHistory.csv',self.timeHistory,delimiter=',')
        
        # save the corrections log
        np.savetxt(self.scanDir + '//' + f'tcorrLog_batch.csv',self.tcorr_log,delimiter=',')
        
        # SCAF, save the run log to a txt file
        with open(self.scanDir + '//' + self.runLogFileName, "w") as text_file:
            text_file.write(self.runLog)
                
        # Disconnect from camera and clean up
        # Disconnect from delay stage and clean up
        
        # Save data plots
        fig, axs = plt.subplots(2,3)
        
        y00 = self.roiScanData.flatten()
        t00 = self.timeHistory.flatten()
        
        sorted_indices = np.argsort(t00)
        
        t00_sorted = t00[sorted_indices]
        y00_sorted = y00[sorted_indices]
        
        mean_y00 = np.mean(y00_sorted)
        std_y00 = np.std(y00_sorted)
        
        fig.suptitle(f't0={self.t0startString}')
        
        axs[0,0].plot(t00_sorted,y00_sorted)
        axs[0,0].set_xlabel('Lab Time (s)')
        # ~ axs[0,0].legend()
        # ~ axs[0,0].set_ylim([mean_y00-6*std_y00],mean_y00+6*std_y00)
        
        deltaT = t00_sorted-np.roll(t00_sorted,1)
        deltaT = deltaT[1:]
        T = np.mean(deltaT)
        deltaS = np.max(self.timeHistory)-np.min(self.timeHistory)
        S = deltaS/(self.s+1)
        N = len(y00_sorted)
        yf = fft.fft(y00_sorted)
        xf = fft.fftfreq(N, T)[:N//2]
        axs[0,1].plot(xf, 2.0/N * np.abs(yf[0:N//2]),alpha=1,c='blue')
        
        print('S',S)
        axs[0,1].axvline(x=1/S,c='red')
        axs[0,1].axvline(x=1/(2*S),c='green')
        
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
        
        scanTimes = (self.tzpos-self.dsPositions)*6.671281904
        
        axs[1,1].errorbar(scanTimes,self.intensityScanTime,self.stdIntensityScanTime,\
                   marker='o',color='green',elinewidth=1,capsize=5,linestyle='') # scan intensity plot with error bars at 1 sigma
        axs[1,1].set_xlabel('t-t0 (ps)')
        
        axs[1,2].imshow(self.camImage)
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
        if self.hardStop == True:
            t00_sorted = t00_sorted[0:len(self.tcorr_log[1:,0])]
            
        if self.translationCorrectionQ:
            print(len(t00_sorted),len(self.tcorr_log[1:,0]))
            axs[0].plot(t00_sorted,self.tcorr_log[1:,0])
            axs[0].set_ylabel('X shift)')
            axs[0].set_xlabel('Lab Time (s)')
            axs[1].plot(t00_sorted,self.tcorr_log[1:,1])
            axs[1].set_ylabel('Y shift)')
            axs[1].set_xlabel('Lab Time (s)')
            fig.suptitle(f't0={self.t0startString}')
            
            plt.tight_layout()
            fig.savefig(self.scanDir + '//' + 'translationCorrectionFigure.pdf')
        
        # Close the scan GUI window
        self.root.quit()
        
        
app = SetupApp()


      
# ~ cam = Camera()
# ~ ds = DelayStage()
# ~ dsPositions = np.array(range(10))

# ~ app = ScanApp(cam=cam,ds=ds,dsPositions=dsPositions,batchSize=5, translationCorrectionQ=True,smartScanQ=True,scanDir=r'C:\Users\thoma\OneDrive - UCLA IT Services\Documents\GitHub\smartScan\overnightScan_test',tzpos=4)

# ~ # Setup the tkinter GUI window
# ~ root = tk.Tk()
# ~ root.geometry("1400x1000")

# ~ dsPosition_label_string = tk.StringVar()
# ~ dsPosition_label_string.set("dsPos = ")
# ~ dsPosition_label = tk.Label(master=root,textvariable=dsPosition_label_string)
# ~ dsPosition_label.pack()
# ~ root.mainloop()
