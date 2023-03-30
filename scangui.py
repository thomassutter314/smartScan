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


class SetupApp():
    def __init__(self, wait = .033):
        self.wait = int(1000*wait) # Convert to ms
        self.root = tk.Tk()
        self.root.geometry("1200x800")
        self.root.configure(bg='white')
        
        self.root.wm_title("Smart Scan")
        self.camWindow = tk.Frame(master=self.root, width=200, height=200, padx=15)
        self.camWindow.pack(side=tk.RIGHT)
        
        # Make the cam plot window for the diffraction image
        
        # SCAF
        self.takeImageTest()
        
        h, w = self.image_0.shape
        gs = gridspec.GridSpec(2, 2,width_ratios=[w,w*.1], height_ratios=[h,h*.1])
        self.camFig = plt.figure()
        self.camAx = [plt.subplot(gs[0]),]
        self.camAx.append(plt.subplot(gs[1],sharey=self.camAx[0]))
        self.camAx.append(plt.subplot(gs[2],sharex=self.camAx[0]))
        
        self.camAx[0].xaxis.set_tick_params(labelbottom=False)
        self.camAx[0].yaxis.set_tick_params(labelleft=False)
        self.camAx[0].set_xticks([])
        self.camAx[0].set_yticks([])
        # ~ self.camAx[0].set_aspect('auto')
        
        self.lineProfile_x = w//2
        self.lineProfile_y = h//2
        self.lineProfile_v = self.camAx[0].axvline(x=self.lineProfile_x,c='black')
        self.lineProfile_h = self.camAx[0].axhline(y=self.lineProfile_y,c='black')
        
        # SCAF
        self.camImage = self.camAx[0].imshow(self.image_0)
        self.lineProfilePlot_v, = self.camAx[1].plot(self.image_0[:,int(self.lineProfile_x)],range(h))
        self.lineProfilePlot_h, = self.camAx[2].plot(range(w),self.image_0[int(self.lineProfile_y),:])
        
        pos0 = self.camAx[0].get_position()
        pos1 = self.camAx[1].get_position()
        pos2 = self.camAx[2].get_position()
        self.camAx[1].set_position([pos1.x0,pos0.y0,pos1.width,pos0.height])
        self.camAx[2].set_position([pos0.x0,pos2.y0,pos0.width,pos2.height])
        
        # key press event on cam figure
        def cam_on_click(event):
            # Find the location of the mouse click and update the crosshair variable to this value
            if event.inaxes == self.camAx[1]:
                self.lineProfile_y = event.ydata
            if event.inaxes == self.camAx[2]:
                self.lineProfile_x = event.xdata

        # Connect the cam figure to key press events
        self.camFig.canvas.mpl_connect('button_press_event', cam_on_click)

        # Attach the matplotlib figures to a tkinter gui window
        self.canvas = FigureCanvasTkAgg(self.camFig, master=self.root)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.root )
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.label = tk.Label(text="")
        self.label.pack()
        
        def on_closing():
            # ~ if tk.messagebox.askokcancel("Quit", "Do you want to quit?"):
            self.root.quit()
            
            
        thr = threading.Thread(target=self.takeImageTestRepeat, args=(), kwargs={})
        thr.start() # Will run takeImageTestRepeat
        
        self.update()
        self.root.protocol("WM_DELETE_WINDOW", on_closing)
        self.root.mainloop()
       
    def takeImageTest(self):
        x, y = np.linspace(-3,3,300), np.linspace(-3,3,300)
        X, Y = np.meshgrid(x,y)
        time.sleep(1)
        self.image_0 = np.exp(-(X**2/4+Y**2)/.2)*np.cos((X**2+Y**2)*10)**2
        print('ok')
        
    def takeImageTestRepeat(self):
        while True:
            time.sleep(.5)
            self.image_0 += 0.01*np.random.random(np.shape(self.image_0))
            print('ok')
    
    def update(self):
        self.lineProfile_v.set_xdata(self.lineProfile_x)
        self.lineProfile_h.set_ydata(self.lineProfile_y)
        
        self.lineProfilePlot_v.set_xdata(self.image_0[:,int(self.lineProfile_x)])
        self.camAx[1].set_xlim(min(self.image_0[:,int(self.lineProfile_x)]),max(self.image_0[:,int(self.lineProfile_x)]))
        self.lineProfilePlot_h.set_ydata(self.image_0[int(self.lineProfile_y),:])
        self.camAx[2].set_ylim(min(self.image_0[int(self.lineProfile_y),:]),max(self.image_0[int(self.lineProfile_y),:]))
        
        self.camImage.set_data(self.image_0)
        
        self.canvas.draw()
        self.root.after(self.wait, self.update)
        
        
        
app = SetupApp()
