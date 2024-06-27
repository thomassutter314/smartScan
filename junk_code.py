
"""
print('________________________________')
print(self.camera.getExposure(),self.camera.getGain())
print('temp',self.camera.cam.DeviceTemperature())
print('power supply current',self.camera.cam.PowerSupplyCurrent())
print('power supply voltage',self.camera.cam.PowerSupplyVoltage())


gs = self.camFig.add_gridspec(2, 2,  width_ratios=(10, 1), height_ratios=(10, 1),
            left=0.1, right=0.9, bottom=0.1, top=0.9,
            wspace=0.05, hspace=0.05)

self.camAx = [plt.subplot(gs[0]),]
self.camAx.append(plt.subplot(gs[1],sharey=self.camAx[0]))
self.camAx.append(plt.subplot(gs[2],sharex=self.camAx[0]))

self.camAx = [self.camFig.add_subplot(gs[1, 0])]
self.camAx.append(self.camFig.add_subplot(gs[0, 0], sharex=self.camAx[0]))
self.camAx.append(self.camFig.add_subplot(gs[1, 1], sharey=self.camAx[0]))


# ~ Configure some aspects of the plot window
gs = gridspec.GridSpec(2, 2, width_ratios=[self.w,self.w*.1], height_ratios=[self.h,self.h*.1])
gs = gridspec.GridSpec(2, 2, width_ratios=[4,1], height_ratios=[1,4])



self.lineProfile_v = self.camAx[0].axvline(x=self.lineProfile_x,c='black',linestyle = 'dashed')




#  Not sure what this code as supposed to be doing...... was on line 816 in scangui.py
pos0 = self.camAx[0].get_position()
pos1 = self.camAx[1].get_position()
pos2 = self.camAx[2].get_position()
self.camAx[1].set_position([pos1.x0,pos0.y0,pos1.width,pos0.height])
self.camAx[2].set_position([pos0.x0,pos2.y0,pos0.width,pos2.height])



# was on line 1197 of scangui.py
self.lineProfile_v.set_xdata([self.lineProfile_x])
self.lineProfile_h.set_ydata([self.lineProfile_y])
self.lineProfile_h.set_xy1(x=self.lineProfile_x, y=self.lineProfile_y)
self.lineProfile_h.set_xy2(x=self.lineProfile_x+1, y=self.lineProfile_y)

self.lineProfilePlot_v.set_xdata(self.displayImage[:,int(self.lineProfile_x)])
self.camAx[2].set_xlim(min(self.displayImage[:,int(self.lineProfile_x)]),max(self.displayImage[:,int(self.lineProfile_x)]))
self.lineProfilePlot_h.set_ydata(self.displayImage[int(self.lineProfile_y1),:])



                
line_profile_data = self.lineProfilePlot_h.get_ydata()

self.camAx[1].set_ylim(min(line_profile_data),max(line_profile_data)) # set the scales for us to see the data
self.camAx[1].set_yticks([min(line_profile_data),0.5*(max(line_profile_data)+min(line_profile_data)),max(line_profile_data)])
                


# ~ # Always set the vertical line at a perpendicular to the horizontal line
self.lineProfile_v.set_xy2(x=self.lineProfile_x1-(self.lineProfile_y2 - self.lineProfile_y1), \
                           y=self.lineProfile_y1 + (self.lineProfile_x2 - self.lineProfile_x1))
                           
                           
"""

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

labels = ['151491', '151493 + 151494', '151495', '151496', '151497', '151498', '151499', '151500', '151501']
I = np.array([0, 1 ,2, 3, 4, 5, 6, 7, 8]) # in mA
R = np.array([2.46, 1.9, 1.5, 1.275, 1.025, 0.860, 0.745, 0.650, 0.560]) # in kOhm
V = I*R


tau = [10249, 5392, 11807, 12103, 4872, 4615, 3963, 2775, 2543]
beta = [1.21, 1.71, 1.69, 0.49, 1.8, 1.64, 1.69, 0.91, 1.72]


#taus
# ~ 0, 10249
# ~ 1, 3128 (5392)
# ~ 2, 11807
# ~ 3, inf (12103)
# ~ 4, 4872
# ~ 5, 4615
# ~ 6, 3963
# ~ 7, 2775 (-)
# ~ 8, 2543

#betas
# ~ 0, 1.21
# ~ 1, 1.71
# ~ 2, 1.69
# ~ 3, 0.49
# ~ 4, 1.8
# ~ 5, 1.64
# ~ 6, 1.69
# ~ 7, 0.91
# ~ 8, 1.72


fig, axs = plt.subplots(1,3)

for i in range(len(I)):
    axs[0].scatter(I[i], V[i], marker = 's', color = cm.jet(i/(len(I)-1)), label = f'uid = {labels[i]}')
    axs[1].scatter(I[i], tau[i], marker = 's', color = cm.jet(i/(len(I)-1)), label = f'uid = {labels[i]}')
    axs[2].scatter(I[i], beta[i], marker = 's', color = cm.jet(i/(len(I)-1)), label = f'uid = {labels[i]}')

fig.suptitle('For constant current scans on 6/24/2024')

# ~ axs[0].set_xlabel('Current (mA)')


# ~ axs[0].legend()

axs[0].set_xlabel('Current (mA)')
axs[0].set_ylabel('Voltage (V) (min)')

axs[1].set_yscale('log')
axs[1].set_ylabel(r'$\tau$ (s)')
axs[1].set_xlabel('Current (mA)')

axs[2].set_ylabel(r'$\beta$')
axs[2].set_xlabel('Current (mA)')

plt.show()
# ~ 2.46Kohms
