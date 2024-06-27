import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import fft
from scipy.ndimage import median_filter
from scipy.ndimage import uniform_filter
from scipy.ndimage import variance
from scipy import optimize as opt


def lee_filter(img, size):
    img_mean = uniform_filter(img, (size, size))
    img_sqr_mean = uniform_filter(img**2, (size, size))
    img_variance = img_sqr_mean - img_mean**2

    overall_variance = variance(img)

    img_weights = img_variance / (img_variance + overall_variance)
    img_output = img_mean + img_weights * (img - img_mean)
    # ~ img_output = img + img_weights*(img_mean - img)
    return img_output

def func_erfexpdecay(x, x0, sigma, A, gamma):
    # ~ decay = np.ones(len(x))
    # ~ mask = x < x0
    # ~ decay[mask] = np.exp(gamma*(x[mask]-x0))
    decay = 0.5*(1 + np.tanh(gamma*x))
    return 1 - A/2*(1+scipy.special.erf((x0-x)/(np.sqrt(2)*sigma))) * decay
                               
def findtzpos(xdata, ydata):
    yrmserr = []
    ydata = (ydata - min(ydata))/(max(ydata)-min(ydata)) # Normalize and scale data
    
    plt.scatter(xdata,ydata)
    plt.figure()
    
    for i in range(len(ydata)):
        yrmserr.append(np.sqrt(np.mean((ydata[:i+1] - 1)**2)))
    
    yrmserr_filtered = scipy.signal.medfilt(yrmserr, kernel_size=3)
    
    plt.scatter(xdata, yrmserr, c = 'blue')
    plt.scatter(xdata, yrmserr_filtered, c = 'red')
    plt.show()

def malusFunc(x,A,phi,offset,invert=True):
    if invert == False:
        return A*(np.cos(np.pi/180*(2*x-phi)))**2 + offset
    else:
        if x < offset:
            return None
        acos_arg = np.sqrt((x-offset)/A)
        if acos_arg**2 > 1:
            return None
        return 0.5*(phi + 180/np.pi * np.arccos(acos_arg))

def getLineProfile(image, p1, p2):
    x1, y1 = p1
    x2, y2 = p2
    
    # avoids 90 degree lines
    if x2 == x1:
        x2 += 1
    # x2 should always be greater than x1, if not, we need to perform a reflection that preserves the line but moves x2, y1
    if x2 < x1:
        x2 += 2*(x1 - x2) 
        y2 += 2*(y1 - y2)
    
    theta = np.arctan2(y2-y1,x2-x1)
    x = np.arange(0, image.shape[1], 1)
    y = (y2 - y1)/(x2 - x1) * (x - x1) + y1
    
    
    # ~ if abs(theta) < np.pi/4:

    # ~ else:
        # ~ y = np.arange(0, image.shape[0], 1)[::-1]
        # ~ x = (x2 - x1)/(y2 - y1) * (y - y1) + x1
        
    x, y = np.array(x, dtype = int), np.array(y, dtype = int)
    
    # ~ x[x < 0] = 0
    # ~ y[y < 0] = 0
    
    line_profile = image[y%image.shape[0], x%image.shape[1]]
    line_profile[y < 0] = 0
    line_profile[x < 0] = 0
    line_profile[y >= len(line_profile)] = 0
    line_profile[x >= len(line_profile)] = 0
    
    return line_profile
