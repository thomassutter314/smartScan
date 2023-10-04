import keyboard
import lab_instruments
import numpy as np
import time
import datetime
import os
import tifffile

def fluenceScan(dsPosList, hwpPosList, exposure, gain, scanDir, dsWait = 0.5, hwpWait = 60, loopsPerHwpPos = 10):
    scanLive = True
    # Initialize camera
    camera = lab_instruments.Camera()
    time.sleep(1)
    camera.setExposure(exposure)
    camera.setGain(gain)
            
    # Connect to the delay stage
    ds = lab_instruments.DelayStage()
    
    # Connect to the hwp
    hwp = lab_instruments.HalfWavePlate()
    
    time.sleep(1)
    
    # take a sample image
    camera.trigger()
    camImage = camera.grabImage()
    h, w = camImage.shape
    image_dtype = np.float32
    zeroImage = np.zeros(np.shape(camImage),dtype=image_dtype) # A zero image that can be used later
    
    print('Starting Fluence Scan')
    print(f'scanDir = {scanDir}')
    
    # Record the metadata of the scan
    now = datetime.datetime.now()
    nowString = now.strftime("%Y_%m_%d-%H_%M_%S")
    scanParamsDict = { \
                    'scanStartTime':nowString,\
                    'gain':str(gain),\
                    'exposure':str(exposure),\
                    'dsPositions':str(dsPosList),\
                    'hwpPositions':str(hwpPosList),\
                    'loopsPerHwpPos':str(loopsPerHwpPos),\
                    'hwpWait (s)':str(hwpWait),\
                    'directory':scanDir}
    with open(f'{scanDir}//scanMetaData_{nowString}.txt','w') as smdf:
        for key, value in scanParamsDict.items():
            smdf.write('%s:%s\n' % (key, value))
            
    original_umask = os.umask(0)

    #pi will be position index, fi will be fluence index, li will be loop index (loops at the same fluence)
    L = -1 # This counts which loop we are on
    while scanLive:
        # Increment the loop number, create a directory to store the images of this loop, populate with empty images
        L += 1
        os.makedirs(f'{scanDir}//loop{L}')
        # Create the files of this loop
        for f in range(len(hwpPosList)):
            for p in range(len(dsPosList)):
                tifffile.imwrite(f'{scanDir}//loop{L}//fi={f}_pi={p}.tiff',zeroImage)
        
        for f in range(len(hwpPosList)):
            if L%2 == 0:
                fi = f
            else:
                fi = len(hwpPosList) - (1+f)
                
            # move to next hwp pos
            hwp.moveAbsolute(hwpPosList[fi])
            print(f'fi = {fi}, hwpPos = {hwpPosList[fi]}')
            print(f'Waiting for {hwpWait} seconds at this hwp position so that temp is steady')
            print(f'hwp read back position = {hwp.getPos()}')
            time.sleep(hwpWait)
            for li in range(loopsPerHwpPos):
                for pi in range(len(dsPosList)):
                    # check for user attempting to exit
                    if keyboard.is_pressed('ESC'):
                        scanLive = False
                        print('Gentle Stop Activated')
                    
                    ds.setPos(dsPosList[pi])
                    print(f'pi = {pi}, dsPos = {dsPosList[pi]}')
                    time.sleep(dsWait)
                    
                    # Trigger the camera
                    camera.trigger()
                    # code waits here until new image is available
                    newImage = camera.grabImage()
                    
                    # In this loop directory, loads the image associated with this fi and pi
                    imAvg = tifffile.imread(f'{scanDir}//loop{L}//fi={fi}_pi={pi}.tiff')
                    imAvg = np.array((newImage + li*imAvg)/(li+1),dtype=image_dtype) # moving average
                    tifffile.imwrite(f'{scanDir}//loop{L}//fi={fi}_pi={pi}.tiff',imAvg) # save the updated image
  





fluenceScan([1,2,3],[4,5,6],1,30,r'C:\Users\thoma\OneDrive\Documents\GitHub\smartScan\New folder')
