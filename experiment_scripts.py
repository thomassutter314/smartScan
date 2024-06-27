import keyboard
import lab_instruments
import numpy as np
import time
import datetime
import os
import tifffile

def fluenceScan(dsPosList, hwpPosList, exposure, gain, scanDir, dsWait = 0.5, hwpWait = 1, loopsPerHwpPos = 3):
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
    nowString = datetime.datetime.now().strftime("%Y_%m_%d-%H_%M_%S")
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
    
    # Wait here for user to begin scan
    print('________________________________')
    beginCode = input('PRESS ENTER TO START FLUENCE SCAN \n')
    
    if beginCode != '':
        scanLive = False
        
    t_start = time.time()
    while scanLive:
        # Increment the loop number, create a directory to store the images of this loop, populate with empty images
        L += 1
        os.makedirs(f'{scanDir}//loop{L}')
        # Create the files of this loop
        for f in range(len(hwpPosList)):
            for p in range(len(dsPosList)):
                tifffile.imwrite(f'{scanDir}//loop{L}//fi={f}_pi={p}.tiff',zeroImage)
        
        # Run log for this loop
        runLog = f'Loop {L} started at {datetime.datetime.now().strftime("%Y_%m_%d-%H_%M_%S")} \n'
        runLogFileName = f'runLog_{L}.txt'  
        
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
                    
                    readBackDsPos = ds.getPos()
                    readBackHwpPos = hwp.getPos()
                    # Record this info in the run log
                    runLog += f'\t fi = {fi}, pi = {pi}, li = {li}: dsSetPos = {dsPosList[pi]}, dsReadBack = {readBackDsPos}, hwpSetPos = {hwpPosList[fi]}, hwpReadBack = {readBackHwpPos}, t-t_start = {round(time.time() - t_start,2)}\n'
                    
                    # In this loop directory, loads the image associated with this fi and pi
                    imAvg = tifffile.imread(f'{scanDir}//loop{L}//fi={fi}_pi={pi}.tiff')
                    imAvg = np.array((newImage + li*imAvg)/(li+1),dtype=image_dtype) # moving average
                    tifffile.imwrite(f'{scanDir}//loop{L}//fi={fi}_pi={pi}.tiff',imAvg) # save the updated image
        
        # Save the run log for this loop
        with open(scanDir + '//' + runLogFileName, "w") as text_file:
            text_file.write(runLog)
        
    # Disconnect from devices    
    ds.disconnect() # disconnect from delay stage
    #  Ending acquisition appropriately helps ensure that devices clean up
    #  properly and do not need to be power-cycled to maintain integrity.
    camera.setExposure(1) # Reset exposure to 1.0 s and disconnect from the camera.
    camera.disconnect() # disconnect from camera
    # disconnect from hwp
    hwp.disconnect()

def heaterlineScan(hwpPosList, exposure, gain, N_per_pos, scanDir, hwp_wait = 1):
    # Initialize camera
    camera = lab_instruments.Camera()
    time.sleep(1)
    camera.setExposure(exposure)
    camera.setGain(gain)
            
    # Connect to the hwp
    hwp = lab_instruments.HalfWavePlate(com = 'COM8')
    
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
    nowString = datetime.datetime.now().strftime("%Y_%m_%d-%H_%M_%S")
    scanParamsDict = { \
                    'scanStartTime':nowString,\
                    'gain':str(gain),\
                    'exposure':str(exposure),\
                    'hwpPositions':str(hwpPosList),\
                    'N_per_pos':str(N_per_pos),\
                    'hwp_wait (s)':str(hwp_wait),\
                    'directory':scanDir}
                    
    with open(f'{scanDir}//scanMetaData_{nowString}.txt','w') as smdf:
        for key, value in scanParamsDict.items():
            smdf.write('%s:%s\n' % (key, value))
    
    original_umask = os.umask(0)

    #pi will be position index, fi will be fluence index, li will be loop index (loops at the same fluence)
    L = -1 # This counts which loop we are on
    
    # Wait here for user to begin scan
    print('________________________________')
    beginCode = input('PRESS ENTER TO START FLUENCE SCAN \n')
    

    t_start = time.time()
    
    for hwp_i in range(len(hwpPosList)):
        hwp_pos = hwpPosList[hwp_i]
        print(f'hwp_i = {hwp_i}, hwp_pos = {hwp_pos}')
        
        # Move to the new hwp position
        hwp.moveAbsolute(hwp_pos) # Sends command to the HWP
        hwp.getPos()
        time.sleep(hwp_wait) # Wait for hwp_wait many seconds
        
        image = np.zeros(np.shape(zeroImage), dtype = np.float32)
        
        for i in range(N_per_pos):   
            print(f'\t image index = {i} out of {N_per_pos}, hwp.getPos() -> {hwp.getPos()}')         
            # Sends the software trigger to the camera
            camera.trigger()
            
            if i > 0:
                # moving average from the previously taken image (doing it this way speeds things up because the camera is capturing while we are doing the math)
                image = (new_image + (i-1)*image)/i 
                
            new_image = camera.grabImage()
        
        # moving average on the final image
        image = (new_image + (N_per_pos-1)*image)/N_per_pos 
        
            
        tifffile.imwrite(f'{scanDir}//hwp_i={hwp_i}.tiff',image)
            
    # Disconnect from devices    
    #  Ending acquisition appropriately helps ensure that devices clean up
    #  properly and do not need to be power-cycled to maintain integrity.
    camera.setExposure(1) # Reset exposure to 1.0 s and disconnect from the camera.
    camera.disconnect() # disconnect from camera
    # disconnect from hwp
    hwp.disconnect()
    

if __name__ == '__main__':
    # ~ fluenceScan([61,62,63],[71,72,73],1,30,r'C:\Users\thoma\OneDrive\Documents\GitHub\smartScan\testDir')
    # ~ hwpPosList = np.linspace(52.43,97.43,200)
    hwpPosList = np.linspace(60,97.43,50,endpoint=True)
    heaterlineScan(hwpPosList = hwpPosList, exposure = 16, gain = 30, N_per_pos = 54, scanDir = r'C:\Users\Kogar\Documents\electron_beam_photos\scans\2024_4_1_heaterlineScan')
