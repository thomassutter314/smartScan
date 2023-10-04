import time
import numpy as np

TESTING = True

if TESTING == False:
    import pyvisa
    import PySpin
    import clr
    from clr import System
    full_filename = r'C:\Windows\Microsoft.NET\assembly\GAC_64\Newport.DLS.CommandInterface\v4.0_1.0.1.0__90ac4f829985d2bf\Newport.DLS.CommandInterface.dll'
    clr.AddReference(full_filename)
    from CommandInterfaceDLS import *
else:
    print('TESTING = True')


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
    def setPos(self, newPos):
        if (newPos > 0) and (newPos < 225):
            self.myDLS.PA_Set(newPos,self.errorStatus)
            #self.pos = newPos
        else:
            print('Invalid Input Position; DS motion denied')

class Camera():
    def __init__(self):        
        # Total number of buffers
        NUM_BUFFERS = 3
        # Number of triggers
        NUM_TRIGGERS = 6
        # Total number of loops
        NUM_LOOPS = 9
        
        result = True
        
        # Retrieve singleton reference to system object
        self.system = PySpin.System.GetInstance()
    
        # Retrieve list of cameras from the system
        self.cam_list = self.system.GetCameras()
    
        num_cameras = self.cam_list.GetSize()
    
        # Finish if there are no cameras
        if num_cameras == 0:
            # Clear camera list before releasing system
            self.cam_list.Clear()
    
            # Release system
            self.system.ReleaseInstance()
    
            print("No Camera Found")
            return False
            
        print(f"Selecting 1st Camera of {num_cameras} camera/s")
        self.cam = self.cam_list.GetByIndex(0)
        
        # Retrieve TL device nodemap and print device information
        self.nodemap_tldevice = self.cam.GetTLDeviceNodeMap()
        result &= self.print_device_info()
        
        # Initialize camera
        self.cam.Init()
        
        # Setup Newest Only Buffer Mode
        # Retrieve Stream Parameters device nodemap
        s_node_map = self.cam.GetTLStreamNodeMap()
        
        # Retrieve Buffer Handling Mode Information
        handling_mode = PySpin.CEnumerationPtr(s_node_map.GetNode('StreamBufferHandlingMode'))
        if not PySpin.IsAvailable(handling_mode) or not PySpin.IsWritable(handling_mode):
            print('Unable to set Buffer Handling mode (node retrieval). Aborting...\n')
            return False

        handling_mode_entry = PySpin.CEnumEntryPtr(handling_mode.GetCurrentEntry())
        if not PySpin.IsAvailable(handling_mode_entry) or not PySpin.IsReadable(handling_mode_entry):
            print('Unable to set Buffer Handling mode (Entry retrieval). Aborting...\n')
            return False

        # Set stream buffer Count Mode to manual
        stream_buffer_count_mode = PySpin.CEnumerationPtr(s_node_map.GetNode('StreamBufferCountMode'))
        if not PySpin.IsAvailable(stream_buffer_count_mode) or not PySpin.IsWritable(stream_buffer_count_mode):
            print('Unable to set Buffer Count Mode (node retrieval). Aborting...\n')
            return False

        stream_buffer_count_mode_manual = PySpin.CEnumEntryPtr(stream_buffer_count_mode.GetEntryByName('Manual'))
        if not PySpin.IsAvailable(stream_buffer_count_mode_manual) or not PySpin.IsReadable(stream_buffer_count_mode_manual):
            print('Unable to set Buffer Count Mode entry (Entry retrieval). Aborting...\n')
            return False

        stream_buffer_count_mode.SetIntValue(stream_buffer_count_mode_manual.GetValue())
        print('Stream Buffer Count Mode set to manual...')

        # Retrieve and modify Stream Buffer Count
        buffer_count = PySpin.CIntegerPtr(s_node_map.GetNode('StreamBufferCountManual'))
        if not PySpin.IsAvailable(buffer_count) or not PySpin.IsWritable(buffer_count):
            print('Unable to set Buffer Count (Integer node retrieval). Aborting...\n')
            return False
        
        # Display Buffer Info
        print('\nDefault Buffer Handling Mode: %s' % handling_mode_entry.GetDisplayName())
        print('Default Buffer Count: %d' % buffer_count.GetValue())
        print('Maximum Buffer Count: %d' % buffer_count.GetMax())
        
        buffer_count.SetValue(NUM_BUFFERS)
        print('Buffer count now set to: %d' % buffer_count.GetValue())
        print('\nCamera will be triggered %d times in a row before %d images will be retrieved' % (NUM_TRIGGERS,(NUM_LOOPS-NUM_TRIGGERS)))

        handling_mode_entry = handling_mode.GetEntryByName('NewestOnly')
        handling_mode.SetIntValue(handling_mode_entry.GetValue())
        handling_mode_entry = PySpin.CEnumEntryPtr(handling_mode.GetCurrentEntry())
        print('\nCurrent Buffer Handling Mode: %s' % handling_mode_entry.GetDisplayName())

        # Configure trigger
        print("*** CONFIGURING SOFTWARE TRIGGER ***\n")
        # Ensure trigger mode off
        # The trigger must be initially disabled in order to configure whether the source is software or hardware.
        self.nodemap = self.cam.GetNodeMap()
        node_trigger_mode = PySpin.CEnumerationPtr(self.nodemap.GetNode("TriggerMode"))
        if not PySpin.IsAvailable(node_trigger_mode) or not PySpin.IsReadable(node_trigger_mode):
            print("Unable to disable trigger mode (node retrieval). Aborting...")
            return False
        node_trigger_mode_off = node_trigger_mode.GetEntryByName("Off")
        if not PySpin.IsAvailable(node_trigger_mode_off) or not PySpin.IsReadable(node_trigger_mode_off):
            print("Unable to disable trigger mode (enum entry retrieval). Aborting...")
            return False
    
        node_trigger_mode.SetIntValue(node_trigger_mode_off.GetValue())
        print("Trigger mode disabled...")
        
        # Select trigger source as software trigger
        node_trigger_source = PySpin.CEnumerationPtr(self.nodemap.GetNode("TriggerSource"))
        if not PySpin.IsAvailable(node_trigger_source) or not PySpin.IsWritable(node_trigger_source):
            print("Unable to get trigger source (node retrieval). Aborting...")
            return False
            
        node_trigger_source_software = node_trigger_source.GetEntryByName("Software")
        if not PySpin.IsAvailable(node_trigger_source_software) or not PySpin.IsReadable(node_trigger_source_software):
                print("Unable to set trigger source (enum entry retrieval). Aborting...")
                return False
        node_trigger_source.SetIntValue(node_trigger_source_software.GetValue())
        
        # Turn trigger mode on
        # Once the appropriate trigger source has been set, turn trigger mode
        # on in order to retrieve images using the trigger.
        node_trigger_mode_on = node_trigger_mode.GetEntryByName("On")
        if not PySpin.IsAvailable(node_trigger_mode_on) or not PySpin.IsReadable(node_trigger_mode_on):
            print("Unable to enable trigger mode (enum entry retrieval). Aborting...")
            return False

        node_trigger_mode.SetIntValue(node_trigger_mode_on.GetValue())
        print("Trigger mode turned back on...")
        
        # Set acquisition mode to continuous
        # In order to access the node entries, they have to be casted to a pointer type (CEnumerationPtr here)
        node_acquisition_mode = PySpin.CEnumerationPtr(self.nodemap.GetNode("AcquisitionMode"))
        if not PySpin.IsAvailable(node_acquisition_mode) or not PySpin.IsWritable(node_acquisition_mode):
            print("Unable to set acquisition mode to continuous (enum retrieval). Aborting...")
            return False
            
        # Retrieve entry node from enumeration node
        node_acquisition_mode_continuous = node_acquisition_mode.GetEntryByName("Continuous")
        if not PySpin.IsAvailable(node_acquisition_mode_continuous) or not PySpin.IsReadable(node_acquisition_mode_continuous):
            print("Unable to set acquisition mode to continuous (entry retrieval). Aborting...")
            return False
            
        # Retrieve integer value from entry node
        acquisition_mode_continuous = node_acquisition_mode_continuous.GetValue()
        
        # Set integer value from entry node as new value of enumeration node
        node_acquisition_mode.SetIntValue(acquisition_mode_continuous)
        
        print("Acquisition mode set to continuous...")
        
        #  Begin acquiring images
        self.cam.BeginAcquisition()

    def print_device_info(self):
        
        print("*** DEVICE INFORMATION ***\n")
    
        try:
            result = True
            node_device_information = PySpin.CCategoryPtr(self.nodemap_tldevice.GetNode("DeviceInformation"))
    
            if PySpin.IsAvailable(node_device_information) and PySpin.IsReadable(node_device_information):
                features = node_device_information.GetFeatures()
                for feature in features:
                    node_feature = PySpin.CValuePtr(feature)
                    print("%s: %s" % (node_feature.GetName(),
                                      node_feature.ToString() if PySpin.IsReadable(node_feature) else "Node not readable"))
    
            else:
                print("Device control information not available.")
    
        except PySpin.SpinnakerException as ex:
            print("Error: %s" % ex)
            return False
    
        return result
        
    def disconnect(self):
        self.cam.EndAcquisition()
        del self.cam
        self.cam_list.Clear()
        self.system.ReleaseInstance()
        print('Camera Disconnected')
    
    def trigger(self):
        # Execute software trigger
        node_softwaretrigger_cmd = PySpin.CCommandPtr(self.nodemap.GetNode("TriggerSoftware"))
        if not PySpin.IsAvailable(node_softwaretrigger_cmd) or not PySpin.IsWritable(node_softwaretrigger_cmd):
            print("Unable to execute trigger. Aborting...")
            return False
        node_softwaretrigger_cmd.Execute()
        
    def grabImage(self):
        image = self.cam.GetNextImage()
        #  Ensure image completion
        if image.IsIncomplete():
            print('Image incomplete with image status %d ...' % image.GetImageStatus())
        else:
            # Getting the image data as a numpy array
            imageArray = image.GetNDArray()
        
        image.Release()
        
        return imageArray
    
    def getExposure(self):
        self.exposure = self.cam.ExposureTime.GetValue()/1e6 # Division by 1e6 to convert from us to s
        return self.exposure
        
    def getGain(self):
        self.gain = self.cam.Gain.GetValue()
        return self.gain
        
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

class HalfWavePlate():
    def __init__(self):
        self.conv = 0.00251175 # Conversion factor between device output and degrees
        self.rm = pyvisa.ResourceManager()
        self.rm.list_resources()
        self.rotationMount = self.rm.open_resource('COM10')
        #COM 10 for pump
        self.rotationMount.read_termination = '\r\n'
        self.rotationMount.write_termination = '\r\n'
        print('Device status: ' + str(self.rotationMount.query('0gs')))
        print('Device settings: ' + str(self.rotationMount.query('0i1')))
        self.pos = self.rotationMount.query('0gp')
        print('raw output = ' + str(self.pos),' ,pos = ' + str(self.conv*self.twos_complement(self.pos[3:],32)) + ' deg')
        
    def twos_complement(self, hexstr, bits):
        value = int(hexstr,16)
        if value & (1 << (bits-1)):
            value -= 1 << bits
        return value

    def hexStringFromDec(self, dec, byteNumber = 32):
        hexString = format(dec,'0X')
        #print(len(hexString))
        if 4*len(hexString) < byteNumber:
            hexString = int(byteNumber/4 - len(hexString))*'0' + hexString
        #print(hexString)
        return hexString
        
    def disconnect(self):
        # Close the instrument
        self.rotationMount.close()
        print('Instrument Closed')
        
    def getPos(self):
        self.pos = self.rotationMount.query('0gp')
        #print('raw output = ' + str(self.pos),' ,pos = ' + str(self.conv*twos_complement(self.pos[3:],32)) + ' deg')
        return self.conv*self.twos_complement(self.pos[3:],32)
        
    def moveAbsolute(self, newPos):
        # The user should enter newPos in degrees; it will be converted to the machine scale.
        machineScalePos = int(newPos/self.conv)
        commandString = "0ma" + self.hexStringFromDec(machineScalePos)
        #print('commandString',commandString)
        self.rotationMount.write(commandString)

class pseudoDelayStage():
    def __init__(self):
        print('This is a pseudo Delay Stage for testing')
        self.pos = 0
    def disconnect(self):
        print('This is a pseudo Delay Stage for testing')
    def getPos(self):
        return self.pos
    def setPos(self, newPos):
        if (newPos > 0) and (newPos < 225):
            self.pos = newPos
        else:
            print('Invalid Input Position; DS motion denied')

class pseudoCamera():
    def __init__(self):
        self.triggerTimer = time.time()
        self.triggered = False
        self.pseudoExposure = 1
        self.pseudoGain = 1
        self.pseudoWait = 1/30
        self.print_device_info()
        self.counter = 0
    def print_device_info(self):
        print('This is a fake camera for testing the program')
    def disconnect(self):
        print('This is a fake camera for testing the program')
    def trigger(self):
        # Record the time of this trigger so that we can make the pseudo camera behave realistically
        self.triggerTimer = time.time()
        if self.triggered == False:
            self.triggered = True
        else:
            print('Attempting to trigger the camera while it is already triggered')
    def grabImage(self):
        while (self.triggered == False) or ((time.time()-self.triggerTimer) < self.pseudoExposure):
            time.sleep(self.pseudoWait)
        self.triggered = False
        self.counter += 1
        # ~ return  np.array((self.counter%11)*1000*np.random.random([500,500]), dtype = 'uint16')
        # ~ return  np.array(9000*np.random.random([500,500]), dtype = 'uint16')
        # ~ return np.array(5000*(np.zeros([500,500]) + 1), dtype = 'uint16')
        return np.array(np.random.normal(loc = 50000,scale = 1000,size = (500,500)), dtype = 'uint16')
    def getExposure(self):
        return self.pseudoExposure
    def getGain(self):
        return self.pseudoGain
    def setExposure(self,setValueSecond):
        self.pseudoExposure = setValueSecond
    def setGain(self,setValue):
        self.pseudoGain = setValue

class pseudoHalfWavePlate():
    def __init__(self):
        print('This is a pseudo half wave plate for testing')
    def moveAbsolute(self, val):
        print(f'Setting HWP Pos to {val} (pseudo hwp)')
    def getPos(self):
        return 0
    def disconnect(self):
        print('disconnect from pseudo half wave plate')

if TESTING == True:
    # In test mode we redefine these classes as test objects that don't connect to anything
    Camera = pseudoCamera
    DelayStage = pseudoDelayStage
    HalfWavePlate = pseudoHalfWavePlate
    
