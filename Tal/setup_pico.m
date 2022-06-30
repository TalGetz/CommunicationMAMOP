%% Clear command window and close any figures

clc;
close all;

%% Load configuration information

PS2000aConfig;

%% Device connection

% Check if an Instrument session using the device object |ps2000aDeviceObj|
% is still open, and if so, disconnect if the User chooses 'Yes' when prompted.
if (exist('ps2000aDeviceObj', 'var') && ps2000aDeviceObj.isvalid && strcmp(ps2000aDeviceObj.status, 'open'))
    
    openDevice = questionDialog(['Device object ps2000aDeviceObj has an open connection. ' ...
        'Do you wish to close the connection and continue?'], ...
        'Device Object Connection Open');
    
    if (openDevice == PicoConstants.TRUE)
        
        % Close connection to device.
        disconnect(ps2000aDeviceObj);
        delete(ps2000aDeviceObj);
        
    else

        % Exit script if User selects 'No'.
        return;
        
    end
    
end

% Create a device object. 
% The serial number can be specified as a second input parameter.
ps2000aDeviceObj = icdevice('picotech_ps2000a_generic.mdd');

% Connect device object to hardware.
connect(ps2000aDeviceObj);

%% Obtain Signalgenerator group object
% Signal Generator properties and functions are located in the Instrument
% Driver's Signalgenerator group.

sigGenGroupObj = get(ps2000aDeviceObj, 'Signalgenerator');
sigGenGroupObj = sigGenGroupObj(1);

%% Turn off signal generator
% Sets the output to 0 V DC.

[status.setSigGenOff] = invoke(sigGenGroupObj, 'setSigGenOff');

%% Arbitrary waveform generator - set parameters
plaster_freq = 1;
hz = 20;
set(ps2000aDeviceObj.Signalgenerator(1), 'startFrequency', 1.0 / plaster_freq * hz);
set(ps2000aDeviceObj.Signalgenerator(1), 'stopFrequency', 1.0 / plaster_freq * hz);
set(ps2000aDeviceObj.Signalgenerator(1), 'offsetVoltage', 0.0);
set(ps2000aDeviceObj.Signalgenerator(1), 'peakToPeakVoltage', 4000.0);

%% physical layer constants
awgBufferSize = get(sigGenGroupObj, 'awgBufferSize');
M = 3;
Max_V = 1;
dV = Max_V / (M-1);
Sample_Size = awgBufferSize;
f = 174300; % HZ
f = 110000; % HZ
f_pico = hz; % HZ
occurences_to_collect = 2;

