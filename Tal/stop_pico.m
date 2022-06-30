%% Turn off signal generator
% Sets the output to 0 V DC.
[status.setSigGenOff] = invoke(sigGenGroupObj, 'setSigGenOff');

%% Disconnect device
% Disconnect device object from hardware.
disconnect(ps2000aDeviceObj);
delete(ps2000aDeviceObj);