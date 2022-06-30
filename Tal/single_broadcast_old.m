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

%% generate data
string = '123456781234567812345678123456'; % (string size + 1) needs to divide 32768
binary = get_bits_from_ascii(string) + 1;
binary = [0 0 0 0 0 0 0 0 binary 0 0 0 0 0 0 0 0];

y = get_AM_Data_of_bits(binary, Sample_Size, f, f_pico, dV);
% plot(y);
%% Arbitrary waveform generator - simple
% Output an arbitrary waveform with constant frequency (defined above).

% Arb. Waveform : y (defined above)
[status.setSigGenArbitrarySimple] = invoke(sigGenGroupObj, 'setSigGenArbitrarySimple', y);

% Increment      : 0 (Hz)
% Dwell Time     : 1 (s)
% Arb. Waveform  : y (defined above)
% Sweep Type     : 0 (ps3000aEnuminfo.enPS3000ASweepType.PS3000A_UP)
% Operation      : 0 (ps3000aEnuminfo.enPS3000AExtraOperations.PS3000A_ES_OFF)
% Shots          : 2 
% Sweeps         : 0
% Trigger Type   : 0 (ps3000aEnuminfo.enPS3000ASigGenTrigType.PS3000A_SIGGEN_RISING)
% Trigger Source : 4 (ps3000aEnuminfo.enPS3000ASigGenTrigSource.PS3000A_SIGGEN_SOFT_TRIG)
% Ext. Threshold : 0
%[status.setSigGenArbitrary] = invoke(sigGenGroupObj, 'setSigGenArbitrary', 0.000001, 1, y, 0, 0, 0, 2, 0, 0, 4, 0);

%[status.sigGenSoftwareControl] = invoke(sigGenGroupObj, 'ps2000aSigGenSoftwareControl', 1);

%%
% 
% <<../images/ps2000a_arbitrary_waveform.PNG>>
% 
% Block data acquisition properties and functions are located in the 
% Instrument Driver's Block group.

blockGroupObj = get(ps2000aDeviceObj, 'Block');
blockGroupObj = blockGroupObj(1);

% Set pre-trigger and post-trigger samples as required - the total of this
% should not exceed the value of |maxSamples| returned from the call to
% |ps2000aGetTimebase2()|. The default of 0 pre-trigger and 8192 post-trigger
% samples is used in this example.


set(ps2000aDeviceObj, 'numPreTriggerSamples', 0);
set(ps2000aDeviceObj, 'numPostTriggerSamples', 8192 * 120 / hz * occurences_to_collect); % 2 occurences of the sequence

%% max detected voltage range
%turn Channel A on with 10v range
invoke(ps2000aDeviceObj, 'ps2000aSetChannel', 0, 1, 1, 5, 0.0); % Channel A

%turn Channel A off
% invoke(ps2000aDeviceObj, 'ps2000aSetChannel', 0, 0, 1, 9, 0.0); % Channel A

% Coupling : 1 (PS2000A_DC)
% Range    : 9 (PS2000A_10V)
% Offset   : 0.00

%% capture data
% This example uses the |runBlock()| function in order to collect a block of
% data - if other code needs to be executed while waiting for the device to
% indicate that it is ready, use the |ps2000aRunBlock()| function and poll
% the |ps2000aIsReady()| function.

% Capture a block of data:
%
% segment index: 0 (The buffer memory is not segmented in this example)

[status.runBlock] = invoke(blockGroupObj, 'runBlock', 0);

% Retrieve data values:

startIndex              = 0;
segmentIndex            = 0;
downsamplingRatio       = 1;
downsamplingRatioMode   = ps2000aEnuminfo.enPS2000ARatioMode.PS2000A_RATIO_MODE_DECIMATE;

% Provide additional output arguments for other channels e.g. chC for
% channel C if using a 4-channel PicoScope.

[numSamples, overflow, chA, chB] = invoke(blockGroupObj, 'getBlockData', startIndex, segmentIndex, ...
                                            downsamplingRatio, downsamplingRatioMode);
                                        
                                        
%% time
timeIntervalNanoseconds = 1;
timeNs = double(timeIntervalNanoseconds) * downsamplingRatio * double(0:numSamples - 1);
timeMs = timeNs / 1e6;


%% Turn off signal generator
% Sets the output to 0 V DC.
[status.setSigGenOff] = invoke(sigGenGroupObj, 'setSigGenOff');

%% Disconnect device
% Disconnect device object from hardware.
disconnect(ps2000aDeviceObj);
delete(ps2000aDeviceObj);

%% INFO:
T_bit = 1 / hz / length(binary) * 1.0085;
% plot(chA);

%% Decode Bit Array from Amplitude Array
time = timeMs;
amplitude_array = chA / 1000;
% find first half bit value, start time there, split in T_BIT times,
% do AM peak to peak in order to find out value then set into bit array
V_max = max(amplitude_array);
H_at_transmittance = V_max;
amp_normalized = abs(amplitude_array) / H_at_transmittance * (M-1);
% plot(time, amp_normalized);
first_level = 1;
below_first_level_threshold = 1/M/2;
%% first zero
amp_normalized_mean = movmean(amp_normalized, 100)*sqrt(2);
bool_arr = (amp_normalized_mean < below_first_level_threshold);
sum_bool = movsum(bool_arr, 100);
%plot(time, sum_bool);
first_ind = find(sum_bool >= max(sum_bool) / M * (M-1), 1, 'first');
% xline(time(first_ind), 'b');
amp_normalized = amp_normalized(first_ind:end);
time = time(first_ind:end);
%% first non-zero
amp_normalized_mean = movmean(amp_normalized, 100)*sqrt(2);

bool_arr = (amp_normalized_mean > below_first_level_threshold);
sum_bool = movsum(bool_arr, 100);
%plot(time, amp_normalized);
second_ind = find(sum_bool >= max(sum_bool) / 2, 1, 'first');
% xline(time(second_ind), 'g');
amp_normalized = amp_normalized(second_ind:end);
time = time(second_ind:end);

%% second zero
amp_normalized_mean = movmean(amp_normalized, 100)*sqrt(2);
bool_arr = (amp_normalized_mean < below_first_level_threshold);
sum_bool = movsum(bool_arr, 100);
%plot(time, amp_normalized);
third_ind = find(sum_bool >= max(sum_bool) / M * (M-1), 1, 'first');
% xline(time(third_ind), 'r');
%% reset arrays
time = time(1:third_ind);
time = time - min(time);
amp_normalized = amp_normalized(1:third_ind);
amp_normalized = amp_normalized / max(amp_normalized) * (M-1) * sqrt(2);
amplitude_array = amplitude_array(1:third_ind);
%% plotting
% plot(time, amp_normalized);
% title("Recieved Data Normalized [At SymbolRate=200symbols/sec]");
% xlabel("time [s]");
% ylabel("Normalized Amplitude");
bits_pred = [];
for i = 0:round(max(time) / T_bit)
    time_mask = (time > i*T_bit) & (time < (i+1)*T_bit);
    bits_pred(end+1) = nearest(prctile(amp_normalized(logical(time_mask)), 50));
end
for i = 0:length(binary)*2
    % xline(i/hz/length(binary) * 1.0085); % 1.0085 is a magic number for phase gained over time (approx.)
end
%bar(bits_pred);
%title("Symbol Predictions [At SymbolRate=200symbols/sec]");
%xlabel("Symbol Index");
%ylabel("Symbol Value");

%%
prediction_correct_length = bits_pred(1: floor(length(bits_pred) / 8)*8);
bit_errors = sum(nonzeros(binary).' ~= prediction_correct_length);
ber = bit_errors / length(bits_pred);
%%
pred_no_zeros = nonzeros(prediction_correct_length) - 1;
val = get_ascii_from_bits(pred_no_zeros.').';
disp(char(val));
%% Channel A
%axisHandleChA = subplot(3,1,1); 
%hold on
% % plot(axisHandleChA, timeMs, chA, 'b');
% plot(time, amplitude_array, 'b');
% scatter(axisHandleChA, timeMs(save_indexes(1, :) * slice_size), zeros(1, length(save_indexes(1, :))), 'g*');%hold off

% ylim(axisHandleChA, [-10 10]); % Adjust vertical axis for signal.
% title(axisHandleChA, 'Channel A');
% xlabel(axisHandleChA, 'Time (s)');
% ylabel(axisHandleChA, 'Voltage (mV)');
% for i = 0:length(binary)*2
%    xline(i/hz/length(binary) * 1.0085); % 1.0085 is a magic number for phase gained over time (approx.)
%end
%grid(axisHandleChA);


%% Channel B
%axisHandleChB = subplot(3,1,2); 
%hold on
% plot(axisHandleChB, timeMs, chB, 'r');
% plot(axisHandleChB, timeMs, moving_data, 'r');
% scatter(axisHandleChB, timeMs(save_indexes(1, :) * slice_size), moving_data(save_indexes(1, :)), 'g*');
% ylim(axisHandleChB, [-1500 1500]); % Adjust vertical axis for signal.
%hold off

%title(axisHandleChB, 'Channel B');
%xlabel(axisHandleChB, 'Time (s)');
%ylabel(axisHandleChB, 'Voltage (mV)');
%grid(axisHandleChB);

%axisHandleBits = subplot(3,1,3);
%hold on
%plot(axisHandleBits, pred_data, 1, 1));
%hold off

%title(axisHandleBits, 'Decoding');
%xlabel(axisHandleBits, 'Bit slices');
%ylabel(axisHandleBits, 'Bit');
%grid(axisHandleBits);


%% generate AM data from bits
function AM_Data = get_AM_Data_of_bits(bits, Sample_Size, f, f_pico, dV)
% constants
N_Bits = length(bits);
Bit_Sample_Size = Sample_Size / N_Bits;
sin_cycles_per_bit = f / f_pico / N_Bits;
% data
time = 0:Sample_Size-1;
mask = repelem(bits, Bit_Sample_Size);
sin_wave = sin(time / Bit_Sample_Size * sin_cycles_per_bit * 2 * pi);
AM_Data = mask .* sin_wave * dV;
%plot(time / Sample_Size / f_pico, AM_Data * 2);
%title("transmitted data [at SymbolRate=2000Symbols/sec]");
%xlabel("time [s]");
%ylabel("Amplitude [V]");
end

%% generate bits from ascii
function bits = get_bits_from_ascii(ascii)
bits = reshape( dec2bin(ascii, 8).', [], 1);
bits(bits=='1') = 1;
bits(bits=='0') = 0;
bits = double(bits).';
end

function ascii = get_ascii_from_bits(bits)
bits_mat = reshape(bits, 8, []).';
ascii = zeros(length(bits)/8, 1);
for i = 1 : length(ascii)
    st = num2str(bits_mat(i, :));
    st(isspace(st)) = '';
    ascii(i) = bin2dec(st);
end
end