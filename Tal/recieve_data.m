%% constants
constants

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

%% INFO:
T_bit = 1 / hz / bytes_length * 1.0085;
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
for i = 0:bytes_length*2
    % xline(i/hz/length(binary) * 1.0085); % 1.0085 is a magic number for phase gained over time (approx.)
end
%bar(bits_pred);
%title("Symbol Predictions [At SymbolRate=200symbols/sec]");
%xlabel("Symbol Index");
%ylabel("Symbol Value");

%%
prediction_correct_length = bits_pred(1: floor(length(bits_pred) / 8)*8);
%bit_errors = sum(nonzeros(binary).' ~= prediction_correct_length);
%ber = bit_errors / length(bits_pred);

%% print data
pred_no_zeros = nonzeros(prediction_correct_length) - 1;
val = get_ascii_from_bits(pred_no_zeros.').';
disp(char(val));

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