constants

%% generate data
string = 'hello johhny'; % (string size + 2) needs to divide 32768 (be a power of 2)
string = pad(string, bytes_length - 2);
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

%% generate signal for <SECONDS> seconds
SECONDS = 0.5;
java.lang.Thread.sleep(SECONDS*1000)

%% Turn off signal generator
% Sets the output to 0 V DC.
%[status.setSigGenOff] = invoke(sigGenGroupObj, 'setSigGenOff');

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