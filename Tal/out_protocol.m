constants
%% input

msg = input(">%> ", "s");

%% generate data
binary = get_AM_Data_of_string(msg, bytes_length, 255); % (string size + 2) needs to divide 32768 (be a power of 2)
y = get_AM_Data_of_bits(binary, Sample_Size, f, f_pico, dV);
% plot(y);
%% Arbitrary waveform generator - simple
[status.setSigGenArbitrarySimple] = invoke(sigGenGroupObj, 'setSigGenArbitrarySimple', y);

%% generate signal for <SECONDS> seconds
SECONDS = 0.5;
java.lang.Thread.sleep(SECONDS*1000)

%% Turn off signal generator
% Sets the output to 0 V DC.
% [status.setSigGenOff] = invoke(sigGenGroupObj, 'setSigGenOff');

%% generate bits from string with index and nulls before&after
function bits = get_AM_Data_of_string(string, bytes_length, message_index)
    string = pad(string, bytes_length - 3);
    binary = get_bits_from_ascii(string) + 1;
    msg_index = de2bi(message_index, 8, 'left-msb') + 1;
    bits = [0 0 0 0 0 0 0 0 msg_index binary 0 0 0 0 0 0 0 0];
end

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