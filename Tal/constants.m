%% physical layer constants
awgBufferSize = 32768;% get(sigGenGroupObj, 'awgBufferSize');
M = 3;
Max_V = 1;
dV = Max_V / (M-1);
Sample_Size = awgBufferSize;
f = 174300; % HZ
f = 110000; % HZ
hz = 15;
f_pico = hz; % HZ
occurences_to_collect = 2;

%% bytes
bytes_length = 32;