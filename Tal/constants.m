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

%% bytes_length
bytes_length = 32;