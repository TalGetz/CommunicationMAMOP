%% constants
constants
delay = 1;


%% flags
load("flags.mat");

%% buffer
buffer = ' ';

%% wait state
while(~out_flag && ~in_flag)
   load("flags.mat");
   java.lang.Thread.sleep(delay*1000);
   [x, ret_val] = evalc('recieve_data_function(blockGroupObj, ps2000aDeviceObj, ps2000aEnuminfo)');
   index_byte = ret_val(1:8);
   index = bin2dec(num2str(index_byte.'));
   ret_val = ret_val(9:end);
   string = char(get_ascii_from_bits(ret_val).');
   if (index == 0)
       in_flag = 1;
       buffer = strcat(buffer, string);
   end
end

% TODO: need to deal with outflag
% TODO: need to save in_flag to flags.mat so that out_protocol doesnt send

last_index = index;
%% read state
while(index ~= 255)
   java.lang.Thread.sleep(delay*1000);
   [x, ret_val] = evalc('recieve_data_function(blockGroupObj, ps2000aDeviceObj, ps2000aEnuminfo)');
   index_byte = ret_val(1:8);
   index = bin2dec(num2str(index_byte.'));
   if (last_index == index)
       continue;
   end
   last_index = index;
   ret_val = ret_val(9:end);
   string = char(get_ascii_from_bits(ret_val).');
   buffer = strcat(buffer, string);
end

%% 
disp(buffer);

%% generate bits from ascii
function ascii = get_ascii_from_bits(bits)
bits_mat = reshape(bits, 8, []).';
ascii = zeros(length(bits)/8, 1);
for i = 1 : length(ascii)
    st = num2str(bits_mat(i, :));
    st(isspace(st)) = '';
    ascii(i) = bin2dec(st);
end
end