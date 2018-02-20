function bitstream = bit_assign(data)

bitstream=[];
data = cell2mat(data);

for i = 1:length(data)
    bitstream = [bitstream str2double(data(i))];
end
end