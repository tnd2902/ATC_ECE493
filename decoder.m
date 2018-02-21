function[rec_sig, r_dct] = decoder(Abit, bitstream,num,txt,g)

%[num,txt,raw] = xlsread('MaxTable2.xlsx');
Abit = Abit;
bitstream = bitstream;


% src2 = num2str(cell2mat(txt(1:2,1)));
% src4 = num2str(cell2mat(txt(1:4,4)));
% src8 = num2str(cell2mat(txt(1:8,7)));
% src16 = num2str(cell2mat(txt(1:16,10)));
% src32 = num2str(cell2mat(txt(1:32,13)));
% src36 = num2str(cell2mat(txt(1:36,16)));

src2 = (cell2mat(txt(1:2,1)));
src4 = (cell2mat(txt(1:4,4)));
src8 = (cell2mat(txt(1:8,7)));
src16 = (cell2mat(txt(1:16,10)));
src32 = (cell2mat(txt(1:32,13)));
src36 = (cell2mat(txt(1:36,16)));

startIdx = 1;
r_dct = [];
rec_sig = [];

for j = 1:size(Abit,2)
    %startIdx = 1;
    for i = 1:size(Abit,1)
    bitNum = Abit(i,j);
    endIdx = startIdx + bitNum - 1;
    partBit = num2str(bitstream(startIdx:endIdx));
    partBit = partBit(find(~isspace(partBit)));
    
    if bitNum == 0
        sample = 0;
    elseif bitNum == 1
        for k = 1:2
            if strcmp(partBit, src2(k,1:end)) == 1
                sample = num(k,2);
            end
        end
    elseif bitNum == 2
        for k = 1:4
            if strcmp(partBit, src4(k,1:end)) == 1
                sample = num(k,5);
            end
        end        
    elseif bitNum == 3
        for k = 1:8
            if strcmp(partBit, src8(k,1:end)) == 1
                sample = num(k,8);
                
            end
        end  
    elseif bitNum == 4
        for k = 1:16
            if strcmp(partBit, src16(k,1:end)) == 1
                sample = num(k,11);
            end
        end  
    elseif bitNum == 5
        for k = 1:32
            if strcmp(partBit, src32(k,1:end)) == 1
                sample = num(k,14);
                %disp('check');
            end
        end  
    elseif bitNum == 6
        for k = 1:36
            if strcmp(partBit, src36(k,1:end)) == 1
                sample = num(k,17);
            end
        end  
    else
        sample = bi2de(partBit{i}).*si_stepSize + si_minDCT;
    end
    
    %r_dct = [r_dct; sample];
    r_dct(i,j) = sample;
    %rec_sig = [rec_sig; idct(r_dct)./sqrt(size(Abit,1))];
    
    startIdx = endIdx +1;
    
    end
    %r_dct = [r_dct sample];
     %inv_dct = idct(r_dct./log10(abs(g(j))));
     rec_sig = [rec_sig; idct(r_dct(1:end,j)./log10(abs(g(j))))];
end


end
