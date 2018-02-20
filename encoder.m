function [quant,bitstream,totalBit] = encoder(Abit,x_dct,num,txt,raw)
Abit = Abit;
x_dct = x_dct;
max_bit = max(Abit);

si_maxDCT = max(x_dct); %Get max value
si_minDCT = min(x_dct); % get min
si_difDCT = si_maxDCT-si_minDCT;
si_stepSize = si_difDCT./(2.^Abit); %Get step size (don't have to send in real system, can calculate from max/min)

%[num,txt,raw] = xlsread('MaxTable.xlsx');
bitstream = [];
for j = 1:length(x_dct)
    level = 2.^Abit(j);
    if level == 1 % 0 bits
        quant(j,1) = 0;
        %bits = [0];
    elseif level == 2 % 1 bits
%         if abs(x_dct(j)) <0.025 
%             quant(j,1) = 0;
%             bits = [0];
%         else
        if x_dct(j) >= 0
            quant(j,1) = num(2,2);
            bits = bit_assign(raw(2,3));
            %bitstream(j,1) = raw(2,3);
        else
            quant(j,1) = num(1,2);
            bits = bit_assign(raw(1,3));
            %bitstream(j,1) = raw(1,3);
        end
%         end
    elseif level == 4 % 2 bits
%         if abs(x_dct(j)) <0.06 
%             quant(j,1) = 0;
%             bits = [0 0];
%         else
        if x_dct(j) <=num(1,4)
            quant(j,1) = num(1,5);
            bits = bit_assign(raw(1,6));
            %bits(j,1) = raw(1,6);
        end
        
        for k = 2:level
            if x_dct(j) > num(k-1,4) && x_dct(j)<=num(k,4)
                if x_dct(j) >=0 
                    quant(j,1) = num(k-1,5);
                    bits = bit_assign(raw(k-1,6));
                    %bitstream(j,1) = raw(k-1,6);
                else
                    quant(j,1) = num(k,5);
                    bits = bit_assign(raw(k,6));
                    %bitstream(j,1) = raw(k,6);
                end
            end
        end
        if x_dct(j) > num(level,4)
            quant(j,1) = num(level,5);
            bits = bit_assign(raw(level,6));
            %bitstream(j,1) = raw(level,6);
        end
%         end
        
    elseif level == 8 % 3 bits
%         if abs(x_dct(j)) <0.06
%             quant(j,1) = 0;
%             bits = [0 0 0];
%         else
        if x_dct(j) <=num(1,7)
            quant(j,1) = num(1,8);
            bits = bit_assign(raw(1,9));
            %bitstream(j,1) = raw(1,9);
        end
        for k = 2:level
            if x_dct(j) > num(k-1,7) && x_dct(j)<=num(k,7)
                if x_dct(j) >=0 
                    quant(j,1) = num(k-1,8);
                    bits = bit_assign(raw(k-1,9));
                    %bitstream(j,1) = raw(k-1,9);
                else
                    quant(j,1) = num(k,8);
                    bits = bit_assign(raw(k,9));
                    %bitstream(j,1) = raw(k,9);
                end
            end
        end
        if x_dct(j) > num(level,7)
            quant(j,1) = num(level,8);
            bits = bit_assign(raw(level,9));
            %bitstream(j,1) = raw(level,9);
        end
%         end
            
    elseif level == 16 % 4 bits
%         if abs(x_dct(j)) <0.06  
%             quant(j,1) = 0;
%             bits = [0 0 0 0];
%         else        
        if x_dct(j) <=num(1,10)
            quant(j,1) = num(1,11);
            bits = bit_assign(raw(1,12));
            %bitstream(j,1) = raw(1,12);
        end
        for k = 2:level
            if x_dct(j) > num(k-1,10) && x_dct(j)<=num(k,10)
                if x_dct(j) >=0 
                    quant(j,1) = num(k-1,11);
                    bits = bit_assign(raw(k-1,12));
                    %bitstream(j,1) = raw(k-1,12);
                else
                    quant(j,1) = num(k,11);
                    bits = bit_assign(raw(k,12));
                    %bitstream(j,1) = raw(k,12);
                end
            end
        end
        if x_dct(j) > num(level,10)
            quant(j,1) = num(level,11);
            bits = bit_assign(raw(level,12));
            %bitstream(j,1) = raw(level,12);
        end
%         end
            
    elseif level == 32 % 5 bits
%         if abs(x_dct(j)) <0.006  
%             quant(j,1) = 0;
%             bits = [0 0 0 0 0];
%         else
        if x_dct(j) <=num(1,13)
            quant(j,1) = num(1,14);
            bits = bit_assign(raw(1,15));
            %bitstream(j,1) = raw(1,15);
        end
        for k = 2:level
            if x_dct(j) > num(k-1,13) && x_dct(j)<=num(k,13)
                if x_dct(j) >=0 
                    quant(j,1) = num(k-1,14);
                    bits = bit_assign(raw(k-1,15));
                    %bitstream(j,1) = raw(k-1,15);
                else
                    quant(j,1) = num(k,14);
                    bits = bit_assign(raw(k,15));
                    %bitstream(j,1) = raw(k,15);
                end
            end
        end
        if x_dct(j) > num(level,13)
            quant(j,1) = num(level,14);
            bits = bit_assign(raw(level,15));
            %bitstream(j,1) = raw(level,15);
        end
%         end
            
    elseif level == 64 % 6 bits
        level = 36;
%         if abs(x_dct(j)) <0.004 
%             quant(j,1) = 0;
%             bits = [0 0 0 0 0 0];
%         else
            if x_dct(j) <=num(1,16)
                quant(j,1) = num(1,17);
                bits = bit_assign(raw(1,18));
                %bitstream(j,1) = num(1,18);
            end
            for k = 2:level
                if x_dct(j) > num(k-1,16) && x_dct(j)<=num(k,16)
                    if x_dct(j) >=0 
                        quant(j,1) = num(k-1,17);
                        bits = bit_assign(raw(k-1,18));
                        %bitstream(j,1) = raw(k-1,18);
                    else
                        quant(j,1) = num(k,17);
                        bits = bit_assign(raw(k,18));
                        %bitstream(j,1) = raw(k,18);
                    end
                end
            end
            if x_dct(j) > num(level,16)
                quant(j,1) = num(level,17);
                bits = bit_assign(raw(level,18));
                %bitstream(j,1) = raw(level,18);
            end
%         end
            
    else % 7, or 8 bits
        %si_stepSize = si_difDCT./(2.^Abit(j)); %Get step size (don't have to send in real system, can calculate from max/min)
        quant(j,1) = round((x_dct(j)-si_minDCT).*(2.^Abit(j)-1)./si_difDCT);
        bits = de2bi(quant(j,1));
        %disp('checkpoint');
    end
bitstream = [bitstream bits];

end
totalBit = length(bitstream);


% quant = mat2cell(quant);
% bitstream = cellstr(bitstream);



end