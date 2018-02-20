% Tien Dang
% Senior Project: ATC
% This code will find the number of bit for each sample in each block
% r: butgeted bitrate
% n = length of block
% v = variance of estimated spectral
% Abit = bit allocation vector
% Tbit = total bits in one block

function [Abit,Tbit] = bitAllocate (r,n,v,limit,mse)


D = v.^(1/(n^2));
%Qnoise = prod(D); % geometric mean
Qnoise = mse;
r = r;
for i = 1:n
    b = r+ 0.5*-log2(v(i)/Qnoise); % from paper
    if b<=0
        Abit(i,1) = 0;
    else
        Abit(i,1) = round(b);
        if round(b) >= 6
            Abit(i,1) = 6;
        end
    end
end
Tbit = sum(Abit);
if Tbit >limit
    r = r-0.20;
    [Abit,Tbit] = bitAllocate (r,n,v,limit,mse);
    
   Tbit = sum(Abit); 
end

end
