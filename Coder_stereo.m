% Tien Dang
% Senior Project: ATC
% This code will test finding formant of STEREO speech
% Plot spectral estimation
% Find: Predicton error e(using levinson)
%       prediction error variances g (using lpc)
% Compare values of both


clear all
close all

% Close all dialog box
if exist('box', 'var')
	delete(box);
	clear('box');
end

load('num.mat');
load('txt.mat');
load('raw.mat');
%[num,txt,raw] = xlsread('MaxTable2.xlsx');

    % input audio file 44.1 Khz
[xorig,Fs] = audioread('ENG_M.wav'); % read the audio file
xMono = sum(xorig, 2) / size(xorig, 2); % convert stereo to mono
 
 
% Downsample from 44.1 Khz to 8 kHz
[P,Q] = rat(8e3/Fs);
abs(P/Q*Fs-8000);

Fs = Fs*P/Q; % New Sampling Frequency Fs

x = resample(xMono,P,Q); % downsample original audio

figure
plot(x);
title('Downsampled audio 8 kHz');



% ==== Create a Buffer that break audio to 20 msec blocks
n = 160; % number of sample size in 1 block
yorig = buffer(x(:,1), n,0); % create <size(x)/160> frames, each frame has 160 samples
% ==== Find the spectrum of audio using lpc
x_idct = [];
ncoeff = 11;
rec_sig = [];

      % ====================================%
      % For transmitter
for i = 1:size(yorig,2) % buffer
    acq{i} = yorig(:,i);      % Assume that this is what our data
                            % acquisition board returns
                         
                            
        % ============== Get linear prediction of whole signal
    [a{i}, g{i}] = lpc(acq{i},ncoeff); % Finding LPC formant
                                     % a: LPC coefficients   
    
        % Find the Gain
    r = xcorr(acq{i},ncoeff,'biased'); % Running Auto-corraltion
    r(1:ncoeff) = [];
    g{i} = sum(a{i}.*r'); % LPC Gain squared of current block
    
    g_db{i} = 20*log10(g{i}); % gain in db
    
        % Perform Frequency Response
                % for formant
    rts{i} = roots(a{i});
    rts{i} = rts{i}(imag(rts{i})>0.01);
    angz{i} = atan2(imag(rts{i}),real(rts{i}));

                % bandwidths of the formants
    [frqs{i},indices{i}] = sort(angz{i}.*(Fs/(2*pi)));
    bw{i} = -1/2*(Fs/(2*pi))*log(abs(rts{i}(indices{i})));    
    
    [h{i},f{i}] = freqz(1,a{i},length(acq{i}),Fs); % for plotting formant
    [h1{i},fv{i}] = freqz(sqrt(g{i}),a{i},length(acq{i}),Fs); % for plotting variance of estimated spectral

    
            % Find Variance for finding Number of Bits
    var1{i} = g{i}./abs(h1{i}).^2; % variance using freqz of a{i}
    
        % Find the DCT
    x_dct{i} = log10(abs(g{i})).*dct(acq{i},n);
    %x_dct{i} = -log2(var1{i}).*dct(acq{i},n);
    %x_dct{i} = dct(acq{i},n);
    maxDCT(i) = max(x_dct{i});
    

    
       
        % Find number of bits
        % 1.75: Avergae number bit/sample
        % 280: Maximum total of bits in 1 block
        temp_mse = prod(var1{i}.^(1/(n^2)));
    	[Abit{i},Tbit(i)] = bitAllocate (1.75,n,var1{i},280,temp_mse);
    

       %[quant{i},bitstream{i},totalBit(i)] = encoder4(Abit{i},x_dct{i},num,txt,raw,var1{i},g{i});
        [quant{i},bitstream{i},totalBit(i)] = encoder2(Abit{i},x_dct{i},num,txt,raw,var1{i});
        mse(i) = immse(x_dct{i},quant{i});
       %[Abit{i},Tbit(i)] = bitAllocate (1.75,n,var1{i},280,mse(i));
       %[quant{i},bitstream{i},totalBit(i)] = encoder2(Abit{i},quant{i},num,txt,raw,var1{i});
    
		%Side info (We still haven't quantized this)
		
		si_maxDCT = max(x_dct{i}); %Get max value
        si_minDCT = min(x_dct{i}); % get min
        si_difDCT = si_maxDCT-si_minDCT;
	si_stepSize = si_difDCT./(2.^Abit{i}); %Get step size (don't have to send in real system, can calculate from max/min)
        si_max_sample = 2.^(max(Abit{i}))./2;
		
		%In a real situation, we will also need to communicate the number of bits per sample (or more efficiently, at which samples it changes)
		%That will also be stored in the side information, since we won't have a 2d vector separating each binary vector neatly like we do in MATLAB
		%--Ben

				%ADC of Ben, tien cmt out feb 8/2018
 		%quantized{i} = round((x_dct{i}-si_minDCT).*(2.^Abit{i}-1)./si_difDCT); %Normalize it from 0 to 2^n-1 + round
%		normalized{i} = (x_dct{i}./si_difDCT) + abs(si_minDCT);
        %quantized{i} = round((x_dct{i}-si_minDCT).*(2.^Abit{i}-1)./si_difDCT); %Normalize it from 0 to 2^n-1 + round
        %quantized{i} = round((normalized{i}-si_minDCT)).*si_stepSize + si_minDCT; 
        %bit_encode{i} = de2bi(quantized{i}); %Convert it to a binary bitstream

        % For Receiver
		% Convert received bitstream back to proper values (convert to decimal and de-normalize)
	%decoded{i} = bi2de(bit_encode{i}).*si_stepSize + si_minDCT;


        % Take the inverse DCT of the dquantized & reconstruct
%     rec_sig = [rec_sig; idct(decoded{i})./sqrt(n)];
%     inv_DCT{i} = idct(decoded{i})./sqrt(n);
        
%rec_sig = [rec_sig; idct(quant{i}./log10(abs(g{i})))];
%rec_sig = [rec_sig; idct(quant{i},n)./(-log2(abs(var1{i})))];
%rec_sig = [rec_sig; idct(quant{i}./(-log2(var1{i})))];
%rec_sig = [rec_sig; idct(quant{i})];
end

b_stream = cell2mat(bitstream);
g = cell2mat(g);
figure; plot(x,'b');hold on; plot(rec_sig,'r');

[rec_sig, r_dct] = decoder(cell2mat(Abit), b_stream,num,txt,g);
%[rec_sig, r_dct] = decoder3(cell2mat(Abit), b_stream,num,txt,g,var1);

err = (x - rec_sig(1:length(x),1));
errDCT = cell2mat(x_dct) - r_dct;
snr_last = snr(x,err);
plotting(215,h,h1,x_dct,var1,Abit,quant,rec_sig,x,f,acq);

bitRate = sum(Tbit)/(n*i)*8000/1000;





