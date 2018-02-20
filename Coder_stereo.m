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

    % input audio file 44.1 Khz
[xorig,Fs] = audioread('testMicrofone.m4a'); % read the audio file
xMono = sum(xorig, 2) / size(xorig, 2); % convert stereo to mono
 
 
% Downsample from 44.1 Khz to 8 kHz
[P,Q] = rat(8e3/Fs);
abs(P/Q*Fs-8000);

Fs = Fs*P/Q; % New Sampling Frequency Fs

x = resample(xMono,P,Q); % downsample original audio

% P44_1 = audioplayer(xMono,44100);
% P8 = audioplayer(x,8000);
% %play(P44_1)
% %play(P8)


% figure
% plot(xorig);
% title('Original audio 44 kHz');

figure
plot(x);
title('Downsampled audio 8 kHz');



% ==== Create a Buffer that break audio to 20 msec blocks
n = 160; % number of sample size in 1 block
yorig = buffer(x(:,1), n); % create <size(x)/160> frames, each frame has 160 samples
% yorig2 = buffer(x(:,2), n);
% yorig = [yorig1 yorig2];
% ==== Find the spectrum of audio using lpc
x_idct = [];
ncoeff = 11;
rec_sig = [];

      % ====================================%
      % For transmitter
for i = 1:  size(yorig,2) % buffer
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

        % Find the DCT
    x_dct{i} = sqrt(n).*dct(acq{i},n);
    
        % Find Variance for finding Number of Bits
    var1{i} = g{i}./abs(h1{i}).^2; % variance using freqz of a{i}
    
       
        % Find number of bits
        % 1.75: Avergae number bit/sample
        % 280: Maximum total of bits in 1 block
    [Abit{i},Tbit(i)] = bitAllocate (1.75,n,var1{i},280);
    
        % zero insertion to DCT vector when zero bits is allocated
%        v_zeros = find(Abit{i}==0);
%        if isempty(v_zeros) == 0
%             si_first_zero = v_zeros(1);
%             si_last_zero = v_zeros(end);
%             x_dct{i}(si_first_zero:si_last_zero) = 0; 
%        end 
       
       
        % Uniform Quantizer
        % Quantize the DCT coefficients
    %quantized{i} = round(x_dct{i}.*2.^Abit{i})./2.^Abit{i}; 
    
        % Image quantizing technique
%         level = flipud(2.^Abit{i});
%     im_quantz_dct{i} = imquantize(x_dct{i},level);






    
    
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
		normalized{i} = (x_dct{i}./si_difDCT) + abs(si_minDCT);
        quantized{i} = round((x_dct{i}-si_minDCT).*(2.^Abit{i}-1)./si_difDCT); %Normalize it from 0 to 2^n-1 + round
        %quantized{i} = round((normalized{i}-si_minDCT)).*si_stepSize + si_minDCT; 
        bit_encode{i} = de2bi(quantized{i}); %Convert it to a binary bitstream

        % For Receiver
		% Convert received bitstream back to proper values (convert to decimal and de-normalize)
	decoded{i} = bi2de(bit_encode{i}).*si_stepSize + si_minDCT;
    %decoded{i} = bi2de(bit_encode{i}).*si_stepSize - abs(si_minDCT);
    
%     if exist('si_first_zero')
%         decoded{i}(si_first_zero-20:si_last_zero) = 0; 
%     end

        % Take the inverse DCT of the dquantized & reconstruct
    rec_sig = [rec_sig; idct(decoded{i})./sqrt(n)];
    inv_DCT{i} = idct(decoded{i})./sqrt(n);

end

    % Residual = Xorig at 8kHz - Reconstructed Signal
%resid = x - rec_sig(1:length(x));
      


% ============== Plot waveform
% t = (0:length(x) - 1)/Fs;
% figure
% %subplot(2,1,1);
% plot(t,x);
% title('Downsampled Waveform');
% legend('Downsampled Waveform');
% xlabel('Time (s)');
% ylabel('Amplitude');



% ========== Block ========================
block = 75;

figure % for plotting estimated spectral
plot(f{block},20*log10(abs(h{block})+eps),'r',f{block},20*log10(abs(h1{block})+eps),'--');
title(strcat('DCT & Formant of Block:   ',num2str(block)));
hold on
plot(f{block},20*log10(abs(x_dct{block})),'b'); % for plotting DCT
legend('LPC Filter','Gain','DCT');
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
hold off



    % Find var using computed gain above inside loop
figure % for testing var using freqz
plot(-log2((var1{block})+eps),'r');
title('Variance log_2(\sigma^2) using frqz of a{i}');
hold on
stairs(Abit{block},'b');
grid on
legend('variance','Number of Bit');


mess = [];
for kk = 1:length(frqs{block})
        mess=[mess sprintf('Formant %d Frequency %f\n',kk,frqs{block}(kk))];
end
mess = [mess sprintf('Gain is: %f\n',sqrt(g{block}))];
box = msgbox(mess,strcat('DCT & Formant of Block (Hz):',num2str(block)));
set(box, 'position', [650 300 300 105]); %makes box bigger



% ========= Find the number of bits of each block 
figure
stairs(Abit{block});
grid on
title('Bit Allocation');

% 
figure
plot(f{block},20*log10(abs(decoded{block})),'b');
hold on
plot(f{block},20*log10(abs(x_dct{block})),'r');
grid on
legend('Decoded Quantized DCT','DCT');
title(strcat('Decoded Quantized DCT vs DCT in Freq domain of Block: ',num2str(block)));

% 
figure
plot(decoded{block},'b');
hold on
plot(x_dct{block},'r');
grid on
legend('Decoded Quantized DCT','DCT');
title('Decoded Quantized DCT vs DCT in Time domain');

% 
figure
plot(x,'b');
hold on
plot(rec_sig,'r');
grid on
legend('Initial Signal','Reconstructed DCT');
title('Initial Signal vs Reconstructed Time domain');

% figure
% plot(resid)
% title('Residual = Original Signal at 8kHz - Reconstructed Signal');
% grid on

figure
subplot(2,1,1)
plot(x,'b');
subplot(2,1,2)
plot(rec_sig,'r');
%plot(x-rec_sig(1:length(x)),'k')
legend('Initial Signal','Reconstructed DCT');


figure
plot(x_dct{block},'r');
hold on; plot(quantized{block},'b');
title('DCT vs normalized DCT (real quantized)');
legend('DCT','Quantized DCT');

figure
plot(acq{block},'r');
hold on; plot(inv_DCT{block},'b');
title('origin signal vs inv DCT (real quantized)');
legend('OG','INV DCT');


figure
plot(idct(decoded{block}),'b');
hold on
plot(idct(x_dct{block}),'r');
grid on
legend('INV of Decoded Quantized DCT','INV of DCT');
title('INV Decoded Quantized DCT vs INV DCT in Time domain');



