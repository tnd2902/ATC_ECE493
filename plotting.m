function plotting(block,h,h1,x_dct,var1,Abit,quant,rec_sig,x,f,acq)


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
plot(-log2((abs(var1{block}))+eps),'r');
title('Variance log_2(\sigma^2) using frqz of a{i}');
hold on
stairs(Abit{block},'b');
grid on
legend('variance','Number of Bit');


% mess = [];
% for kk = 1:length(frqs{block})
%         mess=[mess sprintf('Formant %d Frequency %f\n',kk,frqs{block}(kk))];
% end
% mess = [mess sprintf('Gain is: %f\n',sqrt(g{block}))];
% box = msgbox(mess,strcat('DCT & Formant of Block (Hz):',num2str(block)));
% set(box, 'position', [650 300 300 105]); %makes box bigger



% ========= Find the number of bits of each block 
figure
stairs(Abit{block});
grid on
title('Bit Allocation');


figure
plot(f{block},20*log10(abs(quant{block})),'b');
hold on
plot(f{block},20*log10(abs(x_dct{block})),'r');
grid on
legend('Decoded Quantized DCT','DCT');
title(strcat('Decoded Quantized DCT vs DCT in Freq domain of Block: ',num2str(block)));


figure
plot(quant{block},'b');
hold on
plot(x_dct{block},'r');
grid on
legend('Decoded Quantized DCT','DCT');
title('Decoded Quantized DCT vs DCT in Time domain');

% % 
% figure
% plot(x,'b');
% hold on
% plot(rec_sig,'r');
% grid on
% legend('Initial Signal','Reconstructed DCT');
% title('Initial Signal vs Reconstructed Time domain');


% figure
% subplot(2,1,1)
% plot(x,'b');
% subplot(2,1,2)
% plot(rec_sig,'r');
% %plot(x-rec_sig(1:length(x)),'k')
% legend('Initial Signal','Reconstructed DCT');


figure
plot(x_dct{block},'r');
hold on; plot(quant{block},'b');
title('DCT vs normalized DCT (real quantized)');
legend('DCT','Quantized DCT');

% figure
% plot(acq{block},'r');
% hold on; plot(inv_DCT{block},'b');
% title('origin signal vs inv DCT (real quantized)');
% legend('OG','INV DCT');


figure
plot(idct(quant{block}),'b');
hold on
plot(idct(x_dct{block}),'r');
grid on
legend('INV of Decoded Quantized DCT','INV of DCT');
title('INV Decoded Quantized DCT vs INV DCT in Time domain');

end