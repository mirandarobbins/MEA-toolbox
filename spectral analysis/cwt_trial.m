% temp2=squeeze(data.filtered_lfp(:,spectro.channeltoanalyse,9));
temp2=data.filtered_lfp(:,spectro.channeltoanalyse);
%%
N=size(temp2,1);
t = ( 0:(N-1) )* (1/params.Fs);  % times
f = ( (0:(N-1))/N) * (params.Fs);  % frequencies


freqs = [1:100];
factor = 5/(2*pi);
scale = factor * params.Fs./ freqs;
wave='morl';
freq2=scal2frq(scale,wave,params.Fs^-1);
coef=cwt(temp2, scale, wave);
%
figure; 
subplot(2,1,1)
plot(t,temp2)
subplot(2,1,2)
   imagesc(t, freq2, (coef)); 
   axis('xy'); % flip the vertical axis over
   xlabel('time');
   ylabel('frequency');
   title('scalogram');
   
   %% try MP
     %Build a basis of Gabors 
% sig = temp2;%load('gong.mat'); %included in matlab 
t=0:0.00005:1;t(1)=[];
sig=chirp(t,10,0.9, 100)';
   rg = (-500:500)'; 
  sigmas = exp(log(2):.3:log(200)); 
  gabors = exp(-.5*(rg.^2)*sigmas.^(-2)).*cos(rg*sigmas.^(-1)*2*pi*2); 
   %Express signal as a sparse sum of Gabors 
  [ws,r] = temporalMP(sig,gabors,false,10000); 
  
  subplot(3,1,1);
plot(gabors);
title('Basis');
subplot(3,1,2);
plot([sig,r]);
legend('signal','approximation');
subplot(3,1,3);
%the convolution here is to make the spikes visible
imagesc(conv2(ws,exp(-rg.^2/2/20^2),'same')');
ylabel('increasing time ->');
ylabel('increasing frequency ->');
title(sprintf('weights (%d non-zero weights)',nnz(ws(:))));
  
  