Nfft=2048;
%Load real signals
pss0;
pss1;
pss2;

pss0_t = ifft(pss_0); 
pss0_t = pss0_t/norm(pss0_t);

pss1_t = ifft(pss_1);
pss1_t = pss1_t / norm(pss1_t);

pss2_t = ifft(pss_2);
pss2_t = pss2_t / norm(pss2_t);

%Adding the cyclic prefix
pss0_t = cat(1,pss0_t(((2048-143):2048)), pss0_t).';
pss1_t = cat(1,pss1_t(((2048-143):2048)), pss1_t).';
pss2_t = cat(1,pss2_t(((2048-143):2048)), pss2_t).';

% Reading the rx-signals
load('rxsignal_justnoise.mat');
load('rxsignal_withchannel.mat');
load('rxsignal_withchannelandfreqoff');
