%% Lab 2
%Importing the given signals with cyclic prefix
pss;

%Importing the r_pss
load('rpss.mat');

%% Exercise 2 - OFDM demodulation
%a)

symb_len = 2192;
prefix_len = 144;
rss_len = symb_len-prefix_len;

%To easier change what received signal to use through the code
rxs = rxs0;

%Use NF from lab 1.
NF = 4437;

%Find the r_pss
r_pss = rxs(NF+prefix_len+(0:rss_len-1));

%Find r_sss without prefix and do fft.
r_sss = rxs(NF+2*symb_len+prefix_len+(0:rss_len-1));
R_sss = fftshift(fft(r_sss,2048));
start_point = find(R_sss>10000);
R_sss_relevant = R_sss(start_point(1):start_point(length(start_point)));
figure
stem(abs(R_sss_relevant))
title('Absolute value of SSS')

figure
stem(abs(R_sss))
title('Complex modulus of entire OFDM symbol')

%% b)
R_sss_RE = R_sss(961:961+126);
figure
plot(real(R_sss_RE), imag(R_sss_RE),'bo')
title('Complex modulus of resource elements')

figure
stem(abs(R_sss_RE))

%Can use slide 16 from the OFDM slides to try and explain this.

%% c)
R_pss = fftshift(fft(r_pss,2048));
R_pss_rel = R_pss(961:1087);
pss_2_shift = fftshift(pss_2);
pss_2_rel = pss_2_shift(961:1087);
figure
stem(abs(R_pss_rel))
figure
stem(abs(pss_2_rel))


H = fftshift(conj(pss_2_rel).*R_pss_rel);
h = ifftshift(ifft(H));
H_space = linspace(-length(h)/2, length(h)/2, length(h));
figure
stem(H_space, abs(h))
title('h_PSS')
% See the peak in -20

%% d) 

%% e)

%% f)
SSS = conj(H).*fftshift(R_sss_RE);
figure
scatter(real(SSS),imag(SSS))

%% Exercise 3
res = conj(d(3:3:size(d,1), : )) * SSS;
[seq_val, seq_id] = max(abs(res));
figure(7);
stem(real(res));