function freq_offset = find_freq_offset(m, r, pssi, Fs, n, dfmin)
    exponent = pssi.*exp(-2*1i*pi*n*m*dfmin/Fs);
    
    freq_offset = abs(conj(exponent)*r).^2;
end