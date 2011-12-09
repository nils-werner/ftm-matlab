%% Sinus Plotten

clear
sf = 44000;
x = [0.000001:1/sf:4.0];
sig = sin(2*pi*440*x);
%plot(x,sig);
%sound(sig,sf)

n = length(sig); 
p = fft(sig); % take the fourier transform 
nUniquePts = ceil((n+1)/2); 
p = p(1:nUniquePts); % select just the first half since the second half 
                      % is a mirror image of the first
p = abs(p); % take the absolute value, or the magnitude 

p = p/n; % scale by the number of points so that
           % the magnitude does not depend on the length 
           % of the signal or on its sampling frequency  
p = p.^2;  % square it to get the power 

% multiply by two (see technical document for details)
if rem(n, 2) % odd nfft excludes Nyquist point 
	p(2:end) = p(2:end)*2; 
else 
	p(2:end -1) = p(2:end -1)*2; 
end 
  
freqArray = (0:nUniquePts-1) * (sf / n); % create the frequency array 
plot(freqArray/1000, 10*log10(p), 'k') 
xlabel('Frequency (kHz)') 
ylabel('Power (dB)') 

%%

sound(sig,sf)

%%

clear
load handel
sound(y,Fs)