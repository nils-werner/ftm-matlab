%% Sinus

clear
sf = 44100;
x = [0:1/sf:4.0];
sig = sin(2*pi*440*x);
sig = sig+sin(2*pi*440*2*x);
sig = sig+sin(2*pi*440*4*x);
sig = sig+sin(2*pi*440*8*x);
sig = sig ./ max(abs(sig));
plot(x(1:1000),sig(1:1000));

%% Wiedergabe

sound(sig,sf)

%% Leistungsspektrum

[freq,p] = spektrum(sig,sf);
plot(freq/1000, 10*log10(p), 'k') 
xlabel('Frequency (kHz)') 
ylabel('Power (dB)') 

%% Testmusik

clear
load handel
sound(y,Fs)