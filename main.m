%% Sinus

clear
sf = 44100;
x = [0:1/sf:4.0];
sig = sin(2*pi*440*x);
%plot(x,sig);

%% Wiedergabe

sound(sig,sf)

%% Leistungsspektrum

freq = spektrum(sig);
plot(freq/1000, 10*log10(p), 'k') 
xlabel('Frequency (kHz)') 
ylabel('Power (dB)') 

%% Testmusik

clear
load handel
sound(y,Fs)