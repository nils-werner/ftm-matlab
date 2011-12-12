%% Sinus

clear
sf = 44100;
x = [0:1/sf:4.0];
sig = sin(2*pi*440*x);
sig = sig+sin(2*pi*440*2*x)./2;
sig = sig+sin(2*pi*440*4*x)./4;
sig = sig+sin(2*pi*440*8*x)./8;
sig = sig+sin(2*pi*440*16*x)./16;
%sig = sig+sin(2*pi*440*32*x)./32;
%sig = sig+sin(2*pi*440*64*x)./64;
%sig = sig+sin(2*pi*440*128*x)./128;
%sig = sig+sin(2*pi*440*256*x)./256;
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

%% Floete

[flsig,fs,nbits] = wavread('wav/floete.wav');
[flspe,p] = spektrum(flsig,fs);
plot(flspe/1000, 10*log10(p), 'k') 
xlabel('Frequency (kHz)') 
ylabel('Power (dB)')

%% Wiedergabe

sound(flsig,fs);