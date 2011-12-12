%% Floete
clear
[flsig,fs,nbits] = wavread('wav/floetesoft.wav');
[flspe,p] = spektrum(flsig,fs);
plot(flspe/1000, 10*log10(p), 'k') 
xlabel('Frequency (kHz)') 
ylabel('Power (dB)')

%% Kurzzeit-Spektrum der Floete

frame_overlap = 10; % ms
frame_length  = 20;
window        = 'hamming';

nfft = round(frame_length  * fs / 1000); % convert ms to points
noverlap = round(frame_overlap * fs / 1000); % convert ms to points
window   = eval(sprintf('%s(nfft)', window)); % e.g., hamming(nfft)

spectrogram(flsig(:,1), window, noverlap, nfft, fs, 'yaxis');

%% Wiedergabe

sound(flsig,fs);

%% Synthese-Spektrum

sspek = zeros(1,16000);
sspek(942) = 0.5;
sspek(945) = 0.8;
sspek(947) = 1.1;
sspek(945) = 1;

sspek(974) = 0.4;
sspek(978) = 3;
sspek(980) = 1;

sspek(1891) = 0.2;

sspek = sspek./max(abs(sspek));
%% Sinus

sf = 44100;
x = [0:1/sf:4.0];

%sspek=zeros(0,16000);
%sspek(440) = 1;
%sspek(880) = 1;

f=0;

sig = sin(2*pi*0*2*x).*0;

for i=sspek
	f = f+1;
	if i > 0
		sig = sig+sin(2*pi*f*2*x).*i;
	end
end

%sig = sin(2*pi*440*x);
%sig = sig+sin(2*pi*440*2*x)./2;
%sig = sig+sin(2*pi*440*4*x)./4;
%sig = sig+sin(2*pi*440*8*x)./8;
%sig = sig+sin(2*pi*440*16*x)./16;
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