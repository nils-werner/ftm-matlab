%% Floete
disp('Loading File')
[flsig,sf,nbits] = wavread('wav/floetesoft.wav');
[flspe,p] = spektrum(flsig,sf);
plot(flspe/1000, 10*log10(p), 'k') 
xlabel('Frequency (kHz)') 
ylabel('Power (dB)')


disp('Analyzing Spectrum')

% Kurzzeit-Spektrum der Floete
frame_overlap = 10; % ms
frame_length  = 20;
window        = 'hamming';

nfft = round(frame_length  * sf / 50); % convert ms to points
noverlap = round(frame_overlap * sf / 50); % convert ms to points
window   = eval(sprintf('%s(nfft)', window)); % e.g., hamming(nfft)

[S,F,T,P] = spectrogram(flsig(:,1), window, noverlap, nfft, sf);
surf(T,F,10*log10(P),'edgecolor','none'); axis tight;

disp('Extracting Spectrum')

% 75te Spalte = Spektrum zum Zeitpunkt 0:00:00:75
sspek = P(:,3)';

% Sinus
%sf = 44100;
x = [0:1/sf:4.0];

%sspek=zeros(0,16000);
%sspek(440) = 1;
%sspek(880) = 1;

disp('Normalizing')

sspek = sspek.*max(abs(sspek));


fprintf('Synthesizing')

f=0;
count=0;
sig = sin(2*pi*0*2*x).*0;

for i=1:length(sspek)
	if sspek(i) > 0
		sig = sig+sin(2*pi*i*F(2)*x).*sspek(i);
		count = count+1;
	end
	if mod(i,200) == 0
		fprintf('.')
	end
end

disp([num2str(count) ' Samples'])

disp('Plotting')

sig = sig ./ max(abs(sig));
plot(x(1:1000),sig(1:1000));

% Leistungsspektrum
subplot(2,1,1)
[freq,p] = spektrum(sig,sf);
plot(freq/1000, 10*log10(p), 'k') 
subplot(2,1,2)
[freq,p] = spektrum(flsig,sf);
plot(freq/1000, 10*log10(p), 'k') 
xlabel('Frequency (kHz)') 
ylabel('Power (dB)')

disp('Done')

%% Wiedergabe Synthese

sound(sig,sf)

%% Synthese festhalten

oldsig = sig;

%% Wiedergabe alte Synthese

sound(oldsig,sf)

%% Wiedergabe Original

sound(flsig,sf);

%% Testmusik

clear
load handel
sound(y,Fs)