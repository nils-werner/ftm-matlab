%% Floete

%
% Parameter
%

play = 1;
samples = 8000;


%
% Datei laden
%

disp('Loading File')
[flsig,sf,nbits] = wavread('wav/floetesoft.wav');
[flfrq,flspe] = spektrum(flsig,sf);
plot(flfrq/1000, 10*log10(flspe), 'k') 
xlabel('Frequency (kHz)') 
ylabel('Power (dB)')



%
% Kurzzeit-Spektrum der Floete
%

disp('Analyzing Spectrum')
frame_overlap = 10; % ms
frame_length  = 20;
window        = 'hamming';

nfft = round(frame_length  * sf / 50); % convert ms to points
noverlap = round(frame_overlap * sf / 50); % convert ms to points
window   = eval(sprintf('%s(nfft)', window)); % e.g., hamming(nfft)

[S,F,T,P] = spectrogram(flsig(:,1), window, noverlap, nfft, sf);
surf(T,F,10*log10(P),'edgecolor','none'); axis tight;

disp('Extracting Spectrum')



%
% Spektrum extrahieren und normalisieren
%

% 75te Spalte = Spektrum zum Zeitpunkt 0:00:00:75
sspek = flspe;
ffreq = flfrq;

subsampled = 1:length(sspek)/samples:length(sspek);

sspek = interp1(1:length(sspek),sspek,subsampled);
ffreq = interp1(1:length(ffreq),ffreq,subsampled);

ffreq = ffreq + rand(size(ffreq));

% Normalisierung
disp('Normalizing')
sspek = sspek.*max(abs(sspek));




%
% Synthese
%

fprintf('Synthesizing')

% X-Achse
x = [0:1/sf:2.0];

f=0;
count=0;

% Leeren Signalvektor erzeugen.
sig = sin(x).*0;

% Über Koeffizienten des Spektrums iterieren
for i=1:length(sspek)
	if sspek(i) > 0
		sig = sig+sin(2*pi*ffreq(i)*x).*sspek(i);
		count = count+1;
	end
	if mod(i,floor(length(sspek)/30)) == 0
        % Fortschrittsanzeige
		fprintf('.')
	end
end

disp([num2str(count) ' Samples'])




%
% Spektren plotten
%

disp('Plotting')
sig = sig ./ max(abs(sig));
plot(x(1:1000),sig(1:1000));

subplot(2,1,1)
[freq,p] = spektrum(sig,sf);
plot(freq/1000, 10*log10(p), 'k') 
subplot(2,1,2)
[freq,p] = spektrum(flsig,sf);
plot(freq/1000, 10*log10(p), 'k') 
xlabel('Frequency (kHz)') 
ylabel('Power (dB)')

disp('Done')

%
% Wiedergabe Synthese
%

if play == 1
    sound(sig,sf)
end

%% Synthese festhalten

oldsig = sig;

%% Wiedergabe alte Synthese

sound(oldsig,sf)

%% Wiedergabe Original

sound(flsig,sf);

%% Speichern Synthese

wavwrite(sig,sf,'wav/synthese.wav')

%% Testmusik

clear
load handel
sound(y,Fs)