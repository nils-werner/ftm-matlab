%% Floete

%
% Parameter
%

play = 1;
%file = 'floetesoft';
file = 'gitarre';
samples = 10000;


%
% Datei laden
%

disp('Loading File')
[flsig,sf,nbits] = wavread(sprintf('wav/%s.wav',file));
[flfrq,flspe] = spektrum(flsig,sf);
plot(flfrq/1000, 10*log10(flspe), 'k') 
xlabel('Frequency (kHz)') 
ylabel('Power (dB)')



%
% Spektrum extrahieren und normalisieren
%

disp('Extracting Spectrum')
sspek = flspe;
ffreq = flfrq;

% Subsampling
subsampled = 1:length(sspek)/samples:length(sspek);
sspek = interp1(1:length(sspek),sspek,subsampled);
ffreq = interp1(1:length(ffreq),ffreq,subsampled);

% Frequenzrauschen gegen Schwebungen
ffreq = ffreq + rand(size(ffreq));

% Normalisierung
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
sig = sig .* 0.999;
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

%% Wiedergabe Synthese

sound(sig,sf)

%% Wiedergabe alte Synthese

sound(oldsig,sf)

%% Wiedergabe Original

sound(flsig,sf);

%% Speichern Synthese

wavwrite(sig,sf,sprintf('wav/%s.synth.wav',file))

%% Testmusik

clear
load handel
sound(y,Fs)