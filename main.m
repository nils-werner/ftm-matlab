%% Floete

%
% Parameter
%

play = 1;
samples = 10000;
blocks = 2;


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
% Spektrum extrahieren und normalisieren
%

disp('Extracting Spectrum')
sspek = flspe;
sfreq = flfrq;

% Subsampling
subsampled = 1:length(sspek)/samples:length(sspek);
sspek = interp1(1:length(sspek),sspek,subsampled);
sfreq = interp1(1:length(sfreq),sfreq,subsampled);

tspek = sspek;
tfreq = sfreq;

sspek = zeros(blocks,length(tspek)/blocks);
sfreq = zeros(blocks,length(tspek)/blocks);


for i=1:blocks
    %size(tspek((1:length(tspek)/blocks).*i))
    %size(sspek(1,:))
    sspek(i,:) = tspek((1:length(tspek)/blocks).*i);
    sfreq(i,:) = tfreq((1:length(tspek)/blocks).*i);
end

% Frequenzrauschen gegen Schwebungen
sfreq = sfreq + rand(size(sfreq));

% Normalisierung
sspek = sspek.*max(max(abs(sspek)));




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
tic
parfor j=1:blocks
    for i=1:length(sspek(j,:))
        if sspek(j,i) > 0
            sig = sig+sin(2*pi*sfreq(j,i)*x).*sspek(j,i);
            count = count+1;
        end
        if mod(i,floor(length(sspek)*blocks/30)) == 0
            % Fortschrittsanzeige
            fprintf('.')
        end
    end
end
toc

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

%% Wiedergabe Synthese

sound(sig,sf)

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