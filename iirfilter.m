%%

%            b*z
% G = -----------------
%     z^2 + c_1*z + c_0

clear;

do_play = 1;
do_write = 0;

% Koeffizienten und Daten

% Anzahl Filter
filters = 30;
m = 1:filters;

% Saiten-Koeffizienten
l = 0.65;
Ts = 60.97;
rho = 1140;
A = 0.5188*10^-6;
E = 5.4*10^9;
I = 0.171*10^-12;
d1 = 8*10^-5;
d3 = -1.4*10^-5;

% Abtastpunkt
xa = 0.1;

% Abtastrate und Samplelänge
T = 44100;
seconds = 10;
samples = seconds*T;

% Blockverarbeitungs-Länge
blocksize = 100;

% Ausgangssignal
y = zeros(1,samples);

% Ausgangs- und Übergangsmatrix, Zustandsvektor
block_C = [];
block_A = [];
block_CA = [];
block_state = [];

tic
for i = m;
    % Pol aufstellen
	gamma = i*(pi/l);
	sigma = (1/(2*rho*A)) * (d3*gamma^2 - d1);
	omega = sqrt( ( (E*I)/(rho*A) - (d3^2)/((2*rho*A)^2) )* gamma^4 + (Ts/(rho*A)) * gamma^2 + (d1/(2*rho*A))^2);
    
    % Ausgangsgewichtung
    a = sin(i*pi*xa/l);
    
    % Übertrangungsfunktions-Koeffizienten
	b = T*sin(omega*1/T)/(omega*1/T);
	c1 = -2*exp(sigma*1/T)*cos(omega*1/T);
	c0 = exp(2*sigma*1/T);
		
    % Zustandsraum-Matrizen
	fA = [0 -c0; 1 -c1];
	fC = [0 a];
	state = [1 0]';
	
    block_C = [block_C fC];
    block_A = blkdiag(block_A, fA);
    block_state = [block_state; state];
end

j = 1;
while j <= blocksize
    block_CA = [block_CA; block_C * block_A^j];
    j = j + 1;
end

block_A = block_A^blocksize;

j = 1;
while j <= samples
    y(j:j+blocksize-1) = block_CA * block_state;
    block_state = block_A * block_state;
    j = j + blocksize;
end
toc

y = y./max(abs(y));

if do_play == 1
	sound(y,T);
end

if do_write == 1
	wavwrite(y,T,'wav/filter.synth.wav');
end

%%

wavwrite(y,T,'wav/filter.synth.wav');

%%

sound(y,T);