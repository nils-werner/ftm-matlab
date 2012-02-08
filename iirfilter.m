%%

%            b*z
% G = -----------------
%     z^2 + c_1*z + c_0

clear;
axis auto;

play = 1;
write = 0;
plotting = 1;
customfilter = 1;

% Figures initialisieren/wiederfinden

freqs = findfigure('freqs');
sebene = findfigure('sebene');
signals = findfigure('signals');
result = findfigure('result');



% Koeffizienten und Daten

filters = 30;

l = 0.65;
Ts = 60.97;
rho = 1140;
A = 0.5188*10^-6;
m = 1:filters;
xa = 0.1;

T = 44100;
blocksize = 100;
seconds = 1;
samples = seconds*T;
inputdata = [1 zeros(1,samples-1)];
x = 1:samples;
y = zeros(1,samples);

H = [];
sigmas = [];
omegas = [];
nums = [];
dens = [];
ws = [];
hs = [];
sigs = [];
as = [];

cc=hsv(filters);

tic
for i = m;
	a = sin(i*pi*xa/l);
	sigma = -0.2*i^2;
	omega = i*(pi/l)*sqrt(Ts/(rho*A));
	
	b = T*sin(omega*1/T)/(omega*1/T);
	c1 = -2*exp(sigma*1/T)*cos(omega*1/T);
	c0 = exp(2*sigma*1/T);
	
	num = [0 b 0];
	den = [1 c1 c0];
	
	fA = [0 -c0; 1 -c1];
	fC = [0 1];
	CA = [];
	
	[h,w] = freqz(num,den,[], T);
	
	state = [0 1]';
	sig = zeros(1,samples);
	
	if customfilter == 1;
		j = 1;
		while j <= blocksize
			CA = [CA; fC * fA^j];
			j = j + 1;
		end
		
		Ap = fA^blocksize;
		
		j = 1;
		while j <= samples
			sig(j:j+blocksize-1) = CA * state;
			state = Ap * state;
			j = j + blocksize;
		end

		y = y + a*sig;
	else
		sig = filter(num,den,inputdata);
	
		y = y + a*sig;
	end

	if plotting == 1
		sigmas = [sigmas sigma];
		omegas = [omegas omega];

		nums = [nums num'];
		dens = [dens den'];
		
		as = [as a];
		
		sigs = [sigs sig'];

		hs = [hs h];
		ws = [ws w];
	end
end
toc

y = y./max(abs(y));

hold off

if plotting == 1
	figure(freqs);
	for i = m;
		plot(ws(:,i)',20*log10(abs(as(i))*abs(hs(:,i)')),'color',cc(i,:));
		hold on
	end
	axis([0 11000 70 180]);

	hold off

	figure(signals);
	for i = m;
		plot(x,sigs(:,i),'color',cc(i,:));
		hold on
	end

	hold off

	figure(sebene);
	for i = m;
		plot(sigmas(i),omegas(i), 'x-','color',cc(i,:));
		hold on
	end;

	hold off

	hold on
	figure(result);
	plot(x,y);
	axis([1000 1500 -1.2 1.2]);

	hold off
end

if play == 1
	sound(y,T);
end

if write == 1
	wavwrite(y,T,'wav/filter.synth.wav');
end

%%

wavwrite(y,T,'wav/filter.synth.wav');

%%

sound(y,T);