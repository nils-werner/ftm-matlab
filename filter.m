%%

[b,a] = butter(4,0.3,'s'); 
[bz,az] = impinvar(b,a,10);
sys = tf(b,a)
impulse(sys);
hold on;
impz(10*bz,az,[],10);


%%

%            b*z
% G = -----------------
%     z^2 + c_1*z + c_0

clear;


% Figures initialisieren/wiederfinden

if(findobj('type','figure','name','freqs'))
    freqs = figure(findobj('type','figure','name','freqs'));
else
    freqs = figure('name','freqs');
end

if(findobj('type','figure','name','sebene'))
    sebene = figure(findobj('type','figure','name','sebene'));
else
    sebene = figure('name','sebene');
end

if(findobj('type','figure','name','signals'))
    signals = figure(findobj('type','figure','name','signals'));
else
    signals = figure('name','signals');
end

if(findobj('type','figure','name','result'))
    result = figure(findobj('type','figure','name','result'));
else
    result = figure('name','result');
end



% Koeffizienten und Daten

filters = 30;

l = 0.65;
Ts = 60.97;
rho = 1140;
A = 0.5188*10^-6;
m = 1:filters;

T = 44100;
seconds = 2;
samples = seconds*T;
inputdata = [1 zeros(1,samples-1)];
x = 1:samples;
y = zeros(1,samples);

H = [];
sigmas = [];
omegas = [];

cc=hsv(filters);

for i = m;
	sigma = -0.2*i^2;
	omega = i*(pi/l)*sqrt(Ts/(rho*A));
	
	sigmas = [sigmas sigma];
	omegas = [omegas omega];
	
	b = T*sin(omega*1/T)/(omega*1/T);
	c1 = -2*exp(sigma*1/T)*cos(omega*1/T);
	c0 = exp(2*sigma*1/T);
	
	num = [0 b 0];
	den = [1 c1 c0];
	
	%H = [H tf(num, den, 0.1)]

	figure(freqs);
	[h,w] = freqz(num,den,[], T);
	plot(w,20*log10(abs(h)),'color',cc(i,:));
	hold on
	
	figure(signals);
	plot(x,filter(num,den,inputdata),'color',cc(i,:));
	hold on
	pause(0.1);
	
	y = y + filter(num,den,inputdata);
end

hold off

figure(sebene);
for i = m;
	plot(sigmas(i),omegas(i), 'x-','color',cc(i,:));
	hold on
end;

hold off

figure(result);
plot(x,y);

y = y*min(y)

sound(y,T);

%%

wavwrite(y,T,'wav/filter.synth.wav');

%%

sound(y,T);