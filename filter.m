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

l = 0.65;
Ts = 60.97;
rho = 1140;
A = 0.5188*10^-6;
m = 1:30;

T = 44100;

H = [];
sigmas = [];
omegas = [];

if(findobj('type','figure','name','freqs'))
    freqs = figure(findobj('type','figure','name','freqs'));
else
    freqs = figure('name','freqs');
end

for i = m;
	sigma = -0.00000002*i^2;
	omega = i*(pi/l)*sqrt(Ts/(rho*A));
	
	sigmas = [sigmas sigma];
	omegas = [omegas omega];
	
	b = T*sin(omega*T)/(omega*T);
	c1 = -2*exp(sigma*T)*cos(omega*T);
	c0 = exp(2*sigma*T);
	
	num = [0 b 0];
	den = [1 c1 c0];
	
	%H = [H tf(num, den, 0.1)]

	hold on
	
	freqz(num,den,[], T);

end

hold off

if(findobj('type','figure','name','sebene'))
    sebene = figure(findobj('type','figure','name','sebene'));
else
    sebene = figure('name','sebene');
end

plot(sigmas,omegas, 'x-');

%sound(impulse(H,2),T);