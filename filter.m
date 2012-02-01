%%

[b,a] = butter(4,0.3,'s'); 
[bz,az] = impinvar(b,a,10);
sys = tf(b,a);
impulse(sys);
hold on;
impz(10*bz,az,[],10);


%%

%            b*z
% G = -----------------
%     z^2 + c_1*z + c_0

l = 0.65;
Ts = 60.97;
rho = 1140;
A = 0.5188*10^-6;
m = 1;

sigma = 0;
omega = m*(pi/l)*sqrt(Ts/(rho*A));

num = [0 1 0];
den = [1 1 1];

H = tf(num, den, 0.1)