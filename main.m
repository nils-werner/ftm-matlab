%% Sinus Plotten

clear
x = [0:0.000125:1.0];
sig = sin(2*pi*466.16*x);
%plot(x,sig);
sound(sig)

%%

clear
load handel
sound(y,Fs)