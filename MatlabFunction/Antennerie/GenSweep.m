function [ x, t, ct, N ] = GenSweep( f1,f2,T,ct )
% Generate a sweep signal 
% F1 = bigenning frequency
% F2 = end frequency
% T = number of secong of the sweep


%% Parameter of the sweep signal
w1 = f1*2*pi ;%begin   frequency
w2 = f2*2*pi ;%end frequency

K = T*w1/log(w2/w1) ;% coef 1  
L = T/log(w2/w1) ;% coef 2

%% Parameter of the signal 
t.T_sweep = 0:1/ct.Fs:T ;
x = sin(K.*(exp(t.T_sweep./L)-1))*0.90 ;
N.N_sweep = length(x) ; ct.dfe_sweep=ct.Fs/N.N_sweep;
t.F_sweep = 0:ct.dfe_sweep:(N.N_sweep-1)*ct.dfe_sweep ;% Frequency axis signal



end

