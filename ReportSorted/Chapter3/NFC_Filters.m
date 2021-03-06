% Plot the shape of the function Near Field Compensation filter
% Auteur : Dupont Samuel
% Version : 1.0 June 2016

clear all;close all;clc


M = 5;% order of the different bessel func
ct.r_hp = 1.07;% radius of the simulated struc

%% Initiate default variables
f=linspace(0,20000,1000) ;
ct.k = f.*2*pi/340 ; 
ct.kr_hp = f.*2*pi/340*ct.r_hp;

%% Filters of the speaker structure

for ii=0:M
    for jj=1:length(ct.kr_hp)
        a(jj,ii+1)=1./(ct.k(jj)/(4*pi)*Hankel_sph(ii,2,ct.kr_hp(jj)));
    end
end
figure(1)
% subplot(211)
semilogx(f,20*log10(abs(a)))
grid on
xlim([f(1) f(end)])
xlabel('Freq [Hz]');ylabel('F_{m}(kr)[dB] ref 1')
legend('m=0','m=1','m=2','m=3','m=4','m=5');
% subplot(212)
% plot(f,20*unwrap(angle((a))))
% xlabel('Freq [Hz]');ylabel('Phase [rad]')
% xlim([f(1) f(end)])

%%Filters for source encoding 

ct.r_s = 2.14;% radius of the simulated source
ct.kr_s = f.*2*pi/340*ct.r_s;
for ii=0:M
    for jj=1:length(ct.kr_hp)
        a(jj,ii+1)=Hankel_sph(ii,2,ct.kr_s(jj))./Hankel_sph(ii,2,ct.kr_hp(jj));
    end
end
figure(2)
semilogx(f,20*log10(abs(a)))
grid on
xlim([f(1) f(end)])
xlabel('Freq [Hz]');ylabel('(Hm_{m}(kr) [dB]')
legend('M=0','M=1','M=2','M=3','M=4','M=5');
