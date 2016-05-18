clear all;close all;clc
% Plot the shape of the function 1/(i^m J(Kr)), used in order to obtain the
% Bmn factor for encoding

m = 2;% order of the different bessel func
ct.r_micsph = 0.07;
f=linspace(0,20000,1000) ;
ct.kr = f*2*pi/340*ct.r_micsph ; % radius of the simulated microphone

coef=Bessel_sph_all(m,2,ct.kr);
subplot(211)
plot(f,(coef))
xlabel('Frequency [Hz]');ylabel('J_m(kr)');legend('m=0','m=1','m=2');xlim([f(1) f(end)])
grid on
subplot(212)
plot(f,db(1./real(coef)))
xlabel('Frequency [Hz]');ylabel('(J_m(kr)^{-1} [dB]')
grid on

%% comp rigid, open

for ii=0:m
    for jj=1:length(ct.kr)
        a(jj,ii+1)=((ct.kr(jj).*ct.r_micsph).^2*Hankel_sph_1_deriv(ii,2,ct.kr(jj).*ct.r_micsph));
    end
end
figure

co = get(gca,'ColorOrder') % Initial
% Change to new colors.
set(gca, 'ColorOrder', [ 0.2 0.2 1;1 0.2 0; 1 0.6 0.2 ], 'NextPlot', 'replacechildren');
co = get(gca,'ColorOrder') % Verify it changed
subplot(211)
semilogx(ct.kr,db(1./real(coef)))
grid on
xlabel('kr');ylabel('(J_m(kr)^{-1} [dB]')
ylim([0 150]);xlim([ct.kr(1) ct.kr(end)]);legend('m=0','m=1','m=2');
subplot(212)
semilogx(ct.kr,db(a ))
grid on
ylim([0 150]);xlim([ct.kr(1) ct.kr(end)])
xlabel('kr');ylabel('(Fh_{m}(kr) [dB]')

figure

co = get(gca,'ColorOrder') % Initial
% Change to new colors.
set(gca, 'ColorOrder', [ 0.2 0.3 1;1 0.4 0; 1 0.8 0.2 ], 'NextPlot', 'replacechildren');
co = get(gca,'ColorOrder') % Verify it changed
h1=semilogx(f,db(1./real(coef)),'--');
grid on;hold on
% xlabel('kr');ylabel(' [dB]')

h2=semilogx(f,db(a ));
grid on;legend([h1(1) h2(1)],'(J_m(kr)^{-1}', 'Fh_{m}(kr)');
ylim([0 150]);xlim([f(1) f(end)])
xlabel('Frequency [Hz]');ylabel(' [dB]')

%% Hankel

for ii=0:5
    for jj=1:length(ct.kr)
        a(jj,ii+1)=((ct.kr(jj).*ct.r_micsph).^2*Hankel_sph_1_deriv(ii,2,ct.kr(jj).*ct.r_micsph));
    end
end
figure(3)
semilogx(f,db(a ))
grid on
ylim([0 150]);xlim([f(1) f(end)])
xlabel('Freq [Hz]');ylabel('(Fh_{m}(kr) [dB]')
legend('0','1','2','3','4','5')