clc;clear variables; close all
%% regularisation au sens Tikhonov des filtres de champs proche (encodage)

%% Initialisation
ct.Fs=48000;
ct.M=5;
ct.r_micsph=0.07;
ct.c_air=340;
ct.hankel_order=2;
ct.N_mic=50;
ct.r_target=0.2;
ct.N=2^16;

var.f=FreqVect(ct.Fs,ct.N);

var.k(:,1)=2*pi.*var.f/ct.c_air;
var.m_vect=0:ct.M;
var.m_sum_vect=(var.m_vect+1).^2;
var.nbr_m=(2.*var.m_vect)+1;

N.k=length(var.k);

%% Filtres theoriques
var.Hprim=zeros(N.k,(ct.M+1)^2);% H_prim init

for ii=0:ct.M
    var.Hprim(:,(ii)^2+1:(ii+1)^2) = repmat(1i^(ii-1)./((var.k.*ct.r_micsph).^2.* ...
    Hankel_sph_1_deriv(ii,ct.hankel_order,var.k.*ct.r_micsph)),1,var.nbr_m(ii+1)) ;
end
var.EqFilt=1./ var.Hprim;



%% Regularisation
ct.n=6;
ct.fc=(1:ct.M).*ct.c_air/(2*pi*ct.r_target);
ct.wn=ct.fc*2/ct.Fs;
var.h=ones(ct.N,var.m_sum_vect(ct.M+1));
for ii=1:ct.M
[b,a] = butter(ct.n,ct.wn(ii)/2,'high');
var.h(:,(ii)^2+1)=freqz(b,a,ct.N,ct.Fs*2);
var.h(:,(ii)^2+1:(ii+1)^2)= repmat(var.h(:,(ii)^2+1),1,var.nbr_m(ii+1)) ;
% semilogx(var.f,db(abs(var.h(:,(ii-1)^2+1))))
% hold on
end
var.EqFilt_reg=bsxfun(@times,var.EqFilt,var.h);


%% Affichage
set(gca,'ColorOrderIndex',1)
semilogx(var.f,db(var.EqFilt_reg));
legend('M=0','M=1','M=2','M=3','M=4','M=5');
%back to default
set(gca,'ColorOrderIndex',1);hold on
semilogx(var.f,db(var.EqFilt),'--');

xlabel('Frequency [Hz]');ylabel('Amplitude [dB]');hold off;
ct.pos=closest(80,var.f);
axis([var.f(ct.pos), var.f(end/2) -100 200])
