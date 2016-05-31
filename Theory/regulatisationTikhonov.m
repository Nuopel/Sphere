clc;clear variables; close all
%% regularisation au sens Tikhonov des filtres de champs proche (encodage)

%% Initialisation
ct.Fs=48000;
ct.M=5;
ct.r_micsph=0.07;
ct.c_air=340;
ct.hankel_order=2;
ct.N_mic=50;

var.f=0:24000;
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
ct.ac=20;% maximum noise amplification 

ct.a=10^(ct.ac/20);% maximal amplification for filter
ct.a2=sqrt(ct.N_mic)*10^(ct.ac/20);% maximal amplification for filter

ct.lambda=(1-sqrt(1-1/ct.a^2))/(1+sqrt(1-1/ct.a^2));

%Calculation of the regularised filter with tikhonov formula
var.EqFilt_reg=conj(var.Hprim)./(abs(var.Hprim).^2+ct.lambda^2);


% Affichage
semilogx(var.f,db(var.EqFilt_reg));
%back to default
set(gca,'ColorOrderIndex',1);hold on
semilogx(var.f,db(var.EqFilt),'--');

xlabel('Frequency [Hz]');ylabel('Amplitude [dB]');hold off;
axis([var.f(1), var.f(end) -100 200])
legend('M=0','M=1','M=2','M=3','M=4','M=5');

