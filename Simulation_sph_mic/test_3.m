clear variables; close all;clc

% Simulate the behavior of a pherical microphone
% recording a monopole source positioned on the
% ambisonics set up

%% Define the constant
ct.r_hp_sca = 1.07 ;%rayon de la sphere
ct.r_micsph = 0.1 ;
ct.hankel_order = 2;
ct.M_th =3;
ct.M=1;
ct.nbr_M_th=(ct.M_th+1).^1;

ct.Fs=48000;
ct.c_air=340;

var.m_vect=0:ct.M_th;
var.m_sum_vect=(var.m_vect+1).^2;
var.nbr_m=(2.*var.m_vect)+1;
%% Define ambisonics set up
[ ArraySpeaker, ct.N_speaker ] = CreateSpeakerSystem(ct.r_hp_sca);% create the sphere set up, sort in struc Array



%% Choix de la source (Ae^j(wt-kx))

k = 2*pi*500/340;
N.N_sweep=1;
%% POSITION DE LA SOURCE : select a speaker
speaker = 12 ;
source.x = ArraySpeaker.x(speaker) ;
source.y = ArraySpeaker.y(speaker) ;
source.z = ArraySpeaker.z(speaker) ;
[ source.theta, source.phi, source.r ] = cart2sph(source.x, source.y, source.z) ;
clear speaker ;

%% Define Spherical microphone set up (same as speaker but smaller)
[ Sphmic, ct.N_mic ] = CreateSpeakerSystem(ct.r_micsph);% create the sphere set up, sort in struc Array


%% Encoding source
% Decomposition en harmonique spherique

for ii = 0:ct.M_th
    Fm(:,(ii)^2+1:(ii+1)^2) = repmat(FM_sph( ii, 1, k, ct.r_hp_sca ),1,var.nbr_m(ii+1)) ;
end
% Fm(1,:)=0;% does it right ?
Ymn.source = sph_harmonic( ct.M_th,1,source.theta,source.phi ) ;
Bmn.source = bsxfun(@times,Fm,permute(Ymn.source,[2 1])) ;
% a=load('B2.mat');
%  Bmn.source=a.Expression1(1:var.m_sum_vect(ct.M_th+1))';
%% Fabrication Matrice "C" spherical harmonic
[ Sphmic.theta, Sphmic.phi, Sphmic.r ] = cart2sph( Sphmic.x, Sphmic.y, Sphmic.z ) ;
% conversion  to pherical coord possible problem with ref of cart2sph
Ymn.Mic = sph_harmonic( ct.M_th, ct.N_mic, Sphmic.theta, Sphmic.phi ) ;% % --> calculate Bmn.*Ymn
% % --> change  hfm pour Hankel spherique
Pressure.difract=zeros(N.N_sweep,ct.N_mic);
Pressure.direct=zeros(N.N_sweep,ct.N_mic);

%%M=0
for ii=1:ct.N_mic
p_dir(ii,1)=1i^(0-1)/((k.*ct.r_micsph).^2*Hankel_sph_1_deriv(0,ct.hankel_order,k.*ct.r_micsph))*Bmn.source(1)*Ymn.Mic(1,ii);
p_dir(ii,2)=1i^(1-1)/((k.*ct.r_micsph).^2*Hankel_sph_1_deriv(1,ct.hankel_order,k.*ct.r_micsph))*sum(bsxfun(@times,permute(Bmn.source(1,1:var.m_sum_vect(2)),[2 1 3]),Ymn.Mic(1:var.m_sum_vect(2),ii))); 
p_dir(ii,3)=1i^(2-1)/((k.*ct.r_micsph).^2*Hankel_sph_1_deriv(2,ct.hankel_order,k.*ct.r_micsph))*sum(bsxfun(@times,permute(Bmn.source(1,1:var.m_sum_vect(3)),[2 1 3]),Ymn.Mic(1:var.m_sum_vect(3),ii)));
p_dir(ii,4)=1i^(3-1)/((k.*ct.r_micsph).^2*Hankel_sph_1_deriv(3,ct.hankel_order,k.*ct.r_micsph))*sum(bsxfun(@times,permute(Bmn.source(1,1:var.m_sum_vect(4)),[2 1 3]),Ymn.Mic(1:var.m_sum_vect(4),ii)));
% p_dir(ii,5)=1i^(4-1)/((k.*ct.r_micsph).^2*Hankel_sph_1_deriv(4,ct.hankel_order,k.*ct.r_micsph))*sum(bsxfun(@times,permute(Bmn.source(1,1:var.m_sum_vect(5)),[2 1 3]),Ymn.Mic(1:var.m_sum_vect(5),ii)));
% p_dir(ii,6)=1i^5*Bessel_sph(5,k*Sphmic.r(ii))*sum(bsxfun(@times,permute(Bmn.source(1,1:var.m_sum_vect(6)),[2 1 3]),Ymn.Mic(1:var.m_sum_vect(6),ii)));
% p_dir(ii,7)=1i^6*Bessel_sph(6,k*Sphmic.r(ii))*sum(bsxfun(@times,permute(Bmn.source(1,1:var.m_sum_vect(7)),[2 1 3]),Ymn.Mic(1:var.m_sum_vect(7),ii)));
% p_dir(ii,8)=1i^7*Bessel_sph(7,k*Sphmic.r(ii))*sum(bsxfun(@times,permute(Bmn.source(1,1:var.m_sum_vect(8)),[2 1 3]),Ymn.Mic(1:var.m_sum_vect(8),ii)));
% p_dir(ii,9)=1i^7*Bessel_sph(8,k*Sphmic.r(ii))*sum(bsxfun(@times,permute(Bmn.source(1,1:var.m_sum_vect(9)),[2 1 3]),Ymn.Mic(1:var.m_sum_vect(9),ii)));
% p_dir(ii,10)=1i^7*Bessel_sph(9,k*Sphmic.r(ii))*sum(bsxfun(@times,permute(Bmn.source(1,1:var.m_sum_vect(10)),[2 1 3]),Ymn.Mic(1:var.m_sum_vect(10),ii)));
end
p_dir=sum(p_dir,2);

%% Encodage 
var.Hprim2=zeros(N.N_sweep,ct.M);
for ii=0:ct.M
        var.Hprim2(:,(ii)^2+1:(ii+1)^2) = repmat(1i^(ii-1)./((k.*ct.r_micsph).^2.*Hankel_sph_1_deriv(ii,ct.hankel_order,k.*ct.r_micsph)),1,var.nbr_m(ii+1)) ;
end

Sphmic.w=ArraySpeaker.w;
 Bmn.recons = zeros(var.m_sum_vect(ct.M+1),N.N_sweep) ;
 Ymn.Micrecons = Ymn.Mic(1:var.m_sum_vect(ct.M+1),1:ct.N_mic) ;
for ii=1:N.N_sweep
   Bmn.recons(:,ii) = diag(1./var.Hprim2(ii,:))*Ymn.Micrecons*diag(Sphmic.w)*p_dir(:,ii) ;
   disp(ii) ;
end

% Bmn.source=permute(Bmn.source,[2 1]);