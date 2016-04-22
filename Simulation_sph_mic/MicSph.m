clear variables; close all;clc

% Simulate the behavior of a pherical microphone
% recording a monopole source positioned on the
% ambisonics set up

%% Define the constant
ct.r_hp_sca = 1.07 ;%rayon de la sphere
ct.r_micsph = 0.1 ;
ct.hankel_order = 2;
ct.M_th = 15;
ct.nbr_M_th=(ct.M_th+1).^1;
ct.M=5;
ct.Fs=48000;
ct.c_air=340;

var.m_vect=0:ct.M_th;
var.m_sum_vect=(var.m_vect+1).^2;
var.nbr_m=(2.*var.m_vect)+1;
%% Define ambisonics set up
[ ArraySpeaker, ct.N_speaker ] = CreateSpeakerSystem(ct.r_hp_sca);% create the sphere set up, sort in struc Array


%% Define Spherical microphone set up (same as speaker but smaller)
[ Sphmic, ct.N_mic ] = CreateSpeakerSystem(ct.r_hp_sca);% create the sphere set up, sort in struc Array

%% Choix de la source (Ae^j(wt-kx))
[ source.sweep, t, ct, N ] = GenSweep(20, 20000, 4, ct ) ;
k = 2*pi*500/340;
N.N_sweep=1;
%% POSITION DE LA SOURCE : select a speaker
speaker = 12 ;
source.x = ArraySpeaker.x(speaker) ;
source.y = ArraySpeaker.y(speaker) ;
source.z = ArraySpeaker.z(speaker) ;
[ source.theta, source.phi, source.r ] = cart2sph(source.x, source.y, source.z) ;
clear speaker ;
%% Propagate  spherical microphone
% Decomposition en harmonique spherique

for ii = 0:ct.M_th
    Fm(:,(ii)^2+1:(ii+1)^2) = repmat(-HFm( ii, ct.hankel_order, k, ct.r_hp_sca ),1,var.nbr_m(ii+1)) ;
end
% Fm(1,:)=0;% does it right ?
Ymn.source = sph_harmonic( ct.M_th,1,source.theta,source.phi ) ;
Bmn.source = bsxfun(@times,Fm,Ymn.source') ;


%% Fabrication Matrice "C" spherical harmonic

[ Sphmic.theta, Sphmic.phi, Sphmic.r ] = cart2sph( Sphmic.x, Sphmic.y, Sphmic.z ) ;
% conversion  to pherical coord possible problem with ref of cart2sph
Ymn.Mic = sph_harmonic( ct.M_th, ct.N_mic, Sphmic.theta, Sphmic.phi ) ;
% % --> calculate Bmn.*Ymn
% % --> change  hfm pour Hankel spherique
Pressure=zeros(N.N_sweep,ct.N_mic);
for jj=1:ct.N_mic
    var.Bmn_Ymn = bsxfun(@times,Bmn.source',Ymn.Mic(:,jj)) ;
    
    % var.pressure = zeros(N.N_sweep,ct.M_th) ;------» tic toc method choice
    for ii=0:ct.M_th
        var.sum_Bmn_Ymn = sum(var.Bmn_Ymn(1:var.m_sum_vect(ii+1),:),1) ; 
        var.Hprim = 1i^(var.m_vect(ii+1)-1)./((k.*ct.r_micsph).^2.*Hankel_sph_1_deriv(var.m_vect(ii+1),ct.hankel_order,k.*ct.r_micsph));
        var.Hprim2 = 1i^(var.m_vect(ii+1))*(Bessel_sph(var.m_vect(ii+1),k.*ct.r_micsph)-Bessel_sph_1_deriv(var.m_vect(ii+1),k.*ct.r_micsph)/Hankel_sph_1_deriv(var.m_vect(ii+1),ct.hankel_order,k.*ct.r_micsph)*Hankel_sph(var.m_vect(ii+1),ct.hankel_order,k.*ct.r_micsph))
        if ii==0;
            var.pressure = var.sum_Bmn_Ymn'.*var.Hprim ;
        else
            var.pressure=sum([var.pressure, var.sum_Bmn_Ymn'.*var.Hprim],2);
        end
    end
    Pressure(:,jj)=var.pressure;
    
end
% var.pressure(1,:)=0;% does it rightj ?

%% Encodage 
ct.M=5;
var.Hprim=zeros(N.N_sweep,ct.M);
for ii=0:ct.M
        var.Hprim(:,(ii)^2+1:(ii+1)^2) = repmat(1i^(var.m_vect(ii+1)-1)./((k.*ct.r_micsph).^2.*Hankel_sph_1_deriv(var.m_vect(ii+1),ct.hankel_order,k.*ct.r_micsph)),1,var.nbr_m(ii+1)) ;
end

 Bmn.recons = zeros(var.m_sum_vect(ct.M+1),N.N_sweep) ;
 Ymn.Micrecons = Ymn.Mic(1:var.m_sum_vect(ct.M+1),1:ct.N_mic) ;
for ii=1:N.N_sweep
   Bmn.recons(:,ii) = diag(1./var.Hprim(ii,:))*Ymn.Micrecons*diag(Sphmic.w)*Pressure(ii,:)' ;
   disp(ii) ;
end
% ct=rmfield(ct,'Fs2_sca');
