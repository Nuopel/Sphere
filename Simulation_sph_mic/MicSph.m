clear variables; close all;clc

% Simulate the behavior of a pherical microphone
% recording a monopole source positioned on the
% ambisonics set up

%% Define the constant
ct.r_hp_sca = 1.07 ;%rayon de la sphere
ct.r_micsph = 0.1 ;
ct.hankel_order = 1 ;
ct.M_th = 5 ;
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
k = transpose(2*pi.*t.F_sweep/ct.c_air) ;
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
Ymn.source = sph_harmonic( ct.M_th,1,source.theta,source.phi ) ;
Bmn.source = bsxfun(@times,Fm,Ymn.source') ;


%% Fabrication Matrice "C" spherical harmonic

[ Sphmic.theta, Sphmic.phi, Sphmic.r ] = cart2sph( Sphmic.x, Sphmic.y, Sphmic.z ) ;
% conversion  to pherical coord possible problem with ref of cart2sph
Ymn.Mic = sph_harmonic( ct.M_th, ct.N_mic, Sphmic.theta, Sphmic.phi );
% % --> calculate Bmn.*Ymn
% % --> calculate i^(m-1)/((ct.r_micsph.k).^2.*Hm'(ka))
% % --> change  hfm pour Hankel spherique
for ii=1:1 %ct.N_mic
var.Bmn_Ymn(:,:)=bsxfun(@times,Bmn.source',Ymn.Mic(:,ii));

end
var.pressure=zeros(N.N_sweep,ct.M_th);%cgange
for ii=1:ct.M_th
    var.sum_Bmn_Ymn = sum(var.Bmn_Ymn(1:var.m_sum_vect(ii),:),1);
    var.sum_Hprim = 1i^(var.m_vect(ii))./((k.*ct.r_micsph).^2.*HFm(ii,2,k,ct.r_micsph));
    var.pressure(:,ii) = var.sum_Bmn_Ymn'.*var.sum_Hprim;
end
% Pressure.mic_rigid = sum(i^) 

