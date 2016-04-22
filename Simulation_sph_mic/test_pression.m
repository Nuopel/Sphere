clear variables; close all;clc

% Simulate the behavior of a pherical microphone
% recording a monopole source positioned on the
% ambisonics set up

%% Define the constant
ct.r_hp_sca = 1.07 ;%rayon de la sphere
ct.r_micsph = 0.2 ;
ct.hankel_order = 2;
ct.M_th = 4;
ct.nbr_M_th=(ct.M_th+1).^1;
ct.M=3;
ct.Fs=48000;
ct.c_air=340;

var.m_vect=0:ct.M_th;
var.m_sum_vect=(var.m_vect+1).^2;
var.nbr_m=(2.*var.m_vect)+1;
%% Define ambisonics set up
[ ArraySpeaker, ct.N_speaker ] = CreateSpeakerSystem(ct.r_hp_sca);% create the sphere set up, sort in struc Array



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

%% creation antenne
ct.pas_m = 2.5*2.54e-2; % pas de l'antenne
N.nbrx_sca =50; % nombre de micros par ligne
N.nbry_sca = 50; % nombre de micros par ligne
ct.N_mic=N.nbrx_sca*N.nbry_sca;
[ Antenna ] = AntennArray( ct.pas_m,N.nbrx_sca,N.nbry_sca);
% contruction antenne 2.5 po , 8 mic


%% Encoding source
% Decomposition en harmonique spherique

for ii = 0:ct.M_th
    Fm(:,(ii)^2+1:(ii+1)^2) = repmat(-HFm( ii, 2, k, ct.r_hp_sca ),1,var.nbr_m(ii+1)) ;
end
% Fm(1,:)=0;% does it right ?
Ymn.source = sph_harmonic( ct.M_th,1,source.theta,source.phi ) ;
Bmn.source = bsxfun(@times,Fm,Ymn.source') ;


%% Fabrication Matrice "C" spherical harmonic

[ Sphmic.theta, Sphmic.phi, Sphmic.r ] = cart2sph( Antenna.coord_vect(1,:)', Antenna.coord_vect(2,:)',zeros(length(Antenna.coord_vect(2,:)),1) ) ;
% conversion  to pherical coord possible problem with ref of cart2sph
Ymn.Mic = sph_harmonic( ct.M_th, ct.N_mic, Sphmic.theta, Sphmic.phi ) ;
% % --> calculate Bmn.*Ymn
% % --> change  hfm pour Hankel spherique
Pressure.difract=zeros(N.N_sweep,ct.N_mic);
Pressure.direct=zeros(N.N_sweep,ct.N_mic);

for jj=1:ct.N_mic
    var.Bmn_Ymn = bsxfun(@times,Bmn.source',Ymn.Mic(:,jj)) ;
    
    % var.pressure = zeros(N.N_sweep,ct.M_th) ;------» tic toc method choice
    for ii=0:ct.M_th
        var.sum_Bmn_Ymn = sum(var.Bmn_Ymn(1:var.m_sum_vect(ii+1),:),1) ;
        var.Hprim = 1i^(var.m_vect(ii+1))*(Bessel_sph(var.m_vect(ii+1),k.*Antenna.Rmicro(ii+1))-Bessel_sph_1_deriv(var.m_vect(ii+1),k.*ct.r_micsph)/Hankel_sph_1_deriv(var.m_vect(ii+1),ct.hankel_order,k.*ct.r_micsph)*Hankel_sph(var.m_vect(ii+1),ct.hankel_order,k.*Antenna.Rmicro(ii+1)));
        var.Bessel_int = 1i^(var.m_vect(ii+1))*(Bessel_sph(var.m_vect(ii+1),k.*Antenna.Rmicro(ii+1)));
        if ii==0;
            var.pressure_diffr = var.sum_Bmn_Ymn'.*var.Hprim ;
            var.pressure_direc = var.sum_Bmn_Ymn'.*var.Bessel_int ;

        else
            var.pressure_diffr = sum([var.pressure_diffr, var.sum_Bmn_Ymn'.*var.Hprim],2) ;
            var.pressure_direc = sum([var.pressure_direc, var.sum_Bmn_Ymn'.*var.Bessel_int],2) ;

        end
    end
    Pressure.difract(:,jj)=var.pressure_diffr;
    Pressure.direct(:,jj)=var.pressure_direc;

    
end
% var.pressure(1,:)=0;% does it r
%% Reconstruction

Pmes_mat = reshape(Pressure.direct	,size(Antenna.X_mat));
figure
pcolor(Antenna.x,Antenna.y,real(Pmes_mat));
shading interp
