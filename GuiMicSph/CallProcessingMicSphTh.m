function [ Pressure ] = CallProcessingMicSphTh(Antenna,ct  )

%% Constants
ct.hankel_order=2;
ct.N_mic=50;
ct.M_th=30;
ct.r_hp_sca=1.07;
var.m_vect=0:ct.M_th;
var.m_sum_vect=(var.m_vect+1).^2;
var.nbr_m=(2.*var.m_vect)+1;
ct.r_micsph=0.07;

N.N_sweep=1;
%% Spherical wave target
[ source.x, source.y, source.z ] = sph2cart( -ct.Theta*pi/180+45*pi/180, ct.Phi*pi/180, ct.R ) ;
[ source.theta, source.phi, source.r ] = cart2sph( source.x, source.y, source.z );

%% Spherical wave target on sph mic
%% Encoding signal on Bmn coefficient for diffract wave
Bmn.source = Bmn_monopole_encodage( ct.M_th,source,ct,var ) ;
%% Define Spherical microphone set up (same as speaker but smaller)
[ Sphmic, ct.N_mic ] = CreateSpeakerSystem(ct.r_micsph);% create the sphere set up, sort in struc Array

% Decoding to calculate pressure on the spherical microphone
[ Sphmic,Pressure ] = Decoding_pressure_microphone( Bmn,Sphmic,N,ct,var );
%% Encoding from microphone pressure
Bmn.recons = Bmn_encoding_sph( Pressure,Sphmic,ct,N,var );

        %% Calculate pressure from Bmn coefficient
        ct.N_mic=56;
Pressure.TargetAmbisonics = Decoding_pressure_field(ct.M,Bmn.recons,Antenna,ct,var,N ) ;
[Pressure.monopole ] = monopole_pressure(ct.k,source,Antenna);


end

