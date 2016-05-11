clear variables; close all;clc

% Simulate the behavior of a pherical microphone
% recording a monopole source positioned on the
% ambisonics set up

%% Define constants
ct.r_hp_sca = 1.07 ;%rayon de la sphere
ct.r_micsph = 0.02;
ct.hankel_order =2;
ct.M_th = 15;
ct.M=5;
ct.nbr_M_th=(ct.M_th+1).^1;
ct.Fs=48000;
ct.c_air=340;

var.m_vect=0:ct.M_th;
var.m_sum_vect=(var.m_vect+1).^2;
var.nbr_m=(2.*var.m_vect)+1;

%% Choix de la source (Ae^j(wt-kx))
[ source.sweep, t, ct, N ] = GenSweep(20, 20000, 4, ct ) ;
ct.k = 2*pi.*[755 ]/340;
N.N_sweep=length(ct.k);

%% Define ambisonics set up
[ ArraySpeaker, ct.N_speaker ] = CreateSpeakerSystem(ct.r_hp_sca);% create the sphere set up, sort in struc Array

%% Define Spherical microphone set up (same as speaker but smaller)
[ Sphmic, ct.N_mic ] = CreateSpeakerSystem(ct.r_micsph);% create the sphere set up, sort in struc Array

%% POSITION DE LA SOURCE : select a speaker
source.x = 0 ; source.y = -1.07; source.z = 0 ;
[ source.theta, source.phi, source.r ] = cart2sph(source.x, source.y, source.z) ; clear speaker ;

%% Encoding signal on Bmn coefficient
Bmn.source = Bmn_monopole_encodage( ct.M_th,source,ct,var ) ;

%% Decoding to calculate pressure on the spherical microphone
[ Sphmic,Pressure ] = Decoding_pressure_microphone( Bmn,Sphmic,N,ct,var );

%% Encoding from microphone pressure
Bmn.recons = Bmn_encoding_sph( Pressure,Sphmic,ct,N,var );


%%_________________________________________________________________________

%% _______________________ Affichage data__________________________________
%%_________________________________________________________________________

%% Operation on variable
var.k=ct.k;
N.N_sweep=1;

%% Create antenna 
ct.pas_m = 2e-2; % pas de l'antenne
N.nbrx_sca =50; % nombre de micros par ligne
N.nbry_sca = 50; % nombre de micros par ligne
ct.N_mic=N.nbrx_sca*N.nbry_sca;
[ Antenna ] = AntennArray( ct.pas_m,N.nbrx_sca,N.nbry_sca) ;

%% Conditioning signals
% Select frequency
var.select_freq = 100 ;
var.pos = closest(var.select_freq*2*pi/ct.c_air,var.k);ct.k = var.k(var.pos) ; 
Bmn.source_tronc = Bmn_monopole_encodage( ct.M,source,ct,var ) ;
[Pressure.monopole_exp ] = monopole_pressure(ct.k,source,Antenna);

figure(1);
subplot(311)
var=Pressure_map_( Pressure.monopole_exp,ct,Antenna,var );
subplot(312)
[Pressure.p_recons, ~ ] = Pressure_map_SphMic(ct.M,Bmn.recons(:,var.pos),ct,N,var,Antenna);title('Reconstruction sphere mic');
subplot(313)
[field ,~ ] = erreur_n(Pressure.p_recons,Pressure.monopole_exp);
[~]=Pressure_map_(field,ct,Antenna,var,1);
figure(2)
[Pressure.p_target, ~ ] = Pressure_map_SphMic(ct.M,Bmn.source_tronc.',ct,N,var,Antenna);title('Reconstruction troncation');
figure(3)
[Pressure.monopole, ~ ] = Pressure_map_SphMic(ct.M_th,Bmn.source(var.pos,:).',ct,N,var,Antenna);title('Reconstruction full');


% Bmn.source_tronc=permute(Bmn.source_tronc,[2 1]);
% [ erreur ,norm_e ] = erreur_n(Bmn.source_tronc, Bmn.recons);
