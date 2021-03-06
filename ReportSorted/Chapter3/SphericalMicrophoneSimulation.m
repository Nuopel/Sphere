clear variables; close all;clc

% Simulate the behavior of a spherical microphone
% recording a monopole source positioned on the
% ambisonics set up
% Auteur : Dupont Samuel
% Version : 1.0 May 2016
%% Define constants
ct.r_hp_sca = 1.07 ;%rayon de la sphere
ct.r_micsph = 0.02;
ct.hankel_order =2;
ct.M_th = 5;
ct.M=5;
ct.nbr_M_th=(ct.M_th+1).^1;
ct.Fs=48000;
ct.c_air=340;

var.m_vect=0:ct.M_th;
var.m_sum_vect=(var.m_vect+1).^2;
var.nbr_m=(2.*var.m_vect)+1;

%% Choix de la source (Ae^j(wt-kx))
% [ source.sweep, t, ct, N ] = GenSweep(20, 20000, 4, ct ) ;
ct.k = 2*pi.*(1:20000)/340;
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
Bmn.recons = Bmn_encoding_sph( Pressure.difract,Sphmic,ct,var );


%% _______________________ Affichage data__________________________________
%%_________________________________________________________________________

%% Operation on variable
var.k=ct.k;
N.N_sweep=1;

%% Create antenna 
ct.pas_m = 4.5e-2; % pas de l'antenne
N.nbrx_sca =30; % nombre de micros par ligne
N.nbry_sca = 30; % nombre de micros par ligne
ct.N_mic=N.nbrx_sca*N.nbry_sca;
[ Antenna ] = AntennArray_defined_size( 34,1.5 ) ;

%% Conditioning signals
% Select frequency
var.select_freq = 1500 ;
var.pos = closest(var.select_freq*2*pi/ct.c_air,var.k);ct.k = var.k(var.pos) ;
Bmn.source_tronc = Bmn_monopole_encodage( ct.M,source,ct,var ) ;
[Pressure.monopole_exp ] = monopole_pressure(ct.k,source,Antenna);

figure(1);
% subplot(411)
% var=Pressure_map_( Pressure.monopole_exp,ct,var,0,Antenna );
% title('Real monopole')
subplot(223)
[Pressure.p_recons, ~ ] = Pressure_map_SphMic(ct.M,Bmn.recons(:,var.pos),ct,N,var,Antenna);
title('1500 Hz');
subplot(224)
[field ,~ ] = erreur_n(Pressure.p_recons,Pressure.monopole_exp);
[~]=Pressure_map_(field,ct,var,1,Antenna);
title('Error');


grid_mat_erreur=reshape(field,size(Antenna.X_mat));
[~,hfigc] = contour(Antenna.y,Antenna.x,grid_mat_erreur,[0 14]);
set(hfigc, 'LineWidth',1.0,'Color', [1 1 1]);


