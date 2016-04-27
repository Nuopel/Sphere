clear variables; close all;clc

% Simulate the behavior of a pherical microphone
% recording a monopole source positioned on the
% ambisonics set up

%% Define constants
ct.r_hp_sca = 1.07 ;%rayon de la sphere
ct.r_micsph = 0.07;
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
ct.k = 2*pi*1500/340;
N.N_sweep=1;

%% Define ambisonics set up
[ ArraySpeaker, ct.N_speaker ] = CreateSpeakerSystem(ct.r_hp_sca);% create the sphere set up, sort in struc Array


%% Define Spherical microphone set up (same as speaker but smaller)
[ Sphmic, ct.N_mic ] = CreateSpeakerSystem(ct.r_micsph);% create the sphere set up, sort in struc Array


%% POSITION DE LA SOURCE : select a speaker
speaker = 12;
source.x = ArraySpeaker.x(speaker) ; source.y = ArraySpeaker.y(speaker) ; source.z = ArraySpeaker.z(speaker) ;
[ source.theta, source.phi, source.r ] = cart2sph(source.x, source.y, source.z) ; clear speaker ;

%% Encoding signal on Bmn coefficient

Bmn.source = Bmn_monopole_encodage(ct.M_th,source,ct,var ) ;


%% Decoding to calculate pressure on the spherical microphone

[ Sphmic,Pressure ] = Decoding_pressure_microphone( Bmn,Sphmic,N,ct,var );
% a=load('psca.mat');
% Pressure.difract=a.Expression1.';
%% Encoding from microphone pressure
Bmn.recons = Bmn_encoding_sph( Pressure,Sphmic,ct,N,var );

%% affichage data
Bmn.source_tronc = Bmn_monopole_encodage(ct.M,source,ct,var ) ;
Pressure.p_recons = Pressure_map_SphMic(ct.M,Bmn.recons,ct,N,var);title('Reconstruction sphere mic');
Pressure.p_target = Pressure_map_SphMic(ct.M,Bmn.source_tronc.',ct,N,var);title('Reconstruction troncation');
Pressure.monopole = Pressure_map_SphMic(ct.M_th,Bmn.source.',ct,N,var);title('Reconstruction full');


[field ,norm_e ]=erreur_n(Pressure.p_target,Pressure.p_recons);
Pressure_map_(field,ct,1)


%%
Bmn.source_tronc=permute(Bmn.source_tronc,[2 1]);
[ erreur ,norm_e ]=erreur_n(Bmn.source_tronc, Bmn.recons)
% ct=rmfield(ct,'Fs2_sca');
%%
a=load('psca.mat');
a.Expression1
% ct.k=2*pi*[1:20000]./340;
% for ii=0:5
%     for jj=1:length(ct.k)
%         a(jj)=((ct.k(jj).*ct.r_micsph).^2*Hankel_sph_1_deriv(ii,ct.hankel_order,ct.k(jj).*ct.r_micsph));
%     end
%             semilogx(db(a));hold on;
% 
% end