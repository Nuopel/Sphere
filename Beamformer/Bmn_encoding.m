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
ct.N=2;
ct.nbr_M_th=(ct.M_th+1).^1;
ct.Fs=48000;
ct.c_air=340;

var.m_vect=0:ct.M_th;
var.m_sum_vect=(var.m_vect+1).^2;
var.nbr_m=(2.*var.m_vect)+1;

%% Choix de la source (Ae^j(wt-kx))
[ source.sweep, t, ct, N ] = GenSweep(20, 20000, 5, ct ) ;
source.sweep=0.25*source.sweep;
ct.k = 2*pi.*t.F_sweep/340;
N.N_sweep=length(ct.k);

%% Define ambisonics set up
[ ArraySpeaker, ct.N_speaker ] = CreateSpeakerSystem(ct.r_hp_sca);% create the sphere set up, sort in struc Array

%% POSITION DES SOURCES : select speakers
source.hp(1)=13;
source.hp(2)=12;
% source.hp(1)=3;
% source.hp(2)=4;
for ii=1:ct.N
[ source.theta(:,ii), source.phi(:,ii), source.r(:,ii) ] = cart2sph(ArraySpeaker.x(source.hp(ii)),ArraySpeaker.y...
                                         (source.hp(ii)), ArraySpeaker.z(source.hp(ii))) ; clear speaker ;
end
 Bmn.sources = Bmn_encoding_source_plane( source.sweep,source,ct,N,var );
Bmn.Bmn = sum(Bmn.sources,3);

% PathWriter('Bmn_12_13.txt',sprintf('%i',Bmn.Bmn))
audiowrite('Bmn_12_13.wav',Bmn.Bmn.',ct.Fs)