clear all;close all;clc


[data.EntrySig,ct.Fs]=audioread('../../DataMicSph/Sweep_measure_sphwave.wav');
[data.OutSig,ct.Fs_sca]=audioread('../../DataMicSph/MICSPH_HP2_TEST.w64');

[data.outcalib,ct.Fs_sca]=audioread('../../DataMicSph/calib_sinus_100.w64');
[data.incalib,ct.Fs_sca]=audioread('../../DataMicSph/sinus_calib_100.wav');

data.calib=load('../Calibration/calib_memsbedev_mic.mat')
data.calib=data.calib.max;data.calib=data.calib/data.calib(1);
data.outcalib=bsxfun(@times,data.outcalib,data.calib.');

N.N_sweep=1;
N.N_mic=50;
N.N_sweep=length(data.EntrySig);
N.N_meas=length(data.OutSig);

ct.M=5;
ct.N_sweep_avg=10;
ct.c_air=340;
ct.r_micsph=0.07;
ct.hankel_order=2;

var.m_vect=0:ct.M;
var.m_sum_vect=(var.m_vect+1).^2;
var.nbr_m=(2.*var.m_vect)+1;
%% Define Spherical microphone set up (same as speaker but smaller)
[ Sphmic, ct.N_mic ] = CreateSpeakerSystem(0.07);% create the sphere set up, sort in struc Array

%% Look if data size In equal data size Out
if N.N_meas<N.N_sweep
    data.OutSig=[data.OutSig ;zeros(N.N_sweep-N.N_meas,50)];
    fprintf('%i zeros have been added to the measurement data \n', N.N_sweep-N.N_meas)
else if N.N_meas>N.N_sweep
      data.EntrySig=[data.EntrySig ;zeros(N.N_sweep-N.N_meas,50)];
      fprintf('%i zeros have been added to the sweep data \n Be careful measurement longer than  entry signal \n', -N.N_sweep+N.N_meas)
    end
end

%% Average data
% Select frequency
[~,HData,~,~,t ] = FrfSystemFinal_data(data.OutSig,data.EntrySig,N,ct);
ct.k=2*pi*t.Fsweep_avg/ct.c_air;
var.k=ct.k;

%%

    
%%_________________________________________________________________________

%% _______________________ Affichage data__________________________________
%%_________________________________________________________________________


% Create antenna 
ct.pas_m = 2e-2; % pas de l'antenne
N.nbrx_sca =50; % nombre de micros par ligne
N.nbry_sca = 50; % nombre de micros par ligne
ct.N_mic=N.nbrx_sca*N.nbry_sca;
[ Antenna ] = AntennArray( ct.pas_m,N.nbrx_sca,N.nbry_sca) ;

% Conditioning signals

%% Operation on variable
close
var.select_freq = 5000 ;
ct.N_mic=50;
ct.pos = closest(var.select_freq*2*pi/ct.c_air,var.k);ct.k = var.k(ct.pos) ; 
Bmn.recons = Bmn_encoding_sph(HData.h_sig_fft(ct.pos,:),Sphmic,ct,var );

figure(1);
ct.N_mic=Antenna.N_mic;
N.N_sweep=1;
[Pressure.p_recons, ~ ] = Pressure_map_SphMic(ct.M,Bmn.recons(:,1),ct,N,var,Antenna);title('Reconstruction sphere mic');


