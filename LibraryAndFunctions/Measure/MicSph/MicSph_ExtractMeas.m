%% MicSph_ExtractMeas.m Extract sound pressure field measurements at a certain frequency
%
% INPUT: data.EntrySig, load the .wav used for the measurement (sweep signal).
%
% Output: data.OutSig, load the measured calibrated signal, calibrated (using MeasurementCalibration.m)
%
% The program take the IN and OUT data, perform the FRF and plot the sound
% pressure field at the desired frequency.
% The sound pressure field is obtained from the encoding of the
% measurement, considering a diffracted field by a rigid sphere. Then from
% the obtained Bmn, the sound pressure field is calculated
%
% Define in the program :
% N.N_mic: number of mics on the spherical array
% ct.M: Truncation order
% ct.N_sweep_avg: number of cycle used in the sweep
% ct.r_micsph: radius of the spherical array
% ct.nbr_antenna_mic: number of points per row and line on the reconstructed field
% var.select_freq: frequency at which you want the pressure field
%
% Auteur : Dupont Samuel
% Version : 3.0 June 2016

clear all;close all;clc


[data.EntrySig,ct.Fs]=audioread('sweep_signal_mic_sph.wav');% Data IN
[data.OutSig,ct.Fs_sca]=audioread('sph_hp51_calib.wav');% Data OUT


N.N_sweep=1;
N.N_mic=50;
N.N_sweep=length(data.EntrySig);
N.N_meas=length(data.OutSig);

ct.M=5; % truncation order
ct.N_sweep_avg=10; % number of cycle in the sweep
ct.c_air=340;
ct.r_micsph=0.07;% radius of the microphone
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

%% _______________________ Affichage data__________________________________
%%_________________________________________________________________________

% Create antenna 
r=ct.M/(ct.k); % radius of the accurate reproduction
ct.Size=r*2;
ct.nbr_antenna_mic=20;
Antenna = AntennArray_defined_size(ct.nbr_antenna_mic,ct.Size); %antenna on which the field will be calculated

%% Operation on variable
var.select_freq = 3000 ; % select a frequency
ct.N_mic=50;
ct.pos = closest(var.select_freq*2*pi/ct.c_air,var.k);ct.k = var.k(ct.pos) ; 
ct.M=5;


Bmn.recons = Bmn_encoding_sph(HData.h_sig_fft(ct.pos,:),Sphmic,ct,var,'tik' );% encoding 

figure
ct.N_mic=Antenna.N_mic;
N.N_sweep=1;

% Plot of the sound pressure field
[Pressure.p_recons, ~ ] = Pressure_map_SphMic(ct.M,Bmn.recons(1:var.m_sum_vect(ii+2),1),ct,N,var,Antenna);title(sprintf('Recontructed Sound pressue field %i',var.select_freq));
wavefrontarrow( Pressure.p_recons_target,Antenna );% Intensity arrow


