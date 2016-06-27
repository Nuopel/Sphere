%% Extraction_antenne.m Plot the sound pressure field from the measurement of an antenna
%
% INPUT: data.sweep, load the .wav used for the measurement (sweep signal).
%
% Output: data.antenna, load the measured calibrated signal, calibrated (using MeasurementCalibration.m)
%
% Define in the program :
% ct.pas_m: space in between microphones
% ct.N_sweep_avg: number of cycle used in the sweep
% var.select_freq: frequency at which you want the pressure field
%
% The program plot the pressure field at a specified frequency from the
% measurement set in data.antenna using the FRF ref to data.sweep
%
%
% Auteur : Dupont Samuel
% Version : 3.0 June 2016
clear variables; close all;clc

ct.c_air=340;
ct.M=5;


%% Define antenna set up
ct.pas_m =2.5*2.54e-2; % pas de l'antenne
N.nbrx_sca = 7; % nombre de micros par ligne
N.nbry_sca = 8; % nombre de micros par ligne
offset.x=0;
offset.y=-2.5*2.54e-2/2;
[ Antenna ] = AntennArray( ct.pas_m,N.nbrx_sca,N.nbry_sca,offset);N.N_mic=N.nbrx_sca*N.nbry_sca;

%% Extraction signal original
[data.sweep,ct.Fs_sca]=audioread('sweep_signal_mic_sph.wav');
ct.N_sweep_avg=10; % number of cycle in the sweep


%% Extraction  Data source virtuelle, reelle
[data.antenna, ct.Fs_sca ] =audioread('Antenna_sphwavehp3_ambisonique_decode_encoder(18db)_decoder(26db).w64');

[~,H_antenna,~,~,t ] = FrfSystemFinal_data(data.antenna,data.sweep,N,ct,1200);
%% uncomment to see impulse response movie
% [~,pos]=max(H_antenna.h_sig(:,5));
% MovieMaker( H_antenna.h_sig,Antenna,pos-100,pos+200 ,ct)


%% _______________________ Affichage data__________________________________
%%_________________________________________________________________________

%initialisation
var.select_freq=1300;

pos = closest(var.select_freq,t.Fsweep_avg);
grid_mat=(reshape(H_antenna.h_sig_fft(pos,:),size(Antenna.Y_mat)));

figure
pcolor(Antenna.x,Antenna.y,real(grid_mat)); % plot of the frame
shading interp
colorbar

% plot accurate area circle
ct.k = 2*pi*t.Fsweep_avg(pos)/ct.c_air;
r=ones(1,200)*ct.M/(ct.k); theta=linspace(0,2*pi,200);
[x ,y ]=pol2cart(theta,r);
hold on
plot(x,y,'--r')
xlabel('Position x [m]')
ylabel('Position y [m]')
shading interp
axis equal
axis tight
axis([Antenna.x(1) Antenna.x(end) Antenna.y(1) Antenna.y(end)]	)
hold off

% plot intensity arrow
wavefrontarrow(grid_mat,Antenna);
