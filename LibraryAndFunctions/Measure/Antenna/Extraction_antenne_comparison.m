%% Measurement _Comparison.m plot the sound pressure field of two measurements with their error
% INPUT: data.EntrySig, load the .wav used for the measurement (sweep signal).
%
% Output1: data.OutSig_target, load the measured Target calibrated signal, calibrated (using MeasurementCalibration.m)
% Output2: data.OutSig_simu, load the measured calibrated signal to compare, calibrated (using MeasurementCalibration.m)
%
% The program take the IN and OUT data, perform the FRF and plot the sound
% pressure field at the desired frequency.
%
% The phase shift between target and comp is take into account by the
% difference between the inpulse response time and compensated.
%
% The gain shift is compensated using LMS
%
% The comparison is effectued with the relative error plotted on the non target field
% subplot(311): target field, arrow intensity
% subplot(312): field to compare, arrow intensity(target and comp)
% subplot(313): error field, target, comp
%
% Define in the program :
% ct.fc= cut off frequency of the filter
% ct.M: Truncation order
% ct.N_sweep_avg: number of cycle used in the sweep
% var.select_freq: frequency at which you want the pressure field
%
% Auteur : Dupont Samuel
% Version : 3.0 June 2016

clear variables; close all;clc

ct.M=5;
ct.N_sweep_avg=10;
ct.c_air=340;
ct.r_micsph=0.07;
ct.hankel_order=2;

var.m_vect=0:ct.M;
var.m_sum_vect=(var.m_vect+1).^2;
var.nbr_m=(2.*var.m_vect)+1;

%% Extraction signal original
[data.sweep,ct.Fs_sca]=audioread('sweep_signal_mic_sph.wav');
ct.N_sweep_avg=10;


%% Extraction  Data source virtuelle, reelle

[data.antenna_comp, ct.Fs_sca ] =audioread('Antenna_sphwave_ambisonique_simu_hp51.w64');
[data.antenna_target, ct.Fs_sca ] =audioread('Antenna_sphwave_hp51.w64');

%% Calibration relative microphone
a=load('calib_Antenna_mic_31-06.mat');
a.relativ=a.calib./a.calib(1);
a.relativ(56)=1; %% due to 55 record but 56 points in the array
data.antenna_target = bsxfun(@times,data.antenna_target,a.relativ.');
data.antenna_comp = bsxfun(@times,data.antenna_comp,a.relativ.');

%% Define antenna set up
ct.pas_m =2.5*2.54e-2; % pas de l'antenne
N.nbrx_sca = 7; % nombre de micros par ligne
N.nbry_sca = 8; % nombre de micros par ligne
offset.x=0;
offset.y=-2.5*2.54e-2/2;
[ Antenna ] = AntennArray( ct.pas_m,N.nbrx_sca,N.nbry_sca,offset);N.N_mic=N.nbrx_sca*N.nbry_sca;


%% FRF computation
ct.fc=2000;
[~,H_antenna_target,~,~,t ] = FrfSystemFinal_data(data.antenna_target,data.sweep,N,ct,ct.fc);
[~,H_antenna_comp,~,~,t ] = FrfSystemFinal_data(data.antenna_comp,data.sweep,N,ct,ct.fc);

%% Compensation of the delay
[b t.delay]=xcorr(H_antenna_comp.h_sig(:,29),H_antenna_target.h_sig(:,29),'coeff');% 29 being the center point
[~,pos]=max(abs(b));


% Define the zero needed before the filter
D.FullDelay=mean(t.delay(pos));
N.NOrder=64;
D.maxi=floor(D.FullDelay)-N.NOrder/2;
N.N_sca=length(H_antenna_target.h_sig(:,1));
% Define the fractional delay of the filter
D.tabfrac=D.FullDelay-D.maxi;

h = Fractional_delay_lagrange_matrix(N.NOrder, D.tabfrac);
[ H_antenna_target.h_sig] = DelayImplementation2(h,H_antenna_target.h_sig,D);
H_antenna_target.h_sig_fft=fft(H_antenna_target.h_sig);

%% Compute impulse response movie
[~,pos]=max(H_antenna_target.h_sig(:,25));

% MovieMaker( H_antenna.h_sig,Antenna,pos-50,pos+200,ct)
%% initialisation of the setting
var.select_freq=[250 500 1500];


close all
for ii = 1:length(var.select_freq)
    pos = closest( var.select_freq(ii),t.Fsweep_avg );
    
    ct.k = 2*pi*t.Fsweep_avg(pos)/ct.c_air;
    r=ones(1,200)*ct.M/(ct.k); theta=linspace(0,2*pi,200);
    [x ,y ]=pol2cart(theta,r);
    
    figure
    
    subplot(131)
    grid_mat=reshape(H_antenna_target.h_sig_fft(pos,:),size(Antenna.Y_mat));
    pcolor(Antenna.x,Antenna.y,real(grid_mat)); % plot of the frame
    
    hold on
    
    plot(x,y,'--r')
    title('Target')
    xlabel('Position x [m]')
    ylabel('Position y [m]')
    shading interp
    axis equal
    axis tight
    axis([Antenna.x(1) Antenna.x(end) Antenna.y(1) Antenna.y(end)]	)
    hold off
    colorbar
    var.cax=caxis;
    wavefrontarrow(grid_mat,Antenna );
    
    subplot(132)
    H = fspecial('gaussian',56,10);
    a=H(1:55).*H_antenna_target.h_sig_fft(pos,1:55)*pinv(H_antenna_comp.h_sig_fft(pos,1:55).*H(1:55));
    H_antenna_comp.h_sig_fft(pos,1:55)=abs(a)*H_antenna_comp.h_sig_fft(pos,1:55);
    
    grid_mat_comp=reshape(H_antenna_comp.h_sig_fft(pos,:),size(Antenna.Y_mat));
    pcolor(Antenna.x,Antenna.y,real(grid_mat_comp)); % plot of the frame
    
    hold on
    
    plot(x,y,'--r')
    title('Decoded')
    
    xlabel('Position x [m]')
    ylabel('Position y [m]')
    shading interp
    axis equal
    axis tight
    axis([Antenna.x(1) Antenna.x(end) Antenna.y(1) Antenna.y(end)]	)
    hold off
    colorbar
    
    [erreur,~] = erreur_n( grid_mat, grid_mat_comp );
    Pmes_mat = reshape(erreur ,size(Antenna.X_mat)) ;
    hold on
    [~,hfigc] = contour(Antenna.x,Antenna.y,real(Pmes_mat),[0 15]);
    set(hfigc, 'LineWidth',1.0,'Color', [1 1 1 ]);
    [~,hfigc] = contour(Antenna.x,Antenna.y,real(Pmes_mat),[0 30]);
    set(hfigc, 'LineWidth',1.0,'Color', [1 1 0]);
    caxis(var.cax);
    ct.theta_simu=wavefrontarrow( grid_mat_comp,Antenna,'r' );
    wavefrontarrow( grid_mat,Antenna );
    
    subplot(133)
    Pressure_map_(erreur,ct,var,1,Antenna);
    colorbar
end


