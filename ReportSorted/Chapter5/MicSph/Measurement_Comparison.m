%% Measurement _Comparison.m plot the sound pressure field of two measurements with their error
% INPUT: data.EntrySig, load the .wav used for the measurement (sweep signal).
%
% Output1: data.OutSig_target, load the measured Target calibrated signal, calibrated (using MeasurementCalibration.m)
% Output2: data.OutSig_simu, load the measured calibrated signal to compare, calibrated (using MeasurementCalibration.m)
%
% The program take the IN and OUT data, perform the FRF and plot the sound
% pressure field at the desired frequency.
% The sound pressure field is obtained from the encoding of the
% measurement, considering a diffracted field by a rigid sphere. Then from
% the obtained Bmn, the sound pressure field is calculated
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


%%

[data.EntrySig,ct.Fs]=audioread('sweep_signal_mic_sph.wav');

[data.OutSig_target,ct.Fs_sca]=audioread('sph_hp51_calib.wav');
[data.OutSig_comp,~]=audioread('sph_51_ambisonique_simu_calib.wav');

%% Variable init
N.N_sweep=1;
N.N_mic=50;


ct.M=5;
ct.N_sweep_avg=10;
ct.c_air=340;
ct.r_micsph=0.07;
ct.hankel_order=2;

var.m_vect=0:ct.M;
var.m_sum_vect=(var.m_vect+1).^2;
var.nbr_m=(2.*var.m_vect)+1;

%% Calculate frf
[ Sphmic, ct.N_mic ] = CreateSpeakerSystem(ct.r_micsph);% create the sphere set up, sort in struc Array
System=struct('ct',ct,'var',var,'N',N,'Sphmic',Sphmic);
[ Frf.target,t, var.k] = SphmicFfrDataMeas( data.OutSig_target,data.EntrySig,System);
[ Frf.comp,~, ~] = SphmicFfrDataMeas( data.OutSig_comp,data.EntrySig,System);


%% Delay
for ii=1:1
    [b t.delay]=xcorr(Frf.comp.h_sig(:,ii),Frf.target.h_sig(:,ii),'coeff');
    [~,pos(ii)]=max(abs(b));
end

% Define the zero needed before the filter
D.FullDelay=-mean(t.delay(pos))+0;
N.NOrder=64;% point in the filter
D.maxi=floor(D.FullDelay)-N.NOrder/2;
N.N_sca=length(Frf.target.h_sig(:,1));

% Define the fractional delay of the filter
D.tabfrac=D.FullDelay-D.maxi;

h = Fractional_delay_lagrange_matrix(N.NOrder, D.tabfrac);% calculate delay filter
[ Frf.comp.h_sig ] = DelayImplementation2(h,Frf.comp.h_sig,D);
Frf.comp.h_sig_fft=fft(Frf.comp.h_sig);


%% _______________________ Affichage data__________________________________
%%_________________________________________________________________________


%% show map
close all
var.select_freq = [250 500 ]; %% select desirated frequency

for ii=1:length(var.select_freq)
    ct.N_mic=50;
    ct.pos = closest(var.select_freq(ii)*2*pi/ct.c_air,var.k);ct.k = var.k(ct.pos) ;
    ct.M=5;% order of the desired reconstructed field
    
    % set antenna points for reproduced field.
    ct.r=ct.M/(ct.k);
    ct.Size=ct.r*2;
    Antenna = AntennArray_defined_size(16,Size);
    
    
    Bmn.recons_target = Bmn_encoding_sph(Frf.target.h_sig_fft(ct.pos,:),Sphmic,ct,var,'tik');
    Bmn.recons = Bmn_encoding_sph(Frf.comp.h_sig_fft(ct.pos,:),Sphmic,ct,var,'tik');
    
    N.N_sweep=1;
    
    figure
    subplot(131)
    [Pressure.p_recons_target, ~ ] = Pressure_map_SphMic(ct.M,Bmn.recons_target,ct,N,var,Antenna);
    title('Target')
    var.cax2=caxis;
    hold on
    ct.theta_target=wavefrontarrow( Pressure.p_recons_target,Antenna );
    
    
    
    subplot(132)
    [Pressure.p_recons, ~ ] = Pressure_map_SphMic(ct.M,Bmn.recons,ct,N,var,Antenna);
    H = fspecial('gaussian',16,10);
    a(ii)=H(:).'.*Pressure.p_recons_target*pinv(Pressure.p_recons.*H(:).');
    Pressure.test=abs(a(ii))*Pressure.p_recons;
    Pressure_map_(Pressure.test,ct,var,0,Antenna);
    title('Comparison')
    
    [erreur,~] = erreur_n( Pressure.p_recons_target, Pressure.test );
    Pmes_mat = reshape(erreur ,size(Antenna.X_mat)) ;
    
    hold on
    [~,hfigc] = contour(Antenna.x,Antenna.y,real(Pmes_mat),[0 15]);
    set(hfigc, 'LineWidth',1.0,'Color', [1 1 1 ]);
    [~,hfigc] = contour(Antenna.x,Antenna.y,real(Pmes_mat),[0 30]);
    set(hfigc, 'LineWidth',1.0,'Color', [1 1 0]);
    caxis(var.cax2);
    ct.theta_decode= wavefrontarrow( Pressure.test,Antenna,'r');
    wavefrontarrow( Pressure.p_recons_target,Antenna );
    
    subplot(133)
    Pressure_map_(erreur,ct,var,1,Antenna);
    Pmes_mat = reshape(erreur ,size(Antenna.X_mat)) ;
    Pressure_map_(erreur,ct,var,1,Antenna);
    title('Error')
    
    
end

