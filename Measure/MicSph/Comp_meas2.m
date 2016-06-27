clear all;close all;clc


%%

[data.EntrySig,ct.Fs]=audioread('sweep_signal_mic_sph.wav');

% [data.OutSig_target,ct.Fs_sca]=audioread('sph_hp3_calib.wav');
% [data.OutSig_simu,~]=audioread('sph_3_ambisonique_simu_calib.wav');
% [data.OutSig,~]=audioread('sph_3_ambisonique_decode_encoder(18db)_decoder(26db)_calib.wav');
% 
[data.OutSig_target,ct.Fs_sca]=audioread('sph_hp51_calib.wav');
[data.OutSig_simu,~]=audioread('sph_51_ambisonique_simu_calib.wav');
[data.OutSig,~]=audioread('sph_51_ambisonique_decode_encoder(30db)_decoder(26db)_calib.wav');
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
[ Frf.comp,~, ~] = SphmicFfrDataMeas( data.OutSig,data.EntrySig,System);
[ Frf.comp_simu,~, ~] = SphmicFfrDataMeas( data.OutSig_simu,data.EntrySig,System);


%% Delay
for ii=1:1
    [b t.delay]=xcorr(Frf.comp.h_sig(:,ii),Frf.target.h_sig(:,ii),'coeff');
    [~,pos(ii)]=max(abs(b));
end

% Define the zero needed before the filter
D.FullDelay=-mean(t.delay(pos))+0;
N.NOrder=64;
D.maxi=floor(D.FullDelay)-N.NOrder/2;
N.N_sca=length(Frf.target.h_sig(:,1));
% Define the fractional delay of the filter
D.tabfrac=D.FullDelay-D.maxi;

h = Fractional_delay_lagrange_matrix(N.NOrder, D.tabfrac);
[ Frf.comp.h_sig ] = DelayImplementation2(h,Frf.comp.h_sig,D);
Frf.comp.h_sig_fft=fft(Frf.comp.h_sig);

for ii=1:1
    [b t.delay]=xcorr(Frf.comp_simu.h_sig(:,ii),Frf.target.h_sig(:,ii),'coeff');
    [~,pos(ii)]=max(abs(b));
end

% Define the zero needed before the filter
D.FullDelay=-mean(t.delay(pos))+0;
N.NOrder=64;
D.maxi=floor(D.FullDelay)-N.NOrder/2;
N.N_sca=length(Frf.target.h_sig(:,1));
% Define the fractional delay of the filter
D.tabfrac=D.FullDelay-D.maxi;

h = Fractional_delay_lagrange_matrix(N.NOrder, D.tabfrac);
[ Frf.comp_simu.h_sig ] = DelayImplementation2(h,Frf.comp_simu.h_sig,D);
Frf.comp_simu.h_sig_fft=fft(Frf.comp_simu.h_sig);
%% _______________________ Affichage data__________________________________
%%_________________________________________________________________________


%% show map
close all
var.select_freq = [250 500  750 1000 1500 2000 4000];

for ii=1:length(var.select_freq)
    ct.N_mic=50;
    ct.pos = closest(var.select_freq(ii)*2*pi/ct.c_air,var.k);ct.k = var.k(ct.pos) ;
    ct.M=5;
    
    r=ct.M/(ct.k);
    Size=r*2;
    Antenna = AntennArray_defined_size(16,Size);
    
    
    Bmn.recons_target = Bmn_encoding_sph(Frf.target.h_sig_fft(ct.pos,:),Sphmic,ct,var,1 );
    Bmn.recons = Bmn_encoding_sph(Frf.comp.h_sig_fft(ct.pos,:),Sphmic,ct,var,1 );
    Bmn.recons_simu = Bmn_encoding_sph(Frf.comp_simu.h_sig_fft(ct.pos,:),Sphmic,ct,var,1 );
    
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
    Pressure.test=abs(a(ii))*Pressure.p_recons*exp(1i*pi);
    Pressure_map_(Pressure.test,ct,var,0,Antenna);
    title('Decoded ')
    
    [erreur,~] = erreur_n( Pressure.p_recons_target, Pressure.test );
    Pmes_mat = reshape(erreur ,size(Antenna.X_mat)) ;
    % figure
    % Pressure_map_(erreur,ct,var,1,Antenna);
    hold on
    [~,hfigc] = contour(Antenna.x,Antenna.y,real(Pmes_mat),[0 15]);
    set(hfigc, 'LineWidth',1.0,'Color', [1 1 1 ]);
    [~,hfigc] = contour(Antenna.x,Antenna.y,real(Pmes_mat),[0 30]);
    set(hfigc, 'LineWidth',1.0,'Color', [1 1 0]);
    caxis(var.cax2);
    ct.theta_decode= wavefrontarrow( Pressure.test,Antenna,'r' );
    wavefrontarrow( Pressure.p_recons_target,Antenna );
    subplot(133)
    [Pressure.p_recons, ~ ] = Pressure_map_SphMic(ct.M,Bmn.recons_simu,ct,N,var,Antenna);
    a(ii)=H(:).'.*Pressure.p_recons_target*pinv(Pressure.p_recons.*H(:).');
    Pressure.test=(a(ii))*Pressure.p_recons;
    Pressure_map_(Pressure.test,ct,var,0,Antenna);
    title('Simulated ')
    
    [erreur,~] = erreur_n( Pressure.p_recons_target, Pressure.test );
    Pmes_mat = reshape(erreur ,size(Antenna.X_mat)) ;
    % figure
    % Pressure_map_(erreur,ct,var,1,Antenna);
    hold on
    [~,hfigc] = contour(Antenna.x,Antenna.y,real(Pmes_mat),[0 15]);
    set(hfigc, 'LineWidth',1.0,'Color', [1 1 1 ]);
    [~,hfigc] = contour(Antenna.x,Antenna.y,real(Pmes_mat),[0 30]);
    set(hfigc, 'LineWidth',1.0,'Color', [1 1 0]);
    caxis(var.cax2);
    ct.theta_simu=wavefrontarrow( Pressure.test,Antenna,'r' );
    
    wavefrontarrow( Pressure.p_recons_target,Antenna );
    
    Pressure.test=reshape(Pressure.test.',size(Antenna.X_mat));
end

