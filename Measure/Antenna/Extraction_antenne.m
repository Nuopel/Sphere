clear variables; close all;clc

ct.M=5;
ct.N_sweep_avg=10;
ct.c_air=340;
ct.r_micsph=0.07;
ct.hankel_order=2;

var.m_vect=0:ct.M;
var.m_sum_vect=(var.m_vect+1).^2;
var.nbr_m=(2.*var.m_vect)+1;
%% Define antenna set up
ct.pas_m =2.5*2.54e-2; % pas de l'antenne
N.nbrx_sca = 7; % nombre de micros par ligne
N.nbry_sca = 8; % nombre de micros par ligne
offset.x=0;
offset.y=-2.5*2.54e-2/2;
[ Antenna ] = AntennArray( ct.pas_m,N.nbrx_sca,N.nbry_sca,offset);N.N_mic=N.nbrx_sca*N.nbry_sca;
% contruction antenne 2.5 po , 8 mic


%% Extraction signal original
[data.sweep,ct.Fs_sca]=audioread('sweep_signal_mic_sph.wav');
ct.N_sweep_avg=10;


%% Extraction  Data source virtuelle, reelle
[data.antenna_simu, ct.Fs_sca ] =audioread('Antenna_sphwave_ambisonique_simu_90.w64');
[data.antenna_target, ct.Fs_sca ] =audioread('Antenna_sphwave_hp3.w64');
[data.antenna_decode, ct.Fs_sca ] =audioread('Antenna_sphwavehp3_ambisonique_decode_encoder(18db)_decoder(26db).w64');

% 
% [data.antenna_simu, ct.Fs_sca ] =audioread('Antenna_sphwave_ambisonique_simu_hp51.w64');
% [data.antenna_target, ct.Fs_sca ] =audioread('Antenna_sphwave_hp51.w64');
% [data.antenna_decode, ct.Fs_sca ] =audioread('Antenna_51_ambisonique_decode_encoder(30db)_decoder(26db).w64');

%% Calibration relatve microphone
a=load('calib_Antenna_mic_31-06.mat');
a.relativ=a.calib./a.calib(1);
a.relativ(56)=1;

data.antenna_target = bsxfun(@times,data.antenna_target,a.relativ.');
data.antenna_decode = bsxfun(@times,data.antenna_decode,a.relativ.');
data.antenna_simu = bsxfun(@times,data.antenna_simu,a.relativ.');

ct.fc=2300;
[~,H_antenna_target,~,~,~ ] = FrfSystemFinal_data(data.antenna_target,data.sweep,N,ct,ct.fc);
[~,H_antenna_decode,~,~,t ] = FrfSystemFinal_data(data.antenna_decode,data.sweep,N,ct,ct.fc);
[~,H_antenna_simu,~,~,t ] = FrfSystemFinal_data(data.antenna_simu,data.sweep,N,ct,ct.fc);



[b t.delay]=xcorr(H_antenna_decode.h_sig(:,29),H_antenna_target.h_sig(:,29),'coeff');
[~,pos]=max(abs(b));


% Define the zero needed before the filter
D.FullDelay=mean(t.delay(pos));
N.NOrder=64;
D.maxi=floor(D.FullDelay)-N.NOrder/2;
N.N_sca=length(H_antenna_target.h_sig(:,1));
% Define the fractional delay of the filter
D.tabfrac=D.FullDelay-D.maxi;

% h = Fractional_delay_lagrange_matrix(N.NOrder, D.tabfrac);
% [ H_antenna_decode.h_sig] = DelayImplementation2(h,H_antenna_decode.h_sig,D);
% H_antenna_decode.h_sig_fft=fft(H_antenna_decode.h_sig);
h = Fractional_delay_lagrange_matrix(N.NOrder, D.tabfrac);
[ H_antenna_target.h_sig] = DelayImplementation2(h,H_antenna_target.h_sig,D);
H_antenna_target.h_sig_fft=fft(H_antenna_target.h_sig);


[b t.delay]=xcorr(H_antenna_simu.h_sig(:,29),H_antenna_target.h_sig(:,29),'coeff');
[~,pos]=max(abs(b));


% Define the zero needed before the filter
D.FullDelay=-mean(t.delay(pos));
N.NOrder=64;
D.maxi=floor(D.FullDelay)-N.NOrder/2;
N.N_sca=length(H_antenna_target.h_sig(:,1));
% Define the fractional delay of the filter
D.tabfrac=D.FullDelay-D.maxi;
[ H_antenna_simu.h_sig] = DelayImplementation2(h,H_antenna_simu.h_sig,D);
H_antenna_simu.h_sig_fft=fft(H_antenna_simu.h_sig);

%%
[~,pos]=max(H_antenna_target.h_sig(:,25));

% MovieMaker( H_antenna.h_sig,Antenna,pos-50,pos+200,ct)
%% initialisation of the setting
var.select_freq=[250 500 1000 2000];
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

subplot(132)
H = fspecial('gaussian',56,10);
a=H(1:55).*H_antenna_target.h_sig_fft(pos,1:55)*pinv(H_antenna_decode.h_sig_fft(pos,1:55).*H(1:55));
H_antenna_decode.h_sig_fft(pos,1:55)=abs(a)*H_antenna_decode.h_sig_fft(pos,1:55).*exp(1i*pi);

grid_mat_decode=reshape(H_antenna_decode.h_sig_fft(pos,:),size(Antenna.Y_mat));
pcolor(Antenna.x,Antenna.y,real(grid_mat_decode)); % plot of the frame

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

[erreur,~] = erreur_n( grid_mat, grid_mat_decode );
Pmes_mat = reshape(erreur ,size(Antenna.X_mat)) ;
% figure
% Pressure_map_(erreur,ct,var,1,Antenna);
    hold on
    [~,hfigc] = contour(Antenna.x,Antenna.y,real(Pmes_mat),[0 15]);
    set(hfigc, 'LineWidth',1.0,'Color', [1 1 1 ]);
        [~,hfigc] = contour(Antenna.x,Antenna.y,real(Pmes_mat),[0 30]);
    set(hfigc, 'LineWidth',1.0,'Color', [1 1 0]);
caxis(var.cax);

subplot(133)
H = fspecial('gaussian',56,10);
a=H(1:55).*H_antenna_target.h_sig_fft(pos,1:55)*pinv(H_antenna_simu.h_sig_fft(pos,1:55).*H(1:55));
H_antenna_simu.h_sig_fft(pos,1:55)=H_antenna_simu.h_sig_fft(pos,1:55);

grid_mat_simu=reshape(H_antenna_simu.h_sig_fft(pos,:),size(Antenna.Y_mat));
pcolor(Antenna.x,Antenna.y,real(grid_mat_simu)); % plot of the frame

hold on

plot(x,y,'--r')
title('Simulated')

xlabel('Position x [m]')
ylabel('Position y [m]')
shading interp
axis equal
axis tight
axis([Antenna.x(1) Antenna.x(end) Antenna.y(1) Antenna.y(end)]	)
hold off
colorbar

% subplot(313)
[erreur,~] = erreur_n( grid_mat, grid_mat_simu );
Pmes_mat = reshape(erreur ,size(Antenna.X_mat)) ;
    hold on
    [~,hfigc] = contour(Antenna.x,Antenna.y,real(Pmes_mat),[0 15]);
    set(hfigc, 'LineWidth',1.0,'Color', [1 1 1 ]);
        [~,hfigc] = contour(Antenna.x,Antenna.y,real(Pmes_mat),[0 30]);
    set(hfigc, 'LineWidth',1.0,'Color', [1 1 0]);
caxis(var.cax);
end


