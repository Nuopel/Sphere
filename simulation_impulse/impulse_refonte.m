% Simulation de l'impulsion d'une onde spherique sur un array de microphone
% depuis l'enregistrement des signaux envoye aux haut-parleurs
% todo: implement speaker signal calculation
% Auteur : Dupont Samuel
% Version : 1.0 Mars 2016


close all; clear variables; clc

%% Importation donnee source et entree

[sig_hp_mat,ct.Fs_sca]=audioread('source_simu_2m_45d.wav');% recorded speaker data
[sweep_signal_vect,~]=audioread('sweep_signal_antenne.wav');% sweep signal


%% Definition des constantes et des vecteurs de representation
N.N_sca=length(sig_hp_mat);N.Nfft=2056;
N.N_sweep=length(sweep_signal_vect);

ct.c_air=340;
ct.Ts_sca=1/ct.Fs_sca;% sampling time
ct.dfe=ct.Fs_sca/N.Nfft; % delta Fe for Nfft
ct.dfe_sweep=ct.Fs_sca/N.N_sweep;

t.Tsig_hp=0:1*ct.Ts_sca:(length(sig_hp_mat)-1)*ct.Ts_sca;%time axis
t.Fnfft=0:ct.dfe:(N.Nfft-1)*ct.dfe;% frequency axis Nfft
t.F_sweep=0:ct.dfe_sweep:(N.N_sweep-1)*ct.dfe_sweep;% Frequency axis signal


%% Structure Hp
% take position of the unitary system and change it to cr.r_hp_sca
ct.r_hp_sca=1.07;%rayon de la sphere
imp=load('coords.mat');imp.coords(:,1:3)=imp.coords(:,1:3)*ct.r_hp_sca;
% A=rotz(45*pi/180);
% for ii=1:length(imp.coords)
% imp.coords(ii,1:3)=imp.coords(ii,1:3)*A';
% end
[ArraySpeaker] = CreateSpeakerStructure( imp.coords(:,1), imp.coords(:,2), imp.coords(:,3), imp.coords(:,4) );
clear imp;


%% Antenne de mesure
ct.pas_m = 2.5*2.54e-2; % pas de l'antenne
N.nbrx_sca =8; % nombre de micros par ligne
N.nbry_sca =-7; % nombre de micros par ligne
[ Antenna ] = AntennArray( ct.pas_m,N.nbrx_sca,N.nbry_sca);N.N_mic=abs(N.nbrx_sca*N.nbry_sca);
% contruction antenne 2.5 po , 8 mic



%% Calcul des distances micros-hp
D.R_mat=sqrt(sum(cat(3,bsxfun(@minus,Antenna.coord_vect(1,:),ArraySpeaker.x),bsxfun(@minus,Antenna.coord_vect(2,:),ArraySpeaker.y),bsxfun(@minus,zeros(1,numel(Antenna.coord_vect(2,:))),ArraySpeaker.z)).^2,3));
%create a 3D array with 3 dimensions refering to  x,y,z then sum the square
%value in order to obtain the absolute value



%% Propagation des signaux, creation des outils

%order of the lagrange operator for fractional delay implementation
N.NOrder=64;

%calculation of the number of samples for delay
[D,System,N]=CreateFilterFracDelay(D,ArraySpeaker,ct,N);




%verification signaux
t.Tfir=0:1*ct.Ts_sca:N.NOrder*ct.Ts_sca;%time axis
figure(1)
subplot(211);plot(t.Tfir,System.ht(:,:,1));
ct.dfe_Norder=ct.Fs_sca/N.NOrder;
t.HFfft=0:ct.dfe_Norder:(N.NOrder)*ct.dfe_Norder;
subplot(212);semilogx(t.HFfft,20*log10(abs(System.HF(:,:,1))));xlim([20 ct.Fs_sca/2]);



%% Propagation des signaux, application depuis les signaux des hp
tic
sig_mic_mat  = DelayImplementation(System.ht,sig_hp_mat,ArraySpeaker.L,N,D);
toc
% verification signaux
figure(2)
t.Tsig_mic=0:1*ct.Ts_sca:(length(sig_mic_mat)-1)*ct.Ts_sca;%time axis
plot(t.Tsig_mic,sig_mic_mat(:,:,1));

%% Processing microphone data (obtain fft) for the total system
ct.N_sweep_avg=3 ;
[ sig_mic_mat ,System,N,ct,t] = FrfSystemFinal_simu(sig_mic_mat,sweep_signal_vect,System,N,t,ct);

% verification signaux
figure(3)
subplot(212);semilogx(t.Fsweep_avg,20*log10(abs(System.h_sig_fft(:,1))))
xlim([20 20000])
subplot(211);
plot(System.h_sig(:,1))



 %% Affichage data
 System.h_sig=System.h_sig(1:10000,:)
MovieMaker(System.h_sig,Antenna,7900,8000)
