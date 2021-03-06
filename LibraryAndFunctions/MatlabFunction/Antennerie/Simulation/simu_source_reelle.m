function [  sig_mic_mat ,System,ArraySpeaker,N,ct,t ] = simu_source_reelle( sig_hp_mat, sweep_signal_vect,Antenna,ct,N)

% Generate the impulse response from a set of recorded speaker signal of the ambisonic sphere 50 speakers;

%% Definition des constantes et des vecteurs de representation
N.N_sca=length(sig_hp_mat);N.Nfft=2056;
N.N_sweep=length(sweep_signal_vect);

ct.c_air=340;
ct.Ts_sca=1/ct.Fs_sca;% sampling time
ct.dfe=ct.Fs_sca/N.Nfft; % delta Fe for Nfft
ct.dfe_sweep=ct.Fs_sca/N.N_sweep;

%% Structure Hp

% take position of the unitary system and change it to cr.r_hp_sca
ct.r_hp_sca=1.07;%rayon de la sphere
imp=load('data/coords.mat');imp.coords(:,1:3)=imp.coords(:,1:3)*ct.r_hp_sca;imp.coords(:,2)=-imp.coords(:,2);
A=rotz(45*pi/180);
for ii=1:length(imp.coords)
imp.coords(ii,1:3)=imp.coords(ii,1:3)*A';
end
[ArraySpeaker] = CreateSpeakerStructure( imp.coords(ct.pos,1), imp.coords(ct.pos,2), imp.coords(ct.pos,3), imp.coords(ct.pos,4) );
clear imp;



%% Calcul des distances micros-hp

D.R_mat=sqrt(sum(cat(3,bsxfun(@minus,Antenna.coord_vect(1,:),ArraySpeaker.x),bsxfun(@minus,Antenna.coord_vect(2,:),ArraySpeaker.y),bsxfun(@minus,zeros(1,numel(Antenna.coord_vect(2,:))),ArraySpeaker.z)).^2,3));
%create a 3D array with 3 dimensions refering to  x,y,z then sum the square
%value in order to obtain the absolute value

%% Propagation des signaux, creation des outils
%order of the lagrange operator for fractional delay implementation
N.NOrder=128;

%calculation of the number of samples for delay
[D,System,N]=CreateFilterFracDelay(D,ArraySpeaker,ct,N);

%% Propagation des signaux, application depuis les signaux des hp
tic
sig_mic_mat  = DelayImplementation(System.ht,sig_hp_mat,1,N,D);
toc

%% Processing microphone data (obtain fft) for the total system
t.Fsweep_avg=1;
[ sig_mic_mat ,System,~,~,t] = FrfSystemFinal_simu(sig_mic_mat,sweep_signal_vect,System,N,t,ct);

end

