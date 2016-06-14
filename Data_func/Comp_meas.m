clear all;%close all;clc


%% 

[data.EntrySig,ct.Fs]=audioread('sweep_signal_mic_sph.wav');


[data.OutSig_target,ct.Fs_sca]=audioread('sph_hp3_calib.wav');
[sinus.OutSig_target,~]=audioread('sph_hp3_sinus.w64');

[data.OutSig,~]=audioread('sph_3_ambisonique_decode_encoder(18db)_decoder(26db)_calib.wav');
[sinus.OutSig,~]=audioread('sph_3_ambisonique_decode_encoder(18db)_decoder(26db)_sinus.w64');
[Playback.OutSig,~]=audioread('sph_3_ambisonique_decode_encoder(18db)_decoder(26db)_playback.w64');

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

%% Find sinus value
sinus.OutSig_target=find(abs(sinus.OutSig_target)>0.5e-4);Delay.OutSig_target=sinus.OutSig_target(1);
sinus.OutSig=find(abs(sinus.OutSig)>0.5e-4);Delay.OutSig_sinus=sinus.OutSig(1);
clear sinus

%% Find playback value 
Delay.OutSig_playback=find(abs(Playback.OutSig(:,1))>0.5e-4);Delay.OutSig_playback=Delay.OutSig_playback(1);
clear Playback



%% Calculate frf
[ Sphmic, ct.N_mic ] = CreateSpeakerSystem(ct.r_micsph);% create the sphere set up, sort in struc Array
System=struct('ct',ct,'var',var,'N',N,'Sphmic',Sphmic);
[ Frf.target,t, var.k] = SphmicFfrDataMeas( data.OutSig_target,data.EntrySig,System);
[ Frf.comp,~, ~] = SphmicFfrDataMeas( data.OutSig,data.EntrySig,System);


%% _______________________ Affichage data__________________________________
%%_________________________________________________________________________

Delay.tot=Delay.OutSig_playback+Delay.OutSig_sinus;
%% Operation on variable
var.select_freq = logspace(2,log10(5000),150);
for ii=1:length(var.select_freq)
Delay.exp=exp(1i*Delay.tot/ct.Fs*2*pi*var.select_freq(ii)/ct.c_air);
ct.N_mic=50;
ct.pos = closest(var.select_freq(ii)*2*pi/ct.c_air,var.k);ct.k = var.k(ct.pos) ; 
ct.M=5;




Bmn.recons_target = Bmn_encoding_sph(Frf.target.h_sig_fft(ct.pos,:),Sphmic,ct,var,1 );
Bmn.recons = Bmn_encoding_sph(Frf.comp.h_sig_fft(ct.pos,:),Sphmic,ct,var,1 );
Bmn.recons_delay = Bmn_encoding_sph(Frf.comp.h_sig_fft(ct.pos,:)*Delay.exp,Sphmic,ct,var,1 );

% moindre carre
N.N_sweep=1;

figure
subplot(221)
[Pressure.p_recons_target, ~ ] = Pressure_map_SphMic(ct.M,Bmn.recons_target,ct,N,var);
title('target')
subplot(222)
[Pressure.p_recons, ~ ] = Pressure_map_SphMic(ct.M,Bmn.recons,ct,N,var);
% subplot(222)
[Pressure.p_recons_delay, ~ ] = Pressure_map_SphMic(ct.M,Bmn.recons_delay,ct,N,var);
title('recons delay')


subplot(223)
H = fspecial('gaussian',40,10);
a(ii)=H(:).'.*Pressure.p_recons_target*pinv(Pressure.p_recons.*H(:).');
Pressure.test=a(ii)*Pressure.p_recons;
title('recons inv')

% figure

Pressure_map_(Pressure.test,ct,var,0);

subplot(224)
title('erreur target inv')
[erreur,~] = erreur_n( Pressure.p_recons_target, Pressure.test );
Pressure_map_(erreur,ct,var,1);
subtitle(sprintf('%i',var.select_freq(ii)))


end

