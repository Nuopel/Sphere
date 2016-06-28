clear variables; close all;clc

%Regroupe toute les donnees : source simule, source reelle, source virtuelle


%% Define antenna set up
ct.pas_m = 2.5*2.54e-2; % pas de l'antenne
N.nbrx_sca = 8; % nombre de micros par ligne
N.nbry_sca = 7; % nombre de micros par ligne
[ Antenna ] = AntennArray( ct.pas_m,N.nbrx_sca,N.nbry_sca);N.N_mic=N.nbrx_sca*N.nbry_sca;
% contruction antenne 2.5 po , 8 mic


%% Extraction signal original
[data.sweep,ct.Fs_sca]=audioread('../../Data2/Session1/sweep_signal_antenne.wav');
ct.N_sweep_avg=10;



%% Extraction of the Data source simulated
if isfield(data, 'sc_simu')==0
    [data.hp_simu,ct.Fs2_sca]=audioread('../../Data2/Session1/source_simu_ambi_135d.w64');%data.hp_simu=[data.hp_simu ;zeros(1,50)];
    [data.sc_simu_virt ,H_simu_virt,~,N,ct,~] = simu_impulse( data.hp_simu(1:length(data.hp_simu)/10,:), data.sweep(1:length(data.hp_simu)/10,:),Antenna,ct,N);
    %     H_simu_virt.h_sig=[H_simu_virt.h_sig(1:10000,:) ;zeros(length(H_simu_virt.h_sig)-10000,N.N_mic)  ];
end

%% Calibration relatve microphone
a=load('calib_56_mic.mat');
a.relativ=a.max./a.max(1);

%% Extraction  Data source ambisonics
[data.sc_virt, ct.Fs_sca ] = audioread('../../Data2/Session1/source_ambi_135d.w64');data.sc_virt=[data.sc_virt ;zeros(1,56)];
data.sc_virt = bsxfun(@times,data.sc_virt,a.relativ.');

[~,H_virt,N,~,t ] = FrfSystemFinal_data(data.sc_virt,data.sweep,N,ct);
% H_virt.h_sig=[H_virt.h_sig(1:10000,:) ;zeros(length(H_virt.h_sig)-10000,N.N_mic)  ];


%% Video first comp
% text.b='Source ambisonics';text.a='Source ambisonics simulated';
% MovieMaker_double(H_simu_virt.h_sig,H_virt.h_sig,Antenna,0.003*ct.Fs_sca,0.015*ct.Fs_sca,ct,text)

%% Plot temporel
% compare signal virtuel et signal simuler
% mise en evidence delay

t.h_sig=0:1/ct.Fs_sca:(N.N_sweep_avg-1)/ct.Fs_sca;

ct.micn=28;
comp3(H_simu_virt.h_sig(:,ct.micn),H_virt.h_sig(:,ct.micn),H_virt.h_sig(:,ct.micn),t.h_sig,1)


%% Normalisation exploitation donnees temporelles
[maxi.max_simu_virt, delay.simu_virt]=max(H_simu_virt.h_sig(20:end/2,:));
[maxi.max_virt, delay.virt]=max(H_virt.h_sig(10:end/2,:));

H_norm.h_simu_virt=bsxfun(@rdivide,H_simu_virt.h_sig,maxi.max_simu_virt);

H_norm.h_virt=bsxfun(@rdivide,H_virt.h_sig,maxi.max_virt);

%%calcul des delay entre les differentes configurations, mise en evidence
%%mauvais placements, delay general
for ii=1:N.nbry_sca
    delay.simu_virt_moy(ii,1)= mean(delay.simu_virt(1,(ii-1)*N.nbrx_sca+1:ii*N.nbrx_sca));
    delay.virt_moy(ii,1)= mean(delay.virt(1,(ii-1)*N.nbrx_sca+1:ii*N.nbrx_sca));
end

%%calcul delay moyen configuration, repropagation general
delay.vs=mean(delay.simu_virt_moy-delay.virt_moy);
%% Plot signaux normalise,
% mise en evidence ripple
% comp3( H_norm.h_simu_virt(:,ct.micn),H_norm.h_virt(:,ct.micn),H_norm.h_virt(:,ct.micn),t.h_sig,2)
plot_micro( delay.simu_virt,delay.virt,delay.virt,3)

%% Repropagation
% repropagation pour comparaison sig theorique/meas
[H_simu_virt.h_sig_repro]= Repropagation(delay.virt-delay.simu_virt-10,H_simu_virt.h_sig,N);
H_simu_virt.h_sig_fft_repro=fft(H_simu_virt.h_sig_repro(1:N.N_sweep_avg,:));

[~, delay.simu_repro]=max([zeros(90,56); H_simu_virt.h_sig_repro(100:end/2,:)]);
plot_micro(delay.simu_repro,delay.virt,delay.virt,4);figure
plot(H_simu_virt.h_sig_repro(:,1));hold on ;plot(H_virt.h_sig(:,1));xlim([0 300])

%% Mic 24-32 temp plot
figure(5);

subplot(121)
surf(H_simu_virt.h_sig_repro(65:1200,25:32));
shading interp
view([0 90]);xlim([1 8]);xlabel('Microphone');
set(gca,'XTick',1:2:8);set(gca,'XTickLabel',25:2:32);
title('Simulated ambisonics source')
cax=caxis

subplot(122)
surf(H_virt.h_sig(65:1200,25:32));
shading interp
view([0 90]);xlim([1 8]);xlabel('Microphone');
set(gca,'XTick',1:2:8);set(gca,'XTickLabel',25:2:32);
title('Ambisonics source')
caxis(cax)


%% opt  Erreur Reelle-virtuelle
text.freq_opt=200:20:5000;%[500 800 1200 1800 2500 3000];
texts.a='Simulated ambisonics source';texts.b='Measured ambisonics source';texts.c='Relative error';texts.freq=	t.Fsweep_avg;
amp=10.^((-20:0.1:20)./20);

for jj=1:length(text.freq_opt)
    
    for ii=1:length(amp)
        [~, pos_sca] = min(abs(t.Fsweep_avg-text.freq_opt(jj))); %index of closest value
        erreur.grid_vect_virto = (H_virt.h_sig_fft(pos_sca,:))';
        erreur.grid_vect_simu_virt_reproo = (real(H_simu_virt.h_sig_fft_repro(pos_sca,:)).*amp(ii)+imag(H_simu_virt.h_sig_fft_repro(pos_sca,:)).*1i*amp(ii))';
        erreur.opto(ii)=mean(sqrt(abs(erreur.grid_vect_simu_virt_reproo-erreur.grid_vect_virto).^2./abs(erreur.grid_vect_simu_virt_reproo).^2)*100);
    end
    [~ ,pos(jj) ]= min(erreur.opto );
    erreur.grid_vect_virt(:,jj) = H_virt.h_sig_fft(pos_sca,:)';
    erreur.grid_vect_simu_virt_repro_correct(:,jj) =(real(H_simu_virt.h_sig_fft_repro(pos_sca,:)).*amp(pos(jj))+imag(H_simu_virt.h_sig_fft_repro(pos_sca,:)).*1i*amp(pos(jj)))';
    erreur.grid_vect_simu_virt_repro(:,jj) = H_simu_virt.h_sig_fft_repro(pos_sca,:)';
    
    erreur.opt(:,jj)=sqrt(abs(erreur.grid_vect_simu_virt_repro_correct(:,jj)-erreur.grid_vect_virt(:,jj)).^2./abs(erreur.grid_vect_simu_virt_repro_correct(:,jj)).^2)*100;
    erreur.value_trans(:,jj)=sqrt(abs(erreur.grid_vect_simu_virt_repro(:,jj)-erreur.grid_vect_virt(:,jj)).^2./abs(erreur.grid_vect_simu_virt_repro(:,jj)).^2)*100;
clc;disp(jj)
end
% figure(1)
% plot(20*log10(amp),erreur.opt)
text.a='Simu ambisonics sc';text.b='Meas ambisonics';
% MovieMaker_double_erreur( real(erreur.grid_vect_simu_virt_repro),real(erreur.grid_vect_virt),erreur.value_trans,Antenna,1,length(text.freq) ,ct,text)
erreur.correct=amp(pos);
%% Erreur Reelle-virtuelle


%%
erreur.vectminima=round(linspace(1,length(text.freq),10));

texts.a='Simulated ambisonics source';texts.b='Measured ambisonics source';texts.c='Relative error';texts.freq=text.freq;

for ii=1:length(erreur.vectminima)
    ct.fig=gcf;
    plot_reel_virt_erreur_trans(real(erreur.grid_vect_simu_virt_repro(:,ii)),real(erreur.grid_vect_virt(:,ii)),erreur.value_trans(:,ii),Antenna,erreur.vectminima(ii),texts,ct.fig.Number+1);
    ct.fig=gcf;
    plot_reel_virt_erreur_trans(real(erreur.grid_vect_simu_virt_repro_correct(:,ii)),real(erreur.grid_vect_virt(:,ii)),erreur.opt(:,ii),Antenna,erreur.vectminima(ii),texts,ct.fig.Number+1);
    pause(0.5)
end

