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
    
    
    ct.pos=12;
    [data.hp_simu_reel,ct.Fs2_sca]=audioread('../../Data2/Session1/source_simu_cible_135d.w64');%data.hp_simu_reel=[data.hp_simu_reel ;zeros(1,1)];
    [data.sc_simu_reel ,H_simu_reel,ArraySpeaker,N,ct,t] = simu_source_reelle( data.hp_simu_reel(1:length(data.hp_simu)/10,:), data.sweep(1:length(data.hp_simu)/10,:),Antenna,ct,N);

end

if ct.Fs2_sca~=ct.Fs_sca
    disp('Error the sampling frequency is not the same for input and output');
    return;
else
    ct=rmfield(ct,'Fs2_sca');
end


%% Calibration relatve microphone
a=load('calib_56_mic.mat');
a.relativ=a.max./a.max(1);

%% Extraction  Data source virtuelle, reelle
[data.sc_virt, ct.Fs_sca ] =audioread('../../Data2/Session1/source_ambi_135d.w64');data.sc_virt=[data.sc_virt ;zeros(1,56)];
[data.sc_reel, ct.Fs2_sca ] =audioread('../../Data2/Session1/source_cible_135d.w64');data.sc_reel=[data.sc_reel ;zeros(1,56)];
data.sc_virt = bsxfun(@times,data.sc_virt,a.relativ.');
data.sc_reel = bsxfun(@times,data.sc_reel,a.relativ.');


if ct.Fs2_sca~=ct.Fs_sca
    disp('Error the sampling frequency is not the same for input and output');
    return;
else
    ct=rmfield(ct,'Fs2_sca');
end


[~,H_virt,~,~,~ ] = FrfSystemFinal_data(data.sc_virt,data.sweep,N,ct);
[~,H_reel,N,ct,t ] = FrfSystemFinal_data(data.sc_reel,data.sweep,N,ct);

H_virt.h_sig=[H_virt.h_sig(1:10000,:) ;zeros(length(H_virt.h_sig)-10000,N.N_mic)  ];
H_reel.h_sig=[H_reel.h_sig(1:10000,:) ;zeros(length(H_reel.h_sig)-10000,N.N_mic)  ];


%% Video first comp
% text.b='Source virtuelle';text.a='Source reelle';
%  MovieMaker_double(H_reel.h_sig,H_virt.h_sig,Antenna,0.003*ct.Fs_sca,0.015*ct.Fs_sca,ct,text)

%% Plot temporel
% compare signal virtuel et signal simuler
% mise en evidence delay

t.h_sig=0:1/ct.Fs_sca:(N.N_sweep_avg-1)/ct.Fs_sca;

ct.micn=28;
comp3(H_simu_virt.h_sig(:,ct.micn),H_reel.h_sig(:,ct.micn),H_virt.h_sig(:,ct.micn),t.h_sig,1)


%% Normalisation exploitation donnees temporelles
[maxi.max_simu_reel, delay.simu_reel]=max(H_simu_reel.h_sig);
[maxi.max_simu_virt, delay.simu_virt]=max(H_simu_virt.h_sig);

[maxi.max_reel, delay.reel]=max(H_reel.h_sig);
[maxi.max_virt, delay.virt]=max(H_virt.h_sig);

H_norm.h_simu_reel=bsxfun(@rdivide,H_simu_reel.h_sig,maxi.max_simu_reel);
H_norm.h_simu_virt=bsxfun(@rdivide,H_simu_virt.h_sig,maxi.max_simu_virt);

H_norm.h_reel=bsxfun(@rdivide,H_reel.h_sig,maxi.max_reel);
H_norm.h_virt=bsxfun(@rdivide,H_virt.h_sig,maxi.max_virt);

%%calcul des delay entre les differentes configurations, mise en evidence
%%mauvais placements, delay general
for ii=1:N.nbry_sca
    delay.simu_reel_moy(ii,1)= mean(delay.simu_reel(1,(ii-1)*N.nbrx_sca+1:ii*N.nbrx_sca));
    delay.simu_virt_moy(ii,1)= mean(delay.simu_virt(1,(ii-1)*N.nbrx_sca+1:ii*N.nbrx_sca));

    delay.reel_moy(ii,1)= mean(delay.reel(1,(ii-1)*N.nbrx_sca+1:ii*N.nbrx_sca));
    delay.virt_moy(ii,1)= mean(delay.virt(1,(ii-1)*N.nbrx_sca+1:ii*N.nbrx_sca));
end

%%calcul delay moyen configuration, repropagation general
delay.sr=mean(delay.simu_virt_moy-delay.reel_moy);
delay.rv=mean(delay.virt_moy-delay.reel_moy);
delay.vs=mean(delay.simu_virt_moy-delay.virt_moy);
%% Plot signaux normalise, 
% mise en evidence ripple 
comp3( H_norm.h_simu_virt(:,ct.micn),H_norm.h_reel(:,ct.micn),H_norm.h_virt(:,ct.micn),t.h_sig,2)
plot_micro( delay.simu_virt,delay.reel,delay.virt,3)

% %% Repropagation
% % repropagation pour comparaison sig theorique/meas 
% [H_simu_virt.h_sig]= Repropagation(delay.virt-delay.simu_virt,H_simu_virt.h_sig,N);
% H_simu_virt.h_sig_fft=fft(H_simu_virt.h_sig(1:N.N_sweep_avg,:));
% 
% [~, delay.simu_repro]=max(H_simu_virt.h_sig);
% plot_micro(delay.simu_repro,delay.reel,delay.virt,4)

%% Mic 24-32 temp plot
figure(5);
subplot(121)
surf(H_simu_reel.h_sig(10:1200,25:32));
cax=caxis;

shading interp
view([0 90]);xlim([1 8]);xlabel('Microphone');ylabel('Time in sample ');
set(gca,'XTick',1:2:8);set(gca,'XTickLabel',25:2:32);
title('Simulated target source')


subplot(122)
surf(H_simu_virt.h_sig(10:1200,25:32));
shading interp
view([0 90]);xlim([1 8]);xlabel('Microphone');
set(gca,'XTick',1:2:8);set(gca,'XTickLabel',25:2:32);
title('Simulated ambisonics source')
caxis(cax)

figure(6);
subplot(121)
surf(H_reel.h_sig(10:1200,25:32));
shading interp
title('target source')
view([0 90]);xlim([1 8]);xlabel('Microphone');ylabel('Time in sample ');
set(gca,'XTick',1:2:8);set(gca,'XTickLabel',25:2:32);
caxis(cax)

subplot(122)
surf(H_virt.h_sig(10:1200,25:32));
shading interp
view([0 90]);xlim([1 8]);xlabel('Microphone');
set(gca,'XTick',1:2:8);set(gca,'XTickLabel',25:2:32);
title('Ambisonics source')
caxis(cax)

%% Erreur Reelle-virtuelle 
text.freq=100:1:4000;
erreur.grid_vect_reel=zeros(N.N_mic,length(text.freq));
erreur.grid_vect_virt=erreur.grid_vect_reel;

for ii=1:length(text.freq)
    [~, pos_sca] = min(abs(t.Fsweep_avg-text.freq(ii))); %index of closest value
    
    erreur.grid_vect_reel(:,ii) = (H_reel.h_sig_fft(pos_sca,:))';
    erreur.grid_vect_virt(:,ii) = (H_virt.h_sig_fft(pos_sca,:))';
    
    erreur.grid_vect_simu_reel(:,ii) = (H_simu_reel.h_sig_fft(pos_sca,:))';
    erreur.grid_vect_simu_virt(:,ii) = (H_simu_virt.h_sig_fft(pos_sca,:))';
    
    % grid_mat_reel=reshape(grid_vect_reel,size(Antenna.X_mat));
    % grid_mat_virt=reshape(grid_vect_virt,size(Antenna.X_mat));
end
erreur.value=sqrt(abs(erreur.grid_vect_reel-erreur.grid_vect_virt).^2./abs(erreur.grid_vect_reel).^2)*100;
erreur.test=1./abs(erreur.grid_vect_reel).^2*100;
erreur.value_simu=sqrt(abs(erreur.grid_vect_simu_reel-erreur.grid_vect_simu_virt).^2./abs(erreur.grid_vect_simu_reel).^2)*100;

erreur.value_trans=sqrt(abs(erreur.grid_vect_simu_virt-erreur.grid_vect_virt).^2./abs(erreur.grid_vect_simu_virt).^2)*100;

% MovieMaker_double_erreur( real(erreur.grid_vect_reel),real(erreur.grid_vect_virt),erreur.value,Antenna,1,length(text.freq) ,ct,text)
% MovieMaker_double_erreur( real(erreur.grid_vect_simu_reel),real(erreur.grid_vect_simu_virt),erreur.value_simu,Antenna,800,length(text.freq) ,ct,text)
% text.a='Simu ambisonics sc';text.b='Meas target ';
% MovieMaker_double_erreur( real(erreur.grid_vect_simu_virt),real(erreur.grid_vect_virt),erreur.value_trans,Antenna,1,length(text.freq) ,ct,text)

%% Erreur analysis
% figure(7)
% erreur.value_avg=mean(erreur.value');
% for ii=1:1900
% plot(erreur.value(Antenna.index_sc,ii))
% title(text.freq(ii));
% xlabel('Microphone');ylabel('Error');
% ylim([0 500]);
% pause(0.1);
% end
%%
figure(8)
subplot(211)
a=mean(erreur.value);
a2=mean(erreur.value_simu);
semilogx([a',a2'])
axis([100 2000 0 100]);
xlabel('Frequency');ylabel('Mean error frequency');
legend('Meas','Simu')

subplot(212)
a=mean(erreur.value');
a2=mean(erreur.value_simu');
plot([a(Antenna.index_sc)' a2(Antenna.index_sc)'])
set(gca,'XTick',1:10:56);set(gca,'XTickLabel',Antenna.Rmicro_sc(1:10:end));
xlabel('Microphone distance from center');ylabel('Mean error mic');
legend('Meas','Simu')


figure(9)
subplot(311)
plot(t.Fsweep_avg,20*log10([H_reel.h_sig_fft(:,1) H_virt.h_sig_fft(:,1)]));xlim([20 500])
legend('Meas target','Meas ambisonics')
subplot(312)
plot(t.Fsweep_avg,angle([H_reel.h_sig_fft(:,1) H_virt.h_sig_fft(:,1)]));xlim([20 500])
subplot(313)
plot(erreur.value(1,:));xlim([20 500])


figure(10)
subplot(311)
plot(t.Fsweep_avg,20*log10([H_simu_reel.h_sig_fft(:,1) H_simu_virt.h_sig_fft(:,1)]));xlim([20 500])
legend('Simu target','Simu ambisonics')
subplot(312)
plot(t.Fsweep_avg,angle([H_simu_reel.h_sig_fft(:,1) H_simu_virt.h_sig_fft(:,1)]));xlim([20 500])
subplot(313)
plot(erreur.value_simu(1,:));xlim([20 500])


figure(11)
subplot(211)
plot(t.Fsweep_avg,[permute(mean(H_reel.h_sig_fft'),[2 1]) permute(mean(H_virt.h_sig_fft'),[2 1])]);
subplot(212)
plot(mean(erreur.value))
%%


% erreur.vectminima=[44,107,223,454,1043,1463,1702,1886];
% erreur.vectminima=[46,224,457,739,931,1372,2876]; sweep 3 avg 2
% erreur.vectminima=[48,225,421,740,1044,1307,2114];
erreur.vectminima=[100,256,617,1409,1821];

text.b='Measured target source';text.a='Measured ambisonics source';text.c='Relative error';
texts.b='Simulated target source';texts.a='Simulated ambisonics source';texts.c='Relative error';texts.freq=text.freq;

for ii=1:length(erreur.vectminima)
ct.fig=gcf;
plot_reel_virt_erreur(real(erreur.grid_vect_reel),real(erreur.grid_vect_virt),erreur.value,Antenna,erreur.vectminima(ii),text,ct.fig.Number+1);
plot_reel_virt_erreur(real(erreur.grid_vect_simu_reel),real(erreur.grid_vect_simu_virt),erreur.value_simu,Antenna,erreur.vectminima(ii),texts,ct.fig.Number+2);
pause(0.5)
end

