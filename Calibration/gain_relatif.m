clear variables; close all;clc

%% Constantes
ct.N_sweep_avg=10;
N.N_mic=56;

%% define array
ct.pas_m = 2.5*2.54e-2; % pas de l'antenne
N.nbrx_sca = 8; % nombre de micros par ligne
N.nbry_sca = 7; % nombre de micros par ligne
offset=struct('x',-ct.pas_m/2,'y',0);
[Antenna] = AntennArray( ct.pas_m,N.nbrx_sca,N.nbry_sca,offset);N.N_mic=N.nbrx_sca*N.nbry_sca;

%% Extraction signal original

[data.micro,ct.Fs_sca]=audioread('../../Data2/Session1/source_cible_135d.w64');
[data.sweep,ct.Fs2_sca]=audioread('../../Data2/Session1/sweep_signal_antenne.wav');data.sweep=data.sweep(1:end-1);
N.N_sweep=length(data.sweep);

%% calibration

a=load('calib_56_mic.mat');
data.micro_calib=bsxfun(@rdivide,data.micro,a.max.');
data.micro_calib_sensi=bsxfun(@times,data.micro_calib(:,(1:8:49)+4),4*pi*(1.07-Antenna.R((1:8:49)+3).'.*[1 1 1 1 -1 -1 -1]));

%% Calcul frf

if ct.Fs2_sca~=ct.Fs_sca
    disp('Error the sampling frequency is not the same for input and output');
    return;
else
    ct=rmfield(ct,'Fs2_sca');
end
N.N_mic=7;
[~,H_cible,~,~,t ] = FrfSystemFinal_data(data.micro_calib_sensi,data.sweep,N,ct);

%%


semilogx(t.Fsweep_avg, 20*log10(mean(abs(H_cible.h_sig_fft(:,1:7)),2)/2e-5));

% %% calcul delay
% 
% [maxi.max_cible, delay.cible]=max(H_cible.h_sig(30:700,:));maxi.max_cible=maxi.max_cible.';delay.cible=delay.cible.';
% maxi.grid=reshape(delay.cible,[8 7]);
% ct.c_air=0.381/(mean(abs(maxi.grid(:,1)-maxi.grid(:,end)))/ct.Fs_sca);
% 
% 
% 
% for ii = 1:N.N_mic
%     [data.cor(:,ii), x(ii,:)]=xcorr(data.micro(:,ii),data.sweep);
% end
% 
% plot_micro( delay.cible,1)
% %% 
% select=7*2.5*2.54e-2;
% d = linspace(0,10.5e-3,5000) ;
% AC = zeros(length(d),3); 
% [maxi.max_cible_cor, delay.cible_cor]=max(data.cor);delay.cible_cor=x(1,delay.cible_cor);maxi.max_cible=maxi.max_cible.';delay.cible=delay.cible.';
% plot_micro(delay.cible_cor,2)
% 
% for ii = 1:length(d)
% r_mic(:,ii) = ct.c_air*(delay.cible/ct.Fs_sca-d(ii));
% [Distance, AC(ii,1)] = triangle_c( r_mic(1,ii), r_mic(5,ii),  r_mic(8,ii) );
% [Distance, AC(ii,2)] = triangle_c( r_mic(9,ii), r_mic(13,ii),  r_mic(16,ii) );
% [Distance, AC(ii,3)] = triangle_c( r_mic(17,ii), r_mic(21,ii),  r_mic(24,ii) );
% test(ii)=Distance.verif;
% end
% posd(1)=closest( select,AC(:,1) );d_f(1)=d(posd);
% posd(2)=closest( select,AC(:,2) );d_f(2)=d(posd);
% posd(3)=closest( select,AC(:,3) );d_f(3)=d(posd);
% 
% figure(3)
% plotyy(1:5000,real(AC),1:5000,r_mic(5,:))
% hold on
% plot([d(1) d(end)],[0.447 0.447],'g')
% plot([d(1) d(end)],[0.445 0.445],'r')
% plot([d(1) d(end)],[0.443 0.443],'g')
% 
% 
% figure(4)
% AC = zeros(length(d),1); 
% 
% for ii = 1:length(d)
% r_mic(:,ii) = ct.c_air*(delay.cible_cor/ct.Fs_sca-d(ii));
% [Distance, AC(ii,1)] = triangle_c( r_mic(2,ii), r_mic(5,ii),  r_mic(8,ii) );
% test(1,ii)=Distance.verif;
% [Distance, AC(ii,2)] = triangle_c( r_mic(10,ii), r_mic(13,ii),  r_mic(16,ii) );
% test(2,ii)=Distance.verif;
% [Distance, AC(ii,3)] = triangle_c( r_mic(18,ii), r_mic(21,ii),  r_mic(24,ii) );
% test(3,ii)=Distance.verif;
% [Distance, AC(ii,4)] = triangle_c( r_mic(26,ii), r_mic(29,ii),  r_mic(32,ii) );
% test(4,ii)=Distance.verif;
% [Distance, AC(ii,5)] = triangle_c( r_mic(34,ii), r_mic(37,ii),  r_mic(40,ii) );
% test(5,ii)=Distance.verif;
% end
% posd(1)=closest( select,AC(:,1) );d_f(1)=d(posd(1));
% posd(2)=closest( select,AC(:,2) );d_f(2)=d(posd(2));
% posd(3)=closest( select,AC(:,3) );d_f(3)=d(posd(3));
% posd(4)=closest( select,AC(:,4) );d_f(4)=d(posd(4));
% posd(5)=closest( select,AC(:,5) );d_f(5)=d(posd(5));
% 
% plotyy(d,real(AC),d,r_mic(5,:))
% hold on
% plot([d(1) d(end)],[0.447 0.447],'g')
% plot([d(1) d(end)],[0.445 0.445],'r')
% plot([d(1) d(end)],[0.443 0.443],'g')




