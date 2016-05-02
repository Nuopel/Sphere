clear variables; close all;clc

%% Constantes
ct.N_sweep_avg=10;
N.N_mic=56;

%% Extraction signal original

[data.micro,ct.Fs_sca]=audioread('../../Data2/Session1/source_cible_135d.w64');
[data.sweep,ct.Fs2_sca]=audioread('../../Data2/Session1/sweep_signal_antenne.wav');data.sweep=data.sweep(1:end-1);
N.N_sweep=length(data.sweep);

%% Calcul frf

if ct.Fs2_sca~=ct.Fs_sca
    disp('Error the sampling frequency is not the same for input and output');
    return;
else
    ct=rmfield(ct,'Fs2_sca');
end
[~,H_cible,~,~,t ] = FrfSystemFinal_data(data.micro,data.sweep,N,ct);

%% calcul delay

[maxi.max_cible, delay.cible]=max(H_cible.h_sig(30:700,:));maxi.max_cible=maxi.max_cible.';delay.cible=delay.cible.';
maxi.grid=reshape(delay.cible,[8 7]);
ct.c_air=0.381/(mean(abs(maxi.grid(:,1)-maxi.grid(:,end)))/ct.Fs_sca);

plot_micro( delay.cible,1)
d=linspace(3e-3,5e-3,50000);
for ii=1:length(d)
r_mic(:,ii)=ct.c_air*(delay.cible/ct.Fs_sca-d(ii));
    
end