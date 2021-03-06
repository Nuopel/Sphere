% This program is used to calculate the Frequency Response Function  of all
% the speakers in a set up over all microphones.
% INPUT: The record of all the microphones measurements over all speakers for all the  split in several file if needed 
%        The calibration file .mat obtained from CalibrationPiston.m
% OUTPUT: The FRF  split in a desired number of files.
% Use: set the .wav input and output location and name line 22 and 29
%
% Auteur : Dupont Samuel
% Version : 2.0 Fevrier 2016
clear variables; close all;clc

ct.N_record=5;% number of  record
ct.N_hp_record=10;% number of speaker per record
ct.N_sweep_avg=10;% number of sweep
N.N_mic=50;

%% Extraction signal original
a=load('calib_memsbedev_mic_25-05.mat');
data.calib=a.calib/a.calib(1);clear a; % relative calibration

%%
    [data.EntrySig,ct.Fs_sca]=audioread('sweep_signal_mic_sph.wav');
     N.N_sweep=length(data.EntrySig);

    H_cal=zeros(2^16,N.N_mic,ct.N_hp_record*ct.N_record);

for ii=1:ct.N_record
    
    file=sprintf('TF_anecho_hp/TF_%i-%i.w64',(ii-1)*ct.N_hp_record+1,ii*ct.N_hp_record);disp(file) % select name of the file
    [data.OutSigUnsplit,ct.Fs_sca]=audioread(file); % extract mic signals
    
    
    data.OutSigUnsplit=bsxfun(@times,data.OutSigUnsplit,data.calib.'); % apply calibartion
    N.N_meas_tot=length(data.OutSigUnsplit);
   
    if mod(length(data.OutSigUnsplit),ct.N_hp_record)~=0
    data.OutSigUnsplit=[data.OutSigUnsplit; zeros(mod(length(data.OutSigUnsplit),ct.N_hp_record),N.N_mic)];
    disp('Missing sample in the record compared to the number of recorded speaker ')
    a=input('enter');
    end
    
    for jj=1:ct.N_hp_record
        N.N_sweep_hp=N.N_meas_tot/ct.N_hp_record;
        data.OutSigsplit=data.OutSigUnsplit(1+N.N_sweep_hp*(jj-1):N.N_sweep_hp*jj,:);
        [~,H_cal_int,~,~,t ] = FrfSystemFinal_data_p2(data.OutSigsplit,data.EntrySig,N,ct); 
        H_cal(:,:,(ii-1)*ct.N_hp_record+jj)=H_cal_int.h_sig;
        disp((ii-1)*ct.N_hp_record+jj)
    end
end 

%% Save FRF in the program location
n_save=4;% select number of files output
nvar=50/n_save;
spacing=ceil(linspace(1,N.N_mic+1,n_save+1));

for ii=1:50
test=H_cal(:,:,ii);
if ii>9
a=sprintf('TF_anecho_impulse_calibxMICxHP%i.mat',ii);disp(a);
else
a=sprintf('TF_anecho_impulse_calibxMICxHP0%i.mat',ii);disp(a);
end
save(a,'test')
end
