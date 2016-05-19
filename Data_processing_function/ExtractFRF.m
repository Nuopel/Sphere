clear variables; close all;clc

ct.N_record=10;% number of  record
ct.N_hp_record=5;% number of speaker per record
ct.N_sweep_avg=10;% number of sweep
N.N_mic=50;

%% Extraction signal original
[data.EntrySig,ct.Fs_sca]=audioread('../sweep_signal_mic_sph.wav');
a=load('calib_memsbedev_mic_18-05.mat');
data.calib=a.calib/a.calib(1);clear a; % relative calibration
N.N_sweep=length(data.EntrySig);

%%

for ii=1:ct.N_record
 
    clc;file=sprintf('TF_%i-%i.w64',(ii-1)*ct.N_hp_record+1,ii*ct.N_hp_record) % select name of the file
    [data.OutSigUnsplit,ct.Fs_sca]=audioread(file); % extract mic signals
    data.OutSigUnsplit=bsxfun(@times,data.OutSigUnsplit,data.calib.'); % apply calibartion
    N.N_meas_tot=length(data.OutSigUnsplit);
   
    if mod(length(data.OutSigUnsplit),ct.N_hp_record)~=0
    data.OutSigUnsplit=[data.OutSigUnsplit; zeros(mod(length(data.OutSigUnsplit),ct.N_hp_record),N.N_mic)];
    end
    
    N.N_meas=length(data.OutSigUnsplit);
    for jj=1:ct.N_hp_record
        N.N_sweep_hp=floor(N.N_meas_tot/ct.N_hp_record);
        data.OutSigsplit=data.OutSigUnsplit(1+N.N_sweep_hp*(jj-1):N.N_sweep_hp*jj,:);
        
        [~,H_cal_int,~,~,t ] = FrfSystemFinal_data(data.OutSigsplit,data.EntrySig,N,ct); 
        H_cal(:,:,(ii-1)*ct.N_hp_record+jj)=H_cal_int.h_sig_fft;
        disp((ii-1)*ct.N_hp_record+jj)
    end
   data=rmfield(data,'OutSigsplit');
   clear H_cal_int;
end

%% save
n_save=4;
nvar=50/n_save;
spacing=ceil(linspace(1,N.N_mic+1,n_save+1));
for ii=1:n_save
test=H_cal(:,:,spacing(ii):spacing(ii+1)-1);
a=sprintf('TF_WFS_FRF_calibxMICxHP%i-%i.mat',spacing(ii),spacing(ii+1)-1)
save(a,'test','-v7.3')
end