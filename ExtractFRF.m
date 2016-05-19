clear variables; close all;clc

ct.N_record=5;% number of speaker recorded
ct.N_sweep_avg=10;% number of sweep
N.N_mic=50;

%% Extraction signal original
[data.EntrySig,ct.Fs_sca]=audioread('../sweep_signal_mic_sph.wav');
N.N_sweep=length(data.EntrySig);

%%

for ii=1:ct.N_record
 
    file=sprintf('TF_%i-%i.w64',(ii-1)*ct.N_record+1,ii*ct.N_record);
    [data.OutSigUnsplit,ct.Fs_sca]=audioread(file);
    N.N_meas_tot=length(data.OutSigUnsplit);
    data.OutSigUnsplit=[data.OutSigUnsplit; zeros(mod(length(data.OutSigUnsplit),ct.N_record),N.N_mic)];
    
    N.N_meas=length(data.OutSigUnsplit);
    for jj=1:ct.N_record
        N.N_sweep_avg=floor(N.N_meas_tot/ct.N_record);
        data.sweep_calibration=data.sweep_calibration_all(1+N.N_sweep_avg*(ii-1):N.N_sweep_avg*ii,:);
        [~,H_cal_int,~,~,t ] = FrfSystemFinal_data(data.sweep_calibration,data.sweep,N,ct); 
        H_cal(:,:,(ii-1)*ct.N_record+jj)=H_cal_int.h_sig_fft;
        clc;disp((ii-1)*ct.N_record+jj)
    end
   data=rmfield(data,'sweep_calibration_all');
   clear H_cal_int;
end
% save('calib_frf_short.mat','H','-v7.3')

