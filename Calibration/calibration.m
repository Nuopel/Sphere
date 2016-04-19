clear variables; close all;clc


%% Extraction signal original
[data.sweep_delay,ct.Fs_sca]=audioread('../Data2/sweep_signal_delaydirect.wav');
[data.sweep_meas,ct.Fs_sca]=audioread('../Data2/Meas_delaydirect.w64');data.sweep_meas=[data.sweep_meas; 0];
N.N_mic=1;N.N_sweep=length(data.sweep_delay);ct.N_sweep_avg=10;
[~,H_delay,~,~,~ ] = FrfSystemFinal_data(data.sweep_meas,data.sweep_delay,N,ct);
% figure(1); plot(H_delay.h_sig);xlim([0 150])

%%


[data.piston1,ct.Fs_sca]=audioread('../Data2/calibration_1.w64');ct.p1=rms(data.piston1);
[data.piston2,ct.Fs_sca]=audioread('../Data2/calibration_52.w64');ct.p2=rms(data.piston2);
figure(2);plot(data.piston1);hold on;plot(data.piston2);

[data.sweep,ct.Fs_sca]=audioread('../Data2/sweep_signal_antenne_calibration.wav');

ct.N_record=10;ct.N_sweep_avg=4;N.N_mic=56;

for ii=1:5
    file=sprintf('../Data2/calibration_%i_%i.w64',(ii-1)*10+1,ii*10);
    [data.sweep_calibration_all,ct.Fs_sca]=audioread(file);
    data.sweep_calibration_all=[data.sweep_calibration_all; zeros(mod(length(data.sweep_calibration_all),ct.N_record),56)];
    N.N_sweep=length(data.sweep_calibration_all);
    for jj=1:ct.N_record
        N.N_sweep_avg=floor(N.N_sweep/ct.N_record);
        data.sweep_calibration=data.sweep_calibration_all(1+N.N_sweep_avg*(ii-1):N.N_sweep_avg*ii,:);
        [~,H_cal_int,~,~,t ] = FrfSystemFinal_data(data.sweep_calibration,data.sweep,N,ct); 
        H_cal(:,:,(ii-1)*10+jj)=H_cal_int.h_sig_fft;
        clc;disp((ii-1)*10+jj)
    end
   data=rmfield(data,'sweep_calibration_all');
   clear H_cal_int;
end
% save('calib_frf_short.mat','H','-v7.3')
%%
H_cal_m=smooth(H_cal_m,'rloess');

s1=80;s2=1000;
[~, pos_sca1] = min(abs(t.Fsweep_avg-s1)); %index of closest value
[~, pos_sca2] = min(abs(t.Fsweep_avg-s2)); %index of closest value
