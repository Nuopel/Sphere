clear variables; close all;clc

%% Constantes
ct.N_mic=50;

%% Extraction signaux

for ii = 1:ct.N_mic
    file = sprintf('../../CalibrationMicSph/Audio %i.wav',ii) ;
    
    eval('[data.meas(:,ii), Fs] = audioread(file) ;') ;
end
[freq, data.meas_fft] = ISYfourier(data.meas, Fs);
[data.max, pos] = max(abs(data.meas_fft(1:end/2,:)));
data.max=data.max.';

% [data.meas100, Fs] = audioread('../../DataMicSph/calib_sinus_100.w64') ;
% [data.meas150, Fs] = audioread('../../DataMicSph/calib_sinus_150.w64') ;

% [freq data.meas_fft_100] = ISYfourier(data.meas100, Fs);
% data.max_100=max(abs(data.meas_fft_100(1:end/2,:))).';
% 
% [freq data.meas_fft_150] = ISYfourier(data.meas150, Fs);
% data.max_150=max(abs(data.meas_fft_150(1:end/2,:))).';

% save('calib_56_mic.mat','-struct','data', 'max','-v7.3')
save('calib_memsbedev_mic.mat','-struct','data', 'max','-v7.3')