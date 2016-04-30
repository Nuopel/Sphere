clear variables; close all;clc

%% Constantes
ct.N_mic=56;

%% Extraction signaux

for ii = 1:ct.N_mic
    file = sprintf('../../Data2/Calibration/Audio %i-1.wav',ii) ;
    eval('[data.meas(:,ii), Fs] = audioread(file) ;') ;
end

[freq data.meas_fft] = ISYfourier(data.meas, Fs);
data.max=max(abs(data.meas_fft(1:end/2,:))).';

save('calib_56_mic.mat','-struct','data', 'max','-v7.3')