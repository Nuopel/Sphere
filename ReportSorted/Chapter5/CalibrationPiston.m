%% CalibrationPiston.M Take the measured pistonphone .wav and obtain calibration file
%
% INPUT: ct.N_mic wav files (number of microphones). They must have the
%        same name followed by their number or the program need to be
%        adapted.
%
% OUTPUT: .mat file containg the rms value of all the microphones. set the
%         name in the last line

% The program take the input apply a bandpass 600 Hz-1300 Hz and extract
% the rms value

clear variables; close all;clc

%% Constantes
ct.N_mic=50;

%% Extraction signaux

for ii = 1:ct.N_mic
    file = sprintf('Calibration/Audio %i.wav',ii) ;
    eval('[data.meas(:,ii), Fs] = audioread(file) ;') ;
end


%% filter (antialiasing)

n=10;
fc=1300;
wn=fc*2/Fs;
[b,a] = butter(n,wn);
data.filt=filter(b,a,data.meas);

n=10;
fc=600;
wn=fc*2/Fs;
[b,a] = butter(n,wn,'high');
data.filt=filter(b,a,data.filt);
data.filt=data.filt(5000:end,:);
clear a b fc n wn

%% Extract rms value from a zero value to an other zero value (full period sinus)
for ii = 1:ct.N_mic
    ct.marge=48000;
    [~, ct.pos_sca] = min(abs(data.filt(end-ct.marge:end))); %index of closest value
    [~, ct.pos_sca2] = min(abs(data.filt(1:ct.marge))); %index of closest value
    data.filtz=data.filt(ct.pos_sca2:end-ct.marge+ct.pos_sca-1,ii);
    calib(ii,1)=rms(data.filtz);
%  plot(data.meas(:,ii));hold on; plot(data.filtz);hold off ;xlim([0 15000]); title(sprintf('fig %i',ii))
%     a=input('enter');
% % uncomment if you want to see signals
end
save('calib_Micsph_02-06.mat','calib','-v7.3')


