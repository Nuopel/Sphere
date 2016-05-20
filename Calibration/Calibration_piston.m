clear variables; close all;clc

%% Constantes
ct.N_mic=50;

%% Extraction signaux

for ii = 1:ct.N_mic
    file = sprintf('Audio %i.wav',ii) ;
    eval('[data.meas(:,ii), Fs] = audioread(file) ;') ;
end

%% Calibration fft
% N.N_calib=length(data.meas);
% [freq, data.meas_fft] = ISYfourier(data.meas, Fs);
% [data.max, pos] = max(abs(data.meas_fft(1:end/2,:)));
% data.max=data.max.';

%% Calibration filter + zeros tracking + fft
% for ii = 1:ct.N_mic
%     file = sprintf('Audio %i.wav',ii) ;
%     [data.x, Fs] = audioread(file) ;
%     
%     n=10;
%     fc=900;
%     wn=fc*2/Fs;
%     [b,a] = butter(n,wn);
%     data.x=filter(b,a,data.x);
%     
%     n=10;
%     fc=1100;
%     wn=fc*2/Fs;
%     [b,a] = butter(n,wn,'high');
%     data.x=filter(b,a,data.x);
%     
%     marge=48000;
%     [~, pos_sca] = min(abs(data.x(end-marge:end))); %index of closest value
%     [~, pos_sca2] = min(abs(data.x(1:marge))); %index of closest value
%     data.x=data.x(pos_sca2:end-marge+pos_sca-1);
%     [freq2, data.x_fft] = ISYfourier(data.x, Fs);
%     [data.xmax(ii,1), a] = max(abs(data.x_fft(1:round(end/2),:)));
%     pos2(ii,1)=freq(a);
% end

%%
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
for ii = 1:ct.N_mic
    ct.marge=48000;
    [~, ct.pos_sca] = min(abs(data.filt(end-ct.marge:end))); %index of closest value
    [~, ct.pos_sca2] = min(abs(data.filt(1:ct.marge))); %index of closest value
    data.filtz=data.filt(ct.pos_sca2:end-ct.marge+ct.pos_sca-1,ii);
    calib(ii,1)=rms(data.filtz);
%     plot(data.filtz);hold on;plot(data.meas(:,ii));hold off ;xlim([0 15000])
%     a=input('enter')
end
save('calib_memsbedev_mic_18-05.mat','calib','-v7.3')


