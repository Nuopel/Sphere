%% Post-calibration using band-pass filter and RMS

%% 0) General info

% SPL of the calibrator and other data
SPL = 94;
sf = 48000;

% Filter parameter
low_pass_freq = 1200;
low_pass_order = 12;
high_pass_freq = 800;
high_pass_order = 8;

% Filter design
[bhi,ahi] = butter(high_pass_order,high_pass_freq/(sf/2),'high');
[blo,alo] = butter(low_pass_order,low_pass_freq/(sf/2),'low');

%% 1) Check for the files

M = 96;
prefix = 'raw_recording_no';
ext = '.wav';

% Get dir info and index

%% 2) Operate calibration

% Load the wave file
for mm = 1:M
    
    disp(['File:',num2str(mm),'/',num2str(M)]);
    
    if exist([prefix,num2str(mm),ext])
        
        
        disp(['loading:',prefix,num2str(mm),ext]);
        [data] = wavread([prefix,num2str(mm),ext]);
        
        % Some verifications
        [Pxx_data f]=pwelch(data,hanning(8192),4096,8192,sf);
        rms_straight(mm) = sqrt(mean(data.^2));
        figure(1);
        set(1,'color','w');
        hl = semilogy([1000 1000],[10^-20 1]);
        set(hl,'color','k');
        hold on;
        semilogy(f,Pxx_data,'r');
        xlabel('Freq. [Hz]');ylabel('PSD');
        axis tight;
        hold on;
        
        % Filtering
        xt = filtfilt(bhi,ahi,data);
        xt = filtfilt(blo,alo,xt);
        [Pxx_filtered f]=pwelch(xt,hanning(8192),4096,8192,sf);
        semilogy(f,Pxx_filtered,'g');
        rms_filtered(mm) = sqrt(mean(xt.^2));
        axis tight;
        pause(0.1)
        hold off;
        
        % Actual sensitivity
        sensitivity_buffer(mm) = (10^(SPL/20)*0.000020)/rms_filtered(mm);
        sensitivity_straight(mm) = (10^(SPL/20)*0.000020)/rms_straight(mm);
        
    else
        
        sensitivity_buffer(mm) = 0;
        sensitivity_straight(mm) = 0;
        
    end
    
end

% Display
figure;
bar(1:96,sensitivity_straight,0.8,'r');
hold on;
bar(1:96,sensitivity_buffer,0.4,'g');

%% 3) Store data

sensitivity = sensitivity_buffer;
mic_no = 1:M;

save sensitivity sensitivity mic_no
