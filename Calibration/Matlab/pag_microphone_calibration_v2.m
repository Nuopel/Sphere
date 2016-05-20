% This script is a microphone array calibration script. It is based on
% simple command window interaction with the user. Adjust the relevant
% parameters in the "0) General parameters section".

% P.-A. Gauthier, september 2011

%% 0) General parameters

% Reset
clc
close all
clear

M = 96;% Number of microphones
Sf = 48000; % Sampling frequency
N = Sf; % Record duration in samples
SPL_calibration = 94; % SPL as procured by the pistophone
bandwidth_tolerance = 20; % Bandwidth threshold in Hz
bandwidth_amp = 0.1; % Amplitude for the defenition of the bandwidth ... 0.5 was not enough sensitive.
asiodeveiceid = 0;
set(0,'DefaultFigurePosition',[100 100 560 420]);

% Options
FFTFilterOption = 0;
RMSOption = 1;

close(findobj('type','figure','name','Sensitivity'));
close(findobj('type','figure','name','Oscilloscope and spectra'));
close(findobj('type','figure','name','Sensitivity table'));

%% 1) Setup the workspace, figures and help

% Preset the sensitivity to 0
sensitivity = zeros(M,1);
routing = (1:M)';
last_calibration = repmat('xx-xxx-xxxx xx:xx:xx',M,1);
mic_no = (1:M)';

% Preset the record buffer to 0
data = zeros(N,M);
data_fft = zeros(N,1);
sample_index = (1:N)-1;
f_index = linspace(0,Sf/2,N/2+1);
data_fft = 2*abs(data_fft(1:numel(f_index),:));

% Create the sensitivity bar graph
h_sensitivity = figure(1000);
set(h_sensitivity,'name','Sensitivity','windowstyle','docked','menubar','none','toolbar','none');
h_sensitivity_bar = bar(mic_no,sensitivity);
set(gca,'ylim',[0 40],'xlim',[1 96],'nextplot','replacechildren');
title('Microphone sensitivities [Pa/U]');
xlabel('Micro/chn no');ylabel('Sensitivities [Pa/U]');

% Create the oscilloscope and the analyzer
h_oscillo = figure(1001);
set(h_oscillo,'name','Oscilloscope and spectra','windowstyle','docked','menubar','none','toolbar','none');

% Oscillo
h_oscillo_axis = subplot(1,3,1);
set(gca,'nextplot','replacechildren','box','on');
set(gca,'xlim',[0 sample_index(end)],'ylim',[-0.5 0.5]);
grid on;
h_oscillo_plot = plot(sample_index,data);
title('Recording');
xlabel('Samples index');
ylabel('Amplitude');

% Analyser
h_fft_axis = subplot(1,3,2);
set(gca,'nextplot','replacechildren','box','on','yscale','log');
set(gca,'xlim',[0 f_index(end)],'ylim',[0.00001 1]);
grid on;
h_fft_plot = semilogy(f_index,abs(data_fft));
title('FFT');
xlabel('Frequencies [Hz]');
ylabel('Spectra');

% Analyser zoom
h_fftzoom_axis = subplot(1,3,3);
set(gca,'nextplot','replacechildren','box','on','yscale','log');
set(gca,'xlim',[950 1050],'ylim',[0.00001 1]);
grid on;
h_fftzoom_plot = semilogy(f_index,abs(data_fft));
title('FFT');
xlabel('Frequencies [Hz]');
ylabel('Spectra');

% Create the data display
h_sensitivity_table = figure(1002);
set(h_sensitivity_table,'name','Sensitivity table','windowstyle','docked','menubar','none','toolbar','none');
dblArray = javaArray ('java.lang.Object', 96, 4);
dblArray = pag_update_java_aray(M,sensitivity,mic_no,last_calibration,routing,h_sensitivity_table,dblArray);
pause(1);
t_sensitivity_table = uitable(h_sensitivity_table,dblArray);
set(t_sensitivity_table,'units','normalized','position',[0.05 0.05 0.9 0.9],'editable',0);
set(t_sensitivity_table,'data',dblArray);
set(t_sensitivity_table,'ColumnName',{'Chn no','Sensitivity','Last calibration','Routing'});

%% 2) Define or open an existing project, verification of existing data, etc.

project_mode = input('Start a new project [y or enter] or load an existing project [n]? ','s');
% project_mode = questdlg('Start a new project or load an existing project?','Calibration: Project','Start','Load','Load'); % OLD GUI

if strcmp(project_mode,'y') == 1 | isempty(project_mode) == 1
    
    % Look for the project directory
    project_directory = uigetdir(pwd, 'Select the project directory');
    
    % Make two virgin files
    sensitivity = zeros(M,1);
    last_calibration = repmat('xx-xxx-xxxx xx:xx:xx',M,1);
    routing = (1:M)';
    cd_temp = pwd;
    cd(project_directory); % Go to the project dir
    
    % Check if the project is really a new one!
    if exist('sensitivity.mat','file') == 2;
        user = questdlg('The sensitivity file is already existing! Overwrite?','Warning!','Yes','No','No');
    else
        user = 'Yes';
    end
    
    % Save
    if strcmp(user,'Yes') == 1
        save sensitivity sensitivity last_calibration mic_no
    end
    
    % Check if the project is really a new one!
    if exist('routing.mat','file') == 2;
        user = questdlg('The routing file is already existing! Overwrite?','Warning!','Yes','No','No');
    end
    
    % Save
    if strcmp(user,'Yes') == 1
        save routing routing mic_no
    end
    
    cd(cd_temp); % Back to the past directory
    
elseif strcmp(project_mode,'n') == 1
    
    % Look for the project directory
    project_directory = uigetdir(pwd, 'Select the project directory');
    
end

% Look for the calibration files and routing files
cd_temp = pwd;
cd(project_directory); % Go to the project dir
if exist('sensitivity.mat','file') == 2;
    load('sensitivity.mat');
else
    error('The project does not include the sensitivity.mat file. Please retry.')
end
if exist('routing.mat','file') == 2;
    load('routing.mat');
else
    error('The project does not include the routing.mat file. Please retry.')
end
cd(cd_temp); % Back to the past directory

% Update the figures and display
h_sensitivity_bar = pag_update_calibration_figures(sensitivity,mic_no,h_sensitivity);

% Update the java array
dblArray = pag_update_java_aray(M,sensitivity,mic_no,last_calibration,routing,h_sensitivity_table,dblArray);

% Update the table display
pag_update_table_aray(t_sensitivity_table,dblArray);

%% 3) Make the actual calibration loop

% Start?
user = input('Start [y or enter] or quit [n]?','s');
% user = questdlg('Start or quit?','Calibration: Start','Start','Quit','Start');
if strcmp(user,'y') == 1 | isempty(user) == 1
    
    measure_flag = 1;
    
    % While loop
    while measure_flag == 1
        
        % A new microphone
        disp('================NEW MICROPHONE================');
        
        % Make the recording
        tmp = zeros(N,1);
        data = pa_wavplayrecord(tmp,asiodeveiceid,Sf,N,1,M,asiodeveiceid,'asio');
        %         data = zeros(N,M); % TEMP DEBUG
        %         env = 0.1*randn(1,N/4800)+0.9; % TEMP DEBUG
        %         env = resample(env,4800,1); % TEMP DEBUG
        %         env = env(:,1:N); % TEMP DEBUG
        %         data(:,randint(1,1,M)+1) = env.*sin(2*pi*1000*linspace(0,1,N)); % TEMP DEBUG
        
        % Detect the microphone
        energy = sqrt(sum(data.^2,1))/N;
        [useless,chn_index] = max(energy);
        
        % Remove DC
        data(:,chn_index) = detrend(data(:,chn_index));
        
        % Display in the oscillo
        sfigure(h_oscillo);
        axes(h_oscillo_axis);
        h_oscillo_plot = plot(sample_index,data(:,chn_index));
        % axis tight;
        
        % Clear the tisplay the FFT
        sfigure(h_oscillo);
        axes(h_fft_axis);
        h_fft_plot = loglog(f_index,0*f_index);
        axes(h_fftzoom_axis);
        h_fft_plot = loglog(f_index,0*f_index);
        
        % Report to user
        user_detect = input(['Detected calibration chn: ',num2str(chn_index),', accept [y or enter] or force [n]? '],'s');
        % user_detect = questdlg(['Detected calibration chn: ',num2str(chn_index),', accept or force? '],'Calibration: Detection','Accept','Force','Accept');
        
        if strcmp(user_detect,'y') == 1 | isempty(user_detect) == 1 % Accept the detection
            
            % Nothing happen
            force_flag = 0; % Create aflag
            
        else % Force the channel
            
            % Make the change
            user_chn = str2double(input('Select the channel: ','s'));
            chn_index = user_chn;
            force_flag = 1; % Create a flag
            
        end
        
        % Make the actual verifications and computations
        data_fft = fft(data(:,chn_index));
        data_fft = 2*abs(data_fft(1:numel(f_index),:))/N;
        
        % Display the FFT
        sfigure(h_oscillo);
        axes(h_fft_axis);
        h_fft_plot = semilogy(f_index,data_fft);
        axes(h_fftzoom_axis);
        h_fft_plot = semilogy(f_index,data_fft);
        
        % Find the max
        [max_amp,max_f] = max(data_fft);
        disp(['Detected frequency on channel ',num2str(chn_index),': ',num2str(f_index(max_f)),' Hz']);
        if force_flag == 1 % Make some supplementary verifications
            if max_amp == 0 % Error
                warning(['The channel selected by the user is silent. We restart the calibration iteration for this channel.']);
                continue
            end
        end
        
        % Evaluate the bandwidth
        seek = 1;
        bin_offset = 1;
        delta_f = f_index(2) - f_index(1);
        while seek == 1
            
            fft_amp_left = data_fft(max_f-bin_offset);
            fft_amp_right = data_fft(max_f+bin_offset);
            fft_amp_mean_lr = (fft_amp_left + fft_amp_right)/2;
            
            if le(fft_amp_mean_lr,bandwidth_amp*max_amp) == 1 % Stop
                
                seek = 0; % Stop
                bandwidth = 2*bin_offset*delta_f;
                disp(['Detected bandwidth on channel ',num2str(chn_index),': ',num2str(bandwidth),' Hz']);
                
            else % continue
                
                bin_offset = bin_offset + 1;
                
            end
            
        end
        
        % Depending on options

        if RMSOption ~= 1
            
            % Straight max value
            U = max_amp;
            
        elseif RMSOption == 1
            
            if FFTFilterOption ~= 1
            
                U = sqrt(mean(data(:,chn_index).^2));
                
            elseif FFTFilterOption == 1
                
                data_fft_raw = fft(data(:,chn_index));
                mod_fft = 0*data_fft_raw;
                mod_fft(max_f + [-bin_offset:1:bin_offset]) = data_fft_raw(max_f + [-bin_offset:1:bin_offset]);
                data_filtered = ifft(mod_fft,[],'symmetric');
                U = sqrt(mean(data_filtered.^2));
                
            end
            
        end
        
        % Evaluate sensitivity
        sensitivity_buffer = (10^(SPL_calibration/20)*0.000020)/U;
        
        % Accept or reject automatically, according to bandwidth
        disp(['Detected sensibility on channel ',num2str(chn_index),': ',num2str(sensitivity_buffer),' Pa/U']);
        if ge(bandwidth,bandwidth_tolerance) == 1
            
            disp(['This is above the bandwidth threshold. Please restart: press return.']);
            pause
            continue
            
        end
        
        % Ask the user
        user_accept = input(['Accept [y or enter] or reject [n] the calibration (accepting will save the file and update the displays)? '],'s');
        % user_accept = questdlg(['Accept or reject the calibration (accepting will save the file and update the displays)? '],'Calibration: Accept','Accept','Reject','Accept');
        if strcmp(user_accept,'y') == 1 | isempty(user_accept) == 1 % Accept the calibration
            
            % Pass the value to the array
            sensitivity(chn_index) = sensitivity_buffer;
            
            % Update the bar graph
            h_sensitivity_bar = pag_update_calibration_figures(sensitivity,mic_no,h_sensitivity);
            
            % Update the java array
            last_calibration(chn_index,:) = datestr(clock); % Save current time
            dblArray = pag_update_java_aray(M,sensitivity,mic_no,last_calibration,routing,h_sensitivity_table,dblArray);
            
            % Update the table display
            pag_update_table_aray(t_sensitivity_table,dblArray);
            
            % Save the sensitivity file
            cd_temp = pwd;
            cd(project_directory); % Go to the project dir
            save sensitivity sensitivity last_calibration mic_no % Save calibration
            wavwrite(data(:,chn_index),Sf,16,['raw_recording_no',num2str(chn_index)]); % Save raw data
            cd(cd_temp);
            
        else
            
            disp(['Calibration rejected by the user. Please restart.']);
            continue
            
        end
        
        % Next microphone?
        user_next = input(['Next microphone [y or enter] or stop the calibration [n] (The sensitivity will all be saved) ?'],'s');
        if strcmp(user_next,'y') == 1 | isempty(user_next) == 1 % End the calibration
            measure_flag = 1;
        else
            measure_flag = 0;
        end
        
    end
    
else
    
    return
    
end

%% 4) Close the project

% Back to the direction
cd(cd_temp);

disp('================END================')