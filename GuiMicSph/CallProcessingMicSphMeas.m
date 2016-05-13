function [ HData, t,var ] = CallProcessingMicSphMeas( data,ct  )
 

%% Define constants
ct.r_hp_sca = 1.07 ;%rayon de la sphere
ct.r_micsph = 0.02;
ct.hankel_order =2;
ct.M_th = 15;
ct.M=5;
ct.nbr_M_th=(ct.M_th+1).^1;
ct.Fs=48000;
ct.c_air=340;

var.m_vect=0:ct.M_th;
var.m_sum_vect=(var.m_vect+1).^2;
var.nbr_m=(2.*var.m_vect)+1;


% data.OutSig=[data.OutSig;zeros(2,50)];
% data.OutSig=data.OutSig(length(data.OutSig)/2+1:end,:);

N.N_sweep=length(data.EntrySig);
N.N_meas=length(data.OutSig);

%% Look if data size In equal data size Out
if N.N_meas<N.N_sweep
    data.OutSig=[data.OutSig ;zeros(N.N_sweep-N.N_meas,50)];
    fprintf('%i zeros have been added to the measurement data', N.N_sweep-N.N_meas)
else if N.N_meas>N.N_sweep
      data.EntrySig=[data.EntrySig ;zeros(N.N_sweep-N.N_meas,50)];
      fprintf('%i zeros have been added to the sweep data \n Be careful measurement longer than  entry signal \n', -N.N_sweep+N.N_meas)
    end
end

%% Calibration relative microphone
if isfield(data,'CalibMic') == 1
    data.CalibMic = data.CalibMic./data.CalibMic(1); % calibrate in function of mic 1
    data.OutSig = bsxfun(@times,data.OutSig,data.CalibMic.'); % Apply calibration
else
    disp('No calibration data')
end

%% Average data
ct.N_sweep_avg=10;
ct.N_mic=50;
[~,HData,~,~,t ] = FrfSystemFinal_data_gui(data.OutSig,data.EntrySig,N,ct);


end

