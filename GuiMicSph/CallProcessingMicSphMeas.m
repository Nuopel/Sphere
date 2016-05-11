function [ Pressure ] = CallProcessingMicSph( data,Antenna,ct  )


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

N.N_sweep=length(data.EntrySig);
N.N_meas=length(data.OutSig);

%% Look if data size In equal data size Out
if N.N_meas<N.N_sweep
    data.meas=[data.meas ;zeros(N.N_sweep-N.N_meas,56)];
    fprintf('%i zeros have been added to the measurement data', N.N_sweep-N.N_meas)
else if N.N_meas>N.N_sweep
      data.meas=[data.meas ;zeros(N.N_sweep-N.N_meas,56)];
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

%% Spherical wave target
[ source.x, source.y, source.z ] = cart2sph( ct.Theta, ct.Phi, ct.R ) ;
[Pressure.monopole_exp ] = monopole_pressure(ct.k,source,Antenna);


% %% Define vector
% handles.data.EntrySig
% 
% %% Encoding from microphone pressure
% Bmn.recons = Bmn_encoding_sph( Pressure,Sphmic,ct,N,var );
% 
% 
% %%_________________________________________________________________________
% 
% %% _______________________ Affichage data__________________________________
% %%_________________________________________________________________________
% 
% %% Operation on variable
% var.k=ct.k;
% N.N_sweep=1;

end

