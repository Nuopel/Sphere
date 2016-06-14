function [ HData,t, k ] = SphmicFfrDataMeas( OutSig,EntrySig,System)

v2struct(System)
%% Define Spherical microphone set up (same as speaker but smaller)
N.N_sweep=length(EntrySig);
N.N_meas=length(OutSig);
%% Look if data size In equal data size Out
if N.N_meas<N.N_sweep
    data.OutSig=[OutSig ;zeros(N.N_sweep-N.N_meas,50)];
    fprintf('%i zeros have been added to the measurement data \n', N.N_sweep-N.N_meas)
else if N.N_meas>N.N_sweep
      data.EntrySig=[EntrySig ;zeros(N.N_sweep-N.N_meas,50)];
      fprintf('%i zeros have been added to the sweep data \n Be careful measurement longer than  entry signal \n', -N.N_sweep+N.N_meas)
    end
end

%% Average data
% Select frequency
[~,HData,~,~,t ] = FrfSystemFinal_data(OutSig,EntrySig,N,ct);
ct.k=2*pi*t.Fsweep_avg/ct.c_air;
k=ct.k;


end

