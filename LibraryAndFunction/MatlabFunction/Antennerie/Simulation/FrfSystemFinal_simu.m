function [ sig_mic_mat ,System,N,ct,t] = FrfSystemFinal_simu( sig_mic_mat,sweep_signal_vect,System,N,t,ct )
% Do the Frequency response using the method H=out/in for simulation

%% Fft require matrix composed of column vector, check the dimension of the
% data
[a, b ]=size(sig_mic_mat);
if b>a
    sprintf('The matrix is transposed for fft used which work only with column of data')
    sig_mic_mat=permute(sig_mic_mat,[2 1 3]);
end

[a2, b2 ]=size(sweep_signal_vect);
if b2>a2
sweep_signal_vect=permute(sweep_signal_vect,[2 1 3]);
end

%% Convert the mic (out) signal to the same size than the sweep (in) 
if N.N_sweep<N.N_sca
sig_mic_mat=sig_mic_mat(1:N.N_sweep,:,:);
N.N_sca=N.N_sweep;
else
sweep_signal_vect=sweep_signal_vect(1:N.N_sca,:,:);
N.N_sweep=N.N_sca;
end
%% Some constant definition and useful vector reprensatation
ct.N_sweep_avg=1;
N.N_sweep_avg=floor(N.N_sweep/ct.N_sweep_avg);%  number of point inside the average
ct.dfe_sweep_avg=ct.Fs_sca/N.N_sweep_avg; % frequency padding
t.Fsweep_avg=0:ct.dfe_sweep_avg:(N.N_sweep_avg-1)*ct.dfe_sweep_avg; % Frequency vector for fft from sweep averaged



%% filter
n=6;
fc=2000;
wn=fc*2/ct.Fs_sca;
[b,a] = butter(n,wn);
sig_mic_mat=filter(b,a,sig_mic_mat);
%% Processing of the impulse response and FRF

h_fft_mat=zeros(N.N_sweep_avg,N.N_mic,ct.N_sweep_avg);% initialisation

var=N.N_sweep_avg;%avoid "over communication variable"

% Divide each colomn of fourier vector out by fourier vector in
for ii=1:1 % simulated data 1 avg enought ct.N_sweep_avg
    h_fft_mat(:,:,ii)=bsxfun(@rdivide,fft_norm(sig_mic_mat(1+var*(ii-1):var*ii,:)),fft_norm(sweep_signal_vect(1+var*(ii-1):var*ii,:)));
end

System.h_sig_fft=mean(h_fft_mat,3);% FRF
System.h_sig=real(ifft(System.h_sig_fft)); % Impulse response


%% Revert the matrix if needed to be in the shape of the input 
if b>a
    sig_mic_mat=permute(sig_mic_mat,[2 1 3]);
    System.h_sig_fft=permute(System.h_sig_fft,[2 1 3]);
    System.h_sig=permute(System.h_sig,[2 1 3]);
end


end

