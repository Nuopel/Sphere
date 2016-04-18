function [sig_mic_mat ,System,N,ct,t ] = FrfFarina_simu( sig_mic_mat,sweep_signal_vect,System,N,t,ct )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


sig_mic_mat=sig_mic_mat(:,1:N.N_sweep,:);
N.N_sweep_avg=N.N_sweep/ct.av_sca;
ct.dfe_sweep_avg=ct.Fs_sca/N.N_sweep_avg;
t.Fsweep_avg=0:ct.dfe_sweep_avg:(N.N_sweep_avg-1)*ct.dfe_sweep_avg;


parfor ii=1:ct.N_sweep_avg

%% Preporocessing and variable definition
f_vect=filter([1 -1]*ct.Fs_sca,1,flipud(sweep_signal_vect(1+N_av_sca*(ii-1):N_av_sca*ii)));%time reversal + 3dB compensation

%% Processing of the data
% Use fourier domain in order to do the convolution faster
h_fft_mat(:,:,ii)=bsxfun(@times,fft_norm(f_vect)',fft_norm(data_mat(1+N_av_sca*(ii-1):N_av_sca*ii)));

end
System.h_sig_fft=mean(h_fft_mat,3);
System.h_sig=real(ifft(System.h_sig_fft));

end
