function [ sig_mic_mat ] = DelayImplementation2(ht,sig_hp_mat,D)
% DelayImplementation apply the the delay using the matrix ht to the signal sig_hp_mat

%  HT : Matrix of the fir filter delay.
%  SIG_HP_MAT : signal to convert.
%  D.maxi= full part of the delay


%% Verification size of entry signal
ResizeColumn(sig_hp_mat);

%  initialisation

% select fft or filter way
sig_mic_mat=circshift(filter(ht,1,sig_hp_mat),D.maxi);


    
end



