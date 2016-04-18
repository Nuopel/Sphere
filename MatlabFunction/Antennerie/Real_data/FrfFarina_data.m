function [data_mat ,System,N,ct,t ] = FrfFarina_data( data_mat,sweep_signal_vect,N,ct)

% FrfFarina_data do the Frequency response using Farina sweep method 
% (error in the scale due the the filter +3db compensation)     

%% Fft require matrix composed of column vector, check the dimension of the
% data

[a, b ]=size(data_mat);
if b>a
    sprintf('The matrix is transposed for fft used which work only with column of data')
    data_mat=permute(data_mat,[2 1 3]);
end

[a2, b2 ]=size(sweep_signal_vect);
if b2>a2
    sweep_signal_vect=permute(sweep_signal_vect,[2 1 3]);
end


%% Convert the mic (out) signal to the same size than the sweep (in) 

data_mat=data_mat(1:N.N_sweep,:,:);


%% Some constant definition and useful vector reprensatation

N.N_sweep_avg=N.N_sweep/ct.N_sweep_avg;
ct.dfe_sweep_avg=ct.Fs_sca/N.N_sweep_avg;
t.Fsweep_avg=0:ct.dfe_sweep_avg:(N.N_sweep_avg-1)*ct.dfe_sweep_avg;

%% Processing of the impulse response and FRF

h_fft_mat=zeros(N.N_sweep_avg,N.N_mic,ct.N_sweep_avg);% initialisation
for ii=1:ct.N_sweep_avg
    
    % Preprocessing and variable definition
    f_vect=filter([1 -1]*ct.Fs_sca,1,flipud(sweep_signal_vect(1+N.N_sweep_avg*(ii-1):N.N_sweep_avg*ii)));%time reversal + 3dB compensation
    
    % Processing of the data
    % Use fourier domain in order to do the convolution faster
    h_fft_mat(:,:,ii)=bsxfun(@times,fft_norm(data_mat(1+N.N_sweep_avg*(ii-1):N.N_sweep_avg*ii,:)),fft_norm(f_vect));
    
end
System.h_sig_fft=mean(h_fft_mat,3);% FRF
System.h_sig=real(ifft(System.h_sig_fft)); % Impulse response


%% Revert the matrix if needed to be in the shape of the input 

if b>a
    data_mat=permute(data_mat,[2 1 3]);
    System.h_sig_fft=permute(System.h_sig_fft,[2 1 3]);
    System.h_sig=permute(System.h_sig,[2 1 3]);
end
end
