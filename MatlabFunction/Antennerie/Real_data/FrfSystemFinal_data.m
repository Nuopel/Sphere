function [ data_mat ,System,N,ct,t] = FrfSystemFinal_data( data_mat,sweep_signal_vect,N,ct )
% Do the Frequency response using the method H=out/in
%%TO DO implement the frequency cut to remove low frequency

% Do the Frequency response using the method H=out/in

%% Fft require matrix composed of column vector, check the dimension of the
% data and revert them if needed
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
N.N_sweep=length(sweep_signal_vect);
data_mat=[data_mat ;zeros(abs(mod(length(data_mat),ct.N_sweep_avg)-ct.N_sweep_avg),N.N_mic)];
% data_mat=[data_mat ; zeros(-mod(length(data_mat),ct.N_sweep_avg)+ct.N_sweep_avg,N.N_mic);
% data_mat=data_mat(1:N.N_sweep,:,:);

%% Some constant definition and useful vector reprensatation

N.N_sweep_avg = floor(N.N_sweep/ct.N_sweep_avg);
ct.dfe_sweep_avg=ct.Fs_sca/N.N_sweep_avg;
t.Fsweep_avg=0:ct.dfe_sweep_avg:(N.N_sweep_avg-1)*ct.dfe_sweep_avg;


%% filter (antialiasing)
n=16;
fc=2000;
wn=fc*2/ct.Fs_sca;
[b,a] = butter(n,wn);
% data_mat=filter(b,a,data_mat);
% sweep_signal_vect=filter(b,a,sweep_signal_vect);

%% Processing of the impulse response and FRF


h_fft_mat=zeros(N.N_sweep_avg,N.N_mic,ct.N_sweep_avg);%initialisation
var=N.N_sweep_avg;%avoid "over communication variable"
% Divide each colomn of fourier vector out by fourier vector in
for ii=1:ct.N_sweep_avg
    h_fft_mat(:,:,ii)=bsxfun(@rdivide,fft_norm(data_mat(1+var*(ii-1):var*ii,:)),fft_norm(sweep_signal_vect(1+var*(ii-1):var*ii,:)));
end

System.h_sig_fft=mean(h_fft_mat,3);%frf
System.h_sig=real(ifft(System.h_sig_fft));% impulse response



%% Revert the matrix if needed to be in the shape of the input 

if b>a
    data_mat=permute(data_mat,[2 1 3]);
    System.h_sig_fft=permute(System.h_sig_fft,[2 1 3]);
    System.h_sig=permute(System.h_sig,[2 1 3]);
end

end
