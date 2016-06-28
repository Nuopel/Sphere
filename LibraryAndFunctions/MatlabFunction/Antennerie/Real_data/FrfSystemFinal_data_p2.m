function [ data_mat ,System,N,ct,t] = FrfSystemFinal_data_p2( data_mat,sweep_signal_vect,N,ct )
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
N.sweep=length(sweep_signal_vect);
N.meas=length(data_mat);

%% Some constant definition and useful vector reprensatation

N.N_sweep_avg = 2^16;
ct.dfe_sweep_avg=ct.Fs_sca/(N.N_sweep_avg*2+1);
t.Fsweep_avg=0:ct.dfe_sweep_avg:(N.N_sweep_avg-1)*ct.dfe_sweep_avg;


%% filter (antialiasing)
n=16;
fc=2000;
wn=fc*2/ct.Fs_sca;
[b,a] = butter(n,wn);
% data_mat=filter(b,a,data_mat);
% sweep_signal_vect=filter(b,a,sweep_signal_vect);

%% Processing of the impulse response and FRF

var=N.N_sweep_avg;%avoid "over communication variable"
h_fft_mat=zeros(2*var+1,N.N_mic,ct.N_sweep_avg);%initialisation
% Divide each colomn of fourier vector out by fourier vector in
for ii=1:ct.N_sweep_avg
    h_fft_mat(:,:,ii)=bsxfun(@rdivide,fft(data_mat(1+var*(ii-1):var*ii,:),2*var+1),fft(sweep_signal_vect(1+var*(ii-1):var*ii,:),2*var+1));
end

System.h_sig_fft=mean(h_fft_mat,3);%frf

System.h_sig=real(ifft(System.h_sig_fft));% impulse response
w=tukeywin(2^16,0.25);
w=[ones(2^15,1); w(end/2:end-1,1)];
System.h_sig=bsxfun(@times,System.h_sig(1:2^16,:),w);
System.h_sig_fft=fft(System.h_sig);


%% Revert the matrix if needed to be in the shape of the input 

if b>a
    data_mat=permute(data_mat,[2 1 3]);
    System.h_sig_fft=permute(System.h_sig_fft,[2 1 3]);
    System.h_sig=permute(System.h_sig,[2 1 3]);
end

end
