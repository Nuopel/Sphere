function [ data_mat ,System,N,ct,t] = FrfSystemFinal_data( data_mat,sweep_signal_vect,N,ct,fc )
% Do the Frequency response using the method H=out/in
% In: sweep_signal_vect
% Out: data_mat
% N and ct: structure of points and constants evolving in the main program

% Can set a filter at the cut off frequency fc in option

if nargin>4
    opt=1; 
    sprintf('\n Out data filtered at %i',fc)
else 
    opt=0;
end
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

%% Some constant definition and useful vector reprensatation

N.N_sweep_avg = floor(N.N_sweep/ct.N_sweep_avg);
ct.dfe_sweep_avg=ct.Fs_sca/N.N_sweep_avg;
t.Fsweep_avg=0:ct.dfe_sweep_avg:(N.N_sweep_avg-1)*ct.dfe_sweep_avg;


%% filter (antialiasing)
if opt==1
n=8;
wn=fc*2/ct.Fs_sca;
[b,a] = butter(n,wn);
data_mat=filtfilt(b,a,data_mat);
end
%% Processing of the impulse response and FRF


h_fft_mat=zeros(N.N_sweep_avg,N.N_mic,ct.N_sweep_avg);%initialisation
var=N.N_sweep_avg;%avoid "over communication variable"
% Divide each colomn of fourier vector out by fourier vector in
for ii=1:ct.N_sweep_avg
    h_fft_mat(:,:,ii)=bsxfun(@rdivide,fft_norm(data_mat(1+var*(ii-1):var*ii,:)),fft_norm(sweep_signal_vect(1+var*(ii-1):var*ii,:)));
end

System.h_sig_fft=mean(h_fft_mat,3);%frf
System.h_sig=real(ifft(System.h_sig_fft));% impulse response

w=tukeywin(2^15,0.25);
w=[ones(2^14,1); w(2^14:end-1,1)];

System.h_sig=bsxfun(@times,System.h_sig(1:2^15,:),w);
System.h_sig_fft=fft(System.h_sig);

N.N_sweep_avg = length(System.h_sig_fft);
ct.dfe_sweep_avg=ct.Fs_sca/N.N_sweep_avg;
t.Fsweep_avg=0:ct.dfe_sweep_avg:(N.N_sweep_avg-1)*ct.dfe_sweep_avg;

%% Revert the matrix if needed to be in the shape of the input 

if b>a
    data_mat=permute(data_mat,[2 1 3]);
    System.h_sig_fft=permute(System.h_sig_fft,[2 1 3]);
    System.h_sig=permute(System.h_sig,[2 1 3]);
end

end
