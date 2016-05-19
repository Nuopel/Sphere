function [ freq ,fourier] = ISYfourier(vector,Fs,p)
%Isummonyoufourier
%   Take a vector, and is sampling frequency and sends back the fourier
%   transform and is corresponding frequencies
%  p option to plot
if nargin == 2;
    p=0;
end;

N=length(vector);
dfe=Fs/N;
freq=0:dfe:(N-1)*dfe;
fourier=fft_norm(vector);

if p==1
    plot(freq,20*log10(abs(fourier)/2e-5))
%     xlim([20 20000])
    xlabel('Frequency [Hz]')
    ylabel('Amplitude [dB_{spl}]')

end

end

