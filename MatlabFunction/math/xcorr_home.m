function [ out ] = xcorr_home(a,b,corrLength)
% homemade xcorr funstion
if nargin == 2;
    Na=length(a);Nb=length(b);
    if Na==Nb
        corrLength=Na+Nb;
    else
        corrLength=Na+Nb-1;
    end
end;

out=fftshift(ifft(bsxfun(@times,fft(a,corrLength),conj(fft(b,corrLength)))));
end

