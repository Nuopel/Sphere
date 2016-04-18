function [ out ] = fft_norm(a)
%% do the fft_function but normalise the output by the length of the vector sent by 2
N_sca=length(a);
out=fft(a)*2/N_sca;

end

