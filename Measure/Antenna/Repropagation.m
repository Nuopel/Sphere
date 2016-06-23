function [ var] = Repropagation(delay,sig_mic_mat,N)


%% Verification size of entry signal
[a, b ]=size(sig_mic_mat);
if a>b
    sig_mic_mat=permute(sig_mic_mat,[2 1 3]);
end

%% repropagation

var=zeros(N.N_mic,length(sig_mic_mat));

for ii=1:N.N_mic
    var(ii,:) = circshift(sig_mic_mat(ii,:).',delay).';
end


var=permute(var,[2 1 3]);


end
