function [ var] = Repropagation(delay,sig_mic_mat,N)


%% Verification size of entry signal
[a, b ]=size(sig_mic_mat);
if a>b
    sig_mic_mat=permute(sig_mic_mat,[2 1 3]);
end

%% repropagation
    D.maxi=max(delay);

    var=zeros(N.N_mic,length(sig_mic_mat)+D.maxi);
    for ii=1:N.N_mic
           var(ii,:)=[zeros(1,delay(ii)) sig_mic_mat(ii,:) zeros(1,D.maxi-delay(ii))];
    end
   
if a>b
    var=permute(var,[2 1 3]);
end    
 
end
