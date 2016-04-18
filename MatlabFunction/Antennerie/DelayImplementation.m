function [ sig_mic_mat ] = DelayImplementation(ht,sig_hp_mat,L,N,D,opt )
% DelayImplementation apply the the delay using the matrix ht to the signal sig_hp_mat

%  HT : Matrix of the fir filter delay.
%  SIG_HP_MAT : signal to convert.
%  L : number of speakers.
%  N : structure containing informations regarding length of vectors...
%   N.N_sca= length of the signal
%   N.N_mic=number of microphone on the antenna
%  D : structure containing informations regarding the delay
%   D.tabmax= matrix containing the zeros to add
%   D.maxi= maximum value of D.tabmax

%Implement fft, check if faster

%% Check if fft or filter option (filter default)
if ~exist('opt','var')
    opt= 0;
else
    opt=1;
end

%% Verification size of entry signal
[a, b ]=size(sig_hp_mat);
if a>b
    sig_hp_mat=permute(sig_hp_mat,[2 1 3]);
end
[~, b ]=size(sig_hp_mat);

%  initialisation

% select fft or filter way

if opt==0
    sig_mic_mat=zeros(N.N_mic,N.N_sca+D.maxi);
    var=zeros(L,N.N_sca+D.maxi);
    for ii=1:N.N_mic
        for jj=1:L
            var(jj,:)=[zeros(1,D.tabmax(jj,ii)) filter(ht(jj,:,ii),1,sig_hp_mat(jj,:)) zeros(1,D.maxi-D.tabmax(jj,ii))]/(4*pi*D.R_mat(jj,ii));
        end
        if L~=1
        sig_mic_mat(ii,:)=sum(var);
        else
           sig_mic_mat(ii,:)=var;
        end
    end
    
else
%     disp('fft')
%     h=zeros(L,N.N_sca);
%     sig_mic_mat=zeros(N.N_mic,N.N_sca);
%     for ii=1:N.N_mic
%         for jj=1:L
%             h(jj,1:D.tabmax(jj,ii)+N.NOrder+1)=[zeros(1,D.tabmax(jj,ii)) ht(jj,:,ii)];
%         end
%         h=real(ifft(bsxfun(@times,fft(h(:,:)')',fft(sig_hp_mat(:,:)))));
%         sig_mic_mat(ii,:)=sum(h)';
%     end
    
end



