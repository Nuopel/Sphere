function [ D, System,N ] = CreateFilterFracDelay(D,ArraySpeaker,ct,N )
%CreateFilterFracDelay is used to create fractional delay filter for
%several location (basicaly in between a structure of speaker and a
%structure of microphone). 
%It take the delay in for all the microphone positions refered to all speaker (matrix) .

%[IN]
%R_mat the delay matrix mic-speaker
    %ArraySpeaker structure containing :
        %ArraySpeaker.L number of speaker
    %ct structure containing :
        %Fs_sca sampling frequency
    %N structure containing
        %NOrder order of the lagrange interpolator  
        %N_mic number of microphone
        
%Calculation of the number of samples for delay

D.FullDelay=D.R_mat/ct.c_air*ct.Fs_sca;

% Verification to seee if the minimun delay is possible with the lagrange
% order set ( the less it is the less the approximation is)

if mod(N.NOrder,2)==1
    N.NOrder=N.NOrder-1;
end
while min(D.FullDelay(:))-N.NOrder/2<0
    % change NOrder until it fit delay
    N.NOrder=N.NOrder-2;
end
sprintf('The lagrange order is %i',N.NOrder)


% Define the zero needed before the filter
D.maxi=floor(max(D.FullDelay(:)))-N.NOrder/2;
D.tabmax=floor(D.FullDelay)-N.NOrder/2;

% Define the fractional delay of the filter
D.tabfrac=D.FullDelay-D.tabmax;
System.ht=zeros(ArraySpeaker.L,N.NOrder+1,N.N_mic);

%set the different filter for the speaker and microphone
for ii=1:N.N_mic
    
        System.ht(:,:,ii)=Fractional_delay_lagrange_matrix(N.NOrder, D.tabfrac(:,ii));
   
end

%fourier fir filter 
System.HF=permute(fft(permute(System.ht,[2 1 3])),[2 1 3]); % permute the dimension for fft 



end

