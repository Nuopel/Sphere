function [Bmn]= Bmn_encoding_sph(Pressure,Sphmic,ct,var,opt)
% Encoding the Bmn coefficient from pressure of the spherical microphone
% PRESSURE : vector containing the pressure as function of the frequence
% SPHMIC : structure containing the positions of the microphone
% CT : structure containing the k vector, the distance of the microphones,
%       order of the hankel funtion (2), order of truncature M...
% VAR : structure containing vectors of the sum of harmonics by order and
%        their number
% N : Contains the length of Bmn
%
% Options : you can set regularisation filter by setting opt to :
%             -'tik', tikhonov filters limited at 50 dB 
%             -'rtarget', filter depending on a desired radius
%             - nothing, default theoretical filters
% Examples:
%           Bmn.recons = Bmn_encoding_sph(HData.h_sig_fft(ct.pos,:),Sphmic,ct,var,'tik' );
%           Bmn.recons = Bmn_encoding_sph(HData.h_sig_fft(ct.pos,:),Sphmic,ct,var);  
% Auteur : Dupont Samuel
% Version : 2.0 Fevrier 2016

if nargin<5
    opt=0;
elseif strcmp(opt,'tik')
        opt=1;

elseif strcmp(opt,'rtarget')
        opt=2;
else
    opt=0;
end


%% Initialisation
var.m_vect=0:ct.M;
var.m_sum_vect=(var.m_vect+1).^2;
var.nbr_m=(2.*var.m_vect)+1;

[ ct.k ] = ResizeColumn( ct.k ) ; % check dimension

[a,~]=size(Pressure);
if a~=1
[ Pressure ] = ResizeColumn( Pressure ) ; % check dimension
end
[N.N_sweep,~]=size(Pressure);

Bmn = zeros(var.m_sum_vect(ct.M+1),N.N_sweep) ;% init

var.k=ct.k;
N.k=length(var.k);

var.Hprim=zeros(N.N_sweep,(ct.M+1)^2);% H_prim init

for ii=0:ct.M
    var.Hprim(:,(ii)^2+1:(ii+1)^2) = repmat(1i^(ii-1)./((ct.k.*ct.r_micsph).^2.* ...
        Hankel_sph_1_deriv(ii,ct.hankel_order,ct.k.*ct.r_micsph)),1,var.nbr_m(ii+1)) ;
end


Ymn.Micrecons =  sph_harmonic( ct.M, ct.N_mic, Sphmic.theta, Sphmic.phi ) ; %construction of the spherical harmonics at the mic position

%  Calculation of the Bmn coefficient
% fprintf('\n Bmn Encoding \n') ;
switch opt
    case 0
        disp('Theoretical filters')
        for ii=1:N.N_sweep
            Bmn(:,ii) = diag(1./var.Hprim(ii,:))*Ymn.Micrecons*diag(Sphmic.w)*Pressure(ii,:).' ;
        end
    case 1 %% Tikhonov
        disp('Tikhonov')
        
        %% Filtres theoriques
        var.EqFilt=1./ var.Hprim;
        % Regularisation
        ct.ac=20;% maximum noise amplification
        ct.a=10^(ct.ac/20);% maximal amplification for filter
        ct.a2=sqrt(ct.N_mic)*10^(ct.ac/20);% maximal amplification for filter
        
        ct.lambda=(1-sqrt(1-1/ct.a^2))/(1+sqrt(1-1/ct.a^2));
        
        %Calculation of the regularised filter with tikhonov formula
        var.EqFilt_reg=conj(var.Hprim)./(abs(var.Hprim).^2+ct.lambda^2);
        for ii=1:N.N_sweep
            Bmn(:,ii) = diag(var.EqFilt_reg(ii,:))*Ymn.Micrecons*diag(Sphmic.w)*Pressure(ii,:).' ;
        end
    case 2 % Rtarget
        %% Regularisation
        disp('Rtarget')
        ct.N=2^16;
        ct.r_target=0.2;
        ct.k=FreqVect(ct.Fs,ct.N)*2*pi/ct.c_air;
        ct.f_select=closest(ct.k,var.k);
        ct.n=6;
        ct.fc=(1:ct.M).*ct.c_air/(2*pi*ct.r_target);
        ct.wn=ct.fc*2/ct.Fs;
        var.h=ones(ct.N,var.m_sum_vect(ct.M+1));
        var.EqFilt=1./ var.Hprim;
        
        for ii=1:ct.M
            [b,a] = butter(ct.n,ct.wn(ii)/2,'high');
            var.h(:,(ii)^2+1)=freqz(b,a,ct.N,ct.Fs*2);
            var.h(:,(ii)^2+1:(ii+1)^2)= repmat(var.h(:,(ii)^2+1),1,var.nbr_m(ii+1)) ;
        end
        
        var.EqFilt_reg=bsxfun(@times,var.EqFilt,abs(var.h(ct.f_select,:)));
        
        for ii=1:N.N_sweep
            Bmn(:,ii) = diag(var.EqFilt_reg(ii,:))*Ymn.Micrecons*diag(Sphmic.w)*Pressure(ii,:).' ;
        end
        
end

end



