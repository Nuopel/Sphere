function [Bmn]= Bmn_encoding_sph(Pressure,Sphmic,ct,var,opt)
% Encoding the Bmn coefficient from pressure of the spherical microphone
% PRESSURE : vector containing the pressure as function of the frequence
% SPHMIC : structure containing the positions of the microphone
% CT : structure containing the k vector, the distance of the microphones,
%       order of the hankel funtion (2), order of truncature M...
% VAR : structure containing vectors of the sum of harmonics by order and
%        their number
% N : Contains the length of Bmn
% Samuel Dupont  may 2016

if opt ~=0
    
    opt=1;
    
end

%% Initialisation

[ ct.k ] = ResizeColumn( ct.k ) ; % check dimension
[N.N_sweep,aa]=size(Pressure);
[ Pressure ] = ResizeColumn( Pressure ) ; % check dimension
Bmn = zeros(var.m_sum_vect(ct.M+1),N.N_sweep) ;% init

var.k=ct.k;
var.Hprim=zeros(N.N_sweep,(ct.M+1)^2);% H_prim init

for ii=0:ct.M
    var.Hprim(:,(ii)^2+1:(ii+1)^2) = repmat(1i^(ii-1)./((ct.k.*ct.r_micsph).^2.* ...
        Hankel_sph_1_deriv(ii,ct.hankel_order,ct.k.*ct.r_micsph)),1,var.nbr_m(ii+1)) ;
end


Ymn.Micrecons =  sph_harmonic( ct.M, ct.N_mic, Sphmic.theta, Sphmic.phi ) ; %construction of the spherical harmonics at the mic position

%  Calculation of the Bmn coefficient
disp('Bmn Encoding') ;
if opt==0
    for ii=1:N.N_sweep
        Bmn(:,ii) = diag(1./var.Hprim(ii,:))*Ymn.Micrecons*diag(Sphmic.w)*Pressure(:,ii) ;
    end
else
    N.k=length(var.k);
    
    %% Filtres theoriques
    var.Hprim=zeros(N.k,(ct.M+1)^2);% H_prim init
    
    for ii=0:ct.M
        var.Hprim(:,(ii)^2+1:(ii+1)^2) = repmat(1i^(ii-1)./((var.k.*ct.r_micsph).^2.* ...
            Hankel_sph_1_deriv(ii,ct.hankel_order,var.k.*ct.r_micsph)),1,var.nbr_m(ii+1)) ;
    end
    var.EqFilt=1./ var.Hprim;
    % Regularisation
    ct.ac=20;% maximum noise amplification
    
    ct.a=10^(ct.ac/20);% maximal amplification for filter
    ct.a2=sqrt(ct.N_mic)*10^(ct.ac/20);% maximal amplification for filter
    
    ct.lambda=(1-sqrt(1-1/ct.a^2))/(1+sqrt(1-1/ct.a^2));
    
    %Calculation of the regularised filter with tikhonov formula
    var.EqFilt_reg=conj(var.Hprim)./(abs(var.Hprim).^2+ct.lambda^2);
    
end

end

