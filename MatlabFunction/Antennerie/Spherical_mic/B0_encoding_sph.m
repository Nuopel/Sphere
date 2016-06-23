function [Bmn]= B0_encoding_sph(Pressure,Sphmic,ct,var)
% Encoding the B0 coefficient from pressure of the spherical microphone
% PRESSURE : vector containing the pressure as function of the frequence
% SPHMIC : structure containing the positions of the microphone
% CT : structure containing the k vector, the distance of the microphones,
%       order of the hankel funtion (2), order of truncature M...
% VAR : structure containing vectors of the sum of harmonics by order and
%        their number
% N : Contains the length of Bmn
% Samuel Dupont  may 2016



%% Initialisation

[ ct.k ] = ResizeColumn( ct.k ) ; % check dimension
[N.N_sweep, ~]=size(Pressure);
[ Pressure ] = ResizeColumn( Pressure ) ; % check dimension

var.k=ct.k;


for ii=0:0
    var.Hprim(:,(ii)^2+1:(ii+1)^2) = repmat(1i^(ii-1)./((ct.k.*ct.r_micsph).^2.* ...
        Hankel_sph_1_deriv(ii,ct.hankel_order,ct.k.*ct.r_micsph)),1,var.nbr_m(ii+1)) ;
end


Ymn.Micrecons =  sph_harmonic( 0, ct.N_mic, Sphmic.theta, Sphmic.phi ) ; %construction of the spherical harmonics at the mic position

%  Calculation of the Bmn coefficient
disp('Bmn Encoding') ;
    for ii=1:N.N_sweep
        Bmn(ii) = diag(1./var.Hprim(ii,:))*Ymn.Micrecons*diag(Sphmic.w)*Pressure(ii,:).' ;
    end
    Bmn(1)=0;
end

