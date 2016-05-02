function [ Bmn ] = Bmn_encoding_sph( Pressure,Sphmic,ct,N,var )
% /Encoding the Bmn coefficient from pressure of the  spherical microphone


var.Hprim=zeros(N.N_sweep,ct.M);% H_prim init
% Function taking into account the pressure of thge direct field and the
% diffracted of the rigid sphere
for ii=0:ct.M
    var.Hprim(:,(ii)^2+1:(ii+1)^2) = repmat(1i^(ii-1)./((ct.k.*ct.r_micsph).^2.*Hankel_sph_1_deriv(ii,ct.hankel_order,ct.k.*ct.r_micsph)),1,var.nbr_m(ii+1)) ;
end
% test

Bmn = zeros(var.m_sum_vect(ct.M+1),N.N_sweep) ;% init
Ymn.Micrecons =  sph_harmonic( ct.M, ct.N_mic, Sphmic.theta, Sphmic.phi ) ; %construction of the spherical harmonics at the mic position

%  Calculation of the Bmn coefficient 
for ii=1:N.N_sweep
    Bmn(:,ii) = diag(1./var.Hprim(ii,:))*Ymn.Micrecons*diag(Sphmic.w)*Pressure.difract(:,ii) ;
    disp(ii) ;
end

end

