function [ Sphmic,Pressure ] = Decoding_pressure_microphone( Bmn,Sphmic,N,ct,var )
% Decode the pressure along a spherical array form Bmn coefficient using
% spherical harmonics
% BMN : Coefficients containing information about space
% SPHMIC : struct containg position x, y, z of the microphone
% N : struct
% CT : struct 
% var : struct 
% Samuel Dupont  may 2016

%% Change vector to desired dimension if needed
[a, b ] = size(Bmn.source);
if a>b
    Bmn.source = Bmn.source.';
end

[a, b ] = size(ct.k);
if b>a
    ct.k = ct.k.';
end

%% Sherical harmonics of the microphone
Ymn.Mic = sph_harmonic( ct.M_th, ct.N_mic, Sphmic.theta', Sphmic.phi' ) ;
Pressure.difract=zeros(N.N_sweep,ct.N_mic);%init

%% Carculate Near field filter
for ii=0:ct.M_th
    var.Hprim(:,(ii)^2+1:(ii+1)^2) = repmat(1i^(ii-1)./((ct.k.*ct.r_micsph).^2.* ...
    Hankel_sph_1_deriv(ii,ct.hankel_order,ct.k.*ct.r_micsph)),1,var.nbr_m(ii+1)) ;
end

if sum(isnan(var.Hprim(:)))>0
    var.Hprim(isnan(var.Hprim(:))) = 0+1i*0 ;
    disp('Becareful Nand in near field filter' );
end

%% Calculate pressure
for jj=1:ct.N_mic
    clc ;    fprintf('Microphone %i',jj) ;
    var.Bmn_Ymn = bsxfun(@times,Bmn.source,Ymn.Mic(:,jj)) ;
    Pressure.difract(:,jj) = sum(bsxfun(@times,var.Bmn_Ymn.',var.Hprim),2) ;
end


end
