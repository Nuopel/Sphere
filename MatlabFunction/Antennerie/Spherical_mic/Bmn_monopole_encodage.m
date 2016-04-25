function [ Bmn ] = Bmn_monopole_encodage(order, source,ct,var )
%% Encode Bmn coefficient using spherical harmonics

for ii = 0:order
    Fm(:,(ii)^2+1:(ii+1)^2) = repmat(FM_sph( ii, 1, ct.k, ct.r_hp_sca ),1,var.nbr_m(ii+1)) ;
end

if sum(isnan(Fm))>0
Fm(isnan(Fm))=0;
disp('Careful Nand value change to 0')
end

Ymn.source = sph_harmonic( order,1,source.theta,source.phi ) ;
Bmn = bsxfun(@times,Fm,permute(Ymn.source,[2 1])) ;

end

