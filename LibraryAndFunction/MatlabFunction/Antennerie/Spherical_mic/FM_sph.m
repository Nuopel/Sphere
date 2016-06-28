function [ a ] = FM_sph( m,kind,k,r)
%[a] = FM_sph( m,kind,k,r);
% examplefm=FM_sph( 5,2,2*pi*1000/340,0.07)
% Equalisation filter used in ambisonics for monopole virtualisation
a=-besselh(m+0.5,kind,k.*r).*sqrt(pi./(2.*k*r)).*(1i^(-m+1).*k/(4*pi));

end