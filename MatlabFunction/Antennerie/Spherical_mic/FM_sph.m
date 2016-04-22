function [ a ] = FM_sph( m,kind,k,r)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
a=-besselh(m+0.5,kind,k.*r).*sqrt(pi./(2.*k*r)).*(1i^(-m+1).*k/(4*pi));

end