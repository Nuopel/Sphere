function [ B_sph ] = Bessel_sph(m,kr)
% Send back spherical bessel function of kr, a matrix containing the column
% of m
% http://mathworld.wolfram.com/SphericalBesselFunctionoftheFirstKind.html

% B_sph=sqrt(pi/2)*1./sqrt(kr).*besselh(m+0.5,kind,kr);
B_sph=sqrt(pi/2)/sqrt(kr).*besselj(m+0.5,kr);
end

