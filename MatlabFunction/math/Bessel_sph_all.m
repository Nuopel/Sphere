function [ B_sph ] = Bessel_sph_all(m,kind,kr)
% Send back spherical bessel function of kr, a matrix containing the column
% 0:m spherical bessel function
% http://mathworld.wolfram.com/SphericalBesselFunctionoftheFirstKind.html
B_sph=zeros(length(kr),m(end)+1);
for ii=0:m
B_sph(:,ii+1)=sqrt(pi/2)*1./sqrt(kr).*besselh(ii+0.5,kind,kr);
end

end