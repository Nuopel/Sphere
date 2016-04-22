function [ a ] = Bessel_sph_1_deriv( m,kr)
% Function de bessel spherique
% m = order of the function
% kind = kind of the function (1 or 2)
% kr = evaluated value 
% http://mathworld.wolfram.com/SphericalHankelFunctionoftheFirstKind.html

% a=0.5*(Bessel_sph(m-1,kr)-(m+1)/kr*Bessel_sph(m,kr));
a=(-Bessel_sph(m+1,kr)-(m)/kr*Bessel_sph(m,kr));
a=(Bessel_sph(m-1,kr)-(m+1)/kr*Bessel_sph(m,kr));

end