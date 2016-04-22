function [ a ] = Hankel_sph( m,kind,kr)
% Function de hankel spherique
% m = order of the function
% kind = kind of the function (1 or 2)
% kr = evaluated value 
% http://mathworld.wolfram.com/SphericalHankelFunctionoftheFirstKind.html

a=besselh(m+0.5,kind,kr).*sqrt(pi./(2.*kr)) ;

end
