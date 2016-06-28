function [ a ] = Hankel_sph_1_deriv( m,kind,kr)
% Derive de la fonction de hankel spherique
% m = order of the function
% kind = kind of the function (1 or 2)
% kr = evaluated value 
% http://mathworld.wolfram.com/SphericalHankelFunctionoftheFirstKind.html

a=Hankel_sph(m-1,kind,kr)-(m+1)./kr.*Hankel_sph(m,kind,kr);
end