function R = rotz(alpha)
%ROTX  rotate around Z by ALPHA
%
%	R = ROTZ(ALPHA)
%
% See also: ROTX, ROTY, ROT, POS.

R = [cos(alpha) -sin(alpha) 0; ...
     sin(alpha)  cos(alpha) 0; ...
              0           0 1];

end