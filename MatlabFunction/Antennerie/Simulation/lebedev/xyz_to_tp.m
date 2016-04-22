function [ t, p ] = xyz_to_tp ( x, y, z )

%*****************************************************************************80
%
%% XYZ_TO_TP converts (X,Y,Z) to (Theta,Phi) coordinates on the unit sphere.
%
%  Modified:
%
%    09 September 2010
%
%  Author:
%
%    Dmitri Laikov
%
%  Reference:
%
%    Vyacheslav Lebedev, Dmitri Laikov,
%    A quadrature formula for the sphere of the 131st
%    algebraic order of accuracy,
%    Russian Academy of Sciences Doklady Mathematics,
%    Volume 59, Number 3, 1999, pages 477-481.
%
%  Parameters:
%
%    Input, real X, Y, Z, the Cartesian coordinates of a point
%    on the unit sphere.
%
%    Output, real T, P, the Theta and Phi coordinates of
%    the point.
%

for ii=1:length(x)
  p(ii) = acos ( z(ii) );

  fact = sqrt ( x(ii) * x(ii) + y(ii) * y(ii) );

  if ( 0.0 < fact )
    t(ii) = acos ( x(ii) / fact );
  else
    t(ii) = acos ( x(ii) );
  end

  if ( y < 0.0 )
    t(ii) = - t(ii);
  end
%
%  Convert to degrees.
%
  t(ii) = t(ii) * 180.0 / pi;
  p(ii) = p(ii) * 180.0 / pi;
end
  return
end
