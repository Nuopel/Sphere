function [ Pressure ] = monopole_pressure(k,source,Antenna)
% Create a vector containing the pressure of a monopole along an antenna for :
% K : the wave number
% Source : struct containing x, y, z position of the source
% Antenna : struc containing information about the coordinate of the
%           microphone, ref function: AntennArray
% Samuel Dupont  may 2016
%% creation antenne
[a,b]=size(Antenna.coord_vect);
if a==2
    r_rhp2=sqrt(sum(cat(3,bsxfun(@minus,Antenna.coord_vect(1,:),source.x),bsxfun(@minus,Antenna.coord_vect(2,:),source.y),bsxfun(@minus,zeros(1,length(Antenna.coord_vect(1,:))),source.z)).^2,3));
else
    r_rhp2=sqrt(sum(cat(3,bsxfun(@minus,Antenna.coord_vect(1,:),source.x),bsxfun(@minus,Antenna.coord_vect(2,:),source.y),bsxfun(@minus,zeros(1,length(Antenna.coord_vect(3,:))),source.z)).^2,3));
    
end
Pressure = exp(-1j*k*(r_rhp2))./(4*pi*r_rhp2);

end

