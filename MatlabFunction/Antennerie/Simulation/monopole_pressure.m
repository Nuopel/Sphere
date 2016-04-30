function [ Pressure ] = monopole_pressure(k,source)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% creation antenne
ct.pas_m = 3e-2; % pas de l'antenne
N.nbrx_sca =50; % nombre de micros par ligne
N.nbry_sca = 50; % nombre de micros par ligne
ct.N_mic=N.nbrx_sca*N.nbry_sca;
[ Antenna ] = AntennArray( ct.pas_m,N.nbrx_sca,N.nbry_sca) ;

r_rhp2=sqrt(sum(cat(3,bsxfun(@minus,Antenna.coord_vect(1,:),source.x),bsxfun(@minus,Antenna.coord_vect(2,:),source.y),bsxfun(@minus,zeros(1,length(Antenna.coord_vect(1,:))),source.z)).^2,3));

Pressure = exp(-1j*k*(r_rhp2))./(4*pi*r_rhp2);

end

