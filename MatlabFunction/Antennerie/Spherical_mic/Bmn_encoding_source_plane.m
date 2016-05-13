function [ Bmn ] = Bmn_encoding_source_plane( Pressure,structure,ct,N,var )
% Encoding the Bmn coefficient from pressure of the spherical structure for
% plane wave 
% PRESSURE : vector containing the pressure as function of the frequence
% STRUCTURE : structure containing the positions of the transducer
% CT : structure containing the k vector, the distance of the microphones,
%       order of the hankel funtion (2), order of truncature M...
% VAR : structure containing vectors of the sum of harmonics by order and
%        their number
% Samuel Dupont  may 2016

%% Initialisation
Bmn = zeros(var.m_sum_vect(ct.M+1),N.N_sweep,ct.N) ;% init

Ymn.Micrecons =  sph_harmonic( ct.M, ct.N, structure.theta, structure.phi ) ; %construction of the spherical harmonics at the mic position

%  Calculation of the Bmn coefficient 
disp('Bmn Encoding') ;
for ii=1:ct.N

    Bmn(:,:,ii) = Ymn.Micrecons(:,ii)*Pressure ;

end

end

