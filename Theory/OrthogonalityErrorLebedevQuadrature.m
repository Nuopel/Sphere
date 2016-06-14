% Compute the orthogonality error affilied to the lebedev quadrature
% Auteur : Dupont Samuel
% Version : 1.0 June 2016
clear variables; close all; clc ;

% This program aim to show the orthogonality error affilied to the lebedev quadrature
% D= Id - YtWY

%% Initialisation of the variables
ct.M=8;
ct.r_hp=1.07;
ct.N_hp=50;


%% Calculate  speakers positions
[ LoudspArray.x, LoudspArray.y, LoudspArray.z, LoudspArray.w ] = ld_by_order(ct.N_hp);
LoudspArray.pos_hp = [ LoudspArray.x, LoudspArray.y, LoudspArray.z ]*ct.r_hp;
[ LoudspArray.theta, LoudspArray.phi ,LoudspArray.r ] = cart2sph( LoudspArray.pos_hp(:,1),LoudspArray.pos_hp(:,2),LoudspArray.pos_hp(:,3));

%% Calculate Spherical harmonics
Ymn = sph_harmonic( ct.M,ct.N_hp,LoudspArray.theta,LoudspArray.phi );

%% Orthogonality matrix
D=diag(ones((ct.M+1)^2,1))-(Ymn*diag(LoudspArray.w)*Ymn.');


ph=color(D)
set(h, 'EdgeColor', 'none');
colorbar
% shading interp
