function [ Array ] = CreateSpeakerStructure( x, y, z, w )
%ARRAYCONDITIONING Summary of this function goes here
%   Detailed explanation goes here
% x,y and z in radians

%% Get L

L = numel(x);

%% Create the structure 

Array.L = L; % How may loudspeakers
Array.x = x;
Array.y = y;
Array.z = z;
Array.N = [ones(L,1) zeros(L,2)];
Array.w = w;
[ Array.theta, Array.phi, Array.r ] = cart2sph( Array.x, Array.y, Array.z ) ;

Array.coord_vect = [Array.x';Array.y';Array.z']; % transform the grid in a single row of pair coordinate of the mic

end