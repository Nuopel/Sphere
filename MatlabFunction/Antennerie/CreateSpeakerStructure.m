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

end