clear variables; close all;clc

% Simulate the behavior of a pherical microphone 
% recording a monopole source positioned on the 
% ambisonics set up

%% Define the constant
ct.r_hp_sca = 1.07;%rayon de la sphere
ct.r_micsph = 0.1;

%% Define ambisonics set up
ArraySpeaker = CreateSpeakerSystem(ct.r_hp_sca);% create the sphere set up, sort in struc Array


%% Define Spherical microphone set up (same as speaker but smaller)
Sphmic = CreateSpeakerSystem(ct.r_hp_sca);% create the sphere set up, sort in struc Array

%% Propagate  spherical microphone

