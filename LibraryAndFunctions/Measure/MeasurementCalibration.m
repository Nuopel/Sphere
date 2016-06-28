%% MeasurementCalibration.m
% Used to calibrate measurements from a plane array or spherical microphone
% Use the program CalibrationPiston.m to obtain the calibration file .mat
% from the pistonphone measurement.
%
% It loads the vector containing the rms value of each microphones and do the
% relative ratio compared to microphone 1.
%
% INPUT:the path and name of the measure file is set in 'name'
%
% OUPUT: the file output is the relative calibration of the measurement in .wav
% with '_calib' added to is name
% Auteur : Dupont Samuel
% Version : 1.0 June 2016

clear variables;close all;clc
%% Extraction signal original
a=load('calib_Antenna_mic_31-06.mat');
data.calib=a.calib/a.calib(1);clear a; % relative calibration

name=sprintf('Beamformer/PlaneWave/Antenna_planewave_ambisonique_beamformer0_0_90'); % select name of the file
file=sprintf('%s.w64',name);disp(file)
[data.Out,ct.Fs_sca]=audioread(file); % extract mic signals

%% Calibration 
data.Out=bsxfun(@rdivide,data.Out,[data.calib.', 1]); % apply calibration

%% Write file
file=sprintf('%s_calib.wav',name); disp(file)
audiowrite(file,data.Out,ct.Fs_sca);
