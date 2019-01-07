% MASTER SCRIPT
% Case Studies 2018-19: River Water Quality Model
% by Veronica, Alireza and Freya

%%
clear all
close all
clc

%% global settings for plots
set(groot, 'DefaultLineLineWidth', 3)  % sets line width for all plots in the current session!
set(groot, 'DefaultAxesFontSize', 24)

%% Basic parameters
run('RiverGeometry_SpatialDiscretisation')

%% Temperature Model
% run the Temperature Model
%run('TemperatureModel')
% OR load calculated temperature variables
load('RiverQual_Temperature_variablesCalculations')
run('Temperature_Plots')




