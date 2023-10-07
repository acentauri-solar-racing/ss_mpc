clear all
close all
clc

import casadi.*

warning('off');
addpath('..\..\ss_offline_data\route\BWSC');

addpath('..\..\ss_offline_data\parameters');
addpath('..\..\ss_offline_data\route');

addpath("Functions_Scripts\");
addpath("Functions_Scripts\Model");
addpath("Functions_Scripts\load_function");
addpath("Functions_Scripts\MPC");
addpath("Functions_Scripts\NLP");

addpath("OfflineData\");
addpath("OnlineData\");