%% SetTurningDataPaths.m
% This script sets the paths to the .mat files containing the data and adds
% the required code directories to the Matlab path

% Top-level path; will also be used to store cache files for UMAP
% basepath = '';
basepath = '/Users/bdeangelis/Desktop/Datasets/OptoMethodDatasets/';

% Path to each of the files
% Gr5a
gr5a_datapath = fullfile(basepath, '20190703_Gr5a-Chrimson_20180314-20190702_Merged(Large).mat');
% Bristle
bristle_datapath = fullfile(basepath, '20190121_R38B08(Bristle)-Chrimson_20171103-20171112.mat');
% Control
ctrl_datapath = fullfile(basepath, '20190617_Chrimson_Control_20190607_Merged.mat');

% % Add code directories to Matlab path
% addpath('Utilities');
% addpath('AnalysisUtilities');

