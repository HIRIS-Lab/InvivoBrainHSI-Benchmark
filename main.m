% MAIN Entry point for executing an example of the spectralSpatialMV function.
%
%   This script demonstrates the usage of the spectralSpatialMV function.
%   This function takes a hyperspectral data cube (hsCube) and a model for
%   classification (model), and performs spectral-spatial framework that
%   include classification (classification), PCA reduction (pca_r),
%   K-nearest neighbors (knn) filtering, HKM segmentation (segmentation),
%   majority voting value (mv), and density maps (tmd).
%
%   Authors:  Raquel Leon, Himar Fabelo, Samuel Ortega
%
%   Email address:
%   slmartin@iuma.ulpgc.es, hfabelo@iuma.ulpgc.es, sortega@iuma.ulpgc.es

%   See also: spectralSpatialMV

%% Example using a pre-trained model 

% Clear workspace and command window
clear all;
clc;
addpath('utils/');

%Load the HS Cube
load("data\Op42C2.mat", "hsCube");

%Load the model 
%First, include the model in the folder. Change the variable name
%if is necessary

load("models\model_RF.mat", "model");

%Compute the spatial-spectral framework
[classification, pca_r, knn, clusterMap, mv, tmdMap ] = spectralSpatialMV(model, hsCube);

%Plot results
plotExample(classification, pca_r, knn, clusterMap, mv, tmdMap);

rmpath('utils/');
%% Example using predicted labels and probabilities

% Clear workspace and command window
clear all;
clc;

addpath('utils/');
%Load the HS Cube
load("data\Op42C2.mat", "hsCube");

%Load the predicted labels and probabilities
%First, include the probabilities in the folder. Change the variable name
%if is necessary

load("models\prob_RF.mat", "model");
%Employ the following struct:
%     model.probAllImage: probabilities
%     model.classMap: predicted labels 
%     model.listOfClass: probabilities order

%Compute the spatial-spectral framework
[classification, pca_r, knn, clusterMap, mv, tmdMap ] = spectralSpatialMV(model, hsCube);

%Plot results
plotExample(classification, pca_r, knn, clusterMap, mv, tmdMap);

rmpath('utils/');
% End of main script
