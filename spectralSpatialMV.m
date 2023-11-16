function [classification, pca_r, knn, clusterMap, mv, tmdMap ] = spectralSpatialMV(model, hsCube)
% SPECTRALSPATIALMV Performs spectral-spatial classification.
%
%   [classification, pca_r, knn, segmentation, mv, tmd] =
%   spectralSpatialMV(model, hsCube)
%
%   This function takes a hyperspectral data cube (hsCube) and a model for
%   classification (model), and performs spectral-spatial framework that
%   include classification (classification), PCA reduction (pca_r),
%   K-nearest neighbors (knn) filtering, HKM segmentation (segmentation),
%   majority voting value (mv), and density maps (tmd).
%
%   INPUTS:
%   - model: (obj or struct) 'TreeBagger', 'ClassificationKNN',
%   'ClassificationECOC', 'SeriesNetwork' classification model (e.g., SVM,
%   Random Forest) or Struct with 4-class probability and predicted label for hyperspectral data.    
%       - probAllImage: 4-class probability
%       - classMap: predicted labels
%       - listOfClass: probability order
%   - hsCube: (double) 3D hyperspectral cube to be analyzed.
%
%   OUTPUTS:
%   - classification: (struct) results of the classification usign the model.
%       - classMap: predicted labels
%       - probAllImage:  4-class probability
%       - colorMap: 4-class RGB color map
%   - pca_r: (struct) Principal Component Analysis reduction results.
%       - pcaCube: Principal component scores
%       - pcaComp1: Principal Component
%   - knn: (struct) K-nearest neighbors filtering results.
%       - mapProb: 4-class probability
%       - classMap: predicted labels
%       - colorMap: 4-class RGB color map
%   - clusterMap: Image segmentation results.
%   - mv: (struct) majority voting results.
%       - mapProb: 4-class RGB probability
%       - classMap: predicted labels
%       - colorMap: 4-class RGB color map
%   - tmdMap: density maps results.
%
%   See also: pca, predict, classify
%
%   Authors:  Raquel Leon, Himar Fabelo, Samuel Ortega
%
%   Email address:
%   slmartin@iuma.ulpgc.es, hfabelo@iuma.ulpgc.es, sortega@iuma.ulpgc.es


rng('default') % For reproducibility
addpath(genpath("utils\"));

% Check toolbox dependency
if sum("Statistics and Machine Learning Toolbox" == matlab.addons.installedAddons().Name)
    disp('All Toolbox required installed');
else
    error ('Statistics and Machine Learning Toolbox Missing, install it before continuing.');
end


% Check the dimension of hsCube
if ndims(hsCube) == 3
    % Reshape the 3D matrix to 2D
    hcSize = size(hsCube);
    vectorizedCube = reshape(hsCube, [hcSize(1)*hcSize(2) hcSize(3)]);
else
    error('Input hsCube variable must be a 3D matrix');
end


%% Probabilities generation

if isstruct(model) %Check if model only contains the probabilities
    disp('***********Reading Classification Data...***********');
    listOfClass = model.listOfClass;
    classification.probAllImage = model.probAllImage;
    classification.classMap = model.classMap;

elseif isobject(model) %Check if model is a pre-trained model
    disp('***********Computing Classification...***********');
    algorithmName = class(model);
    % Obtain the probability scores  amd predicted label for each class in each pixel.
    switch algorithmName
        case {'TreeBagger', 'ClassificationKNN', 'ClassificationECOC'}
            [classification.classMap, classification.probAllImage] = predict(model, vectorizedCube);

            if isnumeric(model.ClassNames)
                listOfClass = model.ClassNames;
            else
                listOfClass = str2double(model.ClassNames); %For TreeBagger Models
            end

            if ~isnumeric(classification.classMap) %For TreeBagger Models
                classification.classMap = str2double(classification.classMap);
            end

        case 'SeriesNetwork'
            [classification.classMap, classification.probAllImage] = classify(model, vectorizedCube);
            classification.classMap = single(classification.classMap);

            listOfClass = [1,2,3,4];

        otherwise
            error(['Unexpected model type. Model must be TreeBagger, ClassificationKNN,' ...
                ' ClassificationECOC or SeriesNetwork.'])
    end
end

% Reshape the 3D matrix to 2D
classification.probAllImage = reshape(classification.probAllImage, [hcSize(1), hcSize(2),4]);
classification.classMap = reshape(classification.classMap, [hcSize(1), hcSize(2),1]);

clear model algorithmName
%% PCA
% PCA Generation
disp('***********Computing PCA...***********');

warning('off', 'stats:pca:ColRankDefX')
[~, pcavector, ~] = pca(vectorizedCube); %CHECK OUTPUT
pca_r.pcaCube = reshape(pcavector, [hcSize(1), hcSize(2), hcSize(3)]);
%Get the first component
pca_r.pcaComp1 = pca_r.pcaCube(:,:,1);

clear pcavector
%% KNN
% KNN Filtering Algorithm Parameters:
disp('***********Computing KNN...***********');
lambdaValue = 1;
kValue = 40;

knn.mapProb = knnFilter_window(classification.probAllImage, pca_r.pcaComp1, lambdaValue, kValue);

%Generate a 4 class map
knn.classMap = zeros(hcSize(1), hcSize(2));
[~, currentClassIndex] = max(knn.mapProb, [], 3);
knn.classMap = listOfClass(currentClassIndex);

clear lambdaValue kValue currentClassIndex
%% HKM ----- EVALUAR
% Clustering Algorithm Parameters:
disp('***********Computing HKM...***********');
numberOfClusters = 24;

[clusterMap] = hierclust2nmfMulti(hsCube,numberOfClusters);

clear clusteringAlgoritmName numberOfClusters
%% MV
disp('***********Computing Majority Voting...***********');
[mv.classMap, mv.mapProb] = majorityVoting(knn.classMap, clusterMap);

%% TMD
disp('***********Computing TMD...***********');
[tmdMap] = generateColorMap(mv, 'density');

%% Generate color maps
disp('***********Generating color maps...***********');
knn.colorMap = generateColorMap(knn.classMap, 'fourClass');
mv.colorMap = generateColorMap(mv.classMap, 'fourClass' );
classification.colorMap = generateColorMap(classification.classMap, 'fourClass');

rmpath(genpath("utils\"));
end
