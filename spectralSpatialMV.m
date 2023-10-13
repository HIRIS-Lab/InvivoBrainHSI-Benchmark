% Main

function [pca, knn, segmentation, mv, tmd] = spectralSpatialMV(model, hsCube)
rng('default') % For reproducibility

%Transform the HSCube to 2D
hcSize = size(hsCube);
vectorizedCube = reshape(hsCube, [hcSize(1)*hcSize(2) hcSize(3)]);

% PCA Generation
disp('***********Computing PCA...***********');
[pca.coeff, pca.pcavector, pca.latent, pca.tsquared] = pca(vectorizedCube);
pca.pcaMap = reshape(pca.pcavector, [hcSize(1), hcSize(2), hcSize(3)]);
pca.pca1 = pca.pcaMap(:,:,1);

% Load supervised model and probabilities
probMapR = model.ProbAllImage;

%% KNN
% KNN Filtering Algorithm Parameters:
disp('***********Computing KNN...***********');
lambdaValue = 1;
kValue = 40;

Knn.mapProb = knnFilter_window(probMapR, pca.pca1, lambdaValue, kValue);

listOfClass = model.Label;
knn.map = zeros(hcSize(1), hcSize(2));

for z = 1 : hcSize(1)
    for j = 1 : hcSize(2)
        [~, currentClassIndex] = max( knn.mapProb(z, j, :) );
        if isnumeric( listOfClass(currentClassIndex) )
            knn.map(z,j) = listOfClass(currentClassIndex);
        else
            knn.map(z,j) = str2num( cell2mat( listOfClass(currentClassIndex) ) );
        end
    end
end

%% HKM
% Clustering Algorithm Parameters:
disp('***********Computing HKM...***********');
clusteringAlgoritmName = 'hkm';
numberOfClusters = 24;

segmentation.map = hypercubeClustering(hsCube, clusteringAlgoritmName,...
    numberOfClusters);


%% MV
disp('***********Computing Majority Voting...***********');
[mv.map, mv.mapProb] = majorityVoting(knn.map, segmentation.map);

%% TMD
disp('***********Computing TMD...***********');
[tmd.map,tmd.omdMap] = generateDensityMaps(mvMap, mvMapProb);

%% Generate color maps

knn.ColorMap = generateColorMap(knn.map );
mv.ColorMap = generateColorMap( mv.map );
end