function plotExample(classification, pca_r, knn, clusterMap, mv, tmdMap)
% PLOTEXAMPLE Plots various maps for visualization.
%
%   plotExample(classification, pca_r, knn, clusterMap, mv, tmdMap) 
% 
%   Plots different maps for visualization, including Classification Map, PCA
%   First Component, KNN Map, HKM Map, Majority Voting Map, and TMD Map.
%
%   INPUT:
%   - classification: Structure containing the color map for the
%     Classification Map.
%   - pca_r: Structure containing the first PCA component.
%   - knn: Structure containing the color map for the KNN Map.
%   - clusterMap: Matrix representing the HKM Map.
%   - mv: Structure containing the color map for the Majority Voting Map.
%   - tmdMap: Matrix representing the TMD Map.
%
%   Authors:  Raquel Leon, Himar Fabelo, Samuel Ortega
%
%   Email address:
%   slmartin@iuma.ulpgc.es, hfabelo@iuma.ulpgc.es, sortega@iuma.ulpgc.es

figure;

% Subplot 1: classification
subplot(2, 3, 1);
imshow(classification.colorMap);
title('Classification Map');

% Subplot 2: PCA
subplot(2, 3, 2);
imagesc(pca_r.pcaComp1); 
colormap('jet'); % Colormap opcional
axis image
axis off
title('First Component');

% Subplot 3: KNN
subplot(2, 3, 3);
imshow(knn.colorMap);
title('KNN Map');

% Subplot 4: HKM
subplot(2, 3, 4);
imagesc(clusterMap); 
colormap('jet'); 
axis image
axis off
title('HKM Map');

% Subplot 5: MV
subplot(2, 3, 5);
imshow(mv.colorMap);
title('Majority Voting');

% Subplot 6: TMD
subplot(2, 3, 6);
imshow(tmdMap);
title('TMD Map');
end