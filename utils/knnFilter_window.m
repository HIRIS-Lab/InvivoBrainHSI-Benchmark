function [ filteredData ] = knnFilter_window( classificationResult, guideImage, lambda, knn )
% KNNFILTER_WINDOW Applies a k-Nearest Neighbors (kNN) filter to enhance an image.
%   filteredData = KNNFILTER_WINDOW(classificationResult, guideImage, lambda, knn)
%
%   INPUTS:
%       classificationResult - The probability result of the image classification.
%       guideImage - The image used for filtering.
%       lambda - Regularization parameter controlling the strength of filtering.
%       knn - Number of nearest neighbors to consider in the kNN filter.
%
%   OUTPUT:
%       filteredData - The filtered image data.
%
%   Authors:  Raquel Leon, Himar Fabelo, Samuel Ortega
%
%   Email address:
%   slmartin@iuma.ulpgc.es, hfabelo@iuma.ulpgc.es, sortega@iuma.ulpgc.es

inputSize = size(classificationResult);
numberOfLines   = inputSize(1);
numberOfSamples = inputSize(2);
numberOfClasses = inputSize(3);
windowSize = 14;

% 1) Find de K nearest neighbours of each pixel in the feature space:
vectorGuideImage = reshape(guideImage', [numberOfLines*numberOfSamples 1]);
r = 1:numberOfLines;
c = 1:numberOfSamples;
rc = 1:numberOfLines*numberOfSamples;
[rr, cc] = meshgrid(r,c);
featureMatrixV = [vectorGuideImage lambda*rr(:) lambda*cc(:)];

idx = [rc(:) rr(:) cc(:)]; %1D-2D (vector-matrix) index correspondence

SafeBorderSize = floor(windowSize/2); % Number of lines to safely deal with image borders
SafeBorderSizeSamples = SafeBorderSize*numberOfSamples;
knnMat = zeros(numberOfLines*numberOfSamples, knn);

% Upper lines: from 1 to SafeBorderSizeSamples
ii = ones(SafeBorderSizeSamples,1);
windowSamplesIdx = 1:SafeBorderSizeSamples;
%auxDistance = zeros(SafeBorderSizeSamples,1);
for i_up = 1:SafeBorderSizeSamples

    difference = featureMatrixV(ii,:) - featureMatrixV(windowSamplesIdx,:);
    auxDistance = sum(difference.^2,2);
    auxDistance = [auxDistance idx(windowSamplesIdx,1)]; 	%
    % Fill the KNN Matrix:
    [~,similarIndex] = sortrows(auxDistance);			%
    %clear auxDist;
    knnMat(i_up,:) = squeeze(auxDistance(similarIndex(2:knn+1),2)); %
    % Grow window at lower end
    windowSamplesIdx = [windowSamplesIdx i_up+SafeBorderSizeSamples]; % Add 1 col
    ii = [ii+1; ii(end)+1]; % New pixel. Add 1 row
end

% Mid lines: from SafeBorderSizeSamples+1 to (numberOfLines*numberOfSamples)-(SafeBorderSizeSamples)
for i_mid = i_up+1:(numberOfLines*numberOfSamples)-(SafeBorderSizeSamples)
    %     disp(['Calculating KNN for the measuure number: ', num2str(i_mid), ' of ', num2str(numberOfLines*numberOfSamples) ]);
    difference = featureMatrixV(ii,:) - featureMatrixV(windowSamplesIdx,:);
    auxDistance = sum(difference.^2,2);
    auxDistance = [auxDistance idx(windowSamplesIdx,1)]; 	%
    % Fill the KNN Matrix:
    [~,similarIndex] = sortrows(auxDistance);			%
    %clear auxDist;
    knnMat(i_mid,:) = squeeze(auxDistance(similarIndex(2:knn+1),2)); %
    windowSamplesIdx = windowSamplesIdx+1; % Move window
    ii = ii+1; % New pixel
end

% Lower lines
for i_low = i_mid+1:numberOfLines*numberOfSamples
    difference = featureMatrixV(ii,:) - featureMatrixV(windowSamplesIdx,:);
    auxDistance = sum(difference.^2,2);
    auxDistance = [auxDistance idx(windowSamplesIdx,1)]; 	%
    % Fill the KNN Matrix:
    [~,similarIndex] = sortrows(auxDistance);		%
    %clear auxDist;
    knnMat(i_low,:) = squeeze(auxDistance(similarIndex(2:knn+1),2)); %
    % Shrink window from the upper end
    windowSamplesIdx = windowSamplesIdx(2:end); % Del 1st pos
    ii = ii(2:end)+1;
end

% 2) Apply the knn filtering to the classification results:
% From score matrix to score array:
scoreArray = zeros(numberOfLines*numberOfSamples, numberOfClasses);

index = 1;
for i = 1: numberOfLines
    for j = 1:numberOfSamples
        for k = 1:numberOfClasses
            scoreArray(index,k) = classificationResult(i,j,k);
        end
        index = index + 1;
    end
end

% From scores (P) generate the output (O) of the KNN filtering:
output = zeros(numberOfLines*numberOfSamples, numberOfClasses);
for i = 1:numberOfLines*numberOfSamples
    output(i,:) =  mean(scoreArray(knnMat(i,:),:));
end

% From vector output to matrix output:
filteredData = zeros(numberOfLines, numberOfSamples, numberOfClasses);
index = 1;

for i = 1: numberOfLines
    for j = 1:numberOfSamples
        for k = 1:numberOfClasses
            filteredData(i,j,k) = output(index, k);
        end
        index = index + 1;
    end
end

end

