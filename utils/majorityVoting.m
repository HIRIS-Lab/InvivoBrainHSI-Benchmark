function [ outputMap , frecuencyMap] = majorityVoting( classificationMap, clusteredMap )
%MAJORTYVOTING Performs majority voting to generate an output map
%
%   [outputMap, frequencyMap] = majorityVoting(classificationMap, clusteredMap)
%
%   This function takes the map generated by KNN filter and the
%   segmentation map and performs majority voting to generate an output
%   map. The output map represents the majority vote class results, and the
%   frequency map shows the frequency of each class label in the majority
%   vote.
%
%   INPUTS:
%   - classificationMap: A matrix representing a classification map, where each
%     element contains a class label.
%   - clusteredMap: A matrix representing a clustered map, which can be used for
%     aggregating information from neighboring pixels to assist in majority voting.
%
%   OUTPUTS:
%   - outputMap: The result of majority voting, where each pixel contains the
%     majority class label.
%   - frequencyMap: A map showing the frequency of each class label in the
%     majority vote.
%
%   Authors:  Raquel Leon, Himar Fabelo, Samuel Ortega
%
%   Email address:
%   slmartin@iuma.ulpgc.es, hfabelo@iuma.ulpgc.es, sortega@iuma.ulpgc.es

listOfClusters = unique(clusteredMap);

frecuencyMap1 = zeros(size(clusteredMap));
outputMap1 = zeros(size(clusteredMap));

frecuencyMap2 = zeros(size(clusteredMap));
outputMap2 = zeros(size(clusteredMap));

frecuencyMap3 = zeros(size(clusteredMap));
outputMap3 = zeros(size(clusteredMap));

frecuencyMap = zeros([size(clusteredMap) 3]);
outputMap = zeros([size(clusteredMap) 3]);

for currentCluster = listOfClusters'

    % Find the regions of classified pixels in the current cluster region:
    listOfClassesInRegion = unique(classificationMap(clusteredMap == currentCluster));
    % Calculate the frecuency of each class in region:
    frecuencyCounter = zeros(1,length(listOfClassesInRegion));

    for i = 1:length(listOfClassesInRegion)
        frecuencyCounter(i) = sum(classificationMap(clusteredMap == currentCluster) == listOfClassesInRegion(i));
    end

    [ mostVotedFrecuency, mostVotedIndex] = sort(frecuencyCounter, 'descend');

    % Fill the output map with the highets class in the current cluster region:
    outputMap1(clusteredMap == currentCluster) = listOfClassesInRegion(mostVotedIndex(1));
    frecuencyMap1(clusteredMap == currentCluster) =  mostVotedFrecuency(1) / sum(mostVotedFrecuency);

    if (length(listOfClassesInRegion) >= 2)
        frecuencyMap2(clusteredMap == currentCluster) = mostVotedFrecuency(2) / sum(mostVotedFrecuency);
        outputMap2(clusteredMap == currentCluster) =  listOfClassesInRegion(mostVotedIndex(2));
    else
        frecuencyMap2(clusteredMap == currentCluster) = 0;
        outputMap2(clusteredMap == currentCluster) = 0;
    end

    if (length(listOfClassesInRegion) >= 3)
        outputMap3(clusteredMap == currentCluster)= listOfClassesInRegion(mostVotedIndex(3));
        frecuencyMap3(clusteredMap == currentCluster) = mostVotedFrecuency(3) / sum(mostVotedFrecuency);
    else
        frecuencyMap3(clusteredMap == currentCluster) = 0;
        outputMap3(clusteredMap == currentCluster) = 0;
    end
end

outputMap(:,:,1) = outputMap1;
outputMap(:,:,2) = outputMap2;
outputMap(:,:,3) = outputMap3;

frecuencyMap(:,:,1) = frecuencyMap1;
frecuencyMap(:,:,2) = frecuencyMap2;
frecuencyMap(:,:,3) = frecuencyMap3;

end
