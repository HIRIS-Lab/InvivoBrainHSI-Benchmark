function [HoutVec, nodes] = labelConsistentHierarchy(H)
% Input: H: rows \times cols \times numlevels or H: numLevels \times numpix
% Output: HoutVec
% labels classes consistently with the same label across the hierarchy
% iteratively checks between parent and child if they share the same class
% and propogates the label till root.
% B Ravi Kiran
% Sept 2015

if(length(size(H))==3)
    [r,c,numLevels] = size(H);
    numPix = r*c;
    Hvec = reshape(H,numPix,numLevels); 
    Hvec = Hvec';
else
    [numLevels, numPix] = size(H);
    Hvec = H;
end

%output Hierarchy in vector format
HoutVec = zeros(numLevels,numPix);    
%first level is root label 1
HoutVec(1,:) = ones(1,numPix);

for level = 1:numLevels-1
    parent = HoutVec(level,:);
    child = Hvec(level+1,:);
    parentLabels = unique(parent);
    for j=1:length(parentLabels)
        parentClass = parent==parentLabels(j); 
        label = unique(child(parentClass));
%         disp(['Level = ' num2str(level) '--Labels = ' num2str(label)]);
        if(length(label)==1)
            HoutVec(level+1,parentClass) = parentLabels(j);
        else            
            tempMat = HoutVec;
            childClass1 = child==label(1);
            childClass2 = child==label(2);
            HoutVec(level+1,childClass1) = ones(1,sum(childClass1))*(max(tempMat(:))+1);
            HoutVec(level+1,childClass2) = ones(1,sum(childClass2))*(max(tempMat(:))+2);
        end
    end
end

%create tree structure
nodes(1) = 0; %root
for level=1:numLevels-1
    parent = HoutVec(level,:);
    child = HoutVec(level+1,:);
    parentLabels = unique(parent);
    for j=1:length(parentLabels)
        parentClass = parent==parentLabels(j); 
        label = unique(child(parentClass));
        if(length(label)>1)
            nodes(label) = parentLabels(j);
        end
    end
end
    
