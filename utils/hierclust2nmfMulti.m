function [segmentedHierarchiesCube1] = hierclust2nmfMulti(Hcubes,r)
% Hierarchical custering based on rank-two nonnegative matrix factorization
% applied to hyperspectral image clustering
% B Ravi Kiran, Aug 2015

% Given a data matrix M (m-by-n) representing n data points in an
% m-dimensional space, this algorithm computes a set of clusters obtained
% using the hierarchical rank-two NMF method described in  

% Original implementation by Gillis, Kuang, Park, `Hierarchical Clustering of Hyperspectral Images 
% using Rank-Two Nonnegative Matrix Factorization', arXiv. 
%
% 
% ****** Input ******
%  Hcubes     : cell of cubes{}, each cell is a rows \times cols \times lambdas tensor
%  r          : number of clusters to generate, 
%
% ****** Output ******
%   segmentedHierarchiesCube        : segmented Cube

[rvec,cvec,nb] = size(Hcubes);
M = reshape(Hcubes,rvec*cvec,nb);
M = M';
clear Hcubes

%accumulate the approximation errors for each level
[m,n] = size(M); 

%segmentation labels for r different clusters/levels for n samples
loopIDX = zeros(r,n);
loopIDX(1,:) = ones(1,n); %level 1;

if min(M(:)) < 0
    warning('The input matrix contains negative entries which have been set to zero'); 
    M = max(M,0);
end

manualsplit = 0; % Split according to the proposed criterion

if nargin < 4
    % Intialization of the tree structure
    sol.K{1} = (1:n)'; % All clusters; the first node contains all pixels
    sol.allnodes = 1; % nodes numbering
    sol.maxnode = 1; % Last leaf node added
    sol.parents = [0 0]; % Parents of leaf nodes
    sol.childs = []; % Child(i,:) = child of node i
    sol.leafnodes = 1; % Current clustering: set of leaf nodes corresponding to selected clusters
    sol.e = -1; % Criterion to decide which leafnode to split 
    sol.U(:,1) = ones(m,1); % Centroids
    sol.Ke(1) = 1; % index centroid: index of the endmember
    sol.count = 1; % Number of clusters generated so far
    sol.firstsv = 0; 
    sol.alpha = 0;
end

while sol.count < r 
    %***************************************************************
    % Update: split leaf nodes added at previous iteration
    %***************************************************************
    for k = 1 : length(sol.leafnodes)
        
        if sol.e(sol.leafnodes(k)) == -1 && length(sol.K{sol.leafnodes(k)}) > 1
            % Update leaf node ind(k) by splitting it and adding its child nodes
            [Kc,~,~] = splitclust(M(:,sol.K{sol.leafnodes(k)})); 
            
            if ~isempty(Kc{2}) % the second cluster has to be non-empty: this can occur for a rank-one matrix. 
                % Add the two leaf nodes, child of nodes(sol.leafnodes(k))
                sol.allnodes = [sol.allnodes; sol.maxnode+1; sol.maxnode+2]; 

                sol.parents(sol.maxnode+1,:) = [sol.leafnodes(k) 0]; 
                sol.parents(sol.maxnode+2,:) = [sol.leafnodes(k) 0];

                sol.childs(sol.leafnodes(k), : ) = [sol.maxnode+1 sol.maxnode+2]; 
                sol.childs(sol.maxnode+1 , :) = 0; 
                sol.childs(sol.maxnode+2 , :) = 0; 

                sol.K{sol.maxnode+1} = sol.K{sol.leafnodes(k)}(Kc{1});  
                sol.K{sol.maxnode+2} = sol.K{sol.leafnodes(k)}(Kc{2}); 

                [sol.U(:,sol.maxnode+1),sol.firstsv(sol.maxnode+1), sol.Ke(sol.maxnode+1), sol.alpha(sol.maxnode+1)] = reprvec(M(:,sol.K{sol.maxnode+1})); 
                [sol.U(:,sol.maxnode+2),sol.firstsv(sol.maxnode+2), sol.Ke(sol.maxnode+2), sol.alpha(sol.maxnode+2)] = reprvec(M(:,sol.K{sol.maxnode+2})); 

                % Update criterion to choose next cluster to split
                sol.e([sol.maxnode+1 sol.maxnode+2]) = -1; 

                % Compute the reduction in the error if kth cluster is split 
                sol.e(sol.leafnodes(k)) = sol.firstsv(sol.maxnode+1)^2 + sol.firstsv(sol.maxnode+2)^2 - sol.firstsv(sol.leafnodes(k))^2; 

                sol.maxnode = sol.maxnode+2; 
            end
        end
    end
    %***************************************************************
    % Choose the cluster to split, split it, and update leaf nodes
    %***************************************************************
    if sol.count == 1 % Only one leaf node, the root node: split it.  
        b = 1; 
    elseif manualsplit == 0 % Split the node maximizing the critetion e
        [~,b] = max(sol.e(sol.leafnodes)); 
    end
    
    %Update Leaves    
    sol.leafnodes = [sol.leafnodes; sol.childs(sol.leafnodes(b),:)']; % Add its two children
    sol.leafnodes = sol.leafnodes([1:b-1 b+1:end]); % Remove bth leaf node    
    sol.count = sol.count+1;      
    %update labeling
    loopIDX(sol.count,:) = clu2vec(sol.K(sol.leafnodes));

end %while for h2nmf clusteirng

%label hierarchies consistently and return nodes vector for treeplot
[HoutVec, ~] = labelConsistentHierarchy(loopIDX);

segmentedHierarchiesCube = zeros(rvec,cvec,r);
%reshape clusters to cube row col sizes
for level=1:r
    idxLevel = HoutVec(level,:);
    segmentedHierarchiesCube(:,:,level) = reshape(idxLevel, rvec,cvec);
end

%Remapped between [1,Lmax] hierarchy of segs
cc = unique(segmentedHierarchiesCube(:,:,r));

%get segementation for each cube at input level from different hierarchies
Y = segmentedHierarchiesCube(:,:,r);
[mr,mc] = size(Y);
Y = Y(:); Yout = Y;

for i=1:length(cc)
    Yout(Y==cc(i)) = i;
end

segmentedHierarchiesCube1 = reshape(Yout,mr,mc);


function [u,s,b,alpha] = reprvec(M) 
% Extract "most" representative column from a matrix M as follows: 
% 
% First, it computes the best rank-one approximation u v^T of M. 
% Then, it identifies the column of M minimizing the MRSA with the first
% singular vector u of M. 
% 
% See Section 4.4.1 of 
% Gillis, Kuang, Park, `Hierarchical Clustering of Hyperspectral Images 
% using Rank-Two Nonnegative Matrix Factorization', arXiv. 

[u,s,~] = svds(M,1); 
u = abs(u); 
[m,~] = size(M); 
% Exctract the column of M approximating u the best (up to a translation and scaling)
u = u - mean(u); 
Mm = M - repmat(mean(M),m,1); 
err = acos(  (Mm'*u/norm(u))./( sqrt(sum(Mm.^2)) )' ); 
[~,b] = min(err); 
alpha = mean(err);
u = M(:,b); 


function IDX = clu2vec(K,m,r) 
% Transform a cluster cell to a vector 
if nargin < 3
    r = length(K);
end
if nargin < 2
    % Compute max entry in K
    m = 0; 
    for i = 1 : r
        m = max(0, max(K{i})); 
    end
end
IDX = zeros(m,1); 
for i = 1 : r 
    IDX(K{i}) = i; 
end

        