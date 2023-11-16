function [K,U,s] = splitclust(M)
% Given a matrix M, split its columns into two subsets
%
% See Section 3 in
%
% Gillis, Kuang, Park, `Hierarchical Clustering of Hyperspectral Images
% using Rank-Two Nonnegative Matrix Factorization', arXiv.
%
%
% ****** Input ******
%  M     : m-by-n data matrix (or a H-by-L-by-m tensor)
%
% ****** Output ******
%   K    : two clusters
%   U    : corresponding centroids
%   s    : first singular value of M

[u,s,v] = fastsvds(M,2); % Initialization: SVD+SPA
Kf = FastSepNMF(s*v',2,0);
U0 = u*s*v(Kf,:)';

[IDX,U] = kmeans(M', 2, 'EmptyAction','singleton','Start',U0', 'MaxIter', 1000);
U = U';
K{1} = find(IDX==1);
K{2} = find(IDX==2);
s = s(1);

end



function [u,s,v] = fastsvds(M,r)
% "Fast" but less accurate SVD by computing the SVD of MM^T or M^TM
% ***IF*** one of the dimensions of M is much smaller than the other.
% Note. This is numerically less stable, but useful for large hyperspectral
% images.

[m,n] = size(M);
rationmn = 10; % Parameter, should be >= 1

if m < rationmn*n
    MMt = M*M';
    [u,s,~] = svds(MMt,r);
    v = M'*u;
    v = v.*repmat( (sum(v.^2)+1e-16).^(-0.5),n,1);
    s = sqrt(s);
elseif n < rationmn*m
    MtM = M'*M;
    [~,s,v] = svds(MtM,r);
    u = M*v;
    u = u.*repmat( (sum(u.^2)+1e-16).^(-0.5),m,1);
    s = sqrt(s);
else
    [u,s,v] = svds(M,r);
end

end


function [J,normM,U] = FastSepNMF(M,r,normalize)
% FastSepNMF - Fast and robust recursive algorithm for separable NMF
%
% *** Description ***
% At each step of the algorithm, the column of M maximizing ||.||_2 is
% extracted, and M is updated by projecting its columns onto the orthogonal
% complement of the extracted column.
%
% See N. Gillis and S.A. Vavasis, Fast and Robust Recursive Algorithms
% for Separable Nonnegative Matrix Factorization, arXiv:1208.1237.
%
% See also https://sites.google.com/site/nicolasgillis/
%
% [J,normM,U] = FastSepNMF(M,r,normalize)
%
% ****** Input ******
% M = WH + N : a (normalized) noisy separable matrix, that is, W is full rank,
%              H = [I,H']P where I is the identity matrix, H'>= 0 and its
%              columns sum to at most one, P is a permutation matrix, and
%              N is sufficiently small.
% r          : number of columns to be extracted.
% normalize  : normalize=1 will scale the columns of M so that they sum to one,
%              hence matrix H will satisfy the assumption above for any
%              nonnegative separable matrix M.
%              normalize=0 is the default value for which no scaling is
%              performed. For example, in hyperspectral imaging, this
%              assumption is already satisfied and normalization is not
%              necessary.
%
% ****** Output ******
% J        : index set of the extracted columns.
% normM    : the l2-norm of the columns of the last residual matrix.
% U        : normalized extracted columns of the residual.
%
% --> normM and U can be used to continue the recursion later on without
%     recomputing everything from scratch.
%
% This implementation of the algorithm is based on the formula
% ||(I-uu^T)v||^2 = ||v||^2 - (u^T v)^2.

[~,n] = size(M);
J = [];

if nargin <= 2, normalize = 0; end
if normalize == 1
    % Normalization of the columns of M so that they sum to one
    D = spdiags((sum(M).^(-1))', 0, n, n); M = M*D;
end

normM = sum(M.^2);
nM = max(normM);

i = 1;
% Perform r recursion steps (unless the relative approximation error is
% smaller than 10^-9)
while i <= r && max(normM)/nM > 1e-9
    % Select the column of M with largest l2-norm
    [a,~] = max(normM);
    % Norm of the columns of the input matrix M
    if i == 1, normM1 = normM; end
    % Check ties up to 1e-6 precision
    b = find((a-normM)/a <= 1e-6);
    % In case of a tie, select column with largest norm of the input matrix M
    if length(b) > 1, [~,d] = max(normM1(b)); b = b(d); end
    % Update the index set, and extracted column
    J(i) = b; U(:,i) = M(:,b);

    % Compute (I-u_{i-1}u_{i-1}^T)...(I-u_1u_1^T) U(:,i), that is,
    % R^(i)(:,J(i)), where R^(i) is the ith residual (with R^(1) = M).
    for j = 1 : i-1
        U(:,i) = U(:,i) - U(:,j)*(U(:,j)'*U(:,i));
    end
    % Normalize U(:,i)
    U(:,i) = U(:,i)/norm(U(:,i));

    % Compute v = u_i^T(I-u_{i-1}u_{i-1}^T)...(I-u_1u_1^T)
    v = U(:,i);
    for j = i-1 : -1 : 1
        v = v - (v'*U(:,j))*U(:,j);
    end

    % Update the norm of the columns of M after orhogonal projection using
    % the formula ||r^(i)_k||^2 = ||r^(i-1)_k||^2 - ( v^T m_k )^2 for all k.
    normM = normM - (v'*M).^2;

    i = i + 1;
end

end