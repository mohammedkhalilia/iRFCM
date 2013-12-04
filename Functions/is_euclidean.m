function [euc, min_eval, B, ev] = is_euclidean(D)
%%
%
%   Computes the eigenvalues of the matrix PAP, returns an indicator on
%   whether D is Euclidean or not. If not then the smallest eigenvalue of PAP
%   is returned.
%
%   NOTE: D is the squared distnaces.
%
%   This will center the matrix D according to B = PAP, where A = -0.5*D
%       B = (I - 11') * (-0.5*D) * (I - 11')
%
%           OR if D does not contain squared dissimilarities
%
%       B = (I - 11') * (-0.5*D^2) * (I - 11')
%
%      See ref. [1-3]
%
%   Usage [euc, c, B] = is_euclidean(D) 
%
% euc   - flag, 1 = D is Euclidean, 0 = D is NOT Euclidean
% c     - smallest constant that makes D Euclidean
% B     - the inner product table produced from centering the matrix D
% ev    - the n x n diagonal eigenvalue matrix
%
% Refs:
%   [1] T. Cox and M. Cox, Multidimensional scaling. 2010.
%   [2] K. Mardia and J. Kent, “y BIBBY, JM (1979): Multivariate analysis,” Academic Press, London.
%   [3] J. Benasseni, M. B. Dosse, and S. Joly, “On a General Transformation Making a Dissimilarity 
%       Matrix Euclidean,” Journal of Classification, vol. 24, no. 1, pp. 33–51, Jun. 2007.

    n = size(D,1);
    H = eye(n) - ones(n)/n;
    
    %Compute the centering matrix using H*A*H, where A = -1/2 * D
    %This is the most widely cited way for centering the
    %dissimilarty matrix. Here the smallest constant that makes D PSD
    %(positive semi-definite) is smallest eigenvalue of H*A*H. If one used
    %the Beta-Spread approach
    A = -0.5 .* D;
    B = H*A*H;
    B = single(B);
            
    [~,ev] = eig(B);
    min_eval = min(ev(:));

    %sometimes the smallest eigenvalue is very small
    %something like -4e-16. Basically, -0
    %So I am thresholding the negative 0 and converting it to ZERO
    min_eval(min_eval<0 & min_eval>-0.001) = 0;
    euc = min_eval >= 0;
end