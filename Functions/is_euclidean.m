function varargout = is_euclidean(D)
%%
% Usage [euc, c, B] = is_euclidean(D) 
%   Computes the eigen values of the matrix PAP, returns the
%   smallest constant that makes D Euclidean and a flag indicating if the
%   matrix D is Euclidean or not.
%   NOTE: D is the squared distnaces.
%   This will center the matrix D according to
%       B = (I - 11') * (-0.5*D) * (I - 11')
%           OR if D is the squared distances, then
%       B = (I - 11') * (-0.5*D^2) * (I - 11')
%
%   Then check if the matrix positive semi deifnite and return the
%   smallest eigenvalue as the constant c (see [1-3])
%
% euc   - flag, 1 = D is Euclidean, 0 = D is NOT Euclidean
% c     - smallest constant that makes D Euclidean
% B     - the inner product table produced from centering the matrix D
% D     - the n x n dissimilarity data matrix
%
% Refs:
%   [1] T. Cox and M. Cox, Multidimensional scaling. 2010.
%   [2] K. Mardia and J. Kent, “y BIBBY, JM (1979): Multivariate analysis,” Academic Press, London.
%   [3] J. Benasseni, M. B. Dosse, and S. Joly, “On a General Transformation Making a Dissimilarity 
%       Matrix Euclidean,” Journal of Classification, vol. 24, no. 1, pp. 33–51, Jun. 2007.

    n = size(D,1);
    H = eye(n) - ones(n)/n;
    varargout = {0,0,0,0};
    
    %Compute the centering matrix using H*A*H, where A = -1/2 * d(i,j)
    %This is the most widely cited way for centering the
    %dissimilarty matrix. Here the smallest constant that makes D PSD
    %(positive semi-definite) is smallest eigenvalue of H*A*H
    A = -0.5 .* D;
    B = H*A*H;
    B = single(B);
            
    [~,evalues] = eig(B);
    min_eval = min(evalues(:));

    %sometimes the smallest eigenvalue is very small
    %something like -4e-16. Basically, -0
    %So I am thresholding the negative 0 and converting it to ZERO
    min_eval(min_eval<0 & min_eval>-0.001) = 0;

    varargout{1} = min_eval >= 0;
    varargout{2} = min(evalues(:));
    varargout{3} = B;
    varargout{4} = evalues;  
end