function [Dh, gamma] = euclideanize(D, delta, gamma, transformFunc)
%%
%
%   Transforms a dissimilarity matrix D to an Euclidean matrix Dh. 
%   Simplest way to do that is by using the Beta-Spread transformation.
%   Which first involves centering the matrix D according to
%   
%       B = (I - 11') * (-0.5*D) * (I - 11')
%
%   But if we unsquared distances, then B is
%
%       B = (I - 11') * (-0.5*D.^2) * (I - 11')
%
%   Then check if the matrix positive semi deifnite and return the
%   smallest eigenvalue as the constant c (see [1-2])
%   NOTE: D here is expected to contain squared dissimialrities.
%
% Usage [Dh,c] = make_euclidean(D, c, delta, transfrom_func)
%
% Dh                - the Euclidean realization of D
% gamma             - the constant used to tranform D to Dh
% D                 - dissimilarity matrix
% c                 - Optional Additive Constant provided, if not provided
%                       it will be computed
% transformFunc     - this is a function that is used to tranform D to Dh. 
%                     By default this function is:
%
%                       Dh = D + c*delta;
%                       
%                     However, in some cases as some might want the square root of Dh:
%                     Dh = (D + c*delta).^0.5;
%
%                     Remember that RFCM requires squared dissimilarities
%                     and therefore Dh = Dh = D + c*delta;
%
% delta             - delta is the matrix that is used to tranform R to Euclidean.
%                     delta is of size n x n and can be for instance
%                   
%                       delta = 1-eye(n) (DEFAULT)
%                       delta = subdominant ultrametric matrix
%                       delta = D
%                       delta = power of D
%                       delta = parametric function
%                       
%                     delta is expected to have an Euclidean representation,
%                     otherwise an error will be thrown (see [3-4])
%                  
% Refs:
%   [1] T. Cox and M. Cox, Multidimensional scaling. 2010.
%   [2] K. Mardia and J. Kent, “y BIBBY, JM (1979): Multivariate analysis,” Academic Press, London.
%   [3] J. Benasseni, M. B. Dosse, and S. Joly, “On a General Transformation Making a Dissimilarity 
%       Matrix Euclidean,” Journal of Classification, vol. 24, no. 1, pp. 33–51, Jun. 2007.
%   [4] J. Dattorro, Convex optimization and Euclidean distance geometry. 2005.
    
    %if c is not provided, then first check if the matrix is Euclidean 
    %and compute c
    if nargin < 3 || gamma == 0 
        %if delta is not provided, then use the default, Beta-Spread
        if nargin < 2 || isempty(delta)
            [euc, gamma] = find_gamma(D);
        else
            [euc, gamma] = find_gamma(D, delta);
        end
        
        if euc
            fprintf('The dissimilarity matrix is already Euclidean\n');
        end
    end
    
    %if delta is not provided, then assume that delta = 11' - I
    if nargin < 2 || isempty(delta)
        n = size(D,1);
        delta = ones(n) - eye(n);
    end
    
    if nargin == 4 && ~strcmp(transformFunc,'')
        eval(transformFunc);
    else
        Dh = D + gamma.*delta;
    end
end