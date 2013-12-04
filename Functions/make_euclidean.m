function [Dh, c] = make_euclidean(D, delta, c, transform_func)
%**************************************************************************
% Usage Dh = make_euclidean(D, c, delta, transfrom_func)
%
% Transforms a dissimilarity matrix D to an Euclidean matrix Dh. 
% This will center the matrix D according to
%    B = (I - 11') * (-0.5*D) * (I - 11')
% OR if D is the squared distances, then
%    B = (I - 11') * (-0.5*D^2) * (I - 11')
% Then check if the matrix positive semi deifnite and return the
% smallest eigenvalue as the constant c (see [1-2])
%
% Dh                - the Euclidean realization of D
% c                 - the constant used to tranform D to Dh
% D                 - dissimilarity matrix
% c                 - Optional Additive Constant provided, if not provided it will
%                     computed
% transform_func:   - this is a function of squared D, c and delta and return an
%                     Euclidean dissimilarity matrix. By default this function is:
%
%                       Dh = (D + c*delta).^0.5;
%                       
%                     However, in some cases as in the NERFCM 1994 paper the
%                     tranformation function is slightly different, which is:
%                     Dh = (D + c*delta);
%
%                     The only difference is the square root. Because
%                     NERFCM expects the squared distances rather than the
%                     distances themselves.
%
%                     So, you can provider you own function, but make sure the return
%                     variable name is ALWAYS Dh
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
%                     otherwise an error will be thrown (see [3])
%                  
% Refs:
%   [1] T. Cox and M. Cox, Multidimensional scaling. 2010.
%   [2] K. Mardia and J. Kent, “y BIBBY, JM (1979): Multivariate analysis,” Academic Press, London.
%   [3] J. Benasseni, M. B. Dosse, and S. Joly, “On a General Transformation Making a Dissimilarity 
%       Matrix Euclidean,” Journal of Classification, vol. 24, no. 1, pp. 33–51, Jun. 2007.
%**************************************************************************
    
    %if c is not provided, then first check if the matrix is Euclidean 
    %and compute c
    if nargin < 3 || c == 0 
        %if delta is not provided
        if nargin < 2 || isempty(delta)
            [euc, c] = find_constant_c(D);
        else
            [euc, c] = find_constant_c(D, delta);
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
    
    if nargin == 4 && ~strcmp(transform_func,'')
        eval(transform_func);
    else
        Dh = (D + c.*delta).^0.5;
    end
end