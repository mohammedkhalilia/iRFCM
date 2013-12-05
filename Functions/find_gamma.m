function varargout = find_gamma(D, delta)
%%
%   Computes the additive constant c if the provided matrix D is not
%   Euclidean. The simplest approach is using Beta-Spread, where c is 
%   computed using the doubly centered matrix. 
%   Using Beta-Spread we center the matrix D according to
%
%           B = (I - 11') * (-0.5*D) * (I - 11')
%               OR if D is the squared distances, then
%           B = (I - 11') * (-0.5*D.^2) * (I - 11')
%           
%           Then check if the matrix positive semi deifnite (psd). If D is not psd
%           Then c = -2 * lambda (lambda smallest eigenvalue of B) (see [1-2])
%
%   If delta is provided a different computation is used to find c, which
%   documented in ref. [3]
%
%   Usage c = find_gamma(D, delta)
%
% varargout:
% varargout{1} - 0 (D is not Euclidean) or 1 (D is Euclidean)
% varargout{2} - additive constant
%
% D     - n x n dissimilarity matrix
% delta - delta is the matrix that is used to tranform R to Euclidean.
%         delta is of size n x n and can be for instance
%           delta = 1-eye(n)  (DEFAULT)
%           delta = subdominant ultrametric matrix
%           delta = R
%           delta = power of R
%           delta = parametric function
%           delta is expected to have an Euclidean representation,
%           otherwise an error will be thrown (see [3])
%
% Refs:
%   [1] T. Cox and M. Cox, Multidimensional scaling. 2010.
%   [2] K. Mardia and J. Kent, “y BIBBY, JM (1979): Multivariate analysis,” Academic Press, London.
%   [3] J. Benasseni, M. B. Dosse, and S. Joly, “On a General Transformation Making a Dissimilarity 
%       Matrix Euclidean,” Journal of Classification, vol. 24, no. 1, pp. 33–51, Jun. 2007.
    
    %Get c based on the positive semi definite test without delta
    %See [1-2]. In other words, use Beta-Spread
    if nargin == 1
        [euc, lambda, B, eigValMat] = is_euclidean(D);
        varargout{1} = euc;
        varargout{2} = -2*lambda;
        varargout{3} = B;
        varargout{4} = eigValMat;
        
    %Get c based on the positive semi definite test with a provided delta
    %See [3]
    elseif nargin == 2
        %check if delta is Euclidean, otherwise error 
        if size(D) ~= size(delta)
            error('Dissimilarity matrix and delta must have same dimensions.'); 
        end
            
        %Wdelta is the doubly centered matrix such that 
        %Wdelta = (In - 11') * delta * (In - 11'), where In is the
        %identity matrix
        [euc, ~, wDelta] = is_euclidean(delta);
        
        if ~euc
            error('Delta is not Euclidean'); 
        end
            
        %compute the dounbly centered matrix corresponding to D^2
        [~, ~, wD] = is_euclidean(D); 
        
        %get the eignvalues and vectors of Wdelta
        [eigVecMat eigValMat] = eig(wDelta);
        
        %one of the eignvalues must be zero
        idx = find(sum(eigValMat) > 1e-6);

        %now we end up of diagonal eignvalue of size (n-1) x (n-1)
        eigValMat = eigValMat(idx,idx);

        %and eignvectors matrix of size n X (n-1)
        eigVecMat = eigVecMat(:,idx);
        
        %raise the eignevalues to the power -0.5
        eigValMat(logical(eye(size(eigValMat)))) = diag(eigValMat).^-0.5;

        %compute the additive constant c
        B = (eigValMat*eigVecMat') * wD * (eigVecMat*eigValMat);
        eigValArray = eig(B);
        eigValArray = real(eigValArray);
        
        lambda = min(eigValArray);
        lambda(lambda<0 & lambda>-0.000001) = 0;

        varargout{1} = lambda >= 0;
        varargout{2} = -1 * lambda;
    end
end