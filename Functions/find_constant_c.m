function varargout = find_constant_c(D, delta)
%%
% Usage c = find_constant_c(D, delta)
%
% Computes the additive constant c if the provided matrix D is not
% Euclidean. The simple approach is when delta is not provided, where c is 
% computed using the doubly centered matrix. Also, make sure that the D is
% squared. If you suspect that D is not squared then the input should be
% squared using D.^2.
% When delta is provided computed in terms of D and delta.
% This will center the matrix D according to
%           B = (I - 11') * (-0.5*D) * (I - 11')
%               OR if D is the squared distances, then
%           B = (I - 11') * (-0.5*D^2) * (I - 11')
%           
%           Then check if the matrix positive semi deifnite and return the
%           smallest eigenvalue as the constant c (see [1-2])
%
% c     - the smallest constant to make D Euclidean
% D     - n x n dissimilarity matrix
% delta - delta is the matrix that is used to tranform R to Euclidean.
%         delta is of size n x n and can be for instance
%           delta = 1-eye(n)  (DEFAULT)
%           delta = subdominant ultrametric matrix
%           delta = R
%           delta = power of R
%           delta = parametric function
%           delta is expected to have an Euclidean representation,
%           otherwise an error will be thrown (see [4])
%
% Refs:
%   [1] T. Cox and M. Cox, Multidimensional scaling. 2010.
%   [2] K. Mardia and J. Kent, “y BIBBY, JM (1979): Multivariate analysis,” Academic Press, London.
%   [3] J. Benasseni, M. B. Dosse, and S. Joly, “On a General Transformation Making a Dissimilarity 
%       Matrix Euclidean,” Journal of Classification, vol. 24, no. 1, pp. 33–51, Jun. 2007.
    
    %Get c based on the positive semi definite test without delta
    %See [1-2]
    if nargin == 1
        [euc, c, B, evalues] = is_euclidean(D);
        varargout{1} = euc;
        varargout{2} = -2*c;
        varargout{3} = B;
        varargout{4} = evalues;
        
    %Get c based on the positive semi definite test with a provided delta
    %See [4]
    elseif nargin == 2
        %check if delta is Euclidean, otherwise error 
        if size(D) ~= size(delta)
            error('Dissimilarity matrix and delta must have same dimensions.'); 
        end
            
        %Wdelta is the doubly centered matrix such that 
        %Wdelta = (In - 11') * delta * (In - 11'), where In is the
        %identity matrix
        [euc, ~, Wdelta] = is_euclidean(delta);
        
        if ~euc
            error('Delta is not Euclidean'); 
        end
            
        %compute the dounbly centered matrix corresponding to D^2
        [~, ~, Wd] = is_euclidean(D); 

        %get the eignvalues and vectors of Wdelta
        [vec val] = eig(Wdelta);
        val = real(val);
        
        %one of the eignvalues must be zero
        idx = find(sum(val) > 1e-6);

        %now we end up of diagonal eignvalue of size (n-1) x (n-1)
        val = val(idx,idx);

        %and eignvectors matrix of size n X (n-1)
        vec = vec(:,idx);
        
        %raise the eignevalues to the power -0.5
        val(logical(eye(size(val)))) = diag(val).^-0.5;

        %compute the additive constant c
        B = (val*vec') * Wd * (vec*val);
        evalues = eig(B);
        evalues = real(evalues);
        min_eval = min(evalues);
        min_eval(min_eval<0 & min_eval>-0.000001) = 0;

        varargout{1} = min_eval >= 0;
        varargout{2} = -1 * min_eval;
    end
end