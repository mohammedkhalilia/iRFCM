function output=irfcm(R, c, options)
%% 
%   A new approach to handle dissimilarity matrix tranformation to
%   Euclidean
%
%   Allow the algorithm to either start by initializing U or V
%   NERF options are not passed in a struct "options"
%
% Purpose:The function nerfcm attempts to find a useful clustering of the
%		objects represented by the dissimilarity data in R using the initial
% 		partition	in U0.
%       This will center the matrix D according to
%               B = (I - 11') * (-0.5*D) * (I - 11')
%                   OR if D is the squared distances, then
%               B = (I - 11') * (-0.5*D^2) * (I - 11')
%
%               But you are always assumed to have provided a squared
%               distances.
%
%               Then check if the matrix positive semi deifnite and return the
%               smallest eigenvalue as the constant c (see [1-2])
%
% Usage: [U,beta,num_it]=nerfcm( R,c,options)
%   options is a struct with the following default values:
%       Fuzzifier        = 2;
%       Epsilon          = 0.001;
%       MaxIter          = 100; 
%       InitType         = 2;   
%       AdditiveConstant = 0;
%       Delta            = [];
%
% See detailed explanation of every field below
%
% output    - structure containing:
%               U: fuzzy partition
%               V: cluster centers/coefficients
%               TerminationIter: the number of iterations at termination
%               MaxIter: maximum number iterations allowed
%               Beta: the final value of Beta
%               BetaCount: number of times the Beta block got activiated
%             If you decide to Euclideniate R before running NERF, then the
%             following information will be returned as well
%               KruskalStress: the stress value measured the tranformed R
%               eps: the value that minimized the different between R and D
%               c: the smallest constant to add to R to make it Euclidean
%               D: The Euclideaniated matrix
%               TestType: PSD or NSD (positive/negative semi-definite)
%               Delta
%
% R         - the dissimilarity data matrix of size n x n
% c         - number of clusters
% Fuzzifier - fuzzifier, default 2
% Epsilon   - converage, default 0.0001
% Init      - can take two values: 
%               'U' which initializes the membership matrix 
%               'V' which initializes the cluster centers
% InitType  - this field is tied to the Init value as follows:
%               for 'U' InitType can take either 1 or 2
%               1 = assign first 1/c objects to cluster 1
%                   assign second 1/c objects to cluster 2, etc
%               2 = random initialization
%
%               for 'V'
%               1 = random initialization
%               2 = randomly choose c rows from D
% MaxIter     - the maximum number fo iterations, default 100
% Delta     - delta is the matrix that is used to tranform R to Euclidean.
%             delta is of size n x n and can be for instance
%                  delta = 1-eye(n)  this is the default delta
%                  delta = subdominant ultrametric matrix
%                  delta = R
%                  delta = power of R
%                  delta = log2(1+D), see [5], page 485
%                  delta = parametric function
%              delta is expected to have an Euclidean representation,
%              otherwise an error will be thrown (see [4])
%              Also, providing delta will override the tranform option if
%              it set to 0.
% 
% Refs:
%   [1] T. Cox and M. Cox, Multidimensional scaling. 2010.
%   [2] K. Mardia and J. Kent, “y BIBBY, JM (1979): Multivariate analysis,” Academic Press, London.
%   [3] R. J. Hathaway and J. C. Bezdek, “Nerf c-means: Non-Euclidean 
%       relational fuzzy clustering,” Pattern Recognition, vol. 27, no. 3, pp. 429–437, Mar. 1994.
%   [4] J. Benasseni, M. B. Dosse, and S. Joly, “On a General Transformation Making a Dissimilarity 
%       Matrix Euclidean,” Journal of Classification, vol. 24, no. 1, pp. 33–51, Jun. 2007.
%   [5] J. Dattorro, Convex optimization and Euclidean distance geometry. 2005.

    %% NERF default values
    m = 2; epsilon = 0.0001;max_it = 100;
    delta = [];ac = 0;init_type = 2;
    
    %% Overwrite NERF options by the user defined options
    if nargin == 3 && isstruct(options)
        fields = fieldnames(options);
        for i=1:length(fields)
           switch fields{i}
               case 'Fuzzifier', m = options.Fuzzifier;
               case 'Epsilon', epsilon = options.Epsilon; 
               case 'InitType', init_type = options.InitType; 
               case 'Delta', delta = options.Delta; 
               case 'AdditiveConstant', ac = options.AdditiveConstant; 
               case 'MaxIter', max_it = options.MaxIter; 
           end
        end
    end

    %% Initialize variables
    D = R;n=size(D,1);d = zeros(c,n);
    num_it=0;step_size=epsilon;U=Inf(c,n);
    
    %If transform option is provided, but delta is not, then do the
    %transformation based on default option, delta = 11' - I
    if isfield(options,'Delta')
        [D, b] = make_euclidean(R, delta, ac, 'Dh = (D + c*delta);');
        eps1 = stress(R, D, 'eps1');
        s1 = stress(R, D, 'stress1');
        euc = struct('KruskalStress',s1,'eps',eps1,'c',b,'D',D,'Delta',delta);
    end
    
    %assign the first psize points to cluster 1
    % second psize points to cluster 2, etc.
    V = init_centers(init_type, n, c, D);
    
    if min(min(D)) < 0 || any(any(abs(D - D') > 0)) || max(diag(D)) > 0
		error('R is not properly initialized.')
	elseif size(D) ~= [n,n]
		error('Dimensions R are inconsistent.')
    end
    
    %% Begin the main loop:
    while  num_it < max_it && step_size >= epsilon
        U0 = U;
        
        %Get new (quasi-squared-distance values) d:
        % First, get new initial values for d:
        for i=1:c
            d(i,:)=D*V(i,:)'-V(i,:)*D*V(i,:)'/2;
        end
       
        j = find(d(:) < 0);
        if ~isempty(j)
           output.Error = sprintf('RFCM encountered %d negative relational distances in iteration %d. RFCM terminated execuation.\nPlease re-run iRFCM and provide Delta to Euclideanize D before clustering\n\n', length(j), num_it); 
           return;
        end
        
        %Get new partition matrix U:
        d=d.^(1/(m-1));
		work = 1 ./ d;
		work = sum(work);
        U=1./d;
        U = U./(ones(c,1)*work);
		
        %Get new V prototypes assuming we initialized V first
        V=U.^m;  
        V = V./(sum(V,2) * ones(1,n));
    
        %Calculate step_size and return to top of loop:
		step_size=max(max(abs(U-U0)));
        
        num_it = num_it + 1;
    end
    
    %prepare output structure
    output = struct('U',U,...
                    'V',V,...
                    'TerminationIter',num_it,...
                    'MaxIter',max_it);
                
    if exist('euc','var'),output.Euc = euc;end
    if nargin == 3,output.Options = options;end
end