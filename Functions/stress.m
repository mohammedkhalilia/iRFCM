function varargout = stress(D, Dh, varargin)
%**************************************************************************
% Usage: [s1,ns1,e1,e2,t] = stress(D, Dh, 'stress1', 'nstress1', 'eps1', 'eps2', 'tau')
%  
% Arguments:
%   D: original dissimilarity matrix
%   Dh: tranformed dissimilarity matrix. 
%       Note that if you are computing the
%       stress Q1 mentioned in the paper titled "On a General Transformation 
%       Making a Dissimilarity Matrix Euclidean" Dh should be multipled by eps1
%
%   Additional Arguments:
%       At least you should provide one of the following values:
%       stress1, nstress1, eps1, eps2, tau
%   The function return arguments will correspond to the additional
%   arguments provided.
%
% Note: when computing the stress some approaches do the double summations
% go for all i < j (upper tringular part of the matrix) other approaches
% use all values. This code implements the approach where i < j
%**************************************************************************

    numarg = nargin - 2;
    varargout = cell(1,numarg);
    
    for i=1:numarg
        
       switch varargin{i}
           
            case 'nstress1' %normalized stress1, between 0 and 1
                D = triu(D);
                Dh = triu(Dh);
                varargout{i} = sqrt(sum(sum((D-Dh).^2))/sum(sum(D.^2)));
            
            case 'stress1' %raw stress1, not normalized
                D = triu(D);
                Dh = triu(Dh);
                varargout{i} = sum(sum((D-Dh).^2));
           
            case 'eps1' %epsilon 1
                varargout{i} = sum(sum(D.*Dh))./sum(sum(Dh.^2));
                
            case 'eps2'  %epsilon 2
                Dh = Dh.^2;
                D = D.^2;
                varargout{i} =  sqrt(sum(sum(D.*Dh)))./max(sum(Dh));
                
            case 'tau'
                %ignore off diagonal elements
                off_diag = find(eye(size(D)) == 0);
                tmp = Dh(off_diag)./D(off_diag);
                varargout{i} = max(tmp(:))/min(tmp(:));
             
           otherwise
               error('Incorrect option');
       end
       
    end

end