function sdu = subdominant_ultrametric(D,varargout)
%%
%
% Computes the Sub-dominant ultrametric matrix for a given dissimilarity
% matrix D based on its minimum spanning tree MST. If MST argument is not
% provided it will be computed assuming the MATLAB function graphminspantree
% exists. You need to have the Bioinformatics toolbox for that function.
%
% Usage sdu = subdominant_ultrametric(D,MST)
%
% sdu   - the computed Sub-dominant Ultrametric matrix
% D     - n x n dissimilarity matrix
% MST   - minimum spanning tree as returned by the MATLAB function 
%         graphminspantree

    %if MST not provided, compute it
    if nargin == 1
        MST = graphminspantree(sparse(D));
    else
        MST = varargout{1};
    end
    
    MST = MST + MST';
    n = size(D,1);
    sdu = zeros(n);

    for i=1:n
        for j=i+1:n
            %for every two nodes get the path between them
            [~, path, ~] = graphshortestpath(MST,i,j);
            path_len = length(path) - 1;
            tmp = zeros(1,path_len);
            
            %use the path to find the length of the edges
            for k=1:path_len
                node = path(k:k+1);
                tmp(k) = D(node(1),node(2));
            end
            
            %find the longest edge, the largest distance
            sdu(i,j) = max(tmp);
        end
    end
    
    sdu = sdu + sdu';
end