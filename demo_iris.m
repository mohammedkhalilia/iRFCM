%% Mutation dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
X = load('Data/iris.csv');
D = squareform(pdist(X,'chebychev'));
n = size(D,1);

%% image of the original D^2
D01 = D./max(D(:));
f = figure('Visible','off');imagesc(D.^2);colormap('gray');colorbar;
print(f, '-djpeg', 'Results/Iris/Images/Iris.jpg');

%% Different choices of Delta
% Delta is squared because later it will be added to the squared
% dissimilarities to generate a squared Euclidean distances
deltas = {  'delta = 1 - eye(n);',...
            'delta = (D.^(1/5)).^2;',...
            'delta = (1-exp(-11.4.*D)).^2;',...
            'delta = (log2(1+D.^(1/4))).^2;',...
            'delta = subdominant_ultrametric(D).^2;'};

delta_names = {'beta-spread','power-fit','exp-fit','log-fit','subdominant-ultrametric'};
        
%% NERFCM configurations/options
options.Fuzzifier        = 2;
options.Epsilon          = 0.0001;
options.MaxIter          = 100;
options.InitType         = 2;
options.AdditiveConstant = 0;

c= 3;

out = irfcm(D.^2,c,options);
if isfield(out,'Error')
   fprintf('%s',out.Error); 
end

%% Test the RFCM, Beta-Spread and Ultrametric
for i=1:length(deltas)
    eval(deltas{i});
    options.Delta = delta;
    out = irfcm(D.^2,c,options);
    
    %Save partition matrix
    U = out.U;
    dlmwrite(sprintf('Results/Iris/Partitions/U_%s(%d).csv',delta_names{i},c),U, 'delimiter',',');

    %compute the induced dissimilarity matrix
    uu = 1 - ((U'*U)./max(max(U'*U)));
    f = figure('Visible','off');imagesc(uu);colormap('gray');caxis([0 1]);
    print(f, '-djpeg', sprintf('Results/Iris/Images/UU_%s(%d).jpg',delta_names{i},c));
end