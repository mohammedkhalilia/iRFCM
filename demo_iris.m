%clear MATLAB workspace
clear
close all

%load Iris dataset and compute the sup norm squared dissimilarity
X = load('Data/iris.csv');
D = squareform(pdist(X,'chebychev'));
n = size(D,1);

%compute the normalized dissimilarity image from D
D01 = D./max(D(:));
f = figure('Visible','off');imagesc(D.^2);colormap('gray');colorbar;
print(f, '-djpeg', 'Results/Iris/Images/Iris.jpg');

%set the number of clusters to 3
c= 3;

% Assumed ground truth
labels = [ones(1,50) 2*ones(1,50) 3*ones(1,50)];
GT = sparse(labels, 1:length(labels),1,c,length(labels));
    
% Your five choices of delta
deltas = {  'delta = 1 - eye(n);',...
            'delta = (D.^(1/5)).^2;',...
            'delta = (1-exp(-11.5.*D)).^2;',...
            'delta = (log2(1+D.^(1/4))).^2;',...
            'delta = subdominant_ultrametric(D).^2;'};

%labels assigned to every delta
deltaNames = {'beta-spread','power-fit','exp-fit','log-fit','subdominant-ultrametric'};
        
%% iRFCM configurations/options (those are the default values)
options.fuzzifier        = 2;
options.epsilon          = 0.0001;
options.maxIter          = 100;
options.initType         = 2;
options.gamma            = 0;

%Check if RFCM runs on Iris dataset dissimilarity image. You will see that
%it will return error indicating failure. 
out = irfcm(D.^2,c,options);
if isfield(out,'Error')
   fprintf('%s',out.Error); 
end

%% Since RFCM failed we need to run iRFCM
% loop for every delta
for i=1:length(deltas)
    eval(deltas{i});
    
    %set delta and run iRFCM
    options.delta = delta;
    out = irfcm(D.^2,c,options);
    
    %save the partition matrix for this delta
    U = out.U;
    dlmwrite(sprintf('Results/Iris/Partitions/U_%s(%d).csv',deltaNames{i},c),U, 'delimiter',',');

    %save the induced dissimilarity image for this delta
    %Ref. J. Huband and J. Bezdek, “VCV2– Visual cluster validity,” Comput. Intell. Res. Front., 2008.
    uu = 1 - ((U'*U)./max(max(U'*U)));
    f = figure('Visible','off');imagesc(uu);colormap('gray');caxis([0 1]);
    print(f, '-djpeg', sprintf('Results/Iris/Images/UU_%s(%d).jpg',deltaNames{i},c));
    
    %compute the crisp rand index
    [~,labels] = max(U);
    U = sparse(labels, 1:length(labels),1,c,length(labels));
    r = rand_index(U,GT,2)
end