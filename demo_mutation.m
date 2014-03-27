%clear MATLAB workspace
clear
close all

%load animal mutation dataset
load Data/animal_mutation.csv;
D = animal_mutation;
n = size(D,1);

%compute the normalized dissimilarity image from D
D01 = D./max(D(:));
f = figure('Visible','off');imagesc(D.^2);colormap('gray');
print(f, '-djpeg', 'Results/Mutation/Images/animal_mutation.jpg');

% Your five choices of delta
deltas = {  'delta = 1-eye(n);',...
            'delta = (D.^(1/2)).^2;',...
            'delta = (1-exp(-0.2.*D)).^2;',...
            'delta = (log2(1+D.^(1))).^2;',...
            'delta = subdominant_ultrametric(D).^2;'};

%set the number of clusters to 4
c = 4;

% Assumed ground truth
labels = [ones(1,17) 2 3 4];
GT = sparse(labels, 1:length(labels),1,c,length(labels));

%labels assigned to every delta
deltaNames = {'beta-spread','power-fit','exp-fit','log-fit','subdominant-ultrametric'};
        
%% iRFCM configurations/options (those are the default values)
options.fuzzifier        = 2;
options.epsilon          = 0.0001;
options.maxIter          = 100;
options.initType         = 2;
options.gamma            = 0;

%RFCM does not fail on the mutation dataset
out = irfcm(D.^2,c,options);
U = out.U;
dlmwrite(sprintf('Results/Mutation/Partitions/U-RFCM(%d).csv',c),U, 'delimiter',',');

uu = 1 - ((U'*U)./max(max(U'*U)));
f = figure('Visible','off');imagesc(uu);colormap('gray');caxis([0 1]);
print(f, '-djpeg', sprintf('Results/Mutation/Images/UU-RFCM(%d).jpg',c)); 

%compute the crisp rand index
[~,labels] = max(U);
U = sparse(labels, 1:length(labels),1,c,length(labels));
r = rand_index(U,GT,2)

%% loop for every delta
for i=1:length(deltas)
    eval(deltas{i});
    
    %set delta and run iRFCM
    options.delta = delta;
    out = irfcm(D.^2,c,options);

    %save the partition matrix for this delta
    U = out.U;
    dlmwrite(sprintf('Results/Mutation/Partitions/U-%s(%d).csv',deltaNames{i},c),U, 'delimiter',',');

    %save the induced dissimilarity image for this delta
    %Ref. J. Huband and J. Bezdek, “VCV2– Visual cluster validity,” Comput. Intell. Res. Front., 2008.
    uu = 1 - ((U'*U)./max(max(U'*U)));
    f = figure('Visible','off');imagesc(uu);colormap('gray');caxis([0 1]);
    print(f, '-djpeg', sprintf('Results/Mutation/Images/UU-%s(%d).jpg',deltaNames{i},c));
    
    %compute the crisp rand index
    [~,labels] = max(U);
    U = sparse(labels, 1:length(labels),1,c,length(labels));
    r = rand_index(U,GT,2)
end