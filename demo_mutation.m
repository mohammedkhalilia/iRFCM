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
        
%labels assigned to every delta
delta_names = {'beta-spread','power-fit','exp-fit','log-fit','subdominant-ultrametric'};
        
%% iRFCM configurations/options (those are the default values)
options.Fuzzifier        = 2;
options.Epsilon          = 0.0001;
options.MaxIter          = 100;
options.InitType         = 2;
options.AdditiveConstant = 0;
options.Delta            = [];

%set the number of clusters to 4
c = 4;

%% loop for every delta
for i=1:length(deltas)
    eval(deltas{i});
    
    %set delta and run iRFCM
    options.Delta = delta;
    out = irfcm(D.^2,c,options);

    %save the partition matrix for this delta
    U = out.U;
    dlmwrite(sprintf('Results/Mutation/Partitions/U-%s(%d).csv',delta_names{i},c),U, 'delimiter',',');

    %save the induced dissimilarity image for this delta
    %Ref. J. Huband and J. Bezdek, “VCV2– Visual cluster validity,” Comput. Intell. Res. Front., 2008.
    uu = 1 - ((U'*U)./max(max(U'*U)));
    f = figure('Visible','off');imagesc(uu);colormap('gray');caxis([0 1]);
    print(f, '-djpeg', sprintf('Results/Mutation/Images/UU-%s(%d).jpg',delta_names{i},c)); 
end