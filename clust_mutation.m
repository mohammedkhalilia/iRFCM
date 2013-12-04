%% Mutation dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
load Data/animal_mutation.csv;
D = animal_mutation;
n = size(D,1);

D01 = D./max(D(:));
f = figure('Visible','off');imagesc(D.^2);colormap('gray');
print(f, '-djpeg', 'Results/Mutation/Images/animal_mutation.jpg');

% Delta is squared because later it will be added to the squared
% dissimilarities to generate a squared Euclidean distances
deltas = {  'delta = 1-eye(n);',...
            'delta = (D.^(1/2)).^2;',...
            'delta = (1-exp(-0.2.*D)).^2;',...
            'delta = (log2(1+D.^(1))).^2;',...
            'delta = subdominant_ultrametric(D).^2;'};
        
delta_names = {'beta-spread','power-fit','exp-fit','log-fit','subdominant-ultrametric'};
        
%% NERFCM configurations/options
options.Fuzzifier        = 2;
options.Epsilon          = 0.0001;
options.MaxIter          = 100;
options.InitType         = 2;
options.AdditiveConstant = 0;
options.Delta            = [];

c = 4;

for i=1:length(deltas)
    eval(deltas{i});
    options.Delta = delta;
    out = irfcm(D.^2,c,options);
    s = out.Euc.KruskalStress;
    ac = out.Euc.c;
    eps1 = out.Euc.eps;
        
    U = out.U;
    dlmwrite(sprintf('Results/Mutation/Partitions/U-%s(%d).csv',delta_names{i},c),U, 'delimiter',',');

    uu = 1 - ((U'*U)./max(max(U'*U)));
    f = figure('Visible','off');imagesc(uu);colormap('gray');caxis([0 1]);
    print(f, '-djpeg', sprintf('Results/Mutation/Images/UU-%s(%d).jpg',delta_names{i},c)); 
end