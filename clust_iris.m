%% Mutation dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
X = load('Data/iris.csv');
D = squareform(pdist(X,'chebychev'));

%% iVAT of D^2
R = ivat(D.^2);
f = figure('Visible','off');imagesc(rot90(R,2));colormap('gray');
print(f, '-djpeg', 'Results/Iris/Images/Iris_iVAT.jpg');

%% VAT of D^2
R = vat(D.^2);
f = figure('Visible','off');imagesc(rot90(R,2));colormap('gray');
print(f, '-djpeg', 'Results/Iris/Images/Iris_VAT.jpg');

%% image of the original D^2
D01 = D./max(D(:));
f = figure('Visible','off');imagesc(D.^2);colormap('gray');colorbar;
print(f, '-djpeg', 'Results/Iris/Images/Iris.jpg');

%% Collect some statistics about D^2
fh = fopen('Results/Iris/Iris_stats.csv','w');
fprintf(fh,'Max,%f\n',max(D(:).^2));
fprintf(fh,'Mean,%f\n',mean(D(:).^2));
fprintf(fh,'STD,%f\n',std(D(:).^2));
fprintf(fh,'Variance,%f\n',var(D(:).^2));
fclose(fh);

fh = fopen('Results/Iris/Iris_results.csv','w');
fprintf(fh,'Iterations,Beta,BetaCount,Delta,Additive Constant,Kruskal Stress,Pearson Correlation(D;D(U))\n');

%% Different choices of Delta
% Delta is squared because later it will be added to the squared
% dissimilarities to generate a squared Euclidean distances
deltas = {  'RFCM',...
            'delta=[];',...
            'delta = (D.^(1/5)).^2;',...
            'delta = (1-exp(-11.4.*D)).^2;',...
            'delta = (log2(1+D.^(1/4))).^2;',...
            'delta = subdominant_ultrametric(D).^2;'};

%% NERFCM configurations/options
options.Fuzzifier        = 2;
options.Epsilon          = 0.0001;
options.MaxIter          = 100;
options.Init             = 'V';
options.InitType         = 2;
options.Euclideanize     = 0;
options.AdditiveConstant = 0;
options.Delta            = [];
options.EucType          = 'PSD';

%% Test the RFCM, Beta-Spread and Ultrametric
for i=1:length(deltas)
    prefix = strrep(deltas{i},'/','div');
    
    if i==1
        out = nerfcm(D.^2,3, options);
        s = 0;
        ac = 0;
        eps1 = 0;
    else    
        eval(deltas{i});
        euc = 1;
        if i > 2, euc = is_euclidean(delta);end
            
        %Is delta Euclidean?
        if euc
            options.Delta = delta;
            options.Euclideanize = 1;
            out = nerfcm(D.^2,3,options);
            s = out.Euc.KruskalStress;
            ac = out.Euc.c;
            eps1 = out.Euc.eps;
        else
            fprintf('Delta is not Euclidean: %s\n',deltas{i});
            continue
        end
        
        %Collect some statistics about the Euclideanized D
        fh1 = fopen(sprintf('Results/Iris/Iris_stats_%s.csv',prefix),'w');
        fprintf(fh1,'Max,%f\n',max(out.Euc.D(:)));
        fprintf(fh1,'Min,%f\n',min(out.Euc.D(:)));
        fprintf(fh1,'Mean,%f\n',mean(out.Euc.D(:)));
        fprintf(fh1,'STD,%f\n',std(out.Euc.D(:)));
        fprintf(fh1,'Variance,%f\n',var(out.Euc.D(:)));
        fclose(fh1);
    end
    
    if i > 1
        %show the iVAT and image of the Euclideanized D
        R = ivat(out.Euc.D);
        f = figure('Visible','off');imagesc(rot90(R,2));colormap('gray');colorbar;
        print(f, '-djpeg', sprintf('Results/Iris/Images/Iris_iVAT_%s.jpg',prefix));    
        
        R = vat(out.Euc.D);
        f = figure('Visible','off');imagesc(rot90(R,2));colormap('gray');colorbar;
        print(f, '-djpeg', sprintf('Results/Iris/Images/Iris_VAT_%s.jpg',prefix));  
        
        %show the Euclideanized matrix
        f = figure('Visible','off');imagesc(out.Euc.D);colormap('gray');colorbar;
        print(f, '-djpeg', sprintf('Results/Iris/Images/Iris_D_%s.jpg',prefix));    
    end
    
    %Save partition matrix
    U = out.U;
    dlmwrite(sprintf('Results/Iris/Partitions/U_%s.csv',prefix),U, 'delimiter',',');

    %compute the induced dissimilarity matrix
    uu = 1 - ((U'*U)./max(max(U'*U)));
    f = figure('Visible','off');imagesc(uu);colormap('gray');caxis([0 1]);
    print(f, '-djpeg', sprintf('Results/Iris/Images/Iris_UU_%s.jpg',prefix));

    corr_coeff = corr(D01(:),uu(:),'type','pearson');
        
    fprintf(fh,'%d,%f,%d,%s,%f,%f,%f\n',out.TerminationIter,out.Beta,out.BetaCount,deltas{i},ac,s,corr_coeff);
end

fclose(fh);