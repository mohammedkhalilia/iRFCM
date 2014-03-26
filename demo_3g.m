%clear MATLAB workspace
clear
close all

%generate 3 gaussian clouds
a=normrnd(repmat([1,1],5,1), 1);
b=normrnd(repmat([1,5],5,1), 1);
c=normrnd(repmat([4,3],5,1), 1);
data = [a;b;c];

%compute the squared relational matrix
D = squareform(pdist(data,'euclidean').^2);
n = size(D,1);

%visualized the squared matrix
f = figure('Visible','off');imagesc(D);colormap('gray');colorbar;
print(f, '-djpeg', 'Results/3G/Images/3G_square.jpg');

%plot the 3 gaussians
f = figure('Visible','off');scatter(data(:,1),data(:,2));colormap('gray');colorbar;
print(f, '-djpeg', 'Results/3G/Images/3G_scatter.jpg');

                    
%% iRFCM configurations/options (those are the default values)
options.fuzzifier        = 2;
options.epsilon          = 0.0001;
options.maxIter          = 100;
options.initType         = 2;
options.gamma            = 0;

%set the number of clusters to 3
c= 3;

%run RFCM
out = irfcm(D,c,options);
if isfield(out,'Error')
   fprintf('%s',out.Error); 
end

U = out.U;

%visualize the partition matrix
f = figure('Visible','off');imagesc(U);colormap('gray');caxis([0 1]);
print(f, '-djpeg', 'Results/3G/Images/U(3G).jpg');

%visualize the induced dissimilarity matrix
uu = 1 - ((U'*U)./max(max(U'*U)));
f = figure('Visible','off');imagesc(uu);colormap('gray');caxis([0 1]);
print(f, '-djpeg', 'Results/3G/Images/D(U_3G).jpg');