Improved Relation Fuzzy c-Means
==========================================

Overview
------------------------------------------
iRFCM is an extension to the Relational Fuzzy c-Means algorithm first proposed by Hathaway and Bezdek (see [1]). RFCM expects the input D to be an Euclidean dissimilarity matrix. However, it is not always guaranteed that D is Euclidean, and if it is not Euclidean the duality relationalship between RFCM and FCM will be violated and cause RFCM to fail if the relational distances become negative. 

To overcome this problem, iRFCM will Euclideanize D is it not already Euclidean using various types of transformations that are discussed in [2-3]

NOTE: part of iRFCM, which the Subdominant Ultrametric transformation, relies on two MATLAB built-in fucntions which are part of the Bioinformatics toolbox. The first function is graphminspantree() which is to construct the minimum spanning tree (MST) from D. The second function is graphshortestpath(), which is used to traverse the path between two nodes in MST. 

Directories Included in the Toolbox
------------------------------------------
`Data/` - datasets used by the demo scripts

`Functions/` - the MATLAB functions used in iRFCM

`Results/` - the location where iRFCM toolbox stores the results

Setup
------------------------------------------
You can either download a zip file and extract it to your preferred location. Or clone this repository using git

`git clone https://github.com/mohammedkhalilia/iRFCM.git`

Then add the directory to your MATLAB path.

Known Issues
------------------------------------------
In some cases MATLAB produces complex eigenvalues and vectors in situation where it should not. That problem occurred when using double precision. Some matrices in iRFCM had to be converted to single precision to overcome this problem. Even with that sometimes the problem still occurs where the eigenvalues or vectors have zero imaginary part, in such case only the real part of the number is used.

Despite those work arounds, the iRFCM toolbox performs as expected and the results are verified with other published papers.

iRFCM Configurations
------------------------------------------
iRFCM allows the user to define their own configurations using MATLAB struct. Those configurations are explained in the `Functions/irfcm.m` function, but we will explain here as well. 
Example 4 breifly demonstrates how to define options for iRFCM. The iRFCM options are defined in a structure with the following fields/members:

`fuzzifier` - 

`epsilon` -

`maxIter` - 

`initType` - 

`gamma` - 

`delta` - 

Examples (Mutation Dataset)
-----------------------------------------

### Example 1. iRFCM with Beta-Spread

    %load the mutation dataset (for details on the Mutation dataset see ref. [4])
    %NOTE: the dissimilarities here are not squared
    D = load('Data/animal_mutation.csv');
    
    %initialize delta, delta here being the Beta-Spread
    n = size(D,1);
    delta = 1 - eye(n);
    c = 4; %run with 4 clusters

	%attach delta to a structure options that is inputted to iRFCM. If delta is provided as in this example
	%then iRFCM will test if D is Euclidean if not it will use delta to Euclidean D. If D is found to be Euclidean
	%then iRFCM will ignore delta and cluster D directly.
	options.delta = delta;
	
	%notice that the first input is the Hadamard product of D. Because pdist dissimilarities are not squared.
	out = irfcm(D.^2,c,options);

### Example 2. iRFCM with Subdominant Ultrametric

	%load the mutation dataset (for details on the Mutation dataset see ref. [4])
    %NOTE: the dissimilarities here are not squared
    D = load('Data/animal_mutation.csv');
    
    %initialize delta, delta here being the Beta-Spread
    n = size(D,1);
    
    %compute the subdominant ultrametric of D.^2 (NOT D), unless your dissimilarities are already squared. 
    delta = subdominant_ultrametric(D.^2);
    c = 4; %run with 4 clusters

	%attach delta to a structure options that is inputted to iRFCM
	options.delta = delta;
	
	%notice that the first input is the Hadamard product of D. Because pdist dissimilarities are not squared.
	out = irfcm(D.^2,c,options);

### Example 3. iRFCM without Euclideanizing D
The small code snippet below will run the Mutation dataset without the need for Euclideanizing D first. That is OK because we already know in advance that RFCM does not fail to execute on the Mutation dataset.

	D = load('Data/animal_mutation.csv');
	n = size(D,1);
	c = 4;
	out = irfcm(D.^2,c);

### Example 4. iRFCM with User Defined Configurations
iRFCM in the examples above runs with the default configurations. Those configurations can also be defined by the user

	D = load('Data/animal_mutation.csv');
	n = size(D,1);
	c = 4;
	
	%set iRFCM configurations
	options.fuzzifier        = 2;   	%default
	options.epsilon          = 0.0001; 	%tolerence for the termination criteria
	options.maxIter          = 100;	 	%number of iterations before iRFCM terminates
	options.initType         = 2;    	
	options.gamma            = 0;

	%run iRFCM with user defined values
	out = irfcm(D.^2,c, options);

References
------------------------------------------
1. R. J. Hathaway, J. W. Davenport, and J. C. Bezdek, “Relational duals of the c-means clustering algorithms,” Pattern Recognition, vol. 22, no. 2, pp. 205–212, Jan. 1989.

2. J. Benasseni, M. B. Dosse, and S. Joly, “On a General Transformation Making a Dissimilarity Matrix Euclidean,” J. Classif., vol. 24, no. 1, pp. 33–51, Jun. 2007.

3. J. Dattorro, Convex Optimization & Euclidean Distance Geometry. Meboo Publishing, 2005.

4. W. Fitch and E. Margoliash, “Construction of phylogenetic trees,” Science (80-. )., 1967.

5. E. Anderson, “The Irises of the Gaspe Peninsula,” Bull. Am. Iris Soc., vol. 59, pp. 2 – 5, 1935.

