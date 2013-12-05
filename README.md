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

Examples
-----------------------------------------

### Example 1. Iris Dataset

### Example 2. Mutation Dataset

References
------------------------------------------
1. R. J. Hathaway, J. W. Davenport, and J. C. Bezdek, “Relational duals of the c-means clustering algorithms,” Pattern Recognition, vol. 22, no. 2, pp. 205–212, Jan. 1989.

2. J. Benasseni, M. B. Dosse, and S. Joly, “On a General Transformation Making a Dissimilarity Matrix Euclidean,” J. Classif., vol. 24, no. 1, pp. 33–51, Jun. 2007.

3. J. Dattorro, Convex Optimization & Euclidean Distance Geometry. Meboo Publishing, 2005.
