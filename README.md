Improved Relation Fuzzy c-Means (iRFCM)

iRFCM is an extension to the Relational Fuzzy c-Means algorithm first proposed by Hathaway and Bezdek (see [1]). RFCM expects the input D to be an Euclidean dissimilarity matrix. However, it is not always guaranteed that D is Euclidean, and if it is not Euclidean the duality relationalship between RFCM and FCM will be violated and cause RFCM to fail if the relational distances become negative. 

To overcome this problem, iRFCM will Euclideanize D is it not already Euclidean using various types of transformations that are discussed in [2-3]


[1] R. J. Hathaway, J. W. Davenport, and J. C. Bezdek, “Relational duals of the c-means clustering algorithms,” Pattern Recognition, vol. 22, no. 2, pp. 205–212, Jan. 1989.

[2] J. Benasseni, M. B. Dosse, and S. Joly, “On a General Transformation Making a Dissimilarity Matrix Euclidean,” J. Classif., vol. 24, no. 1, pp. 33–51, Jun. 2007.

[3] J. Dattorro, Convex Optimization & Euclidean Distance Geometry. Meboo Publishing, 2005.
