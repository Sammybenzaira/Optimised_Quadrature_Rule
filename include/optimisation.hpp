#ifndef optimisation_header
#define optimisation_header

//////////////////////////////////////////////////////////////
//Module that produces an optimised quadrature
//////////////////////////////////////////////////////////////


#include <vector>
#include "polygon.hpp"

//////////////////////////////////////////////////////////////
//Function Prototypes
//////////////////////////////////////////////////////////////


double* ConstructOptimisation(DoubleMatrix A, double* b, double tol);
DoubleMatrix ConstructOptimisedQuad(GeneralPolygon Omega, int degree);
void ConstructMatrixAVectorB(double** quadPoints, int degree, DoubleMatrix& A, double* b, int numberOfQuadPoints);

#endif