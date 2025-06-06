#ifndef legendre_header
#define legendre_header

//////////////////////////////////////////////////////////////
//Module that computes the Legendre polynomials and finds
//their roots
//////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>

//////////////////////////////////////////////////////////////
//Function Prototypes
//////////////////////////////////////////////////////////////
;

double ComputeLegendre(double x, int p);
double ComputeDiffLegendre(double x, int p);
double FindRootNewton(double x, int p, double tol);
#endif