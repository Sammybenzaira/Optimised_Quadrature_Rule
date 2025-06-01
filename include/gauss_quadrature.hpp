#ifndef gauss_quadrature_header
#define gauss_quadrature_header

//////////////////////////////////////////////////////////////
//Module that constructs Gauss Quadrature on various domains
//////////////////////////////////////////////////////////////

#include "polygon.hpp"

//////////////////////////////////////////////////////////////
//Function Prototypes
//////////////////////////////////////////////////////////////

double** ConstructGaussQuad(int n);
double** ConstructSquareGaussQuad(int n);
double** ConstructTriangleGaussQuad(int n);
double** TransformTriangleQuad(Triangle Tri, int n);
double** ConstructPolygonQuad(GeneralPolygon Omega, int n);

#endif