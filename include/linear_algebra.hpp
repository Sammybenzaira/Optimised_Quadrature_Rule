#ifndef liner_algebra_header
#define liner_algebra_header

//////////////////////////////////////////////////////////////
//Module that implements some linear algebra routines
//////////////////////////////////////////////////////////////

#include <iostream>
#include <math.h>
#include <fstream>

//////////////////////////////////////////////////////////////
//Matrix Struct Definitions
//////////////////////////////////////////////////////////////

//General double valued matrix struct
struct DoubleMatrix{
    double** matrix_entries; //Dynamically allocated matrix entries
    int no_rows = 0; //No Rows
    int no_cols = 0; //No Cols
};

//Symmetric (square) double valued matrix struct
struct SymmetricMatrix{
    double** matrix_entries;
    int no_rows = 0;
};

//////////////////////////////////////////////////////////////
//Function Prototypes
//////////////////////////////////////////////////////////////

///Vector Functions///
double* AllocateVector(int n);
void ZeroVector(double* vector, int n);
double* PrintVector(double* vector,int n);
void CopyVector(double* vector,double* copied_vector,int n);
void SubtractVectors(double* vec1,double* vec2,int n);
void CombineVectors(double* vector1,double* vector2,double scale,int n);
void ScaleVector(double* vec,double scaleFactor,int n);
double NormVector(double* vector,int n);
void FindMaximum(double* vector,int n, double& max_val, int& index);
void FindMinimum(double* vector,int n, double& min_val, int& index);
void DeallocateVector(double* vector);

///Matrix Functions///
DoubleMatrix AllocateDoubleMatrix(int noRows,int noCols);
SymmetricMatrix AllocateSymmetricMatrix(int noRows);
void DeallocateMatrix(DoubleMatrix& matrix);
void DeallocateMatrix(SymmetricMatrix& matrix);
void PrintMatrix(DoubleMatrix);
void PrintMatrix(SymmetricMatrix matrix);
void MultiplyVectorByMatrix(DoubleMatrix& matrix, double* vec, double* product);
void PerformSymmetricBackSubstitution(SymmetricMatrix& matrix, double* rhs,
	double* solution);
void PerformSymmetricForwardSubstitution(SymmetricMatrix& matrix, double* rhs,
	double* solution);
void ComputeLDLFactorisation(SymmetricMatrix& matrix);
void PerformSymmetricSolve(SymmetricMatrix& matrix, double* rhs, double* solution);
void TransposeMatrix(DoubleMatrix& matrix, DoubleMatrix& transposeMatrix);
DoubleMatrix ConvertSymmetricMatrixToDouble(SymmetricMatrix matrix);
void SolveRestrictedProblem(DoubleMatrix& A, double* b, double* s_p);
SymmetricMatrix ConstructSymmetricMatrix(DoubleMatrix matrix);
void ConstructVectorV(DoubleMatrix A, double* b, double* x, double* v);
void PrintNormBMinusAX(DoubleMatrix& A, double* b, double* x);
#endif