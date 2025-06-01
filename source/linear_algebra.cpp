//////////////////////////////////////////////////////////////
//Module that implements some linear algebra routines
//////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <fstream>
#include <cassert>

#include "linear_algebra.hpp"
#include "legendre.hpp"

//////////////////////////////////////////////////////////////
/// Vector Functions
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
double* AllocateVector(int n)
//Allocates a double precision vector of length n
//and sets all entries to zero
{
	double* vector = new double[n];
	for (int i=0;i<n;i++)
	{
		vector[i] = 0.0;
	}

	return vector;
}

//////////////////////////////////////////////////////////////
void ZeroVector(double* vector, int n)
//Sets all entries in a vector of length n to zero
{
	for (int i=0;i<n;i++)
	{
		vector[i] = 0.0;
	}
}

//////////////////////////////////////////////////////////////
double* PrintVector(double* vector,int n)
//Prints a vector of length n to the screen
{
	for (int i=0;i<n;i++)
	{
		std::cout << vector[i] << std::endl;
	}

	return vector;
}

//////////////////////////////////////////////////////////////
void CopyVector(double* vector,double* copied_vector,int n)
//Makes a copy of a vector of length n- assumes that the memory for the
//copy has already been allocated
{
	for (int k=0;k<n;k++)
	{
		copied_vector[k] = vector[k];
	}
}
//////////////////////////////////////////////////////////////
void SubtractVectors(double* vec1,double* vec2,int n)
//Overwrites vec1 with vec1-vec2 (both of length n)
{
    for (int k=0;k<n;k++)
    {
        vec1[k] -= vec2[k];
    }
}
//////////////////////////////////////////////////////////////
void ScaleVector(double* vec,double scaleFactor,int n)
//Overwrites vec of length n with scaleFactor*vec
{
    for (int k=0;k<n;k++)
    {
        vec[k] *= scaleFactor;
    }
}
//////////////////////////////////////////////////////////////
void CombineVectors(double* vector1,double* vector2,double scale,int n)
//Carries out the operation x <- x+scale*y, where x = vector1, y = vector2
//and both of length n
{
	for (int k=0;k<n;k++)
	{
		vector1[k] += scale*vector2[k];
	}
}

//////////////////////////////////////////////////////////////
double NormVector(double* vector,int n)
//Compute the l2-norm of a vector of length n
{
	double norm = 0.0;
	for (int k=0;k<n;k++)
	{
		norm += pow(vector[k],2);
	}

	return sqrt(norm);
}

//////////////////////////////////////////////////////////////
void FindMaximum(double* vector,int n, double& max_val, int& index)
//Function to find maximum value of a vector of length n.
//Outputs:
// max_val - the maximum value
// index - index of the maximum value
{
	max_val = vector[0];
	index = 0;
    for (int k=1;k<n;k++)
    {
        if (vector[k] > max_val)
        {   
            max_val = vector[k];
			index = k;
        }
    }
}

//////////////////////////////////////////////////////////////
void FindMinimum(double* vector,int n, double& min_val, int& index)
//Function to find minimum value of a vector of length n.
//Outputs:
// min_val - the minimum value
// index - index of the minimum value
{
	min_val = vector[0];
	index = 0;
    for (int k=1;k<n;k++)
    {
        if (vector[k] < min_val)
        {   
            min_val = vector[k];
			index = k;
        }
    }
}

//////////////////////////////////////////////////////////////
void DeallocateVector(double* vector)
//Deletes storage for a vector
{
	delete[] vector;
}

//////////////////////////////////////////////////////////////
/// Matrix Functions
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
DoubleMatrix AllocateDoubleMatrix(int noRows,int noCols)
//Allocates a rectangular matrix of type DoubleMatrix
//and sets the entries to 0.
{
	DoubleMatrix matrix;

	matrix.no_rows = noRows;
	matrix.no_cols = noCols;

	matrix.matrix_entries = new double*[noRows];
	for (int i=0;i<noRows;i++)
	{
		matrix.matrix_entries[i] = new double[noCols];
		for (int j=0;j<noCols;j++)
		{
			matrix.matrix_entries[i][j] = 0.0;
		}
	}

	return matrix;
}

//////////////////////////////////////////////////////////////
SymmetricMatrix AllocateSymmetricMatrix(int noRows)
//Allocates a double precision symmetric matrix and sets all entries to zero
//Stored as a lower triangular matrix
{
	SymmetricMatrix matrix;

	matrix.no_rows = noRows;
	matrix.matrix_entries = new double*[noRows];

	for (int k=0;k<noRows;k++)
	{
		matrix.matrix_entries[k] = new double[k+1];
		for (int j=0;j<=k;j++)
		{
			matrix.matrix_entries[k][j] = 0.0;
		}
	}

	return matrix;

}

//////////////////////////////////////////////////////////////
void DeallocateMatrix(SymmetricMatrix& matrix)
//Deletes a matrix of type SymmetrixMatrix
{
	for (int i=0;i<matrix.no_rows;i++)
	{
		delete[] matrix.matrix_entries[i];
	}

	delete[] matrix.matrix_entries;

	matrix.no_rows = 0;
}

//////////////////////////////////////////////////////////////
void DeallocateMatrix(DoubleMatrix& matrix)
//Deletes a matrix of type DoubleMatrix
{

	for (int i=0;i<matrix.no_rows;i++)
	{
		delete[] matrix.matrix_entries[i];
	}

	delete[] matrix.matrix_entries;

	matrix.no_rows = 0;
	matrix.no_cols = 0;
}

//////////////////////////////////////////////////////////////
void PrintMatrix(SymmetricMatrix matrix)
//Prints a double precision matrix to the screen
{
	for (int i=0;i<matrix.no_rows;i++)
	{
		for (int j=0;j<=i;j++)
		{
			std::cout << matrix.matrix_entries[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

//////////////////////////////////////////////////////////////
void PrintMatrix(DoubleMatrix matrix)
//Prints a double type matrix to the screen
{

	for (int i=0;i<matrix.no_rows;i++)
	{
		for (int j=0;j<matrix.no_cols;j++)
		{
			std::cout << matrix.matrix_entries[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

//////////////////////////////////////////////////////////////
void PerformSymmetricForwardSubstitution(SymmetricMatrix& matrix, double* rhs,
	double* solution)
//For a lower triangular matrix, performs forward subsitition to solve Ax = b
//Assumes there are 1s on the main diagonal
{
	solution[0] = rhs[0];
	for (int i=1;i<matrix.no_rows;i++)
	{
		solution[i] = rhs[i];
		for (int j=0;j<i;j++)
		{
			solution[i] -= matrix.matrix_entries[i][j]*solution[j];
		}
	}
}

//////////////////////////////////////////////////////////////
void PerformSymmetricBackSubstitution(SymmetricMatrix& matrix, double* rhs,
	double* solution)
//For a lower triangular matrix, performs back subsitition to solve A^Tx = b
//Assumes there are 1s on the main diagonal
{
	solution[matrix.no_rows-1] = rhs[matrix.no_rows-1];
	for (int i=matrix.no_rows-1;i>=0;i--)
	{
		solution[i] = rhs[i];
		for (int j=i+1;j<matrix.no_rows;j++)
		{
			solution[i] -= matrix.matrix_entries[j][i]*solution[j];
		}
	}
}

//////////////////////////////////////////////////////////////
void ComputeLDLFactorisation(SymmetricMatrix& matrix)
//Finds the LDLT factorisation of a symmetrix matrix and
//overwrites the lower half of the matrix
{
	for (int i=1;i<matrix.no_rows;i++)
	{
		
		for (int j=0;j<i;j++)
		{
			for (int k=0;k<j;k++)
			{
				matrix.matrix_entries[i][j] -= 
					matrix.matrix_entries[i][k]*matrix.matrix_entries[j][k]*
					matrix.matrix_entries[k][k];
			}

			matrix.matrix_entries[i][j] /= matrix.matrix_entries[j][j];

		}

		for (int k=0;k<i;k++) //Find diagonal terms
		{

			matrix.matrix_entries[i][i] -= 
				matrix.matrix_entries[i][k]*matrix.matrix_entries[i][k]*
				matrix.matrix_entries[k][k];

		}

		if (std::fabs(matrix.matrix_entries[i][i]) < 1.0-12)
		{
			std::cout << "Warning: LDL Factorisation not possible" << std::endl;
			break;
		}

	}
}

//////////////////////////////////////////////////////////////
void PerformSymmetricSolve(SymmetricMatrix& matrix, double* rhs, double* solution)
//Solves the problem Ax = b, where A is a general symmetrix matrix
{
	//Compute LDL factorisation - in place

	ComputeLDLFactorisation(matrix);

	
	//Forward substitution

	PerformSymmetricForwardSubstitution(matrix,rhs,solution);

	

	//Diagonal Division

	for (int i=0;i<matrix.no_rows;i++)
	{
		rhs[i] = solution[i]/matrix.matrix_entries[i][i];
	}
	

	//Backward Substitution
	PerformSymmetricBackSubstitution(matrix,rhs,solution);
	


}

//////////////////////////////////////////////////////////////
void MultiplyVectorByMatrix(DoubleMatrix& matrix, double* vec, double* product)
//Computes Ax, where A is a general double precision matrix
{
	int no_rows = matrix.no_rows;
	int no_cols = matrix.no_cols;

	for (int i=0;i<no_rows;i++)
	{
		product[i] = 0.0;
		for (int j=0;j<no_cols;j++)
		{
			product[i] += matrix.matrix_entries[i][j]*vec[j];
		}
	}
}

/////////////////////////////////////////////////////////////
DoubleMatrix ConvertSymmetricMatrixToDouble(SymmetricMatrix matrix)
//Converts a symmetric matrix into a type double matrix
{
	
	//Initialise
	DoubleMatrix new_matrix = AllocateDoubleMatrix(matrix.no_rows, matrix.no_rows);
	
	
	for (int i=0; i<matrix.no_rows; i++)
	{
		for (int j=0; j<matrix.no_rows; j++)
		{   
			if (j<=i)
			//Fill the entries as normal if it is in the lower triangular portion
			{
				new_matrix.matrix_entries[i][j] = matrix.matrix_entries[i][j];
			}
			else
			//Fills the entrys of the upper triangluar portion by finding the transpose
			{
				new_matrix.matrix_entries[i][j]= matrix.matrix_entries[j][i];
			}
		}
	}
	
    return new_matrix;
}

////////////////////////////////////////////////////////////
SymmetricMatrix ConstructSymmetricMatrix(DoubleMatrix matrix)
//Constructs the Symmetric matrix A^tA for a given transpose matrix A^t
//Outputs the matrix as type SymmetricMatrix
{
	//Initialise
	SymmetricMatrix new_matrix = AllocateSymmetricMatrix(matrix.no_rows);
	
	//We loop over the new matrix adding new entries
	for (int i=0; i<matrix.no_rows; i++)
	{
		for (int j=0; j<=i; j++)
		{   
			for (int k=0; k<matrix.no_cols; k++)
			{
				//We only use the transpose to calculate since it is more efficient going through columns first then the rows
				new_matrix.matrix_entries[i][j] += matrix.matrix_entries[i][k]*matrix.matrix_entries[j][k];
			}
			
			
		}
	}
	
    return new_matrix;
}

/////////////////////////////////////////////////////////////
void TransposeMatrix(DoubleMatrix& matrix, DoubleMatrix& transposeMatrix)
//Computes the transpose of a given matrix
//Assumes the memory has already been allocated
{

	for (int i=0; i<matrix.no_cols; i++)
	{
		for (int j=0; j<matrix.no_rows; j++)
		{
			transposeMatrix.matrix_entries[i][j] = matrix.matrix_entries[j][i];
		}
	}
}

///////////////////////////////////////////////////////////////////
void ConstructVectorV(DoubleMatrix A, double* b, double* x, double* v)
//Constructs the vector v=A^t(b-Ax)
//Assumes the memory for v has already been allocated 
{
	//Initialise the values of the variables that we need
	
	double* A_tAx = AllocateVector(A.no_cols);
	DoubleMatrix A_Transpose = AllocateDoubleMatrix(A.no_cols, A.no_rows);  
    TransposeMatrix(A, A_Transpose);
	SymmetricMatrix A_tA_symmetric = ConstructSymmetricMatrix(A_Transpose);
	DoubleMatrix A_tA = ConvertSymmetricMatrixToDouble(A_tA_symmetric);
	

	//Apply linear algebra functions to find V
	MultiplyVectorByMatrix(A_Transpose, b, v);
    MultiplyVectorByMatrix(A_tA, x, A_tAx);
    SubtractVectors(v,A_tAx,A.no_cols);

	//Deallocate
	DeallocateMatrix(A_Transpose);
	DeallocateMatrix(A_tA_symmetric);
	DeallocateMatrix(A_tA);
	DeallocateVector(A_tAx);

}

////////////////////////////////////////////////////////////////////
void SolveRestrictedProblem(DoubleMatrix& A_p, double* b,  double* s_p)
//Solves the restriced problem A^ps^p = b
//Assumes the memory for v has already been allocated
{
	//Initialise
	int p_size = A_p.no_cols;
	double* A_p_tb = AllocateVector(p_size);
    DoubleMatrix A_p_T = AllocateDoubleMatrix(p_size, A_p.no_rows);
	TransposeMatrix(A_p, A_p_T);
	SymmetricMatrix A_tA_symmetric = ConstructSymmetricMatrix(A_p_T);
	
	
	
	
	//We construct the equation (A^p)^t(A^p)s^p = (A^P)^tb
    MultiplyVectorByMatrix(A_p_T, b, A_p_tb); 
	
	//We solve by symmetric factorisation
    PerformSymmetricSolve(A_tA_symmetric, A_p_tb, s_p);

	//Deallocate
	DeallocateVector(A_p_tb);
	DeallocateMatrix(A_p_T);
	DeallocateMatrix(A_tA_symmetric);

}

//////////////////////////////////////////////////////////////
void PrintNormBMinusAX(DoubleMatrix& A, double* b, double* x)
//Prints the norm of the vector (b-Ax) 
{	
	//Initialise
	double* Ax = AllocateVector(A.no_rows);
    double* b_copy = AllocateVector(A.no_rows); 
	
	//Finds the residual
	CopyVector(b, b_copy, A.no_rows);
	MultiplyVectorByMatrix(A, x, Ax);
	SubtractVectors(b_copy, Ax, A.no_rows);
	
	//Prints the residual to the screen
	std::cout << NormVector(b_copy, A.no_rows) << std::endl;

	//Deallocate
	DeallocateVector(Ax);
	DeallocateVector(b_copy);

}


