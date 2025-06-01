#include "optimisation.hpp"
#include "linear_algebra.hpp"
#include "legendre.hpp"
#include "gauss_quadrature.hpp"

double* ConstructOptimisation(DoubleMatrix A, double* b, double tol)
//Constructs optimisation 1 that returns the optimised quadrature where the weights of some points are
//Zero while still being able to integrate polynomials over a given degree
//The function takes as input the matrix A and the vector b where A and b are defined in question 4b and a given tolerance
//Returns the vector x which is the optimised weights with some values being zero  
{    
    //Initialisation of vectors and arrays
    std::vector<int> p;
    double* x = AllocateVector(A.no_cols);
    double* v = AllocateVector(A.no_cols);
    double* v_copy = AllocateVector(A.no_cols);
    double* s = AllocateVector(A.no_cols);
    double* s_copy = AllocateVector(A.no_cols);

    int index; double max_val_v; double min_val_s_p;
    int Iteration_no = 0;
    
    //Constructs the vector v given as v=A^t(b-Ax)  
    ConstructVectorV(A, b, x, v);
    
    //Finds the maximum value of the vector v
    FindMaximum(v, A.no_cols, max_val_v, index);
    
    
    while(p.size() < A.no_cols && max_val_v > tol)
    {   
        
        //Construct a copy of v where for every index of v that is in the set p,
        //We set that value to a very large negative number
        CopyVector(v, v_copy, A.no_cols);
        for (int j=0; j<p.size(); j++ )
        {
            v_copy[p[j]] = -1E+300;
        }
        
        //We then find the maximum value in this copy of v which shouldn't be in P
        //Then we append this value to the set P
        FindMaximum(v_copy, A.no_cols, max_val_v, index);
        p.push_back(index);
        
        //We initialise the matrix A^P as describe in question 4
        DoubleMatrix A_p = AllocateDoubleMatrix(A.no_rows, p.size());
        
        //Fill the A^P matrix with values
        for (int i=0; i<A.no_rows; i++ )
        {
            
            for (int j=0; j<p.size(); j++ )
            {
                A_p.matrix_entries[i][j] = A.matrix_entries[i][p[j]];
            }
        
        }
        
        //Allocate the vector s_p  
        double* s_p = AllocateVector(p.size());
        
        
        //We solve the restricted problem A^p*s^p=b 
        SolveRestrictedProblem(A_p, b, s_p);
        
        //Find the minimum value of the vector s^p
        
        FindMinimum(s_p, p.size(), min_val_s_p, index);
         
        DeallocateMatrix(A_p);
        
        while(min_val_s_p<=0)
        {
            
            //Set the initial Alpha using the minimum value of s_p since this must be less 
            //Than zero if it entered the while loop 
            double alpha = x[p[index]]/(x[p[index]]-s_p[index]);

            //Finds the minimum value and index of x_i/(x_i-s_i) where s_i<0
            //Whilst also filling the vector s with the elements in s_p
            //At the corresponding index
            for (int i=0; i<p.size(); i++)
            {
                if ((alpha > x[p[i]]/(x[p[i]]-s_p[i])) && (s_p[i]<=0))
                {
                    alpha = x[p[i]]/(x[p[i]]-s_p[i]);
                    index = i;
                }
                s[p[i]] = s_p[i];
            }
            
            //Erases the index found since this would be the ith value where
            //x_i = 0
            p.erase(p.begin()+index);
    
            
            //Constructs the new vectorx=x+alpha(s-x)
            CopyVector(s, s_copy, A.no_cols);
            SubtractVectors(s_copy, x, A.no_cols );
            CombineVectors(x, s_copy, alpha, A.no_cols);
            
            //We initialise the matrix A^P again since the size of p has changed
            A_p = AllocateDoubleMatrix(A.no_rows, p.size());
            
            //Fill the A^P matrix with values
            for (int i=0; i<A_p.no_rows; i++ )
            {
                for (int j=0; j<p.size(); j++ )
                {
                    A_p.matrix_entries[i][j] = A.matrix_entries[i][p[j]];
                }
            }
            
            //Create a new s_p vector using the new cardinality of set P.
            DeallocateVector(s_p);
            s_p = AllocateVector(p.size());
            
            //We solve the restricted problem A^p*s^p=b and find the new minimum value of the vector s^p
            SolveRestrictedProblem(A_p, b, s_p);
            FindMinimum(s_p, p.size() , min_val_s_p, index);
            
            
            DeallocateMatrix(A_p);
            
        
        }
        
        //Initialise the vector s and fill it with elements in s^p at the corresponding indices
        ZeroVector(s,A.no_cols);
        for (int i=0; i<p.size(); i++)
        {
            s[p[i]] = s_p[i];
        }
        
        //Set x is equal to s
        CopyVector(s, x, A.no_cols);
        
        //Create a new vecotr v and find the new maximum value
        ConstructVectorV(A, b, x, v);
        FindMaximum(v, A.no_cols, max_val_v, index);

        //Print to the screen the iteration number and the norm of the residual
        Iteration_no++;
        
        std::cout << "Iteration number" << std::endl;
        std::cout << Iteration_no << std::endl;
        
        std::cout << "Norm of residual" << std::endl;
        PrintNormBMinusAX(A, b, x);
        
        DeallocateVector(s_p);
        
    }

    //Deallocate the used vectors
    DeallocateVector(v);
    DeallocateVector(v_copy);
    DeallocateVector(s_copy);
    DeallocateVector(s);

    return x;
}

void ConstructMatrixAVectorB(double** quadPoints, int degree, DoubleMatrix& A, double* b, int numberOfQuadPoints)
//Constructs the matrix A and the vector b as shown in Q4b for a given degree and a given set of quadrature points
//Assumes that A and b have already been allocated
{
    
    for (int j=0; j<degree+1; j++)
    {
        for (int k=0; k<degree+1; k++)
        {
                
            b[j*(degree+1)+k] = 0; //Initialises the b vector
        
            
        
            for (int i=0; i<numberOfQuadPoints; i++)
            {

            double x_val = quadPoints[i][0];
            double y_val = quadPoints[i][1]; //Initialise
            double w_val = quadPoints[i][2];
                
            //Constructs the matrix A and the vector b
            b[j*(degree+1)+k] += w_val*ComputeLegendre(x_val, j)*ComputeLegendre(y_val, k);
            A.matrix_entries[j*(degree+1)+k][i] = ComputeLegendre(x_val, j)*ComputeLegendre(y_val, k);

            }
        }
    }
}

DoubleMatrix ConstructOptimisedQuad(GeneralPolygon Omega, int degree)
//Produces an optimised quadrature over a general polygon and to an exact degree specified by the input.
//The function returns a new matrix pointer containing the new optimised quadrature points
{
    //Initialise
    //We let the number of quad points be the minimum needed to integrate the chosen degree exactly
    int n = int (degree+2)/2;
    DoubleMatrix A = AllocateDoubleMatrix(degree*degree+2*degree+1, n*n*Omega.no_triangles);
    double* b = AllocateVector(degree*degree+2*degree+1);
    std::vector<int> p_set;

    //Construct the unoptimised quadrature points over the chosen shape
    //With n points in each direction
    double** quad_points = ConstructPolygonQuad(Omega, n);

    //Construct the matrix A and the vector b
    ConstructMatrixAVectorB(quad_points, degree, A, b, n*n*Omega.no_triangles);

    //Construct the new optimised weights
    double* new_weights = ConstructOptimisation(A, b, 1E-12);
    
    //Deallocate
    DeallocateMatrix(A);
    DeallocateVector(b);


    //Loop through the new weights and add the index where the weight is not zero into p_set.
    for (int i=0; i<n*n*Omega.no_triangles; i++)
    {
        quad_points[i][3] = new_weights[i];
        if (new_weights[i] != 0)
        {
            p_set.push_back(i);
        }
    }
    
    //Initialise the new optimised quad points with no zero weights
    DoubleMatrix optimised_quad_points;
    optimised_quad_points.no_rows = p_set.size();
    optimised_quad_points.no_cols = 3;

    optimised_quad_points.matrix_entries = new double*[p_set.size()];
    
    //Fill the new quad points array with the optimised weight and the corresponding quad points
    for (int j=0; j<p_set.size(); j++)
    {
        optimised_quad_points.matrix_entries[j] = new double[3];
        optimised_quad_points.matrix_entries[j][0] = quad_points[p_set[j]][0];
        optimised_quad_points.matrix_entries[j][1] = quad_points[p_set[j]][1];
        optimised_quad_points.matrix_entries[j][2] = new_weights[p_set[j]];
        
    }

    return optimised_quad_points;

}

