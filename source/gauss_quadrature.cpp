//////////////////////////////////////////////////////////////
//Module that constructs Gauss Quadrature on various domains
//////////////////////////////////////////////////////////////


#include "gauss_quadrature.hpp"
#include "polygon.hpp"
#include "legendre.hpp"


//////////////////////////////////////////////////////////////
/// Gauss Quadrature Functions
//////////////////////////////////////////////////////////////

double** ConstructGaussQuad(int n)
//Constructs the one dimensional gauss quadrature for a given number of points n.
//Output returns a pointer to a matrix where the first column contains the x values
//and the second column contains the corresponding weights
{
    double** quad_points = new double* [n]; //Initialise
    
    for (int i=1; i<=n; i++)
    {
        quad_points[i-1]=new double[2]; //Creates two columns for the x values and weights

        double x_value=(double) cos((M_PI*(i-0.5))/n); //Initialise the first guess of x as the roots of the Chebyshev polynomials
        x_value = FindRootNewton(x_value, n, 1E-12); //Find the root of the nth lagrange polynomial 
        double w_value = 2/((1-pow(x_value,2))*pow(ComputeDiffLegendre(x_value, n),2)); //Calculate the corresponding weights
        
        quad_points[i-1][0] = x_value;
        quad_points[i-1][1] = w_value; 
    }

    return quad_points;
}

/////////////////////////////////////////////////////////////////////////////
double** ConstructSquareGaussQuad(int n)
//Constructs a 2D quadrature on a square domain S:=(-1,1)^2 for a given number of Quadrature points in each direction n.
//Output returns a pointer to a matrix where the first two columns contain the x and y values
//And the third column contains the corresponding weights with n^2 rows
{
    double** quad_points_square = new double* [n*n];
    double** quad_points = ConstructGaussQuad(n);    //Initialise
    
    for (int i=0; i<n; i++)
    {
        
        for (int j=0; j<n; j++)
        {
            quad_points_square[n*i+j] = new double [3];
            quad_points_square[n*i+j][0] = quad_points[i][0];   //Uses the 1D Quadrature to construct the 2D Quadrature points
            quad_points_square[n*i+j][1] = quad_points[j][0];  
            quad_points_square[n*i+j][2] = quad_points[i][1]*quad_points[j][1];
        }
        
    }
    for (int i=0;i<n;i++)
    {
        delete[] quad_points[i]; //Deallocate the quad_points matrix
    }
    delete[] quad_points;

    return quad_points_square;
}

////////////////////////////////////////////////////////////////////////
double** ConstructTriangleGaussQuad(int n)
//Uses a mapping from S to T to construct the Quadrature points of the reference triangle for n quad points in each direction
//Output returns a pointer to a matrix where the first two columns contain the x and y values
//And the third column contains the corresponding weights with n^2 rows
{
    
    double** quad_points_triangle = ConstructSquareGaussQuad(n); //Initialise
    
    for (int i=0; i<n; i++)
    {
        
        for (int j=0; j<n; j++)
        {
            
            
            double w_value = quad_points_triangle[n*i+j][2]; //Initialise
            double y_value = quad_points_triangle[n*i+j][1];
            double x_value = quad_points_triangle[n*i+j][0];
            
            //Finds the new x values of the Quadrature points of the reference triangle
            //And the corresponding new weights
            //Note that the y values stay the same so they dont need changing
            quad_points_triangle[n*i+j][0] = -1+((1-y_value)*(x_value+1))/2; //Finds the new x values of the Quadrature poins of the reference triangle
            quad_points_triangle[n*i+j][2] = w_value*(1-y_value)/2;
        }
    }
    return quad_points_triangle;

}

/////////////////////////////////////////////////////////////////////////////////////
double** TransformTriangleQuad(Triangle Tri, int n )
//Uses a mapping from the reference triangle to a any general triangle 
//To construct the Quadrature points of a general triangle for n quad points in each direction
//Output returns a pointer to a matrix where the first two columns contain the x and y values
//And the third column contains the corresponding weights with n^2 rows
{
    double** new_quad_points = ConstructTriangleGaussQuad(n);
    
    double T_x1=Tri.vertices[0][0]; double T_y1=Tri.vertices[0][1];
    double T_x2=Tri.vertices[1][0]; double T_y2=Tri.vertices[1][1]; //Initialising the vertices of the triangle
    double T_x3=Tri.vertices[2][0]; double T_y3=Tri.vertices[2][1];


    for (int i=0; i<n; i++)
    {
        
        for (int j=0; j<n; j++)
        {
            
            
            double w_value = new_quad_points[n*i+j][2];
            double y_value = new_quad_points[n*i+j][1];
            double x_value = new_quad_points[n*i+j][0];
            
            //Find the new quad points using an appropriate affine map.
            new_quad_points[n*i+j][0] = ((T_x2-T_x1)*x_value+(T_x3-T_x1)*y_value+(T_x2+T_x3))/2;  
            new_quad_points[n*i+j][1] = ((T_y2-T_y1)*x_value+(T_y3-T_y1)*y_value+(T_y2+T_y3))/2;     
            new_quad_points[n*i+j][2] = w_value*((T_y3-T_y1)*(T_x2-T_x1)-(T_y2-T_y1)*(T_x3-T_x1))/4;
        }
    }

    return new_quad_points;
}

////////////////////////////////////////////////////////////////////////////////////
double** ConstructPolygonQuad(GeneralPolygon Omega, int n)
//Constructs the Quadrature points of any polygon made up of triangles
//Output returns a pointer to a matrix where the first two columns contain the x and y values
//And the third column contains the corresponding weights
//Where the first n^2 rows correspond to the Quadrature points of the first triangle and so on.
{
    int n_squared = pow(n,2);
    double** quad_points_poly = new double* [n_squared*Omega.no_triangles]; //Initialise new matrix
    ;
    for (int i=0; i<Omega.no_triangles; i++)
    {
        Triangle Tri = Omega.component_triangles[i];
        double** quad_points_triangle = TransformTriangleQuad(Tri, n);
        
        for (int j=0; j<n_squared; j++)
        {
            quad_points_poly[n_squared*i+j] = new double [3];
            //For each triangle of the polygon we construct the corresponding quadrature points and add them to our matrix
            quad_points_poly[n_squared*i+j][0] = quad_points_triangle[j][0]; 
            quad_points_poly[n_squared*i+j][1] = quad_points_triangle[j][1];
            quad_points_poly[n_squared*i+j][2] = quad_points_triangle[j][2];
            
            
        } 
        for (int k=0;k<n_squared;k++)
        {
            delete[] quad_points_triangle[k]; //Deallocate the matrix used.
        }
        delete[] quad_points_triangle;

    }
    
    
    return quad_points_poly;

}
