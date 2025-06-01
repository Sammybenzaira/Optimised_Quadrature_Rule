#include "gauss_quadrature.hpp"
#include "linear_algebra.hpp"

//////////////////////////////////////////////////////////////
//Function Prototypes
//////////////////////////////////////////////////////////////

double func1(double x, double y);
double func2(double x, double y);
double func3(double x, double y);

//////////////////////////////////////////////////////////////
int main()
{
     
    int n;
    double sq_val1=0; double sq_val2=0; double sq_val3=0; //initialise

    std::cout << "Enter the number of quadrature points please" << std::endl;
    std::cin >> n;

	DoubleMatrix quad_points_square = AllocateDoubleMatrix(n*n, 3); //Initialise the matrix for the square quad points
    quad_points_square.matrix_entries = ConstructSquareGaussQuad(n);   //Initialise the quadrature points                                                     
    
     
    
	for (int i=0; i<n; i++)
    {  for (int j=0; j<n; j++) 
        {
			double x_i = quad_points_square.matrix_entries[i*n+j][0];
			double y_i = quad_points_square.matrix_entries[i*n+j][1];
            double w_i = quad_points_square.matrix_entries[i*n+j][2];
            sq_val1 += w_i * func1(x_i,y_i);  //Finds the approximation of the integral of the given functions over the square region
            sq_val2 += w_i * func2(x_i,y_i);
			sq_val3 += w_i * func3(x_i,y_i);
		}
    }
     
	double tri_val1=0; double tri_val2=0; double tri_val3=0; //Initialise

	DoubleMatrix quad_points_tri = AllocateDoubleMatrix(n*n, 3); //Initialise the matrix for the square quad points
    quad_points_tri.matrix_entries = ConstructTriangleGaussQuad(n);   //Initialise the quadrature points
	
	for (int i=0; i<n; i++)
    {  for (int j=0; j<n; j++) 
        {
			double x_i = quad_points_tri.matrix_entries[i*n+j][0];
			double y_i = quad_points_tri.matrix_entries[i*n+j][1];
            double w_i = quad_points_tri.matrix_entries[i*n+j][2];
            tri_val1 += w_i * func1(x_i,y_i); //Finds the approximation of the integral of the given functions over the reference triangle
            tri_val2 += w_i * func2(x_i,y_i);
			tri_val3 += w_i * func3(x_i,y_i);
		}
    }
    
	//Prints to the screen the approximations to the integrals
	std::cout << "The square integrals are " << std::endl;
	std::cout << sq_val1 << std::endl;
	std::cout << sq_val2 << std::endl;
	std::cout << sq_val3 << std::endl;
	std::cout << "The Triangle integrals are " << std::endl;
	std::cout << tri_val1 << std::endl;
	std::cout << tri_val2 << std::endl;
	std::cout << tri_val3 << std::endl;

	DeallocateMatrix(quad_points_square);
	DeallocateMatrix(quad_points_tri);   //Deallocate the used matrices
}

/////////////////////////////////////////////////////////////////
double func1(double x, double y)
//Function that returns 3(x^3)(y^2)+2y(x^2)+2x for given x and y
{
	return 3*pow(x,3)*pow(y,2)+2*pow(x,2)*y+2*x;
}

////////////////////////////////////////////////////////////////
double func2(double x, double y)
//Function that returns sin(pi(x+y)) for given x and y
{
	return sin(M_PI*(x+y));
}

/////////////////////////////////////////////////////////////////
double func3(double x, double y)
//Function that returns exp(x^2+y^2) for given x and y
{
	return exp(pow(x,2)+pow(y,2));
}