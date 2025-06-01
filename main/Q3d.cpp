#include "gauss_quadrature.hpp"
#include "legendre.hpp"
#include "polygon.hpp"

//////////////////////////////////////////////////////////////
//Function Prototypes
//////////////////////////////////////////////////////////////

double func1(double x, double y);
double func2(double x, double y);
double func3(double x, double y);

//////////////////////////////////////////////////////////////
int main()
{
     
     int n; int refine; char polygon;
     double val1=0; double val2=0; double val3=0; //Initialise
	 
	 GeneralPolygon Omega; 

     std::cout << "Choose the polygon to construct a quadrature on" << std::endl;
     std::cout << "Enter s for square, l for L shape or h for hexagon. The default is hexagon." << std::endl;
     std::cin >> polygon;

     std::cout << "Enter the number of refinements and quadrature points" << std::endl;
     std::cin >> refine >> n;
     
    switch(polygon)
    {
        //Switch statement to generate the polgon depending on the choice by the user
		case('s'):
     		Omega = GenerateSquareMesh(refine);
			break;
		case('l'):
			Omega = GenerateLShapedMesh(refine);
			break;
		default:
			Omega = GenerateHexagonalMesh(refine); 
     	
	} 
     
	DoubleMatrix quad_points_poly = AllocateDoubleMatrix(n*n*(Omega.no_triangles),3); //Initialise the matrix for the polygon quad points
    quad_points_poly.matrix_entries = ConstructPolygonQuad(Omega, n);   //Initialise the quadrature points
	
	int no_rows = quad_points_poly.no_rows; //Initialise the rows and columns
	int no_cols = quad_points_poly.no_cols;

	

    for (int i=0;i<no_rows;i++)
	{   
	    double x_i = quad_points_poly.matrix_entries[i][0];
	    double y_i = quad_points_poly.matrix_entries[i][1];
        double w_i = quad_points_poly.matrix_entries[i][2];
        val1 += w_i * func1(x_i,y_i);   //Finds the approximation of the integral of the given functions over the given polygon
        val2 += w_i * func2(x_i,y_i);
		val3 += w_i * func3(x_i,y_i);
	}
				
	PrintMatrix(quad_points_poly); //Prints the quadrature to the screen

	//Prints to the screen the approximations to the integrals
    std::cout << "The polygon integrals are " << std::endl;
	std::cout << val1 << std::endl;
	std::cout << val2 << std::endl;
	std::cout << val3 << std::endl;

	DeallocateMatrix(quad_points_poly);
	DeallocateGeneralPolygon(Omega);   //Deallocate the polygon and the quadrature.
	
     
    return 0; 

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