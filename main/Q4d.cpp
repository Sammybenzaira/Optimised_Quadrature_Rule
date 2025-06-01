#include "optimisation.hpp"
#include "gauss_quadrature.hpp"
#include "polygon.hpp"

double func1(double x, double y);
double func2(double x, double y);
double func3(double x, double y);

int main()
{
    //Initialise
     int degree; int refine; char polygon; int n;
     double val1=0; double val2=0; double val3=0;
     
	 GeneralPolygon Omega;

     std::cout << "Choose the polygon to construct a quadrature on" << std::endl;
     std::cout << "Enter s for square, l for L shape or h for hexagon. The default is hexagon." << std::endl;
     std::cin >> polygon;

     std::cout << "Enter the number of refinements and the polynomial degree that should be integrated exactly" << std::endl;
     std::cin >> refine >> degree;

     
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
     
	//Construct the optimised quadrature 
	DoubleMatrix matrix = ConstructOptimisedQuad(Omega, degree);

	
	//Constructs the integral as well as printing the optimised quadrature to the screen
    for (int i=0;i<matrix.no_rows;i++)
	{   
	    double x_i = matrix.matrix_entries[i][0];
	    double y_i = matrix.matrix_entries[i][1];
        double w_i = matrix.matrix_entries[i][2];
        val1 += w_i*func1(x_i,y_i);
        val2 += w_i*func2(x_i,y_i);
		val3 += w_i*func3(x_i,y_i);
				
		for (int j=0;j<3;j++)
		{
			std::cout << matrix.matrix_entries[i][j] << " ";    
		}
		std::cout << std::endl;
	} 

    
    std::cout << "The polygon integrals are " << std::endl;
	std::cout << val1 << std::endl;
	std::cout << val2 << std::endl;
	std::cout << val3 << std::endl;

	//Deallocate
	DeallocateMatrix(matrix);
	DeallocateGeneralPolygon(Omega);
     
    return 0; 

}

double func1(double x, double y)

{
	return 3*pow(x,3)*pow(y,2)+2*pow(x,2)*y+2*x;
}

double func2(double x, double y)

{
	return sin(M_PI*(x+y));
}

double func3(double x, double y)

{
	return exp(pow(x,2)+pow(y,2));
}