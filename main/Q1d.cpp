#include "gauss_quadrature.hpp"
#include "linear_algebra.hpp"

int main()
{

     
    int n;
    double val1=0; double val2=0; //Initialise


    std::cout << "Enter the number of quadrature points" << std::endl;
    std::cin >> n; 

    DoubleMatrix quad_points = AllocateDoubleMatrix(n, 2); //Use functions from linear_algebra.hpp to print the quad points
                                                         
    quad_points.matrix_entries = ConstructGaussQuad(n);   //Initialise the quadrature points                                                  
    
    PrintMatrix(quad_points);   //Prints the quad points
    
     
    for (int i=0; i<n; i++)
    {   
        double x_i = quad_points.matrix_entries[i][0];
        double w_i = quad_points.matrix_entries[i][1]; 
        val1 += w_i * exp(x_i);            //Finds the approximation of the given integrals
        val2 += w_i * cos(M_PI*x_i/2);  
    }
     
    std::cout << "The integral I_1 and I_2 are " << std::endl;
    std::cout << val1 << " and " << val2 << std::endl; //Prints the integral values

}

