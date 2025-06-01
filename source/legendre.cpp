//////////////////////////////////////////////////////////////
//Module that computes the Legendre polynomials and finds
//their roots
//////////////////////////////////////////////////////////////

#include "legendre.hpp"

//////////////////////////////////////////////////////////////


double ComputeLegendre(double x, int p)
//Computes the value of the pth order Legendre polynomial at
//a given point x
{
	if (p==0)
	{
		return 1;
	}
	else if (p==1)
	{
		return x;
	}
	else
	{
        double val2 = 1;
        double val1 = x;
        double val;
        for (int k=2;k<=p;k++)
        {
            val = ((2.0*double(k)-1.0)*x*val1-double(k-1)*val2)/double(k);

            val2 = val1;
            val1 = val;
        }

        return val;
	}
}

double ComputeDiffLegendre(double x, int p)
//Computes the value of the derivative of the pth order Legendre polynomial at
//a given point x
{
	if (p==0)
	{
		return 0;
	}
	else if (p==1)
	{
		return 1;
	}
	else
	{
        double val2 = 0;
        double val1 = 1;
        double val;
        for (int k=2;k<=p;k++) //We use a similar technique to finding a legendre polynomial 
        
		{
            val = ((2.0*double(k)-1.0)*(ComputeLegendre(x, k-1)+x*val1)-double(k-1)*val2)/double(k);
			                     //Initialise the first value as the derivative of the legendre polynomial recursive equation
            val2 = val1;
            val1 = val;
        }

        return val;
	}
}
double FindRootNewton(double x, int p, double tol)
//Computes the root of the pth order legendre polynomial using the Newton Raphson method
// Given an initial guess with an error less than tol. 
{   
	double val;
	bool running = true; //Initialise
    
	while (running)
	{
      val=x-(ComputeLegendre(x, p))/ComputeDiffLegendre(x, p); //Repeatedly calculates value using the Newton Raphson equation
	  
	  if (fabs(val-x)<tol) //Checks if the value is less than tol
	  {
		running=false;
	  }
      else
	  {
		x=val;
	  }
	}

    return val; //Returns the last calculated root

}