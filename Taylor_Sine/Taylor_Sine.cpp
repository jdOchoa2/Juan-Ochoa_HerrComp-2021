#include <iostream>
#include <cmath>

double mysin(double x, int N, double []);                  //Returns the aproximation of sin(x) with N terms
  
int main()
{
  std::cout.setf(std::ios::scientific);
  std::cout.precision(6);

  const double x = M_PI/3;
  const double exact = std::sin(x);
  const int NMAX = 10000;

  double terms[NMAX]={x};
  double diff = std::fabs(mysin(x,1,terms)-exact)/exact;   //Return relative error
  std::cout << "1" << "\t" << diff << "\n";                //Print the relative error with N=1
  
  for(int N  = 2; N <= NMAX; N++)
    {
      terms[N-1]=-terms[N-2]*std::pow(x,2)/(4*N*N-6*N+2);  //Returns the Nth term of the Taylor expansion using the previus one
      diff = std::fabs(mysin(x,N,terms)-exact)/exact;      //Return relative error
      std::cout << N << "\t" << diff << "\n";              //Print the relative error with N terms
    }
  return 0;
}

double mysin(double x, int N, double a[])
{
  double suma = 0.0;
  for(int n = N-1; n >= 0; n--)                            //Sums the terms backwrds
    {
      suma += a[n];
    }
  return suma;
}
