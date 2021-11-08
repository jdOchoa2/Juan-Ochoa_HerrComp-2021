#include <iostream>
#include <cstdlib>
#include <cmath>

int foo(int a, int b);
int bar(int a, int b);
double baz(double x);
void print_array(const double data[], const int & size);

int main (int argc, char **argv)
{
  int a, b;
  a =  0;
  b = -1;

  foo(a, b);
  foo(b, a);
  baz(25.9);

  const int NX = 2, NY = 3, NZ = 5;
  double *x, y[NY]={0}, z[NZ]={0};
  x = new double [NX]();

  print_array(x, NX);
  print_array(y, NY);
  print_array(z, NZ);
  std::cout << std::endl;

  for (int ii = 0; ii < NX; ++ii) {
    x[ii] = ii;
  }

  print_array(x, NX);
  print_array(y, NY);
  print_array(z, NZ);
  std::cout << std::endl;

  for (int jj = 0; jj < NY; ++jj) {
    y[jj] = jj;
  }

  print_array(x, NX);
  print_array(y, NY);
  print_array(z, NZ);
  std::cout << std::endl;

  for (int kk = 0; kk < NZ; ++kk) {
    z[kk] = kk;
  }

  print_array(x, NX);
  print_array(y, NY);
  print_array(z, NZ);
  std::cout << std::endl;

  delete [] x;

   return EXIT_SUCCESS;
}

int foo(int a, int b)
{
  if (a*b*bar(a,b) == 0)
    {
      std::cout << "ATENCIÓN: División por 0"<<std::endl;
      return 0;
    }
  else return a/b + b/bar(a, b) + b/a;
}

int bar(int a, int b)
{
  int c = std::abs(a);
  return c + a - b;
}

double baz(double x)
{
  if (std::abs(x) < 0.0000001) return 1-(x+1);
  if (x<0)
  {
    std::cout<<"ATENCIÓN: Raíz imaginaria"<<std::endl;
  }
  return std::sqrt(x);
}

void print_array(const double data[], const int & size)
{
  std::cout << std::endl;
  for (int ii = 0; ii < size; ++ii) {
    std::cout << data[ii] << "  " ;
  }
}
