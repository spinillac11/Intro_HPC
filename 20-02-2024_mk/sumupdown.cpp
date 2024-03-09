#include <iostream>
#include <cmath>

float sumup(int nterms);
float sumdown(int nterms);

int main(int argc, char **argv)
{
  std::cout.precision(6); std::cout.setf(std::ios::scientific);
  for (int ii = 1; ii <= 10000; ++ii) {
    float sum1 = sumup(ii);
    float sum2 = sumdown(ii);
    std::cout << ii << "\t" << sum1 << "\t" << sum2
              << "\t" << std::fabs(sum1-sum2)/sum2 << "\n";
  }

  return 0;
}

float sumup(int nterms)
{
    float result = 0.0;
    for(int k = 1; k <= nterms; k++) {
        result += 1.0/k;
    }
    return result;
}

float sumdown(int nterms)
{
    float result = 0.0;
    for(int k = nterms; k >= 1; k--) {
        result += 1.0/k;
    }
    return result;
}
