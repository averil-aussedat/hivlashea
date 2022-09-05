#include <cstdlib>
#include <cmath>
double pi = 4*atan(1.0);
double mu = 1./100;

/* Initial datum and boundary */
double fi_0(double x,double v)
{
    return exp(-0.5 * v*v)/sqrt(2*pi);
}

double fe_0(double x,double v)
{
    return sqrt(mu)*exp(-0.5 * mu * v*v)/sqrt(2*pi);
}

double E0(double x)
{
    return 0.;
}

