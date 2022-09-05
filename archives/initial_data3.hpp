#include <cstdlib>
#include <cmath>
double pi = 4*atan(1.0);
double mu =  1./100;  //  1./10

/* Initial datum and boundary */
 double fi_0(double x,double v)
{
    return exp(-0.5 * v*v)/sqrt(2*pi);
} 

double fe_0(double x,double v)
{
    return sqrt(mu)*exp(-0.5 * mu * v*v)/sqrt(2*pi);             // exp( -100*(x-1)*(x-1) - 100*v*v );      //sqrt(mu)*exp(-0.5 * mu * v*v)/sqrt(2*pi);
}

double E0(double x)
{
    return 0;
}

 double fe_L(double v, double tn)
{
   // X = -cos(tn)-v*sin(tn); Y = -sint(tn)+v*cos(tn); 
    return  exp( -100*(-cos(tn)-v*sin(tn)-1)*(-cos(tn)-v*sin(tn)-1) - 100*(-sin(tn)+v*cos(tn))*(-sin(tn)+v*cos(tn)) );
} 

 double fe_R(double v, double tn)
{
    //X = cos(tn)-v*sin(tn); Y = sint(tn)+v*cos(tn); 
    return  exp( -100*(cos(tn)-v*sin(tn)-1)*(cos(tn)-v*sin(tn)-1) - 100*(sin(tn)+v*cos(tn))*(sin(tn)+v*cos(tn)) );
} 
