#include <cassert>
#include <iostream>
#include "Vector3.h"
#include <math.h>
using namespace std;

Vector::Vector(int taille,double* coordinate):size(taille)
{
	if(taille==0)
		Coord=NULL;
	else if(coordinate==NULL)
    {
		Coord=new double [taille];
		for(int i=0;i<size;i++)
		{
			Coord[i]=0;
		}
    }
	else
    {
		Coord=new double [taille];
		for(int i=0;i<taille;i++)
		{
			Coord[i]=coordinate[i];
		}
    }
}

Vector::~Vector()
{
	
	delete [] Coord;
}
int Vector::getSize()const
{
	return size;
}
Vector::Vector(const Vector& V):size(V.size)
{
	if(V.size==0)
    {
		Coord=NULL;
    }
	else
    {
		Coord=new double [V.size];
		for(int i=0;i<V.size;i++)
		{
			Coord[i]=V.Coord[i];
		}
    }
}

Vector& Vector::operator=(const Vector&V)
{
	if(this!=&V)
    {
		delete [] Coord;
		size=V.size;
		Coord=new double[V.size];
		for(int i=0;i<V.size;i++)
		{
			Coord[i]=V.Coord[i];
		}
    }
	else
    {
		cout<<"Affectation Impossible -> Adresse memoire identique"<<endl;
    }
	return *this;
}

double& Vector::operator[](int j)const
{
	return Coord[j];
}

double Vector::operator*(const Vector& V)
{
	if(size !=V.size)
    {
		cout<<"Error Scalar product -> Incompatible dimension"<<endl;
		return ErrorDimension;
    }
	else
    {
		double result=0;
		for(int i=0;i<size;i++)
		{
			result=result+Coord[i]*V.Coord[i];
		}
		return result;
    }
}
Vector Vector::operator+(const Vector& V)
{
	if(size!=V.size)
    {
		cout<<"Error sum -> Incompatible dimension"<<endl;
		return ErrorDimension;
    }
	else
    {
		Vector resultat(size);
		for(int i=0;i<size;i++)
		{
			resultat.Coord[i]=Coord[i]+V.Coord[i];
		}
		return resultat;
    }
}
Vector Vector::operator-(const Vector& V)
{
	if(size!=V.size)
    {
		cout<<"Error sum -> Incompatible dimension"<<endl;
		return ErrorDimension;
    }
	else
    {
		Vector resultat(size);
		for(int i=0;i<size;i++)
		{
			resultat.Coord[i]=Coord[i]-V.Coord[i];
		}
		return resultat;
    }
}
double Vector::Norm2()
{
double Norme=0;
for(int i=0;i<size;i++)
{
Norme+=Coord[i]*Coord[i];
}
return sqrt(Norme);
}

double Vector::linfty()
{
double temp=fabs(Coord[0]);
for(int i=0;i<size;i++)
{
if(temp<fabs(Coord[i]))
{
temp=fabs(Coord[i]);
}
}
return temp;
}

double Vector::min()
{
  double r = Coord[0];
  for(int i= 0 ; i < size;i++)
    {
      if( Coord[i] < r )
	{
	  r = Coord[i];
	}
    }
  return r;
}

double Vector::max()
{
  double r = Coord[0];
  for(int i = 0; i < size;i++)
    {
      if( Coord[i] > r)
	{
	  r = Coord[i];
	}
    }
  return r;
}

Vector operator*(const double& k,const Vector& V)
{
	if(V.size==0)
    {
		cout<<"Error multiplication -> Not a Vector or a RealNumber"<<endl;
		return ErrorDimension;
    }
	else
    {
		Vector resultat(V.size);
		for(int i=0;i<V.size;i++)
		{
			resultat.Coord[i]=k*V.Coord[i];
		}
		return resultat;
    }
}

Vector Vector::operator^(const Vector & V)
{
assert(V.size==3); // Produit Vectoriel de R3

Vector res(3);

res[0]=Coord[1]*V[2]-Coord[2]*V[1];
res[1]=Coord[2]*V[0]-Coord[0]*V[2];
res[2]=Coord[0]*V[1]-Coord[1]*V[0];

return res;

}
ostream& operator<<(ostream& flux,const Vector& V)
{
	flux<<"[";
	for(int i=0;i<V.size;i++)
    {
		flux<<" "<<V[i]<<" ";
    }
	flux<<"];";
	//flux<<"     "<<"dim["<<V.size<<"]";
	return flux;
}




