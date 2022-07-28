#include "Matrix.h"
#include <iostream>
#include <cassert>
#include "math.h"
#include <ctime>
#define zero 1e-12
using namespace std;
Matrix::Matrix(int M,int N, double ** Coefficient):length(M),width(N)
{
  if(M==0 || N==0)
    {
      Coef=NULL;
    }
  else if(Coefficient==NULL)
    {
      Coef=new Vector[M];
      for(int i=0;i<M;i++)
	{
	  Coef[i]=Vector(N);
	}
    }
  else
    {
  Coef=new Vector[M];
  for(int i=0;i<M;i++)
    {
      Coef[i]=Vector(N,Coefficient[i]);
    }
   }
}

Matrix::~Matrix()
{
  delete [] Coef;
}
int Matrix::getLength()const
{
  return length;
}
int Matrix::getWidth()const
{
  return width;
}
Matrix::Matrix(const Matrix& M):length(M.length),width(M.width)
{
  if(M.width==0 || M.length==0)
    {
      Coef=NULL;
    }
  else
    {
  Coef=new Vector[M.length];
  for(int i=0;i<M.length;i++)
    {
      Coef[i]=Vector(M.Coef[i]);
    }
   }
}

Matrix& Matrix::operator=(const Matrix& M)
{
  if(this!=&M)
    {
      for(int i=0;i<length;i++)
	{
	  Coef[i].~Vector();
	}
      // delete [] Coef;
     length=M.length;
      width=M.width;
      Coef=new Vector[M.length];
      for(int i=0;i<M.length;i++)
	{
	  Coef[i]=M.Coef[i];
	  }
    }
  else
    {
      cout<<"Affectation Impossible -> Adresse memoire identique"<<endl;
    }
  return *this;
}

Vector& Matrix::operator[](int j)
{
  return Coef[j];
}

Matrix Matrix::operator+(const Matrix& M)const
{
  if(length==M.length && width==M.width)
    {
      Matrix resultat(length,width);
      for(int i=0;i<length;i++)
	{
	  for(int j=0;j<width;j++)
	    {
	      resultat[i][j]=Coef[i][j]+M.Coef[i][j];
	    }
	}
      return resultat;
    }
  else
    {
      cout<<"Error Addition -> Incompatible dimension"<<endl;
      return ErrorDimension;
    }
}
Matrix Matrix::operator-(const Matrix& M)const
{
if(length==M.length && width==M.width)
    {
      Matrix resultat(length,width);
      for(int i=0;i<length;i++)
	{
	  for(int j=0;j<width;j++)
	    {
	      resultat[i][j]=Coef[i][j]-M.Coef[i][j];
	    }
	}
      return resultat;
    }
  else
    {
      cout<<"Error Addition -> Incompatible dimension"<<endl;
      return ErrorDimension;
    }
}
double Matrix::Max()const
{
double max=0;
for(int i=0;i<length;i++)
{
double temp=Coef[i].linfty();
if(max<temp)
{
max=temp;
}
}

return max;
}
Matrix Matrix::operator*(const Matrix& M)const
{
  if(width==M.length)
    {
      Matrix resultat(length,M.width);
      for(int i=0;i<length;i++)
	{
	  for(int j=0;j<M.width;j++)
	    {
	      double S=0;
	      for(int k=0;k<width;k++)
		{
		  S=S+Coef[i][k]*M.Coef[k][j];
		}
	      resultat.Coef[i][j]=S;
	    }
	}
      return resultat;
    }
  else
    {
cout<<"Error Multiplication -> Incompatible dimension"<<endl;
      return ErrorDimension;
    }
}
Matrix Matrix::transp()const
{
  if(length==0 || width==0)
    {
cout<<"Error Transp -> Incompatible dimension"<<endl;
      return ErrorDimension;
    }
  else
    {
      Matrix resultat(width,length);
      for(int i=0;i<width;i++)
	{
	  for(int j=0;j<length;j++)
	    {
	      resultat.Coef[i][j]=Coef[j][i];
	    }
	}
      return resultat;
    }
}

Matrix operator*(double k,const Matrix& A)
{
	Matrix Resultat(A.length,A.width);
	for(int i=0;i<A.length;i++)
	{
		for(int j=0;j<A.width;j++)
		{
			Resultat.Coef[i][j]=k*A.Coef[i][j];
		}
	}
	return Resultat;
}

Vector operator*(const Matrix&A,const Vector& X)
{
	Vector Resultat(A.length);
	//cout<<"Résultat ====="<<endl;
	//cout<<Resultat<<endl;
	if(A.width != X.getSize())
	{
		cout<<"Matrix Vector Product Errror -> Incompatible Dimension"<<endl;
		return NULL;
	}
	else
	{
		//Vector Resultat(A.length);
		for(int i=0;i<A.length;i++)
		{
			for(int j=0;j<A.width;j++)
			{
			Resultat[i]+=A.Coef[i][j]*X[j];
			}
			//cout<<Resultat[i]<<endl;
		}
		//cout<<Resultat<<endl;
		return Resultat;
	}
}

int Gauss(Matrix& A,Vector& P)
{
  if(A.length!=A.width)
    {
      cout<<"ErrorInversion -> Not a square Matrix"<<endl;
      return ErrorDimension;
    }
  else
    {
  double max,C;
  int stock,line;
  for(int j=0;j<A.width;j++)
    {
      //Test sur le pivot
      if(fabs(A[P[j]][j])<zero)
	{
	  //Ligne du pivot nul
	  line=P[j];
	  //Initialisation du max
	  max=A[P[j]][j];
	  for(int i=j;i<A.length;i++)
	    {
	      if(fabs(max)<fabs(A[P[i]][j]))
		{
		  max=fabs(A[P[i]][j]);
		  line=P[i];
		}
	    }
	  if(fabs(max)<zero)
	    {
	      cout<<"Error Singular Matrix"<<endl;
	      return 0;
	    }
	  stock=P[j];
	  P[j]=P[line];
	  P[line]=stock;
	}
      for(int i=j+1;i<A.length;i++)
	{
	  C=(A[P[i]][j])/(A[P[j]][j]);
	  A[P[i]][j]=C;
	  for(int k=j+1;k<A.length;k++)
	    {
	      A[P[i]][k]=A[P[i]][k]-C*(A[P[j]][k]);
	    }
	}
    }

  
     return 1;
     }
  return 0;
}

Vector Solve(Matrix& A,Vector& P,Vector& B)
{
  Vector resultat=Vector(A.length);
  for(int i=0;i<B.getSize();i++)
    {
    resultat[i]=B[P[i]];
    for(int k=0;k<i;k++)
      {
	resultat[i]=resultat[i]-(A[P[i]][k])*resultat[k];
      }
    }
  for(int i=B.getSize()-1;i>=0;i--)
    {
      for(int k=i+1;k<B.getSize();k++)
	{
	  resultat[i]=resultat[i]-(A[P[i]][k])*resultat[k];
	}
      resultat[i]=(resultat[i])/(A[P[i]][i]);
    }
  return resultat;
}


Vector TestGauss(Matrix& A,Vector& B)
{
  if(A.length != B.getSize())
    {
      cout<<"ErrorTestGauss -> Not a Square Matrix"<<endl;
      return ErrorDimension;
    }
  else
    {
      Vector P(B.getSize());
      for(int i=0;i<B.getSize();i++)
	{
	  P[i]=i;
	}
      Vector X(B.getSize());
      Matrix C(A);
      Vector BB(B);
      if(Gauss(C,P)!=0)
	{
	  X=Solve(C,P,B);
	  cout<<"Matrice du Systeme"<<endl;
	  cout<<A<<endl;
	  cout<<"Second Membre"<<endl;
	  cout<<BB<<endl;
	  cout<<"Solution"<<endl;
	  cout<<X<<endl;
	}
      return X;
    }
}

void ConjuguateGradient(Matrix& A,Vector& b,Vector& x0,Matrix& C,int itermax,double eps)
{
clock_t start,finish;
double duration;
  int n=b.getSize();
  assert(n==x0.getSize());
  Vector g(n),h(n),Ah(n),Cg(n);
  g=A*x0-b;
  Cg=C*g;
  h=-1*Cg;
  double g2=Cg*g;
 //cout<<"g2= "<<g2<<endl;
  double reps=eps*eps*g2;
 //cout<<"reps"<<reps<<endl;
 start=clock();
  for(int i=0;i<=itermax;i++)
    {
      Ah=A*h;
      double ro=-(g*h)/(h*Ah);
      x0=x0+ro*h;
      g=g+ro*Ah;
      Cg=C*g;
      double g2p=g2;
      g2=Cg*g;
     // cout<<i<<"ro ="<<ro<<"    ||g||^2= "<<g2<<endl;
      if(g2<reps)
        {
		finish=clock();
		duration=(double)(finish-start)/CLOCKS_PER_SEC;
			//cout<<"Convergence en "<<i<<" iterations"<<"    ro  = "<<ro<<"   ||g||^2===   "<<g2<<endl;
			//cout<<"Temps d'execution du gradient conjugue  "<<duration<<" s"<<endl;
          return;
        }
      double gamma=g2/g2p;
      h=gamma*h;
      h=h-Cg;
    }
	finish=clock();
	duration=(double)(finish-start)/CLOCKS_PER_SEC;
  cout<<"Non convergence de la méthode du gradient conjugue"<<endl;
  //cout<<"Temps d'execution du gradient conjugue  "<<duration<<" s"<<endl;
  return;

}

ostream& operator<<(ostream& flux, const Matrix&M)
{
  for(int i=0;i<M.length;i++)
    {
      flux<<M.Coef[i]<<endl;
    }
  flux<<endl;
  return flux;
}

