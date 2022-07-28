#ifndef VECTOR_H
#define VECTOR_H
enum ErrorType{ ErrorDimension=NULL};
class Vector
{
private:
  int size;
  double * Coord;

public:
  Vector(int=0,double* =NULL);
  Vector(const Vector&);
  ~Vector();
  int getSize()const;
  //Operation Elementaire sur les Vecteurs
  Vector& operator=(const Vector&); //Operateur d'affectation
  double& operator[](int)const; //Operateur de projection
  double operator*(const Vector&); //Produit Scalaire
  Vector operator+(const Vector&); //Addtion
  Vector operator-(const Vector&); //Soustraction
  Vector operator^(const Vector&); //Produit Vectoriel dans R3
  // Norme L2
  double Norm2();
  
  //Norme L8
  double linfty();
  //Min
  double min();
  //Max
  double max();
  //Operateur scalaire fois vecteur
  friend Vector operator*(const double&,const Vector&);
  //Affichage contenu
  friend std::ostream& operator<<(std::ostream&,const Vector&);
};
#endif
