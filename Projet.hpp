#ifndef Projet.hpp
#define Projet.hpp
#include <iostream>
using namespace std ;




#############################################################
###################### Classe Matrice #######################
#############################################################





#############################################################
###################### Classe Maillage ######################
#############################################################

class Point
{public :
    double x;
    double y;
    Point(double x0=0,double y0=0):x(x0) ,y(y0) {}
    Point& operator +=( const Point &P) {x+=P.x;y+=P.y;return *this;}
    Point& operator -=( const Point &P) {x-=P.x;y-=P.y;return *this;}
    Point& operator *=( double a) {x*=a;y*=a;return *this;}
    Point& operator /=( double a)
    {
        if (a==0)
        {
            exit(-1);
        }
        else
        {
            x/=a;
            y/=a;
        }
        return *this;
    }
};

Point operator +(const Point & P, const Point & Q);
Point operator -(const Point & P, const Point & Q);
Point operator *(const Point & P,double a);
Point operator *(double a,const Point & P);
Point operator /(const Point & P, double a);
bool operator ==(const Point & P, const Point & Q);
bool operator !=(const Point & P, const Point & Q);
ostream & operator<<(ostream & out, const Point & P);


class Numeros : public vector<int>
{public :
    Numeros(int i1 , int i2 , int i3)
    {resize(3);
     vector<int>::iterator it =begin();
     *it++=i1; *it++=i2; *it=i3;}
    int operator()(int i) const {return (*this)[i-1];}
    int& operator()(int i){return (*this)[i-1];}
};

ostream & operator<<(ostream& out, const Numeros & N);


class Maillage
{public :
    vector<Point> sommets;
    list<Numeros> numelts;
    Maillage(int m,int n)
    {maille_carre_unite(m,n);}
    Maillage(double a,double b,double c,double d, int m, int n)
    {maille_rectangle(a,b,c,d,m,n);}
    void maille_carre_unite(int m, int n);
    void maille_rectangle(double a,double b, double c, double d, int m, int n);
    void affiche() const;
    Maillage& tf_affine(const vector<double> &,const vector<double> &);
    Maillage & operator +=(const Maillage &);
    void save(const char *fn) const;
};
Maillage operator +(const Maillage &, const Maillage &);


#############################################################
###################### Classe Vecteur #######################
#############################################################

// utilitaires
void stop(const char * msg);                     // message d'arrêt
void test_dim(int d1, int d2, const char * org); // test dimension

// classe vecteur de réels double précision
class vecteur
{
private :
  int dim_;          // dimension du vecteur
  double * val_;     // tableaux de valeurs

public:
  vecteur(int d=0, double v0=0); // dim et val constante
  vecteur(const vecteur & v);    // constructeur par copie
  ~vecteur();

  // tools
  void init(int d);        // allocation
  void clear();            // désallocation
  int dim() const {return dim_;}                // accès dimension

  // opérateur d'assignation
  vecteur & operator=(const vecteur & v);  // affectation d'un vecteur
  vecteur & operator=(double x);           // affectation d'une valeur

  // opérateurs d'accès (pour les utilisateurs)
  double  operator [](int i) const {return val_[i];}   // valeur    0->dim-1
  double& operator [](int i) {return val_[i];}         // référence 0->dim-1
  double  operator ()(int i) const {return val_[i-1];} // valeur    1->dim
  double& operator ()(int i) {return val_[i-1];}       // référence 1->dim
  vecteur operator ()(int i, int j) const;             // valeurs i à j (0<i<=j<=dim)

  // opérateurs algébriques
  vecteur& operator +=(const vecteur & v);             // u += v
  vecteur& operator -=(const vecteur & v);             // u -= v
  vecteur& operator +=(double x);                      // u += x
  vecteur& operator -=(double x);                      // u -= x
  vecteur& operator *=(double x);                      // u *= x
  vecteur& operator /=(double x);                      // u /= x

}; // fin de définition de la classe

// opérateurs externes
vecteur operator +(const vecteur & u); //+ unaire ne fait rien !
vecteur operator -(const vecteur & u); //- unaire : chgt de signe
vecteur operator +(const vecteur & u,const vecteur & v); // u + v
vecteur operator -(const vecteur & u,const vecteur & v); // u - v
vecteur operator +(const vecteur & u,double x);          // u + x
vecteur operator +(double x, const vecteur & u);         // x + u
vecteur operator -(const vecteur & u,double x);          // u - x
vecteur operator -(double x, const vecteur & u);         // x - u
vecteur operator *(const vecteur & u,double x);          // u * x
vecteur operator /(const vecteur & u,double x);          // u / x
vecteur operator +(double x,const vecteur & u);          // x + u
vecteur operator -(double x,const vecteur & u);          // x - u
vecteur operator *(double x,const vecteur & u);          // x * u
double operator |(const vecteur & u,const vecteur & v);  // u | v
vecteur operator,(const vecteur & u,const vecteur & v);  // [u1,...,um,v1,...,vn]

// opérateurs de comparaison
bool operator == (const vecteur & u,const vecteur & v);  // u == v
bool operator != (const vecteur & u,const vecteur & v);  // u != v

// opérateurs de lecture et d'écriture
istream & operator>>(istream & is, vecteur & u);         // is >> u
ostream & operator<<(ostream & os, const vecteur & u);   // os << u

//opérations composante à composante  u*~v et u/~v
class operande
{
 public:
  const vecteur *vp;
  operande(const vecteur *p=0) : vp(p) {}
};

operande operator~(const vecteur & u);                     // &u (unaire) -> opu
vecteur operator*(const vecteur & u, const operande & ov); // u * opv
vecteur operator/(const vecteur & u, const operande & ov); // u / opv










#endif
