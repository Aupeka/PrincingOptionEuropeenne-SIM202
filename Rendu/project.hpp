#ifndef Project.hpp
#define Project.hpp
#include <iostream>
#include <vector>
#include <list>
#include <map>
using namespace std ;




/*
#############################################################
###################### Classe Maillage ######################
#############################################################
*/
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
    Point& tf_affine(const vector<double> &,const vector<double> &);
};

Point operator +(const Point & P, const Point & Q);
Point operator -(const Point & P, const Point & Q);
Point operator *(const Point & P,double a);
Point operator *(double a,const Point & P);
Point operator /(const Point & P, double a);
bool operator ==(const Point & P, const Point & Q);
bool operator !=(const Point & P, const Point & Q);
bool operator <(const Point & P, const Point & Q);
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
    vector<Point> sommets; //=Coorneu
    list<Numeros> numelts; //=Numtri
    Maillage(int m,int n)
    {maille_carre_unite(m,n);}
    Maillage(double a,double b,double c,double d, int m, int n)
    {maille_rectangle(a,b,c,d,m,n);}
    void maille_carre_unite(int m, int n);
    void maille_rectangle(double a,double b, double c, double d, int m, int n);
    void affiche() const;
    Maillage& affine(const vector<double> &,const vector<double> &);
    Maillage & operator +=(const Maillage &);
    void savecoord(const char *fn) const; //Coordon�es de tous les points du maillage et num�ro associ� de ces sommets
    void savenumtri(const char *fn) const; // Num�ros des sommets constituants chaque triangle et num�rotation des triangles

};
Maillage operator +(const Maillage &, const Maillage &);

/*
#############################################################
###################### Classe Vecteur #######################
#############################################################
*/

// utilitaires
void stop(const char * msg);                     // message d'arr�t
void test_dim(int d1, int d2, const char * org); // test dimension

// classe vecteur de r�els double pr�cision
class vecteur
{
public :
  int dim_;          // dimension du vecteur
  double * val_;     // tableaux de valeurs

  vecteur(int d=0, double v0=0); // dim et val constante
  vecteur(const vecteur & v);    // constructeur par copie
  ~vecteur();

  // tools
  void init(int d);        // allocation
  void clear();
  void resize(int ni);          // d�sallocation
  int dim() const {return dim_;}                // acc�s dimension

  // op�rateur d'assignation
  vecteur & operator=(const vecteur & v);  // affectation d'un vecteur
  vecteur & operator=(double x);           // affectation d'une valeur

  // op�rateurs d'acc�s (pour les utilisateurs)
  double  operator [](int i) const {return val_[i];}   // valeur    0->dim-1
  double& operator [](int i) {return val_[i];}         // r�f�rence 0->dim-1
  double  operator ()(int i) const {return val_[i-1];} // valeur    1->dim
  double& operator ()(int i) {return val_[i-1];}       // r�f�rence 1->dim
  vecteur operator ()(int i, int j) const;             // valeurs i � j (0<i<=j<=dim)

  // op�rateurs alg�briques
  vecteur& operator +=(const vecteur & v);             // u += v
  vecteur& operator -=(const vecteur & v);             // u -= v
  vecteur& operator +=(double x);                      // u += x
  vecteur& operator -=(double x);                      // u -= x
  vecteur& operator *=(double x);                      // u *= x
  vecteur& operator /=(double x);                      // u /= x

}; // fin de d�finition de la classe

// op�rateurs externes
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

// op�rateurs de comparaison
bool operator == (const vecteur & u,const vecteur & v);  // u == v
bool operator != (const vecteur & u,const vecteur & v);  // u != v

// op�rateurs de lecture et d'�criture
istream & operator>>(istream & is, vecteur & u);         // is >> u
ostream & operator<<(ostream & os, const vecteur & u);   // os << u

//op�rations composante � composante  u*~v et u/~v
class operande
{
 public:
  const vecteur *vp;
  operande(const vecteur *p=0) : vp(p) {}
};

operande operator~(const vecteur & u);                     // &u (unaire) -> opu
vecteur operator*(const vecteur & u, const operande & ov); // u * opv
vecteur operator/(const vecteur & u, const operande & ov); // u / opv

/*
#############################################################
###################### Classe En plus #######################
#############################################################
*/

//Numeros get(list<Numeros>, int);

int minimum(int, int);


/*
#############################################################
###################### Classe Matrice #######################
#############################################################
*/

class matrice
{
public:
    vecteur* cols_;
    int m,n;
    matrice(int mi, int ni, double v=0.);
    matrice(const matrice& v);
    ~matrice();
    matrice& operator=(const matrice& A);
    double operator()(int i, int j) const;
    double& operator()(int i, int j);
    matrice& operator*=(double a);
    matrice& operator/=(double a);
    matrice& operator+=(const matrice& A);
    matrice& operator-=(const matrice& A);
};

matrice operator*(const matrice& A, double a);
matrice operator/(const matrice& A, double a);
matrice operator*(double a, const matrice& A);
matrice operator+(const matrice& A, const matrice& B);
matrice operator-(const matrice& A, const matrice& B);

vecteur operator*(const matrice& A, const vecteur& u);
matrice operator*(const matrice& A, const matrice& B);
matrice transpose(const matrice& A);
ostream & operator<<(ostream & os, const matrice& A);
/*
#############################################################
###################### Matrice Profil #######################
#############################################################
*/

class matrice_profil //pour matrices � profil sym�trique
{
public:
    int n; //on ne consid�re que des matrices carr�es car � profil sym�triques
    vecteur Profil;
    vecteur Posdiag;
    //vecteur Lower;
    //vecteur Upper;
    matrice_profil(int ni, const vecteur Pi); //constructeur de la matrice � partir du profil
    matrice_profil(const matrice_profil& A); //constructeur par copie
    //~matrice_profil();
};

/*
#############################################################
################### Matrice Sym�trique ######################
#############################################################
*/

class matrice_sym : public matrice_profil
{
public:
    vecteur Lower;

    matrice_sym(int ni, const vecteur Pi); //constructeur de la matrice
    matrice_sym(const matrice_sym& A); //constructeur par copie
    //~matrice_sym();
    matrice_sym& operator=(const matrice_sym& A);
    double operator()(int i, int j) const;
    double& operator()(int i, int j);
    matrice_sym& operator*=(double a);
    matrice_sym& operator/=(double a);
    matrice_sym& operator+=(const matrice_sym& A); //On suppose que les matrices ont le m�me profil
    matrice_sym& operator-=(const matrice_sym& A); //On suppose que les matrices ont le m�me profil
};

ostream & operator<<(ostream & os, const matrice_sym& A);
void print(const matrice_sym& A);

matrice_sym operator*(const matrice_sym& A, double a);
matrice_sym operator/(const matrice_sym& A, double a);
matrice_sym operator*(double a, const matrice_sym& A);
vecteur operator*(const matrice_sym& A, const vecteur& V);
matrice_sym operator+(const matrice_sym& A, const matrice_sym& B); //On suppose que les matrices ont le m�me profil
matrice_sym operator-(const matrice_sym& A, const matrice_sym& B); //On suppose que les matrices ont le m�me profil

matrice_sym transpose(const matrice_sym& A);

/*
#############################################################
################ Matrice non-sym�trique #####################
#############################################################
*/

class matrice_nonsym : public matrice_profil
{
public:
    vecteur Lower;
    vecteur Upper;

    matrice_nonsym(int ni, const vecteur Pi); //constructeur de la matrice vide
    matrice_nonsym(const matrice_nonsym& A); //constructeur par copie
    matrice_nonsym(const matrice_sym& A);//constructeur par copie � partir d'une matrice sym�trique
    //~matrice_nonsym();
    matrice_nonsym& operator=(const matrice_nonsym& A);
    double operator()(int i, int j) const;
    double& operator()(int i, int j);
    matrice_nonsym& operator*=(double a);
    matrice_nonsym& operator/=(double a);
    matrice_nonsym& operator+=(const matrice_nonsym& A); //On suppose que les matrices ont le m�me profil
    matrice_nonsym& operator-=(const matrice_nonsym& A); //On suppose que les matrices ont le m�me profil
};

ostream & operator<<(ostream & os, const matrice_nonsym& A);
void print(const matrice_nonsym& A);

matrice_nonsym operator*(const matrice_nonsym& A, double a);
matrice_nonsym operator/(const matrice_nonsym& A, double a);
vecteur operator*(const matrice_nonsym& A, const vecteur& V);
matrice_nonsym operator*(double a, const matrice_nonsym& A);
matrice_nonsym operator+(const matrice_nonsym& A, const matrice_nonsym& B); //On suppose que les matrices ont le m�me profil
matrice_nonsym operator-(const matrice_nonsym& A, const matrice_nonsym& B); //On suppose que les matrices ont le m�me profil

matrice_nonsym transpose(const matrice_nonsym& A);

/*
#############################################################
##################### Op�rations Mixtes #####################
#############################################################
*/

//Op�rations entre une amtrice sym�trique et une non-sym�trique
//Retourne n�cessairement une matrice non sym�trique
matrice_nonsym operator+(const matrice_nonsym& A, const matrice_sym& B); //On suppose que les matrices ont le m�me profil
matrice_nonsym operator-(const matrice_nonsym& A, const matrice_sym& B); //On suppose que les matrices ont le m�me profil
matrice_nonsym operator+(const matrice_sym& A, const matrice_nonsym& B); //On suppose que les matrices ont le m�me profil
matrice_nonsym operator-(const matrice_sym& A, const matrice_nonsym& B); //On suppose que les matrices ont le m�me profil
matrice_nonsym LUdecomposition(const matrice_nonsym& A);
matrice_nonsym LUdecomposition(const matrice_sym& A);


/*
#############################################################
################## Matrices Elementaires ####################
#############################################################
*/

matrice matM_elem(const Point& P1,const Point& P2,const Point& P3);
matrice matK_elem(const Point& P1,const Point& P2,const Point& P3);
matrice A(const Point& P);
vecteur V(const Point& P);
matrice matB_elem(const Point& P1,const Point& P2,const Point& P3);
int minimum(int a, int b);
matrice_sym matM(const vector<Point>& Coorneu, const list<Numeros>& Numtri);
matrice_sym matK(const vector<Point>& Coorneu, const list<Numeros>& Numtri);
matrice_nonsym matB(const vector<Point>& Coorneu, const list<Numeros>& Numtri);
matrice_nonsym matD(const vector<Point>& Coorneu, const list<Numeros>& Numtri, double r);

/*
#############################################################
####################### Resolution ##########################
#############################################################
*/

vecteur resolsys(const matrice_nonsym A, const vecteur b);
double pos(double x);
vecteur Q(double K,const vector<Point> & sommets);
vecteur resolution_1(double deltaT, double K, const vector<Point> & sommets, const list<Numeros> & numelts);
vecteur resolution_2(double deltaT, double K, const vector<Point> & sommets, const list<Numeros> & numelts);


#endif
