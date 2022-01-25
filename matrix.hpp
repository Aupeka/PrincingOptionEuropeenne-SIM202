#ifndef matrix.hpp
#define matrix.hpp
#include <iostream>
using namespace std ;

/*
#################################"Classe Vecteur"###################################
####################################################################################
*/

// utilitaires
void stop(const char * msg);                     // message d'arr�t
void test_dim(int d1, int d2, const char * org); // test dimension

// classe vecteur de r�els double pr�cision
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
  void clear();            // d�sallocation
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

  //Fonction suppl�mentaire : allonge le vecteur d'une longueur 1, ins�re l'�l�ment et d�cale les suivants
  vecteur& add(int i, double v);
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



/*
##########################"Classe Matrice"####################################
##############################################################################
*/

class matrice
{
public:
    vecteur* cols_;
    int m,n;
    matrice(int mi=0, int ni=0; double vi=0.);
    matrice(const matrice& v);
    ~matrice();
    matrice& operator=(const matrice& A);
    double operator[](int i, int j) const;
    double& operator[](int i, int j);
};

vecteur produit(const matrice& A, const vecteur& u);

/*
####################################"Matrice profil"##################################
######################################################################################
*/

class matrice_profil //pour matrices � profil sym�trique
{
public:
    int n; //on ne consid�re que des matrices carr�es car � profil sym�triques
    vecteur Profil;
    vecteur Posdiag;
    //vecteur Lower;
    //vecteur Upper;
    matrice_profil(int ni); //constructeur de la matrice vide
    matrice_profil(const matrice_profil& A); //constructeur par copie
    ~matrice_profil();
    /*
    matrice_profil& operator=(const matrice_profil& A);
    double operator[](int i, int j) const;
    double& operator[](int i, int j);
    matrice_profil& operator*(double a);
    matrice_profil& operator/(double a);
    matrice_profil& operator+(const matrice_profil& A);
    matrice_profil& operator-(const matrice_profil& A);
    */
};

/*
matrice_profil& operator*(const matrice_profil& A, double a);
matrice_profil& operator*(const matrice_profil& A, double a);
matrice_profil& operator*(double a, const matrice_profil& A);
matrice_profil& operator+(const matrice_profil& A, const matrice_profil& B);
matrice_profil& operator-(const matrice_profil& A, const matrice_profil& B);
*/

/*
#############################"Matrice sym�trique"#####################################
######################################################################################
*/

class matrice_sym : public matrice_profil
{
public:
    vecteur Lower;

    matrice_sym(int ni); //constructeur de la matrice vide
    matrice_sym(const matrice_sym& A); //constructeur par copie
    ~matrice_sym();
    matrice_sym& operator=(const matrice_sym& A);
    double operator[](int i, int j) const;
    double& operator[](int i, int j);
    matrice_sym& operator*(double a);
    matrice_sym& operator/(double a);
    matrice_sym& operator+(const matrice_sym& A);
    matrice_sym& operator-(const matrice_sym& A);
};

matrice_sym operator*(const matrice_sym& A, double a);
matrice_sym operator*(const matrice_sym& A, double a);
matrice_sym operator*(double a, const matrice_sym& A);
matrice_sym operator+(const matrice_sym& A, const matrice_sym& B);
matrice_sym operator-(const matrice_sym& A, const matrice_sym& B);


/*
###################"Matrice non sym�trique"#######################################
##################################################################################
*/

class matrice_nonsym : public matrice_profil
{
public:
    vecteur Lower;
    vecteur Upper;

    matrice_nonsym(int ni); //constructeur de la matrice vide
    matrice_nonsym(const matrice_nonsym& A); //constructeur par copie
    matrice_nonsym(const matrice_sym& A);//constructeur par copie � partir d'une matrice sym�trique
    ~matrice_nonsym();
    matrice_nonsym& operator=(const matrice_nonsym& A);
    double operator[](int i, int j) const;
    double& operator[](int i, int j);
    matrice_nonsym& operator*(double a);
    matrice_nonsym& operator/(double a);
    matrice_nonsym& operator+(const matrice_nonsym& A);
    matrice_nonsym& operator-(const matrice_nonsym& A);
};

matrice_nonsym operator*(const matrice_nonsym& A, double a);
matrice_nonsym operator*(const matrice_nonsym& A, double a);
matrice_nonsym operator*(double a, const matrice_nonsym& A);
matrice_nonsym operator+(const matrice_nonsym& A, const matrice_nonsym& B);
matrice_nonsym operator-(const matrice_nonsym& A, const matrice_nonsym& B);


/*
###########################"Op�rations mixtes"######################################
####################################################################################
*/

//Op�rations entre une amtrice sym�trique et une non-sym�trique
//Retourne n�cessairement une matrice non sym�trique
matrice_nonsym operator+(const matrice_nonsym& A, const matrice_sym& B);
matrice_nonsym operator-(const matrice_nonsym& A, const matrice_sym& B);
matrice_nonsym operator+(const matrice_sym& A, const matrice_nonsym& B);
matrice_nonsym operator-(const matrice_sym& A, const matrice_nonsym& B);


#endif