#ifndef matrix_hpp
#define matrix_hpp
#include <iostream>
using namespace std ;

/*
#################################"Classe Vecteur"###################################
####################################################################################
*/

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
  void resize(int ni);
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

  //Fonction supplémentaire : allonge le vecteur d'une longueur 1, insère l'élément et décale les suivants
  vecteur& add(int i, double v);
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



/*
##########################"Classe Matrice"####################################
##############################################################################
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
    matrice& operator*(double a);
    matrice& operator/(double a);
    matrice& operator+(const matrice& A);
    matrice& operator-(const matrice& A);
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
####################################"Matrice profil"##################################
######################################################################################
*/

class matrice_profil //pour matrices à profil symétrique
{
public:
    int n; //on ne considère que des matrices carrées car à profil symétriques
    vecteur Profil;
    vecteur Posdiag;
    //vecteur Lower;
    //vecteur Upper;
    matrice_profil(int ni, const vecteur Pi); //constructeur de la matrice à partir du profil
    matrice_profil(const matrice_profil& A); //constructeur par copie
    //~matrice_profil();
};

/*
#############################"Matrice symétrique"#####################################
######################################################################################
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
    matrice_sym& operator*(double a);
    matrice_sym& operator/(double a);
    matrice_sym& operator+(const matrice_sym& A); //On suppose que les matrices ont le même profil
    matrice_sym& operator-(const matrice_sym& A); //On suppose que les matrices ont le même profil
};

ostream & operator<<(ostream & os, const matrice_sym& A);
void print(const matrice_sym& A);

matrice_sym operator*(const matrice_sym& A, double a);
matrice_sym operator/(const matrice_sym& A, double a);
matrice_sym operator*(double a, const matrice_sym& A);
vecteur operator*(const matrice_sym& A, const vecteur& V);
matrice_sym operator+(const matrice_sym& A, const matrice_sym& B); //On suppose que les matrices ont le même profil
matrice_sym operator-(const matrice_sym& A, const matrice_sym& B); //On suppose que les matrices ont le même profil

matrice_sym transpose(const matrice_sym& A);

/*
###################"Matrice non symétrique"#######################################
##################################################################################
*/

class matrice_nonsym : public matrice_profil
{
public:
    vecteur Lower;
    vecteur Upper;

    matrice_nonsym(int ni, const vecteur Pi); //constructeur de la matrice vide
    matrice_nonsym(const matrice_nonsym& A); //constructeur par copie
    matrice_nonsym(const matrice_sym& A);//constructeur par copie à partir d'une matrice symétrique
    //~matrice_nonsym();
    matrice_nonsym& operator=(const matrice_nonsym& A);
    double operator()(int i, int j) const;
    double& operator()(int i, int j);
    matrice_nonsym& operator*(double a);
    matrice_nonsym& operator/(double a);
    matrice_nonsym& operator+(const matrice_nonsym& A); //On suppose que les matrices ont le même profil
    matrice_nonsym& operator-(const matrice_nonsym& A); //On suppose que les matrices ont le même profil
};

ostream & operator<<(ostream & os, const matrice_nonsym& A);
void print(const matrice_nonsym& A);

matrice_nonsym operator*(const matrice_nonsym& A, double a);
matrice_nonsym operator/(const matrice_nonsym& A, double a);
vecteur operator*(const matrice_nonsym& A, const vecteur& V);
matrice_nonsym operator*(double a, const matrice_nonsym& A);
matrice_nonsym operator+(const matrice_nonsym& A, const matrice_nonsym& B); //On suppose que les matrices ont le même profil
matrice_nonsym operator-(const matrice_nonsym& A, const matrice_nonsym& B); //On suppose que les matrices ont le même profil

matrice_nonsym transpose(const matrice_nonsym& A);

/*
###########################"Opérations mixtes"######################################
####################################################################################
*/

//Opérations entre une amtrice symétrique et une non-symétrique
//Retourne nécessairement une matrice non symétrique
matrice_nonsym operator+(const matrice_nonsym& A, const matrice_sym& B); //On suppose que les matrices ont le même profil
matrice_nonsym operator-(const matrice_nonsym& A, const matrice_sym& B); //On suppose que les matrices ont le même profil
matrice_nonsym operator+(const matrice_sym& A, const matrice_nonsym& B); //On suppose que les matrices ont le même profil
matrice_nonsym operator-(const matrice_sym& A, const matrice_nonsym& B); //On suppose que les matrices ont le même profil

/*
################################ancienne version############################################

void LUdecomposition(const matrice& A);
//void LUdecomposition(const matrice_sym& A,matrice& l,matrice& u, int n);
//void LUdecomposition(const matrice_nonsym& A,matrice& l,matrice& u, int n);
*/
matrice_nonsym LUdecomposition(const matrice_nonsym& A);
matrice_nonsym LUdecomposition(const matrice_sym& A);

#endif
