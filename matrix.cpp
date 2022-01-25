#include "matrix.hpp"
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <cstdlib>
using namespace std;

/*
#################"Classe vecteur"#################
#################################################
*/

// utilitaire de messages d'erreur
void stop(const char * msg)
{
  cout << "ERREUR classe vecteur, " << msg;
  exit(-1);
}

void test_dim(int d1, int d2, const char * org)
{
  if (d1==d2)  return;
  cout << org << " (" << d1 << "," << d2 << ") dimensions incompatibles";
  exit(-1);
}


// constructeurs
vecteur::vecteur(int d, double v0) // dim et val constante
{
  init(d);
  for (int i=0;i<dim_;i++) val_[i]=v0;
}

vecteur::vecteur(const vecteur & v) // constructeur par copie
{
  init(v.dim_);
  for (int i=0;i<dim_;i++) val_[i]=v[i];
}

vecteur::~vecteur()
{ clear(); }

void vecteur::init(int d) // initialisation avec allocation dynamique
{
  if (d<=0) stop("init() : dimension <= 0");
  dim_=d;
  val_ = new double[d];
}

void vecteur::clear()    // désallocation
{
  if (val_!=0) delete [] val_;
  dim_=0;
}

// affectation
vecteur & vecteur::operator=(const vecteur & v)  // affectation d'un vecteur
{
  if (dim_!=v.dim_) // redimensionnement
  {
    clear();
    init(v.dim_);
  }
  // recopie
  for (int i=0;i<dim_;i++) val_[i]=v[i];
  return *this;
}
vecteur & vecteur::operator=(double x)  // affectations d'une valeur
{
  for (int i=0;i<dim_;i++) val_[i]=x;
  return *this;
}

// opérateur d'accès à une plage i à j (1<i<=j<=dim)
vecteur vecteur::operator ()(int i, int j) const
{
  if (i<1 || i>j || j>dim_) stop("plage inconsistante");
  vecteur u(j-i+1);
  for (int k=0;k<=j-i;k++) u[k]=val_[i+k-1];
  return u;
}

// opérateurs algébriques
vecteur& vecteur::operator +=(const vecteur & v)
{
  test_dim(dim_,v.dim_,"op +=");
  for (int i=0;i<dim_;i++) val_[i]+=v[i];
  return *this;
}
vecteur& vecteur::operator -=(const vecteur & v)
{
  test_dim(dim_,v.dim_,"op -=");
  for (int i=0;i<dim_;i++) val_[i]-=v[i];
  return *this;
}
vecteur& vecteur::operator +=(double x)
{
  for (int i=0;i<dim_;i++) val_[i]+=x;
  return *this;
}
vecteur& vecteur::operator -=(double x)
{
  for (int i=0;i<dim_;i++) val_[i]-=x;
  return *this;
}
vecteur& vecteur::operator *=(double x)
{
  for (int i=0;i<dim_;i++) val_[i]*=x;
  return *this;
}
vecteur& vecteur::operator /=(double x)
{
  if (x==0) stop("op /= : division par 0");
  for (int i=0;i<dim_;i++) val_[i]/=x;
  return *this;
}

//Fonction supplémentaire : allonge le vecteur d'une longueur 1, insère l'élément et décale les suivants
vecteur& vecteur::add(int i, double v)
{
    dim_+=1;
    vecteur temp(dim_);
    for (int k=0; k<dim_;++k)
    {
        if (k!=i)
        {
            temp[k] = val_[k];
        }
        else
        {
            temp[k] = v;
        }
    }
    delete [] val_;
    val_ = temp.val_;
    return *this;
}


// opérateurs externes
vecteur operator +(const vecteur & u) //+ unaire (ne fait rien)!
{ return u; }
vecteur operator -(const vecteur & u) //- unaire : chgt de signe
{ vecteur w=u; return w*=-1.; }
vecteur operator +(const vecteur & u,const vecteur & v)
{ vecteur w=u; return w+=v; }
vecteur operator -(const vecteur & u,const vecteur & v)
{ vecteur w=u; return w-=v; }
vecteur operator +(const vecteur & u,double x)
{ vecteur w=u; return w+=x; }
vecteur operator -(const vecteur & u,double x)
{ vecteur w=u; return w-=x; }
vecteur operator *(const vecteur & u,double x)
{ vecteur w=u; return w*=x; }
vecteur operator /(const vecteur & u,double x)
{ vecteur w=u; return w/=x; }
vecteur operator +(double x,const vecteur & u)
{ vecteur w=u; return w+=x; }
vecteur operator -(double x,const vecteur & u)
{ vecteur w=-u; return w-=x; }
vecteur operator *(double x,const vecteur & u)
{ vecteur w=u; return w*=x; }

// opérateurs de comparaison
bool operator == (const vecteur & u,const vecteur & v)
{
  if (u.dim()!=v.dim()) return false;
  for (int i=0;i<u.dim();i++)
    if (u[i]!=v[i]) return false;
  return true;
}
bool operator != (const vecteur & u,const vecteur & v)
{
  return !(u==v);
}

// produit scalaire
double operator |(const vecteur & u,const vecteur & v)
{
  test_dim(u.dim(),v.dim(),"operateur |");
  double ps=0.;
  for (int i=0;i<u.dim();i++) ps+=u[i]*v[i];
  return ps;
}

// opérateur de concaténation
vecteur operator,(const vecteur &u, const vecteur& v)
{
  vecteur uv(u.dim()+v.dim());
  int k=0;
  for (int i=0;i<u.dim();i++, k++)uv[k]=u[i];
  for (int i=0;i<v.dim();i++, k++)uv[k]=v[i];
  return uv;
}

// opérateurs de lecture et d'écriture
istream & operator>>(istream & is, vecteur & u)
{
  for (int i=0;i<u.dim();i++) is>>u[i];
  return is;
}
ostream & operator<<(ostream & os, const vecteur & u)
{
  os << "(";
  for (int i=0;i<u.dim()-1;i++) os << u[i] << ",";
  os << u[u.dim()-1] << ")";
  return os;
}

/*
############################"Classe Matrice"###########################
#######################################################################
*/

matrice::matrice(int mi, int ni, double v) : cols_(0), m(mi), n(ni)
{
    if (n>0 || m<0)
    {
        cout<<"dimension négative"<<endl;
        exit(-1);
    }
    cols_ = new vecteur[n];
    vecteur V(m,v);
    for (int j=0; j<n; ++j)
    {
        cols_[j] = V;
    }
}

matrice::matrice(const matrice& A): cols_(0), m(A.m), n(A.n)
{
    if (A.cols_!=0)
    {
        cols_ = new vecteur[n];
    }
    for (int j=0; j<n; ++j)
    {
        cols_[j] = A.cols_[j];
    }
}

matrice& matrice::operator=(const matrice& A)
{
    if ((m!=A.m || n!=A.n) && cols_!=0)
    {
        delete [] cols_;
        cols_ = new vecteur[A.n];
        m = A.m;
        n = A.n;
    }
    for (int j=0; j<n; ++j)
    {
        cols_[j] = A.cols_[j];
    }
    return *this;
}

double matrice::val(int i, int j) const
{
    if (i>m || j>n || i<0 || j<0)
    {
        cout <<"Indice en dehors des bornes"<<endl;
        exit(-1);
    }
    return cols_[j][i]; //surcharge du vecteur
}

double& matrice::val(int i, int j)
{
    if (i>m || j>n || i<0 || j<0)
    {
        cout <<"Indice en dehors des bornes"<<endl;
        exit(-1);
    }
    return cols_[j][i]; //surcharge du vecteur
}

vecteur produit(const matrice& A, const vecteur& u)
{
    if (A.n!=u.dim())
    {
        cout<<"produit matrice vecteur : dimensions incompatibles"<<endl;
        exit(-1);
    }
    vecteur Res(A.m,0.);
    for (int i=1; i<<A.m; ++i)
    {
        for (int j=1; j<A.n;++j)
        {
           Res[i] = A.val(i,j)*u[j];
        }
    }
    return Res;
}

/*
###################################################"Matrice profil"####################################################
#######################################################################################################################
*/

matrice_profil::matrice_profil(int ni) : n(ni), Profil(0), Posdiag(0)
{
    if (n<0)
    {
        cout<<"dimension négative"<<endl;
        exit(-1);
    }
}


 matrice_profil::matrice_profil(const matrice_profil& A):n(A.n)
 {
     int d_prof = A.Profil.dim();
     int d_pos = A.Posdiag.dim();
     if (d_prof!=0)
    {
        Profil = vecteur(d_prof);
        Posdiag = vecteur(d_pos);
    }
    for (int k=0; k<d_prof; ++k)
    {
        Profil[k] = A.Profil[k];
    }
    for (int k=0; k<d_pos; ++k)
    {
        Posdiag[k] = A.Posdiag[k];
    }
 }

  matrice_profil::~matrice_profil()
 {
     if (Profil!=0)
     {
         delete [] Profil;
     }
     if (Posdiag!=0)
     {
         delete [] Posdiag;
     }
     n=0;
 }

/*
 #####################################"Matrice sym"########################################################
 ##########################################################################################################
*/

 matrice_sym::matrice_sym(int ni) : matrice_profil(ni), Lower(0) {}

 matrice_sym::matrice_sym(const matrice_sym& A):matrice_profil(A.n)
 {
     int d_prof = A.Profil.dim();
     int d_pos = A.Posdiag.dim();
     int d_low = A.Lower.dim();
     if (d_prof!=0)
    {
        Profil = vecteur(d_prof);
        Posdiag = vecteur(d_pos);
        Lower = vecteur(d_low);
    }
    for (int k=0; k<d_prof; ++k)
    {
        Profil[k] = A.Profil[k];
    }
    for (int k=0; k<d_pos; ++k)
    {
        Posdiag[k] = A.Posdiag[k];
    }
    for (int k=0; k<d_low; ++k)
    {
        Lower[k] = A.Lower[k];
    }
 }

 matrice_sym::~matrice_sym()
 {
     if (Profil!=0)
     {
         delete [] Profil;
     }
     if (Posdiag!=0)
     {
         delete [] Posdiag;
     }
     if (Lower!=0)
     {
         delete [] Lower;
     }
     n = 0;
 }

matrice_sym& matrice_sym::operator=(const matrice_sym& A)
{
    int d_prof = A.Profil.dim();
    int d_pos = A.Posdiag.dim();
    int d_low = A.Lower.dim();
    delete [] Profil;
    delete [] Posdiag;
    delete [] Lower;
    Profil = vecteur(d_prof);
    Posdiag = vecteur(d_pos);
    Lower = vecteur(d_low);
    if (n!=A.n)
    {
        n = A.n;
    }
    for (int k=0; k<d_prof; ++k)
    {
        Profil[k] = A.Profil[k];
    }
    for (int k=0; k<d_pos; ++k)
    {
        Posdiag[k] = A.Posdiag[k];
    }
    for (int k=0; k<d_low; ++k)
    {
        Lower[k] = A.Lower[k];
    }
    return *this;
}

double matrice_sym::val(int i, int j) const //Opérateur de lecture
{
    if (i==j)
    {
        return Lower[Posdiag[i]];
    }
    if (i>j) //triangle inférieur
    {
        if (j>=Profil[i]) //après le premier terme non-nul de la ligne
        {
            return Lower[Posdiag[i-1]+(j-Profil[i]+1)];
        }
        else //avant le premier terme non-nul de la ligne
        {
            return 0.;
        }
    }
    else //triangle supérieur
    {
        if (j>=Profil[i]) //après le premier terme non-nul de la colonne
        {
            return Lower[Posdiag[i-1]+(j-Profil[i]+1)];
        }
        else//avant le premier terme non-nul de la colonne
        {
            return 0.;
        }
    }
}


double& matrice_sym::val(int i, int j) //Lecture et écriture
{
    if (i==j)
    {
        return Lower[Posdiag[i]]; //n'incrémente que dans le lower
    }
    if (i>j) //triangle inférieur
    {
        if (j>=Profil[i]) //après le premier terme non-nul de la ligne
        {
            return Lower[Posdiag[i-1]+(j-Profil[i]+1)];
        }
        else //avant le premier terme non-nul de la ligne
        {
            //On doit l'ajouter au bon endroit
            for (int k=i; k<Profil[i];++k)//termes à rajouter
            {
                Lower.add(Posdiag[i-1]+k+1,0.); //rajoute des 0 entre le nouveau point à écrire et celui existant //Fonction à implémenter (indice, valeur)
                Posdiag[i]+=1;//décale le terme diagonal de la ligne à chaque ajout d'un terme
            }
            Profil[i]=j; //nouveau premier terme non-nul
            return Lower[Posdiag[i-1]+Profil[i]+1];
        }
    }
    else //triangle supérieur
    {
        if (j>=Profil[i]) //après le premier terme non-nul de la colonne
        {
            return Lower[Posdiag[i-1]+(j-Profil[i]+1)];
        }
        else//avant le premier terme non-nul de la colonne
        {
            //On doit l'ajouter au bon endroit
            for (int k=i; k<Profil[i];++k)//termes à rajouter
            {
                Lower.add(Posdiag[i-1]+k+1,0.); //rajoute des 0 entre le nouveau point à écrire et celui existant //Fonction à implémenter (indice, valeur)
                Posdiag[i]+=1;//décale le terme diagonal de la ligne à chaque ajout d'un terme
            }
            Profil[i]=j; //nouveau premier terme non-nul
            return Lower[Posdiag[i-1]+Profil[i]+1];
        }
    }
}

matrice_sym& matrice_sym::operator*(double a)
{
    Lower*=a;
    return *this;
}

matrice_sym& matrice_sym::operator/(double a)
{
    if (a==0.)
    {
        cout<<"division par zéro"<<endl;
        exit(-1);
    }
    Lower/=a;
    return *this;
}

matrice_sym& matrice_sym::operator+(const matrice_sym& A)
{
    if (n=!A.n)
    {
        cout<<"Les matrices ne sont pas de mêmes dimensions"<<endl;
        exit(-1);
    }
    for (int i=0; i<n; ++i)
    {
        (*this).val(i,i)+=A.val(i,i); //n'incrémente que dans le lower
        //Lower[Posdiag[i]]+=A.Lower[Posdiag[i]];

        for (int j=A.Profil[i]; j<i; ++j) //du premier non-nul jusqu'à la diagonale exclu
        {
            (*this).val(i,j)+=A.val(i,j);
        }

        /*
        if (A.profil[i]>=profil[i])//A n'enlève pas de 0 sur cette ligne
        {
            for (int j=A.Profil[i]; j<i; ++j) //du premier non-nul jusqu'à la diagonale exclu
            {
                (this)[i,j]+=A[i,j];
                (this)[j,i]+=A[j,i];
                //Lower[Posdiag[i-1]+j+1]+=A.Lower[A.Posdiag[i-1]+j+1];
                //Upper[Posdiag[i-1]+j+1]+=A.Upper[A.Posdiag[i-1]+j+1];
            }
        }
        else//A modifie la ligne
        {
            for (int j=A.profil[i]; j<i; ++j)
            {
                (this)[i,j]+=A[i,j];
                (this)[j,i]+=A[j,i];
                //
            }
        }
        */
    }
    return *this;
}

matrice_sym& matrice_sym::operator-(const matrice_sym& A)
{
    if (n=!A.n)
    {
        cout<<"Les matrices ne sont pas de mêmes dimensions"<<endl;
        exit(-1);
    }
    for (int i=0; i<n; ++i)
    {
        (*this).val(i,i)-=A.val(i,i); //n'incrémente que dans le lower
        //Lower[Posdiag[i]]+=A.Lower[Posdiag[i]];

        for (int j=A.Profil[i]; j<i; ++j) //du premier non-nul jusqu'à la diagonale exclu
        {
            (*this).val(i,j)-=A.val(i,j);
        }
    }
    return *this;
}

matrice_sym operator*(const matrice_sym& A, double a)
{
    return matrice_sym(A)*a;
}

matrice_sym operator*(double a, const matrice_sym& A)
{
    return matrice_sym(A)*a;
}

matrice_sym operator/(const matrice_sym& A, double a)
{
    return matrice_sym(A)/a;
}

matrice_sym operator+(const matrice_sym& A, const matrice_sym& B)
{
    return matrice_sym(A)+B;
}

matrice_sym operator-(const matrice_sym& A, const matrice_sym& B)
{
    return matrice_sym(A)-B;
}

/*
#####################################"Matrice non sym"#####################################################
###########################################################################################################
*/

 matrice_nonsym::matrice_nonsym(int ni) : matrice_profil(ni), Lower(0), Upper(0) {}

 matrice_nonsym::matrice_nonsym(const matrice_nonsym& A):matrice_profil(A.n)
 {
     int d_prof = A.Profil.dim();
     int d_pos = A.Posdiag.dim();
     int d_low = A.Lower.dim();
     int d_up = A.Upper.dim();
     if (d_prof!=0)
    {
        Profil = vecteur(d_prof);
        Posdiag = vecteur(d_pos);
        Lower = vecteur(d_low);
        Upper = vecteur(d_up);
    }
    for (int k=0; k<d_prof; ++k)
    {
        Profil[k] = A.Profil[k];
    }
    for (int k=0; k<d_pos; ++k)
    {
        Posdiag[k] = A.Posdiag[k];
    }
    for (int k=0; k<d_low; ++k)
    {
        Lower[k] = A.Lower[k];
    }
    for (int k=0; k<d_up; ++k)
    {
        Upper[k] = A.Upper[k];
    }
 }

 matrice_nonsym::matrice_nonsym(const matrice_sym& A):matrice_profil(A.n)//constructeur par copie à partir d'une matrice symétrique
 {
     int d_prof = A.Profil.dim();
     int d_pos = A.Posdiag.dim();
     int d_low = A.Lower.dim();
     if (d_prof!=0)
    {
        Profil = vecteur(d_prof);
        Posdiag = vecteur(d_pos);
        Lower = vecteur(d_low);
        Upper = vecteur(d_low);
    }
    for (int k=0; k<d_prof; ++k)
    {
        Profil[k] = A.Profil[k];
    }
    for (int k=0; k<d_pos; ++k)
    {
        Posdiag[k] = A.Posdiag[k];
    }
    for (int k=0; k<d_low; ++k)
    {
        Lower[k] = A.Lower[k];
    }
    for (int k=0; k<d_low; ++k)
    {
        Upper[k] = A.Lower[k]; //Upper est identique à Lower
    }
 }

 matrice_nonsym::~matrice_nonsym()
 {
     if (Profil!=0)
     {
         delete [] Profil;
     }
     if (Posdiag!=0)
     {
         delete [] Posdiag;
     }
     if (Lower!=0)
     {
         delete [] Lower;
     }
     if (Upper!=0)
     {
         delete [] Upper;
     }
     n = 0;
 }

matrice_nonsym& matrice_nonsym::operator=(const matrice_nonsym& A)
{
    int d_prof = A.Profil.dim();
    int d_pos = A.Posdiag.dim();
    int d_low = A.Lower.dim();
    int d_up = A.Upper.dim();
    delete [] Profil;
    delete [] Posdiag;
    delete [] Lower;
    delete [] Upper;
    Profil = vecteur(d_prof);
    Posdiag = vecteur(d_pos);
    Lower = vecteur(d_low);
    Upper = vecteur(d_up);
    if (n!=A.n)
    {
        n = A.n;
    }
    for (int k=0; k<d_prof; ++k)
    {
        Profil[k] = A.Profil[k];
    }
    for (int k=0; k<d_pos; ++k)
    {
        Posdiag[k] = A.Posdiag[k];
    }
    for (int k=0; k<d_low; ++k)
    {
        Lower[k] = A.Lower[k];
    }
    for (int k=0; k<d_up; ++k)
    {
        Upper[k] = A.Upper[k];
    }
    return *this;
}

double matrice_nonsym::val(int i, int j) const //Opérateur de lecture
{
    if (i==j)
    {
        return Lower[Posdiag[i]];
    }
    if (i>j) //triangle inférieur
    {
        if (j>=Profil[i]) //après le premier terme non-nul de la ligne
        {
            return Lower[Posdiag[i-1]+(j-Profil[i]+1)];
        }
        else //avant le premier terme non-nul de la ligne
        {
            return 0.;
        }
    }
    else //triangle supérieur
    {
        if (j>=Profil[i]) //après le premier terme non-nul de la colonne
        {
            return Upper[Posdiag[i-1]+(j-Profil[i]+1)];
        }
        else//avant le premier terme non-nul de la colonne
        {
            return 0.;
        }
    }
}


double& matrice_nonsym::val(int i, int j) //Lecture et écriture
{
    if (i==j)
    {
        return Lower[Posdiag[i]]; //n'incrémente que dans le lower
    }
    if (i>j) //triangle inférieur
    {
        if (j>=Profil[i]) //après le premier terme non-nul de la ligne
        {
            return Lower[Posdiag[i-1]+(j-Profil[i]+1)];
        }
        else //avant le premier terme non-nul de la ligne
        {
            //On doit l'ajouter au bon endroit
            for (int k=i; k<Profil[i];++k)//termes à rajouter
            {
                Lower.add(Posdiag[i-1]+k+1,0.); //rajoute des 0 entre le nouveau point à écrire et celui existant //Fonction à implémenter (indice, valeur)
                Posdiag[i]+=1;//décale le terme diagonal de la ligne à chaque ajout d'un terme
            }
            Profil[i]=j; //nouveau premier terme non-nul
            return Lower[Posdiag[i-1]+Profil[i]+1];
        }
    }
    else //triangle supérieur
    {
        if (j>=Profil[i]) //après le premier terme non-nul de la colonne
        {
            return Upper[Posdiag[i-1]+(j-Profil[i]+1)];
        }
        else//avant le premier terme non-nul de la colonne
        {
            //On doit l'ajouter au bon endroit
            for (int k=i; k<Profil[i];++k)//termes à rajouter
            {
                Upper.add(Posdiag[i-1]+k+1,0.); //rajoute des 0 entre le nouveau point à écrire et celui existant //Fonction à implémenter (indice, valeur)
                Posdiag[i]+=1;//décale le terme diagonal de la ligne à chaque ajout d'un terme
            }
            Profil[i]=j; //nouveau premier terme non-nul
            return Upper[Posdiag[i-1]+Profil[i]+1];
        }
    }
}

matrice_nonsym& matrice_nonsym::operator*(double a)
{
    Lower*=a;
    Upper*=a;
    return *this;
}

matrice_nonsym& matrice_nonsym::operator/(double a)
{
    if (a==0.)
    {
        cout<<"division par zéro"<<endl;
        exit(-1);
    }
    Lower/=a;
    Upper/=a;
    return *this;
}

matrice_nonsym& matrice_nonsym::operator+(const matrice_nonsym& A)
{
    if (n=!A.n)
    {
        cout<<"Les matrices ne sont pas de mêmes dimensions"<<endl;
        exit(-1);
    }
    for (int i=0; i<n; ++i)
    {
        (*this).val(i,i)+=A.val(i,i); //n'incrémente que dans le lower
        //Lower[Posdiag[i]]+=A.Lower[Posdiag[i]];
        Upper[Posdiag[i]]+=A.Upper[Posdiag[i]];

        for (int j=A.Profil[i]; j<i; ++j) //du premier non-nul jusqu'à la diagonale exclu
        {
            (*this).val(i,j)+=A.val(i,j); //modification de Lower
            (*this).val(j,i)+=A.val(j,i);//modification de Upper
        }

        /*
        if (A.profil[i]>=profil[i])//A n'enlève pas de 0 sur cette ligne
        {
            for (int j=A.Profil[i]; j<i; ++j) //du premier non-nul jusqu'à la diagonale exclu
            {
                (this)[i,j]+=A[i,j];
                (this)[j,i]+=A[j,i];
                //Lower[Posdiag[i-1]+j+1]+=A.Lower[A.Posdiag[i-1]+j+1];
                //Upper[Posdiag[i-1]+j+1]+=A.Upper[A.Posdiag[i-1]+j+1];
            }
        }
        else//A modifie la ligne
        {
            for (int j=A.profil[i]; j<i; ++j)
            {
                (this)[i,j]+=A[i,j];
                (this)[j,i]+=A[j,i];
                //
            }
        }
        */
    }
    return *this;
}

matrice_nonsym& matrice_nonsym::operator-(const matrice_nonsym& A)
{
    if (n=!A.n)
    {
        cout<<"Les matrices ne sont pas de mêmes dimensions"<<endl;
        exit(-1);
    }
    for (int i=0; i<n; ++i)
    {
        (*this).val(i,i)-=A.val(i,i); //n'incrémente que dans le lower
        //Lower[Posdiag[i]]+=A.Lower[Posdiag[i]];
        Upper[Posdiag[i]]-=A.Upper[Posdiag[i]];

        for (int j=A.Profil[i]; j<i; ++j) //du premier non-nul jusqu'à la diagonale exclu
        {
            (*this).val(i,j)-=A.val(i,j);
            (*this).val(j,i)-=A.val(j,i);
        }
    }
    return *this;
}

matrice_nonsym operator*(const matrice_nonsym& A, double a)
{
    return matrice_nonsym(A)*a;
}

matrice_nonsym operator*(double a, const matrice_nonsym& A)
{
    return matrice_nonsym(A)*a;
}

matrice_nonsym operator/(const matrice_nonsym& A, double a)
{
    return matrice_nonsym(A)/a;
}

matrice_nonsym operator+(const matrice_nonsym& A, const matrice_nonsym& B)
{
    return matrice_nonsym(A)+B;
}

matrice_nonsym operator-(const matrice_nonsym& A, const matrice_nonsym& B)
{
    return matrice_nonsym(A)-B;
}

/*
###########################"Opérations mixtes"############################################
##########################################################################################
*/

matrice_nonsym operator+(const matrice_nonsym& A, const matrice_sym& B)
{
    return matrice_nonsym(B)+A;
}

matrice_nonsym operator-(const matrice_nonsym& A, const matrice_sym& B)
{
    return matrice_nonsym(B)-A;
}

matrice_nonsym operator+(const matrice_sym& A, const matrice_nonsym& B)
{
    matrice_nonsym(A)+B;
}

matrice_nonsym operator-(const matrice_sym& A, const matrice_nonsym& B)
{
    matrice_nonsym(A)-B;
}
