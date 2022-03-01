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
  if (d<0) stop("init() : dimension <= 0");
  dim_=d;
  val_ = new double[d];
}

void vecteur::clear()    // désallocation
{
  if (val_!=0) delete [] val_;
  dim_=0;
}

void vecteur::resize(int ni)
{
    if (ni<=0)
    {
        cout<<"Dimension demandée <=0"<<endl;
        exit(-1);
    }
    if (ni>dim_)
    {
        double * temp_ = new double[ni];
        for (int k=0; k<dim_; ++k)
        {
            temp_[k] = val_[k];
        }

        for (int k=dim_; k<ni; ++k)
        {
            temp_[k] = 0.;
        }
        dim_ = ni;
        delete [] val_;
        val_ = temp_;
    }
    else
    {
        dim_ = ni;
        cout<<"Pas besoin de redimensionnement car la taille actuelle >= à la taille demandée"<<endl;
        exit(-1);
    }
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
    if (n<0 || m<0)
    {
        cout<<"dimension negative"<<endl;
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


matrice::~matrice()
{
    if (cols_!=0)
    {
        delete [] cols_;
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

double matrice::operator()(int i, int j) const
{
    if (i>m || j>n || i<0 || j<0)
    {
        cout <<"Indice en dehors des bornes"<<endl;
        exit(-1);
    }
    return cols_[j][i]; //surcharge du vecteur
}

double& matrice::operator()(int i, int j)
{
    if (i>m || j>n || i<0 || j<0)
    {
        cout <<"Indice en dehors des bornes"<<endl;
        exit(-1);
    }
    return cols_[j][i]; //surcharge du vecteur
}

matrice& matrice::operator*(double a)
{
    for (int i=0; i<m; i++)
    {
        for (int j=0; j<n; ++j)
        {
            (*this)(i,j)*=a;
        }
    }
    return (*this);
}

matrice& matrice::operator/(double a)
{
    if (a==0)
    {
        cout<<"Division par zéro"<<endl;
        exit(-1);
    }
    for (int i=0; i<m; i++)
    {
        for (int j=0; j<n; ++j)
        {
            (*this)(i,j)/=a;
        }
    }
    return (*this);
}

matrice& matrice::operator+(const matrice& A)
{
    if (A.m!=m || A.n!=n)
    {
        cout<<"Les matrices ne sont pas de mêmes dimensions"<<endl;
        exit(-1);
    }
    for (int i=0; i<m; i++)
    {
        for (int j=0; j<n; ++j)
        {
            (*this)(i,j)+=A(i,j);
        }
    }
    return (*this);
}

matrice& matrice::operator-(const matrice& A)
{
    if (A.m!=m || A.n!=n)
    {
        cout<<"Les matrices ne sont pas de mêmes dimensions"<<endl;
        exit(-1);
    }
    for (int i=0; i<m; i++)
    {
        for (int j=0; j<n; ++j)
        {
            (*this)(i,j)-=A(i,j);
        }
    }
    return (*this);
}

matrice operator*(const matrice& A, double a)
{
    return matrice(A)*a;
}

matrice operator/(const matrice& A, double a)
{
    return matrice(A)/a;
}

matrice operator*(double a, const matrice& A)
{
    return matrice(A)*a;
}

matrice operator+(const matrice& A, const matrice& B)
{
    return matrice(A)+B;
}

matrice operator-(const matrice& A, const matrice& B)
{
    return matrice(A)-B;
}

matrice operator*(const matrice& A, const matrice& B)
{
    if (A.n!=B.m)
    {
        cout<<"Problèmes de dimensions, produit impossible"<<endl;
        exit(-1);
    }
    int d = A.n;
    matrice Res(A.m,B.n);

    for (int i=0; i<A.m; ++i)
    {
        for (int j=0; j<B.n; ++j)
        {
            double res = 0.;
            for (int k=0; k<d; ++k)
            {
                res += A(i,k)*B(k,j);
            }
            Res(i,j) = res;
        }
    }
    return Res;
}


vecteur operator*(const matrice& A, const vecteur& u)
{
    if (A.n!=u.dim())
    {
        cout<<"produit matrice vecteur : dimensions incompatibles"<<endl;
        exit(-1);
    }
    vecteur Res(A.m,0.);
    for (int i=1; i<A.m; ++i)
    {
        for (int j=1; j<A.n;++j)
        {
           Res[i] = A(i,j)*u[j];
        }
    }
    return Res;
}


matrice transpose(const matrice& A)
{
    matrice Res(A.m,A.n);
    for (int i=0; i<A.m; ++i)
    {
        for (int j=0; j<A.n; ++j)
        {
            Res(i,j) = A(j,i);
        }
    }
    return Res;
}

ostream & operator<<(ostream & os, const matrice& A)
{
    if (A.cols_!=0)
    {
        for (int i=0; i<A.m; ++i)
        {
            os<<"[";
            for (int j=0; j<A.n-1; ++j)
            {
                os<<A(i,j)<<",";
            }
            os<<A(i,A.n-1)<<"]"<<endl;
        }
    }
    return os;
}

/*
###################################################"Matrice profil"####################################################
#######################################################################################################################
*/

matrice_profil::matrice_profil(int ni, const vecteur Pi) : n(ni), Profil(Pi), Posdiag(0)
{
    if (n<0)
    {
        cout<<"dimension négative"<<endl;
        exit(-1);
    }
    Posdiag.resize(n);
    for (int k=1; k<n; ++k) //1er terme est nécessairement sur la diagonale
    {
        Posdiag[k] = Posdiag[k-1] + (k-Profil[k]+1); //nombre de terme par ligne ajouté
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

/*
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
 */

/*
 #####################################"Matrice sym"########################################################
 ##########################################################################################################
*/

matrice_sym::matrice_sym(int ni, const vecteur Pi) : matrice_profil(ni, Pi), Lower(0)
{
    //Lower = vecteur(int(Posdiag[ni-1])+1,0.); //initialisation de longueur connue : dernier terme de Posdiag + 1
    Lower.resize(int(Posdiag[ni-1])+1);
}

matrice_sym::matrice_sym(const matrice_sym& A):matrice_profil(A.n,A.Profil)
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

 /*
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
*/

matrice_sym& matrice_sym::operator=(const matrice_sym& A)
{
    int d_prof = A.Profil.dim();
    int d_pos = A.Posdiag.dim();
    int d_low = A.Lower.dim();
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

double matrice_sym::operator()(int i, int j) const //Opérateur de lecture
{
    if (i==j)
    {
        return Lower[Posdiag[i]];
    }
    if (i<j) //triangle supérieur
    {
        int temp = j;
        j = i;
        i = temp;
    }
    //code pour le triangle inférieur : i>j mais inversion indice permet de tout traiter
    if (j>=Profil[i]) //après le premier terme non-nul de la ligne
    {
        return Lower[Posdiag[i-1]+(j-Profil[i]+1)];
    }
    else //avant le premier terme non-nul de la ligne
    {
        return 0.;
    }
}


double& matrice_sym::operator()(int i, int j) //Lecture et écriture
{
    if (i==j)
    {
        return Lower[Posdiag[i]]; //n'incrémente que dans le lower
    }
    if (i<j) //triangle supérieur
    {
        int temp = j;
        j=i;
        i=temp;
    }
    //code pour le triangle inférieur : i>j mais inversion indice permet de tout traiter
    if (j>=Profil[i]) //après le premier terme non-nul de la ligne
    {
        return Lower[Posdiag[i-1]+(j-Profil[i]+1)];
    }
    else //avant le premier terme non-nul de la ligne
    {
        cout<<"("<<i<<","<<j<<")"<<": Coordonnees hors profil"<<endl;
        exit(-1);

        /*
        //On doit l'ajouter au bon endroit
        for (int k=j; k<Profil[i];++k)//termes à rajouter
        {
            Lower.add(Posdiag[i-1]+k+1,0.); //rajoute des 0 entre le nouveau point à écrire et celui existant
            Posdiag[i]+=1;//décale le terme diagonal de la ligne à chaque ajout d'un terme
        }
        Profil[i]=j; //nouveau premier terme non-nul
        return Lower[Posdiag[i-1]+Profil[i]+1];
        */
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
    if (n!=A.n)
    {
        cout<<"Les matrices ne sont pas de mêmes dimensions"<<endl;
        exit(-1);
    }
    for (int i=0; i<n; ++i)
    {
        (*this)(i,i)+=A(i,i);

        for (int j=A.Profil[i]; j<Posdiag[i]; ++j) //du premier non-nul jusqu'à la diagonale exclu
        {
            (*this)(i,j)+=A(i,j);
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
    if (n!=A.n)
    {
        cout<<"Les matrices ne sont pas de mêmes dimensions"<<endl;
        exit(-1);
    }
    for (int i=0; i<n; ++i)
    {
        (*this)(i,i)-=A(i,i);

        for (int j=A.Profil[i]; j<Posdiag[i]; ++j) //du premier non-nul jusqu'à la diagonale exclu
        {
            (*this)(i,j)-=A(i,j);
        }
    }
    return *this;
}

ostream & operator<<(ostream & os, const matrice_sym& A)
{
  for (int i=0; i<A.n; ++i)
    {
        os<<"[";
        for (int j=0; j<A.n-1; ++j)
        {
            os<<A(i,j)<<",";
        }
        os<<A(i,A.n-1)<<"]"<<endl;
    }
  return os;
}

void print(const matrice_sym& A)
{
    for (int i=0; i<A.n; ++i)
    {
        cout<<"[";
        for (int j=0; j<A.n-1; ++j)
        {
            cout<<A(i,j)<<",";
        }
        cout<<A(i,A.n-1)<<"]"<<endl;
    }
}

matrice_sym operator*(const matrice_sym& A, double a)
{
    return matrice_sym(A)*a;
}

matrice_sym operator*(double a, const matrice_sym& A)
{
    return matrice_sym(A)*a;
}

vecteur operator*(const matrice_sym& A, const vecteur& V)
{
    if (A.n!=V.dim())
    {
        cout<<"La matrice n'a pas la même dimension que le vecteur"<<endl;
        exit(-1);
    }
    vecteur Res(A.n);
    for (int k=0;k<A.n;++k)
    {
        for (int i=0; i<A.n; ++i)
        {
            Res[k]+=A(k,i)*V[i];
        }
    }
    return Res;
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

matrice_sym transpose(const matrice_sym& A)
{
    return matrice_sym(A);
}

/*
#####################################"Matrice non sym"#####################################################
###########################################################################################################
*/

 matrice_nonsym::matrice_nonsym(int ni, const vecteur Pi) : matrice_profil(ni, Pi), Lower(0), Upper(0)
 {
    //Lower = vecteur(int(Posdiag[ni-1])+1,0.); //initialisation de longueur connue : dernier terme de Posdiag + 1
    //Upper = vecteur(int(Posdiag[ni-1])+1,0.);
    Lower.resize(int(Posdiag[ni-1])+1);
    Upper.resize(int(Posdiag[ni-1])+1);
 }

 matrice_nonsym::matrice_nonsym(const matrice_nonsym& A):matrice_profil(A.n,A.Profil)
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

 matrice_nonsym::matrice_nonsym(const matrice_sym& A):matrice_profil(A.n, A.Profil)//constructeur par copie à partir d'une matrice symétrique
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

 /*
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
 */

matrice_nonsym& matrice_nonsym::operator=(const matrice_nonsym& A)
{
    int d_prof = A.Profil.dim();
    int d_pos = A.Posdiag.dim();
    int d_low = A.Lower.dim();
    int d_up = A.Upper.dim();
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

double matrice_nonsym::operator ()(int i, int j) const //Opérateur de lecture
{
    if (i==j)
    {
        return Upper[Posdiag[i]];
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
        if (i>=Profil[j]) //après le premier terme non-nul de la colonne
        {
            return Upper[Posdiag[j-1]+(i-Profil[j]+1)];
        }
        else//avant le premier terme non-nul de la colonne
        {
            return 0.;
        }
    }
}


double& matrice_nonsym::operator()(int i, int j) //Lecture et écriture
{
    if (i==j)
    {
        return Upper[Posdiag[i]]; //n'incrémente que dans le Upper (mieux que dans Lower à cause de la factorisation LU)
    }
    if (i>j) //triangle inférieur
    {
        if (j>=Profil[i]) //après le premier terme non-nul de la ligne
        {
            return Lower[Posdiag[i-1]+(j-Profil[i]+1)];
        }
        else //avant le premier terme non-nul de la ligne
        {
            cout<<"("<<i<<","<<j<<")"<<": Coordonnees hors profil"<<endl;
            exit(-1);

            /*
            //On doit l'ajouter au bon endroit
            for (int k=j; k<Profil[i];++k)//termes à rajouter
            {
                Lower.add(Posdiag[i-1]+k+1,0.); //rajoute des 0 entre le nouveau point à écrire et celui existant
                Posdiag[i]+=1;//décale le terme diagonal de la ligne à chaque ajout d'un terme
            }
            Profil[i]=j; //nouveau premier terme non-nul
            return Lower[Posdiag[i-1]+Profil[i]+1];
            */
        }
    }
    else //triangle supérieur
    {
        if (i>=Profil[j]) //après le premier terme non-nul de la colonne
        {
            return Upper[Posdiag[j-1]+(i-Profil[j]+1)];
        }
        else//avant le premier terme non-nul de la colonne
        {
            cout<<"("<<i<<","<<j<<")"<<": Coordonnees hors profil"<<endl;
            exit(-1);

            /*
            //On doit l'ajouter au bon endroit
            for (int k=j; k<Profil[i];++k)//termes à rajouter
            {
                Upper.add(Posdiag[i-1]+k+1,0.); //rajoute des 0 entre le nouveau point à écrire et celui existant
                Posdiag[i]+=1;//décale le terme diagonal de la ligne à chaque ajout d'un terme
            }
            Profil[i]=j; //nouveau premier terme non-nul
            return Upper[Posdiag[i-1]+Profil[i]+1];
            */
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
    if (n!=A.n)
    {
        cout<<"Les matrices ne sont pas de mêmes dimensions"<<endl;
        exit(-1);
    }
    for (int i=0; i<n; ++i)
    {
        (*this)(i,i)+=A(i,i); //n'incrémente que dans le lower
        //Lower[Posdiag[i]]+=A.Lower[Posdiag[i]];
        Upper[Posdiag[i]]+=A.Upper[Posdiag[i]];

        for (int j=A.Profil[i]; j<Posdiag[i]; ++j) //du premier non-nul jusqu'à la diagonale exclu
        {
            (*this)(i,j)+=A(i,j); //modification de Lower
            (*this)(j,i)+=A(j,i);//modification de Upper
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
    if (n!=A.n)
    {
        cout<<"Les matrices ne sont pas de mêmes dimensions"<<endl;
        exit(-1);
    }
    for (int i=0; i<n; ++i)
    {
        (*this)(i,i)-=A(i,i); //n'incrémente que dans le lower
        //Lower[Posdiag[i]]+=A.Lower[Posdiag[i]];
        Upper[Posdiag[i]]-=A.Upper[Posdiag[i]];

        for (int j=A.Profil[i]; j<Posdiag[i]; ++j) //du premier non-nul jusqu'à la diagonale exclu
        {
            (*this)(i,j)-=A(i,j);
            (*this)(j,i)-=A(j,i);
        }
    }
    return *this;
}

ostream & operator<<(ostream & os, const matrice_nonsym& A)
{
  for (int i=0; i<A.n; ++i)
    {
        os<<"[";
        for (int j=0; j<A.n-1; ++j)
        {
            os<<A(i,j)<<",";
        }
        os<<A(i,A.n-1)<<"]"<<endl;
    }
  return os;
}

void print(const matrice_nonsym& A)
{
    for (int i=0; i<A.n; ++i)
    {
        cout<<"[";
        for (int j=0; j<A.n-1; ++j)
        {
            cout<<A(i,j)<<",";
        }
        cout<<A(i,A.n-1)<<"]"<<endl;
    }
}

matrice_nonsym operator*(const matrice_nonsym& A, double a)
{
    return matrice_nonsym(A)*a;
}

matrice_nonsym operator*(double a, const matrice_nonsym& A)
{
    return matrice_nonsym(A)*a;
}

vecteur operator*(const matrice_nonsym& A, const vecteur& V)
{
    if (A.n!=V.dim())
    {
        cout<<"La matrice n'a pas la même dimension que le vecteur"<<endl;
        exit(-1);
    }
    vecteur Res(A.n);
    for (int k=0;k<A.n;++k)
    {
        for (int i=0; i<A.n; ++i)
        {
            Res[k]+=A(k,i)*V[i];
        }
    }
    return Res;
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

matrice_nonsym transpose(const matrice_nonsym& A)
{
    matrice_nonsym Res(A.n,A.Profil);
    Res.Upper = A.Lower;
    Res.Lower = A.Upper;
    return Res;
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
    return matrice_nonsym(A)+B;
}

matrice_nonsym operator-(const matrice_sym& A, const matrice_nonsym& B)
{
    return matrice_nonsym(A)-B;
}

/*
matrice_nonsym LUdecomposition(const matrice_sym& A)
{
   int n = A.n;
   matrice_nonsym LU(n,A.Profil);
   //int i = 1, j = 1, k = 1;
   for (int i = 1; i < n+1; i++)
   {
       for (int j = 1; j < n+1; j++)
       {
           if (i>=j)
           {
               LU(j,i) = A(i,j);
                for (int k = 1; k < i+1; k++)
                {
                    LU(j,i) = LU(j,i) - LU(j,k) * LU(k,i);
                }
           }
       }
       for (int j = 1; j < n+1; j++)
       {
           if (j==i)
           {
               LU(i,j)=1;
           }
           else if (j>i)
           {
               LU(i,j) = A(i,j) / LU(i,i);
               for (int k = 1; k < i+1; k++)
               {
                   LU(i,j) = LU(i,j) - ((LU(i,k) * LU(k,j)) / LU(i,i));
               }
           }
       }
   }
}
*/
/*
################################ancienne version############################################

void LUdecomposition(const matrice& A)
{
   int n = A.n;
   matrice l(n,n);
   matrice u(n,n);
   //int i = 1, j = 1, k = 1;
   for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
         if (i < j)
         l(i,j) = 0;
         else if (i==j)
         {
             l(i,i) = 1;
         }
         else {
            l(i,j) = A(i,j);
            for (int k = 0; k < i; k++) {
               l(i,j) -= l(i,k) * u(k,j);
            }
         }
      }
      for (int j = 0; j < n; j++) {
         if (i > j)
         u(i,j) = 0;
         //else if (j == i)
         //u(i,j) = 1;
         else {
            u(i,j) = A(i,j) / l(i,i);
            for (int k = 0; k < i; k++) {
               u(i,j) -= ((l(i,k) * u(k,j)) / l(i,i));
            }
         }
      }
   }
   cout<<u<<endl;
   cout<<l<<endl;
}

void LUdecomposition(const matrice_nonsym& A,matrice& l,matrice& u, int n)
{
   int i = 1, j = 1, k = 1;
   for (i = 1; i < n+1; i++) {
      for (j = 1; j < n+1; j++) {
         if (j < i)
         l(j,i) = 0;
         else {
            l(j,i) = A(i,j);
            for (k = 1; k < i+1; k++) {
               l(j,i) = l(j,i) - l(j,k) * u(k,i);
            }
         }
      }
      for (j = 1; j < n+1; j++) {
         if (j < i)
         u(i,j) = 0;
         else if (j == i)
         u(i,j) = 1;
         else {
            u(i,j) = A(i,j) / l(i,i);
            for (k = 1; k < i+1; k++) {
               u(i,j) = u(i,j) - ((l(i,k) * u(k,j)) / l(i,i));
            }
         }
      }
   }
}

*/

/*
matrice_nonsym LUdecomposition(const matrice_nonsym& A)
{
    int n=A.n;
    matrice_nonsym lu(n,A.Profil); //Utilise la conservation du profil p
    for (int i=0; i<n; i++)
    {
        lu.Lower[lu.Posdiag[i]] = 1;//mets diagonale lower à 0 car accesseur classique rempli Upper et non lower
    }
    for (int p=0; p<n; p++)
    {
        for (int j=p; j<n; j++) //boucle sur la matrice triangulaire supérieure
        {
            if (p>=lu.Profil[j]) //ne modifie que les termes dans le profil
            {
                lu(p,j) = A(p,j);
                for (int k=lu.Profil[p]; k<p-1; k++)
                {
                    lu(p,j) -= lu(p,k)*lu(k,j);
                }
            }
        }
        if (lu(p,p)==0)
        {
            cout<<lu<<endl;
            cout<<"Erreur : Matrice non-factorisable"<<endl;
            exit(-1);
        }
        for (int i=p+1; i<n; i++) //boucle sur la matrice triangulaire inférieure
        {
            if (p>=lu.Profil[i])
            {
                lu(i,p) = A(i,p);
                for (int k=0; k<p-1; k++)
                {
                    lu(i,p) -= lu(i,k)*lu(k,p);
                }
                lu(i,p) /= lu(p,p);
            }
        }
    }
    return lu;
}
*/

matrice_nonsym LUdecomposition(const matrice_nonsym& A)
{
    int n=A.n;
    matrice_nonsym lu(n,A.Profil); //Utilise la conservation du profil p
    for (int i=1; i<=n; i++)
    {
        lu.Lower[lu.Posdiag[i-1]] = 1;//mets diagonale lower à 0 car accesseur classique rempli Upper et non lower
    }
    for (int p=1; p<=n; p++)
    {
        for (int j=p-1; j<n; j++) //boucle sur la matrice triangulaire supérieure
        {
            if (p-1>=lu.Profil[j]) //ne modifie que les termes dans le profil
            {
                lu(p-1,j) = A(p-1,j);
                for (int k=lu.Profil[p-1]; k<=(p-1)-1; k++)
                {
                    lu(p-1,j) -= lu(p-1,k)*lu(k,j);
                }
            }
        }
        if (lu(p-1,p-1)==0)
        {
            cout<<lu<<endl;
            cout<<"Erreur : Matrice non-factorisable"<<endl;
            exit(-1);
        }
        for (int i=(p-1)+1; i<n; i++) //boucle sur la matrice triangulaire inférieure
        {
            if (p-1>=lu.Profil[i])
            {
                lu(i,p-1) = A(i,p-1);
                for (int k=lu.Profil[i]; k<=(p-1)-1; k++)
                {
                    lu(i,p-1) -= lu(i,k)*lu(k,p-1);
                }
                lu(i,p-1) /= lu(p-1,p-1);
            }
        }
    }
    return lu;
}

matrice_nonsym LUdecomposition(const matrice_sym& A)
{
    return LUdecomposition(matrice_nonsym(A));
}
