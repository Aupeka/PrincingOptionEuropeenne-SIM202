#include "Projet.hpp"
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <cstdlib>
using namespace std;

#################Classe vecteur ##############
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

// opérations composante à composante u*~v et u/~v
operande operator~(const vecteur & u)
{
  return operande(&u);
}

vecteur operator*(const vecteur & u, const operande & ov)
{
  const vecteur& v=*ov.vp;
  int d=v.dim();
  test_dim(u.dim(),d,"vector * tensor");
  vecteur uv(d);
  for (int i=0;i<d;i++) uv[i]=u[i]*v[i];
  return uv;
}

vecteur operator/(const vecteur & u, const operande & ov)
{
  const vecteur& v=*ov.vp;
  int d=v.dim();
  test_dim(u.dim(),d,"vector / tensor");
  vecteur uv(d);
  for (int i=0;i<d;i++) uv[i]=u[i]/v[i];
  return uv;
}

#############################################################

######################Classe Maillage #######################

 #####Classe Point########
Point operator +(const Point & P, const Point & Q)
{
    P+=Q;
    return P;
}
Point operator -(const Point & P, const Point & Q)
{
    P-=Q;
    return P;
}
Point operator *(const Point & P,double a)
{
    P*=a;
    return P;
}
Point operator *(double a,const Point & P)
{
    P*=a;
    return P;
}
Point operator /(const Point & P, double a)
{
    P/=a;
    return P;
}
bool operator ==(const Point & P, const Point & Q)
{
    if ((P.x == Q.x) && (P.y==Q.y))
    {
        return true;
    }
    else
    {
        return false;
    }
}
bool operator !=(const Point & P, const Point & Q)
{
     if ((P.x != Q.x) || (P.y != Q.y))
    {
        return true;
    }
    else
    {
        return false;
    }
}
ostream & operator<<(ostream & out, const Point & P)
{
    out<<"("<<P.x<<" , "<<P.y<<")";
    return out;
}

#### Classe Numeros ####

ostream& operator<<(ostream& out, const Numeros & N)
{ int i=1;
  for (; i<int(N.size()); i++) out<<N(i)<<" ";
  out<<N(i);
  return out;
}

#####Classe Maillage ######

void Maillage::maille_carre_unite(int m,int n)
{
    double dx=1./m, dy=1./n;
    sommets.resize((m+1)*(n+1));
    vector<Point>::iterator its=sommets.begin();
    for (int j=0; j<n+1;j++)
    {
        double y=j*dy;
        for (int i=0; i<m+1;i++)
        {
            *its=Point(i+dx,y);
        }
    }
    for (int j=0;j<n;j++)
    {
        for (int i=0;i<m;i++)
        {
            int q=j*(m+1)+i;
            if (rand()%2)
            {
                numelts.push_back(Numeros(q,q+1,q+m+2));
                numelts.push_back(Numeros(q,q+m+2,q+m+1));
            }
            else
            {
                numelts.push_back(Numeros(q,q+1,q+m+1));
                numelts.push_back(Numeros(q+1,q+m+2,q+m+1));
            }
        }
    }
}

void Maillage::affiche() const
{
    cout<<'Liste des sommets ('<<sommets.size()<<' points)\n';
    vector<Point>::const_iterator its=sommets.begin();
    int i=0;
    while(its!=sommets.end())
}
