#include "Projet.hpp"
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <cstdlib>
using namespace std;
/*
#################Classe vecteur ##############
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
/*
#############################################################

######################Classe Maillage #######################

 #####Classe Point########
*/

 Point& Point::tf_affine(const vector<double> & A, const vector<double> &t)
 {
     double xx=x, yy=y;
     x=A[0]*xx+A[1]*yy+t[0];
     y=A[2]*xx+A[3]*yy+t[1];
     return *this;
 }



Point operator +(const Point & P, const Point & Q)
{
    return Point(P)+=Q;
}
Point operator -(const Point & P, const Point & Q)
{
    return Point(P)-=Q;
}
Point operator *(const Point & P,double a)
{
    return Point(P)*=a;
}
Point operator *(double a,const Point & P)
{
    return Point(P)*=a;
}
Point operator /(const Point & P, double a)
{
    return Point(P)/=a;
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

bool operator <(const Point & P, const Point & Q)
{
    return (P.x<Q.x)||((P.x==Q.x)&&(P.y<Q.y));
}

ostream & operator<<(ostream & out, const Point & P)
{
    out<<"("<<P.x<<" , "<<P.y<<")";
    return out;
}
/*
#### Classe Numeros ####
*/

ostream& operator<<(ostream& out, const Numeros & N)
{ int i=1;
  for (; i<int(N.size()); i++) out<<N(i)<<" ";
  out<<N(i);
  return out;
}
/*
#####Classe Maillage ######
*/

void Maillage::maille_carre_unite(int m,int n)
{
    double dx=1./m, dy=1./n;
    sommets.resize((m+1)*(n+1));
    vector<Point>::iterator its=sommets.begin();
    for (int j=0; j<n+1;j++)
    {
        double y=j*dy;
        for (int i=0; i<m+1;i++,its++)
        {
            *its=Point(i*dx,y);
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
    cout<<"Liste des sommets ("<<sommets.size()<<" points)\n";
    vector<Point>::const_iterator its=sommets.begin();
    int i=1;
    while(its!=sommets.end())
    {
         cout<<"sommet "<<i<<" : "<<*its<<"\n";
         its++;
         i++;
    }
    cout<<"Liste des triangles ("<<numelts.size()<<" triangles)\n";
    list<Numeros>::const_iterator itn=numelts.begin();
    i=1;
    while(itn!=numelts.end())
    {
        cout<<"triangle "<<i<<" : "<<*itn<<"\n";
        itn++;
        i++;
    }
}

Maillage& Maillage::affine(const vector<double> &A, const vector<double> &t)
{
    vector<Point>::iterator its=sommets.begin();
    for(; its!=sommets.end();its++) its->tf_affine(A,t);
    return *this;
}

void Maillage::maille_rectangle(double a, double b, double c, double d, int m, int n)
{
    maille_carre_unite(m,n);
    vector<double> A(4,0.); A[0]=b-a; A[3]=d-c;
    vector<double> t(2 ,0.); t[0]=a; t[1]=c;
    affine(A,t);
}

Maillage & Maillage::operator+= (const Maillage &M)
{
    vector<Point>::const_iterator itp;
    map<Point,int> ptrang;
    map<Point,int>::iterator itm;
    vector<int> num2;
    int k=0;
    for(itp=sommets.begin();itp!=sommets.end();++itp,k++)
    {
        ptrang.insert(pair<Point,int>(*itp,k));
    }
    int l=0;
    num2.resize(M.sommets.size());
    for(itp=M.sommets.begin();itp!=M.sommets.end();++itp,l++)
    {
        map<Point,int>::iterator itm=ptrang.find(*itp);
        if(itm!=ptrang.end()) {num2[l]=itm->second;}
        else
        {
            ptrang.insert(pair<Point,int>(*itp,k));
            num2[l]=k++;
        }
    }
    sommets.resize(ptrang.size());
    for(itm=ptrang.begin();itm!=ptrang.end();++itm)
    {
        sommets[itm->second]=itm->first;
    }
    list<Numeros>::const_iterator itn=M.numelts.begin();
    for(;itn!=M.numelts.end();itn++)
    {
        Numeros nums=*itn;
        for(int i=0;i<int(nums.size());i++) nums[i]=num2[nums[i]];
        numelts.push_back(nums);
    }
    return(*this);
}

Maillage operator+(const Maillage &M1, const Maillage &M2)
{
    Maillage M(M1);
    return M+=M2;
}

void Maillage::savecoord(const char *fn) const
{
    ofstream out(fn);
    out<<sommets.size()<<endl;
    vector<Point>::const_iterator itn=sommets.begin();
    int i=1;
    for(;itn!=sommets.end();itn++,i++)
    {
      out<<i<<" "<<itn->x<<" "<<itn->y<<endl;
    }
    out.close();
}

void Maillage::savenumtri(const char *fn) const
{
    ofstream out(fn);
    out<<numelts.size()<<endl;
    list<Numeros>::const_iterator itn=numelts.begin();
    int i=1;
    for(;itn!=numelts.end();itn++,i++)
    {
      out<<i<<" "<<(*itn)<<endl;
    }
    out.close();
}

/*
##############Classe matrice #################
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

matrice::~ matrice(){if(cols_!=0) delete [] cols_ ;}
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

double matrice::val(int i, int j) const {return (cols_[j-1].val_)[i-1];}
double& matrice::val(int i, int j) {return (cols_[j-1].val_)[i-1];}

vecteur produit(const matrice& A, const vecteur& u)
{
    if (A.n!=u.dim_)
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

void LUdecomposition(const matrice& A, int n)
{
   matrice l(n,n,0);
   matrice u(n,n,0);
   int i = 0, j = 0, k = 0;
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         if (j < i)
         l.val(j,i) = 0;
         else {
            l.val(j,i) = A.val(i,j);
            for (k = 0; k < i; k++) {
               l.val(j,i) = l.val(j,i) - l.val(j,k) * u.val(k,i);
            }
         }
      }
      for (j = 0; j < n; j++) {
         if (j < i)
         u.val(i,j) = 0;
         else if (j == i)
         u.val(i,j) = 1;
         else {
            u.val(i,j) = A.val(i,j) / l.val(i,i);
            for (k = 0; k < i; k++) {
               u.val(i,j) = u.val(i,j) - ((l.val(i,k) * u.val(k,j)) / l.val(i,i));
            }
         }
      }
   }
}


