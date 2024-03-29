#include "Projet.hpp"
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <cstdlib>
using namespace std;

/*
#############################################################
###################### Classe Vecteur #######################
#############################################################
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

void vecteur::clear()    // d�sallocation
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

// op�rateur d'acc�s � une plage i � j (1<i<=j<=dim)
vecteur vecteur::operator ()(int i, int j) const
{
  if (i<1 || i>j || j>dim_) stop("plage inconsistante");
  vecteur u(j-i+1);
  for (int k=0;k<=j-i;k++) u[k]=val_[i+k-1];
  return u;
}

// op�rateurs alg�briques
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

void vecteur::resize(int ni)
{
    if (ni<=0)
    {
        cout<<"Dimension demand�e <=0"<<endl;
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
        cout<<"Pas besoin de redimensionnement car la taille actuelle >= � la taille demand�e"<<endl;
        exit(-1);
    }
}


// op�rateurs externes
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

// op�rateurs de comparaison
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

// op�rateur de concat�nation
vecteur operator,(const vecteur &u, const vecteur& v)
{
  vecteur uv(u.dim()+v.dim());
  int k=0;
  for (int i=0;i<u.dim();i++, k++)uv[k]=u[i];
  for (int i=0;i<v.dim();i++, k++)uv[k]=v[i];
  return uv;
}

// op�rateurs de lecture et d'�criture
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

// op�rations composante � composante u*~v et u/~v
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
###################### Classe Point #########################
#############################################################
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
#############################################################
##################### Classe Num�ros ########################
#############################################################
*/

ostream& operator<<(ostream& out, const Numeros & N)
{ int i=1;
  for (; i<int(N.size()); i++) out<<N(i)<<" ";
  out<<N(i);
  return out;
}
/*
#############################################################
##################### Classe Maillage #######################
#############################################################
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
    list<Numeros>::const_iterator itn=numelts.begin();
    int i=1;
    for(;itn!=numelts.end();itn++,i++)
    {
      out<<i<<" "<<(*itn)<<endl;
    }
    out.close();
}

/*
#############################################################
###################### Classe Matrice #######################
#############################################################
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

matrice& matrice::operator*=(double a)
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

matrice& matrice::operator/=(double a)
{
    if (a==0)
    {
        cout<<"Division par z�ro"<<endl;
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

matrice& matrice::operator+=(const matrice& A)
{
    if (A.m!=m || A.n!=n)
    {
        cout<<"Les matrices ne sont pas de m�mes dimensions"<<endl;
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

matrice& matrice::operator-=(const matrice& A)
{
    if (A.m!=m || A.n!=n)
    {
        cout<<"Les matrices ne sont pas de m�mes dimensions"<<endl;
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
    return matrice(A)*=a;
}

matrice operator/(const matrice& A, double a)
{
    return matrice(A)/=a;
}

matrice operator*(double a, const matrice& A)
{
    return matrice(A)*=a;
}

matrice operator+(const matrice& A, const matrice& B)
{
    return matrice(A)+=B;
}

matrice operator-(const matrice& A, const matrice& B)
{
    return matrice(A)-=B;
}

matrice operator*(const matrice& A, const matrice& B)
{
    if (A.n!=B.m)
    {
        cout<<"Probl�mes de dimensions, produit impossible"<<endl;
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
    for (int i=0; i<A.m; ++i)
    {
        for (int j=0; j<A.n;++j)
        {
           Res[i] += A(i,j)*u[j];
        }
    }
    return Res;
}


matrice transpose(const matrice& A)
{
    matrice Res(A.n,A.m);
    for (int i=0; i<A.n; ++i)
    {
        for (int j=0; j<A.m; ++j)
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
#############################################################
###################### Matrice Profil #######################
#############################################################
*/

matrice_profil::matrice_profil(int ni, const vecteur Pi) : n(ni), Profil(Pi), Posdiag(0)
{
    if (n<0)
    {
        cout<<"dimension n�gative"<<endl;
        exit(-1);
    }
    Posdiag.resize(n);
    for (int k=1; k<n; ++k) //1er terme est n�cessairement sur la diagonale
    {
        Posdiag[k] = Posdiag[k-1] + (k-Profil[k]+1); //nombre de terme par ligne ajout�
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
#############################################################
################### Matrice Sym�trique ######################
#############################################################
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

double matrice_sym::operator()(int i, int j) const //Op�rateur de lecture
{
    if (i==j)
    {
        return Lower[Posdiag[i]];
    }
    if (i<j) //triangle sup�rieur
    {
        int temp = j;
        j = i;
        i = temp;
    }
    //code pour le triangle inf�rieur : i>j mais inversion indice permet de tout traiter
    if (j>=Profil[i]) //apr�s le premier terme non-nul de la ligne
    {
        return Lower[Posdiag[i-1]+(j-Profil[i]+1)];
    }
    else //avant le premier terme non-nul de la ligne
    {
        return 0.;
    }
}


double& matrice_sym::operator()(int i, int j) //Lecture et �criture
{
    if (i==j)
    {
        return Lower[Posdiag[i]]; //n'incr�mente que dans le lower
    }
    if (i<j) //triangle sup�rieur
    {
        int temp = j;
        j=i;
        i=temp;
    }
    //code pour le triangle inf�rieur : i>j mais inversion indice permet de tout traiter
    if (j>=Profil[i]) //apr�s le premier terme non-nul de la ligne
    {
        return Lower[Posdiag[i-1]+(j-Profil[i]+1)];
    }
    else //avant le premier terme non-nul de la ligne
    {
        cout<<"Coordonnees hors profil"<<endl;
        exit(-1);

        /*
        //On doit l'ajouter au bon endroit
        for (int k=j; k<Profil[i];++k)//termes � rajouter
        {
            Lower.add(Posdiag[i-1]+k+1,0.); //rajoute des 0 entre le nouveau point � �crire et celui existant
            Posdiag[i]+=1;//d�cale le terme diagonal de la ligne � chaque ajout d'un terme
        }
        Profil[i]=j; //nouveau premier terme non-nul
        return Lower[Posdiag[i-1]+Profil[i]+1];
        */
    }
}

matrice_sym& matrice_sym::operator*=(double a)
{
    Lower*=a;
    return *this;
}

matrice_sym& matrice_sym::operator/=(double a)
{
    if (a==0.)
    {
        cout<<"division par z�ro"<<endl;
        exit(-1);
    }
    Lower/=a;
    return *this;
}


matrice_sym& matrice_sym::operator+=(const matrice_sym& A)
{
    if (n!=A.n)
    {
        cout<<"Les matrices ne sont pas de m�mes dimensions"<<endl;
        exit(-1);
    }
    for (int i=0; i<n; ++i)
    {
        (*this)(i,i)+=A(i,i);

        for (int j=A.Profil[i]; j<Posdiag[i]; ++j) //du premier non-nul jusqu'� la diagonale exclu
        {
            (*this)(i,j)+=A(i,j);
        }

        /*
        if (A.profil[i]>=profil[i])//A n'enl�ve pas de 0 sur cette ligne
        {
            for (int j=A.Profil[i]; j<i; ++j) //du premier non-nul jusqu'� la diagonale exclu
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

matrice_sym& matrice_sym::operator-=(const matrice_sym& A)
{
    if (n!=A.n)
    {
        cout<<"Les matrices ne sont pas de m�mes dimensions"<<endl;
        exit(-1);
    }
    for (int i=0; i<n; ++i)
    {
        (*this)(i,i)-=A(i,i);

        for (int j=A.Profil[i]; j<Posdiag[i]; ++j) //du premier non-nul jusqu'� la diagonale exclu
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
    return matrice_sym(A)*=a;
}

matrice_sym operator*(double a, const matrice_sym& A)
{
    return matrice_sym(A)*=a;
}

vecteur operator*(const matrice_sym& A, const vecteur& V)
{
    if (A.n!=V.dim())
    {
        cout<<"La matrice n'a pas la m�me dimension que le vecteur"<<endl;
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
    return matrice_sym(A)/=a;
}

matrice_sym operator+(const matrice_sym& A, const matrice_sym& B)
{
    return matrice_sym(A)+=B;
}

matrice_sym operator-(const matrice_sym& A, const matrice_sym& B)
{
    return matrice_sym(A)-=B;
}

matrice_sym transpose(const matrice_sym& A)
{
    return matrice_sym(A);
}

/*
#############################################################
################ Matrice non-sym�trique #####################
#############################################################
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

 matrice_nonsym::matrice_nonsym(const matrice_sym& A):matrice_profil(A.n, A.Profil)//constructeur par copie � partir d'une matrice sym�trique
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
        Upper[k] = A.Lower[k]; //Upper est identique � Lower
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

double matrice_nonsym::operator ()(int i, int j) const //Op�rateur de lecture
{
    if (i==j)
    {
        return Lower[Posdiag[i]];
    }
    if (i>j) //triangle inf�rieur
    {
        if (j>=Profil[i]) //apr�s le premier terme non-nul de la ligne
        {
            return Lower[Posdiag[i-1]+(j-Profil[i]+1)];
        }
        else //avant le premier terme non-nul de la ligne
        {
            return 0.;
        }
    }
    else //triangle sup�rieur
    {
        if (i>=Profil[j]) //apr�s le premier terme non-nul de la colonne
        {
            return Upper[Posdiag[j-1]+(i-Profil[j]+1)];
        }
        else//avant le premier terme non-nul de la colonne
        {
            return 0.;
        }
    }
}


double& matrice_nonsym::operator()(int i, int j) //Lecture et �criture
{
    if (i==j)
    {
        return Lower[Posdiag[i]]; //n'incr�mente que dans le lower
    }
    if (i>j) //triangle inf�rieur
    {
        if (j>=Profil[i]) //apr�s le premier terme non-nul de la ligne
        {
            return Lower[Posdiag[i-1]+(j-Profil[i]+1)];
        }
        else //avant le premier terme non-nul de la ligne
        {
            cout<<"Coordonnees hors profil"<<endl;
            exit(-1);

            /*
            //On doit l'ajouter au bon endroit
            for (int k=j; k<Profil[i];++k)//termes � rajouter
            {
                Lower.add(Posdiag[i-1]+k+1,0.); //rajoute des 0 entre le nouveau point � �crire et celui existant
                Posdiag[i]+=1;//d�cale le terme diagonal de la ligne � chaque ajout d'un terme
            }
            Profil[i]=j; //nouveau premier terme non-nul
            return Lower[Posdiag[i-1]+Profil[i]+1];
            */
        }
    }
    else //triangle sup�rieur
    {
        if (i>=Profil[j]) //apr�s le premier terme non-nul de la colonne
        {
            return Upper[Posdiag[j-1]+(i-Profil[j]+1)];
        }
        else//avant le premier terme non-nul de la colonne
        {
            cout<<"Coordonnees hors profil"<<endl;
            exit(-1);

            /*
            //On doit l'ajouter au bon endroit
            for (int k=j; k<Profil[i];++k)//termes � rajouter
            {
                Upper.add(Posdiag[i-1]+k+1,0.); //rajoute des 0 entre le nouveau point � �crire et celui existant
                Posdiag[i]+=1;//d�cale le terme diagonal de la ligne � chaque ajout d'un terme
            }
            Profil[i]=j; //nouveau premier terme non-nul
            return Upper[Posdiag[i-1]+Profil[i]+1];
            */
        }
    }
}

matrice_nonsym& matrice_nonsym::operator*=(double a)
{
    Lower*=a;
    Upper*=a;
    return *this;
}


matrice_nonsym& matrice_nonsym::operator/=(double a)
{
    if (a==0.)
    {
        cout<<"division par z�ro"<<endl;
        exit(-1);
    }
    Lower/=a;
    Upper/=a;
    return *this;
}

matrice_nonsym& matrice_nonsym::operator+=(const matrice_nonsym& A)
{
    if (n!=A.n)
    {
        cout<<"Les matrices ne sont pas de m�mes dimensions"<<endl;
        exit(-1);
    }
    for (int i=0; i<n; ++i)
    {
        (*this)(i,i)+=A(i,i); //n'incr�mente que dans le lower
        //Lower[Posdiag[i]]+=A.Lower[Posdiag[i]];
        Upper[Posdiag[i]]+=A.Upper[Posdiag[i]];

        for (int j=A.Profil[i]; j<Posdiag[i]; ++j) //du premier non-nul jusqu'� la diagonale exclu
        {
            (*this)(i,j)+=A(i,j); //modification de Lower
            (*this)(j,i)+=A(j,i);//modification de Upper
        }

        /*
        if (A.profil[i]>=profil[i])//A n'enl�ve pas de 0 sur cette ligne
        {
            for (int j=A.Profil[i]; j<i; ++j) //du premier non-nul jusqu'� la diagonale exclu
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

matrice_nonsym& matrice_nonsym::operator-=(const matrice_nonsym& A)
{
    if (n!=A.n)
    {
        cout<<"Les matrices ne sont pas de m�mes dimensions"<<endl;
        exit(-1);
    }
    for (int i=0; i<n; ++i)
    {
        (*this)(i,i)-=A(i,i); //n'incr�mente que dans le lower
        //Lower[Posdiag[i]]+=A.Lower[Posdiag[i]];
        Upper[Posdiag[i]]-=A.Upper[Posdiag[i]];

        for (int j=A.Profil[i]; j<Posdiag[i]; ++j) //du premier non-nul jusqu'� la diagonale exclu
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
    return matrice_nonsym(A)*=a;
}

matrice_nonsym operator*(double a, const matrice_nonsym& A)
{
    return matrice_nonsym(A)*=a;
}

vecteur operator*(const matrice_nonsym& A, const vecteur& V)
{
    if (A.n!=V.dim())
    {
        cout<<"La matrice n'a pas la m�me dimension que le vecteur"<<endl;
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
    return matrice_nonsym(A)/=a;
}

matrice_nonsym operator+(const matrice_nonsym& A, const matrice_nonsym& B)
{
    return matrice_nonsym(A)+=B;
}

matrice_nonsym operator-(const matrice_nonsym& A, const matrice_nonsym& B)
{
    return matrice_nonsym(A)-=B;
}

matrice_nonsym transpose(const matrice_nonsym& A)
{
    matrice_nonsym Res(A.n,A.Profil);
    Res.Upper = A.Lower;
    Res.Lower = A.Upper;
    return Res;
}

/*
#############################################################
##################### Op�rations Mixtes #####################
#############################################################
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

matrice_nonsym LUdecomposition(const matrice_nonsym& A)
{
    int n=A.n;
    matrice_nonsym lu(n,A.Profil); //Utilise la conservation du profil p
    for (int i=1; i<=n; i++)
    {
        lu.Lower[lu.Posdiag[i-1]] = 1;//mets diagonale lower � 0 car accesseur classique rempli Upper et non lower
    }
    for (int p=1; p<=n; p++)
    {
        for (int j=p-1; j<n; j++) //boucle sur la matrice triangulaire sup�rieure
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
        for (int i=(p-1)+1; i<n; i++) //boucle sur la matrice triangulaire inf�rieure
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



vecteur resolsys(const matrice_nonsym A, const vecteur b) //Ax=b
{
    int n=A.n;
    matrice_nonsym lu=LUdecomposition(A);
    vecteur y(n);
    for (int i=1;i<n+1;++i)
    {
        double S=b[n-1];
        for (int j=1;j<i;++j)
        {
            S=S-lu(i-1,j-1)*y[j-1];
        }
        if (lu(i-1,i-1)==0)
        {
            exit(-1);
        }
        y[i-1]=S;
    }
    vecteur x(n);
    for (int i=n;i>0;--i)
    {
        double SS=y[i-1];
        for (int j=n;j>i-1;--j)
        {
            SS=SS-lu(i-1,j-1)*x[j-1];
        }
        if (lu(i-1,i-1)==0)
        {
            exit(-1);
        }
        x[i-1]=double(SS)/lu(i-1,i-1);
    }
    return(x);
}

/*
#############################################################
################## Matrices Elementaires ####################
#############################################################
*/

matrice matM_elem(const Point& P1,const Point& P2,const Point& P3)
{
    double x1 = P1.x;
    double y1 = P1.y;
    double x2 = P2.x;
    double y2 = P2.y;
    double x3 = P3.x;
    double y3 = P3.y;

    double D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
    matrice Mel(3,3);
    for (int i=0;i<3;++i)
    {
        for (int j=0;j<3;++j)
        {
            if (i==j)
            {
                Mel(i,j)=abs(D)/12;
            }
            else
            {
                Mel(i,j)=abs(D)/24;
            }
        }
    }
    return(Mel);

}

matrice matK_elem(const Point& P1,const Point& P2,const Point& P3)
{
    double x1 = P1.x;
    double y1 = P1.y;
    double x2 = P2.x;
    double y2 = P2.y;
    double x3 = P3.x;
    double y3 = P3.y;

    double D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
    matrice norm(3,2);
    norm(0,0)=y2-y3;
    norm(0,1)=x3-x2;
    norm(1,0)=y3-y1;
    norm(1,1)=x1-x3;
    norm(2,0)=y1-y2;
    norm(2,1)=x2-x1;
    matrice Kel(3,3);
    Point P(x1+(x2-x1)/6+(x3-x1)/6,y1+(y2-y1)/6+(y3-y1)/6);
    Point Q(x1+(x2-x1)*2/3+(x3-x1)/6,y1+(y2-y1)*2/3+(y3-y1)/6);
    Point R(x1+(x2-x1)/6+(x3-x1)*2/3,y1+(y2-y1)/6+(y3-y1)*2/3);

    Kel=norm*(A(P)/6 + A(Q)/6 + A(R)/6)*transpose(norm)/(3*abs(D));
    return(Kel);
}


matrice A(const Point& P)
{
    double x1=P.x;
    double x2=P.y;
    matrice A(2,2);
    A(0,0)=0.5*0.04*x1*x1;
    A(0,1)=0.5*(-0.024)*x1*x2;
    A(1,0)=0.5*(-0.024)*x1*x2;
    A(1,1)=0.5*0.04*x2*x2;
    return(A);
}

vecteur V(const Point& P)
{
    double x1=P.x;
    double x2=P.y;
    vecteur V(2);
    V[0]=(0.04+0.5*(-0.024)-0.05)*x1;
    V[1]=(0.04+0.5*(-0.024)-0.05)*x2;
    return(V);
}

matrice matB_elem(const Point& P1,const Point& P2,const Point& P3)
{
    double x1 = P1.x;
    double y1 = P1.y;
    double x2 = P2.x;
    double y2 = P2.y;
    double x3 = P3.x;
    double y3 = P3.y;

    double D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
    matrice norm(3,2);
    norm(0,0)=y2-y3;
    norm(0,1)=x3-x2;
    norm(1,0)=y3-y1;
    norm(1,1)=x1-x3;
    norm(2,0)=y1-y2;
    norm(2,1)=x2-x1;


    //omega(j,q)
    matrice omega_chap(3,3);
    omega_chap(0,0) = 1 - x1 - y1;
    omega_chap(1,0) = x1;
    omega_chap(2,0) = y1;

    omega_chap(0,1) = 1 - x2 - y2;
    omega_chap(1,1) = x2;
    omega_chap(2,1) = y2;

    omega_chap(0,2) = 1 - x3 - y3;
    omega_chap(1,2) = x3;
    omega_chap(2,2) = y3;

    //S^q
    matrice S(3,2);
    S(0,0) = 1/6;
    S(0,1) = 1/6;

    S(1,0) = 2/3;
    S(1,1) = 1/6;

    S(2,0) = 1/6;
    S(2,1) = 2/3;

    //Initialisation
    matrice Bel(3,3);

    //Constantes
    double omega_0 = 1/6;
    double k = -0.038;

    //Calculs Moches
    for (int i = 0; i<3; ++i)
    {

        for (int j = 0; j<3; ++j)
        {
            //S^q * norm (puis on ne garde que la valeur qui nous int�resse)
            matrice temp = norm*transpose(S);


            for (int q = 0; q<3; ++q)
            {
                Bel(i,j) += temp(i,q)*omega_chap(q,j);
            }
        }
    }

    Bel = k*omega_0*Bel;
    return(Bel);
}
/*
#############################################################
######################## Assemblage #########################
#############################################################
*/

int minimum(int a, int b){
    if (a <= b){return a;}
    else {return b;}

}

    //Assemblage des Matrices

matrice_sym matM(const vector<Point>& Coorneu, const list<Numeros>& Numtri){

    // Nombre de triangle : Nbtri

int Nbtri = Numtri.size();

    //D�finition du profil

vecteur Prof2(Nbtri);

    //Initialisation du profil
for (int I = 0; I < Nbtri; ++I){
    Prof2[I] = I;
}

    //On identifie le profil via des minimums
/*
for (int l = 0; l < Nbtri; ++l)
{
    for (int i = 0; i<3; ++i){
        int I = get(Numtri,l)(i);
        for (int j = 0; j<3; ++j){
            int J = get(Numtri,l)(j);
            if (J<I) { Prof2(I) = minimum(Prof2(I), J);}
        }
    }
}
*/
matrice Numtri_l(3,Nbtri);

for (int l = 0; l < Nbtri; ++l)
{
    list<Numeros>::const_iterator it = Numtri.begin();
        for(int n=0; n<l; n++){
            ++it;
        }
    for (int i = 0; i<3; ++i){
        int I = (*it)(i+1); //Parentheses =accesseur
        Numtri_l(i,l) = I;
        for (int j = 0; j<3; ++j){
            int J =(*it)(j+1);

            if (J<I)
                {
                    Prof2[I] = minimum(Prof2[I], J);
                }
        }
    }
}



matrice_sym MM(Nbtri, Prof2);


for (int l = 0; l < Nbtri; ++l) //"l" est un "L" + On va jamais aller l� o� il ne faut pas aller !
{
    int indice_1 = Numtri_l(0,l);
    int indice_2 = Numtri_l(1,l);
    int indice_3 = Numtri_l(2,l);
    //cout<<indice_1<<" "<<indice_2<<" "<<indice_3<<endl;
    matrice Mel = matM_elem(Coorneu[indice_1], Coorneu[indice_2], Coorneu[indice_3]);


    for (int i = 0; i < 3; ++i){

       int I=Numtri_l(i,l);

       for (int j = 0; j < 3; ++j){

           int J=Numtri_l(i,l);

           MM(I,J)=MM(I,J) + Mel(i,j);
       }
    }
}

return MM;


}



matrice_sym matK(const vector<Point>& Coorneu, const list<Numeros>& Numtri){

    // Nombre de triangle : Nbtri

int Nbtri = Numtri.size();

    //D�finition du profil

vecteur Prof2(Nbtri);

    //Initialisation du profil
for (int I = 0; I < Nbtri; ++I){
    Prof2[I] = I;
}


matrice Numtri_l(3,Nbtri);

for (int l = 0; l < Nbtri; ++l)
{
    list<Numeros>::const_iterator it = Numtri.begin();
        for(int n=0; n<l; n++){
            ++it;
        }

    for (int i = 0; i<3; ++i){

        int I = (*it)(i+1);
        Numtri_l(i,l) = I;


        for (int j = 0; j<3; ++j){

            int J =(*it)(j+1);

            if (J<I) { Prof2[I] = minimum(Prof2[I], J);}
        }
    }
}



matrice_sym KK(Nbtri, Prof2);


for (int l = 0; l < Nbtri; ++l) //"l" est un "L" + On va jamais aller l� o� il ne faut pas aller !
{
    int indice_1 = Numtri_l(0,l);
    int indice_2 = Numtri_l(1,l);
    int indice_3 = Numtri_l(2,l);

    matrice Kel = matK_elem(Coorneu[indice_1], Coorneu[indice_2], Coorneu[indice_3]);

    for (int i = 0; i < 3; ++i){

       int I=Numtri_l(i,l);

       for (int j = 0; j < 3; ++j){

           int J=Numtri_l(j,l);

           KK(I,J)=KK(I,J) + Kel(i,j);
       }
    }
}

return KK;


}

matrice_nonsym matB(const vector<Point>& Coorneu, const list<Numeros>& Numtri){

    // Nombre de triangle : Nbtri

int Nbtri = Numtri.size();

    //D�finition du profil

vecteur Prof2(Nbtri);

    //Initialisation du profil
for (int I = 0; I < Nbtri; ++I){
    Prof2[I] = I;
}

    //On identifie le profil via des minimums

matrice Numtri_l(3,Nbtri);

for (int l = 0; l < Nbtri; ++l)
{
    list<Numeros>::const_iterator it = Numtri.begin();
        for(int n=0; n<l; n++){
            ++it;
        }

    for (int i = 0; i<3; ++i){

        int I = (*it)(i+1);
        Numtri_l(i,l) = I;


        for (int j = 0; j<3; ++j){

            int J =(*it)(j+1);

            if (J<I) { Prof2[I] = minimum(Prof2[I], J);}
        }
    }
}


matrice_nonsym BB(Nbtri, Prof2);


for (int l = 0; l < Nbtri; ++l) //"l" est un "L" + On va jamais aller l� o� il ne faut pas aller !
{
    int indice_1 = Numtri_l(0,l);
    int indice_2 = Numtri_l(1,l);
    int indice_3 = Numtri_l(2,l);

    matrice Bel = matB_elem(Coorneu[indice_1], Coorneu[indice_2], Coorneu[indice_3]);
    //cout<<Bel<<endl;

    for (int i = 0; i < 3; ++i){

       int I=Numtri_l(i,l);

       for (int j = 0; j < 3; ++j){

           int J=Numtri_l(j,l);

           BB(I,J)=BB(I,J) + Bel(i,j);
       }
    }
}

return BB;

}


matrice_nonsym matD(const vector<Point>& Coorneu, const list<Numeros>& Numtri, double r){

        // Nombre de triangle : Nbtri

    int Nbtri = Numtri.size();

        //D�finition du profil

    vecteur Prof2(Nbtri);

        //Initialisation du profil
    for (int I = 0; I < Nbtri; ++I){
        Prof2[I] = I;
    }

    matrice Numtri_l(3,Nbtri);

    for (int l = 0; l < Nbtri; ++l)
    {
        list<Numeros>::const_iterator it = Numtri.begin();
            for(int n=0; n<l; n++){
                ++it;
            }
        for (int i = 0; i<3; ++i){
            int I = (*it)(i+1); //Parentheses =accesseur
            Numtri_l(i,l) = I;
            for (int j = 0; j<3; ++j){
                int J =(*it)(j+1);

                if (J<I)
                    {
                        Prof2[I] = minimum(Prof2[I], J);
                    }
            }
        }
    }



matrice_nonsym DD(Nbtri, Prof2);

    matrice_sym MM = matM(Coorneu, Numtri);
    matrice_sym KK = matK(Coorneu, Numtri);
    matrice_nonsym BB = matB(Coorneu, Numtri);
    //cout<<BB<<endl;

    DD = KK + BB + MM*r;
    return DD;
}


/*
#############################################################
####################### Resolution ##########################
#############################################################
*/

double pos(double x)
{
    if (x>0)
    {
        return(x);
    }
    return(0);
}

vecteur Q(double K,vecteur V)
{
    vecteur Q(V.dim_);
    for (int i=0;i<V.dim_;++i)
    {
        Point P=V[i];
        Q[i]=pos(P.x+P.y-K);
    }
    return(Q);
}


vecteur resolution_1(double deltaT, double K,const vecteur & V, const vector<Point> & sommets, const list<Numeros> & numelts)
{
  ofstream out("data_1");
  matrice_nonsym D = matD(sommets, numelts, 0.05);
  cout<<D<<endl;
  matrice_sym M = matM(sommets, numelts);
  matrice_nonsym E = M + deltaT*D;
  vecteur P = Q(K,V);
  out<<P<<endl;
  for (int k=0;k<K;++k)
  {
    vecteur temp(P);
    P = resolsys(E,M*temp);
    out<<P<<endl;
  }
  out.close();
  return P;
}

vecteur resolution_2(double deltaT, double K,const vecteur & V, const vector<Point> & sommets, const list<Numeros> & numelts)
{
  ofstream out("data_2");
  matrice_nonsym D = matD(sommets, numelts, 0.05);
  matrice_sym M = matM(sommets, numelts);
  matrice_nonsym E = M + deltaT/2*D;
  matrice_nonsym F = M - deltaT/2*D;
  vecteur P = Q(K,V);
  out<<P<<endl;
  for (int k=0;k<K;++k)
  {
    vecteur temp(P);
    P = resolsys(E,F*temp);
    out<<P<<endl;
  }
  out.close();
  return P;
}
