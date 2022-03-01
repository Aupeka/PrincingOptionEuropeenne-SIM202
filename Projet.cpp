#include "Projet.hpp"
#include "matrix.hpp"
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
  if (d<0) stop("init() : dimension <= 0");
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
##############Classe matrice #################
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
        cout<<"Coordonnees hors profil"<<endl;
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
            cout<<"Coordonnees hors profil"<<endl;
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
            cout<<"Coordonnees hors profil"<<endl;
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

void LUdecomposition(const matrice& A,matrice& l,matrice& u, int n)
{
   int i = 0, j = 0, k = 1;
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


/*#####################################Matrice éléments finis##########################
*/


matrice_sym matM_elem(const Point& P1,const Point& P2,const Point& P3)
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
        for (int j=0;i<3;++j)
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

matrice_sym matK_elem(const Point& P1,const Point& P2,const Point& P3)
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

    Kel=(A(P)/6 + A(Q)/6 + A(R)/6)*norm*transpose(norm)/(3*abs(D));
  
    return Kel;
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
    vecteur V(2);
    double x1=P.x;
    double x2=P.y;
    V[0]=(0.04+0.5*(-0.024)-0.05)*x1;
    V[1]=(0.04+0.5*(-0.024)-0.05)*x2;
    return(V);
}

matrice_nonsym matB_elem(const Point& P1,const Point& P2,const Point& P3)
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
    omega_chap(1,1) = 1 - x1 - y1;
    omega_chap(2,1) = x1;
    omega_chap(3,1) = y1;

    omega_chap(1,2) = 1 - x2 - y2;
    omega_chap(2,2) = x2;
    omega_chap(3,2) = y2;

    omega_chap(1,3) = 1 - x3 - y3;
    omega_chap(2,3) = x3;
    omega_chap(3,3) = y3;

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
            //S^q * norm (puis on ne garde que la valeur qui nous intéresse)
            matrice temp = transpose(S)*norm;

            for (int q = 0; q<3; ++q)
            {
                Bel(i,j) += temp(q,j)*omega_chap(j,q);
            } 
        }
    }

    Bel *= k*omega_0*Bel;
  
    return Bel;
}







