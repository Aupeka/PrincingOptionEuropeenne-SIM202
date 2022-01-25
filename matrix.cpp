
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

double matrice::operator[](int i, int j) const
{
    if (i>m || j>n || i<0 || j<0)
    {
        cout <<"Indice en dehors des bornes"<<endl;
        exit(-1);
    }
    return cols_[j][i]; //surcharge du vecteur
}

double& matrice::operator[](int i, int j)
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
    if (A.n!=u.n)
    {
        cout<<"produit matrice vecteur : dimensions incompatibles"<<endl;
        exit(-1);
    }
    vecteur Res(A.m,0.);
    for (int i=1; i<<A.m; ++i)
    {
        for (int j=1; j<A.n;++j)
        {
           Res[i] = A[i,j]*u[j];
        }
    }
    return Res;
}
