
################Classe Matrice##################
################################################

class matrice
{
public:
    vecteur* cols_;
    int m,n;
    matrice (int mi=0, int ni=0; double vi=0.);
    matrice(const matrice& v);
    ~matrice();
    matrice& operator=(const matrice& A);
    double operator[](int i, int j) const;
    double& operator[](int i, int j);
};

vecteur produit(const matrice& A, const vecteur& u);


class matrice_profil
{
public:
    int n,m;
    vecteur Profil;
    vecteur Posdiag;
    vecteur Lower;
    vecteur Upper;
    matrice_profil(int ni, int mi); //constructeur de la matrice vide
    matrice_profil(const matrice_profil& A); //constructeur par copie
    ~matrice_profil();
    matrice_profil& operator=(const matrice_profil& A);
    double operator[](int i, int j) const;
    double& operator[](int i, int j);
};


class matrice_sym : public matrice_profil
{
public:

};

class matrice_nonsym :public matrice_profil
{
public:

};
