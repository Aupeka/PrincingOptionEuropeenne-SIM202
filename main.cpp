#include <iostream>
#include "Projet.hpp"

using namespace std;

int main()
{

    Maillage M(0,14.,0,14.,10,10);
    Point P(0,2);
    vecteur VV =V(P);
    //matrice_nonsym BB = matB(M.sommets, M.numelts);
    //cout<<BB<<endl;
    vecteur PP = resolution_1(0.1,5,VV,M.sommets,M.numelts);
    cout<<PP<<endl;
    return 0;
}
