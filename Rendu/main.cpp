#include <iostream>
#include "Projet.hpp"

using namespace std;

int main()
{

    Maillage M(0,14.,0,14.,30,30);
    M.savecoord("Coorneu.txt");
    M.savenumtri("Numtri.txt");
    vecteur PP = resolution_1(0.01,10,M.sommets,M.numelts);
    cout<<"résultat d'ordre 1:"<<PP<<endl;
    vecteur PP2 = resolution_2(0.001,25,M.sommets,M.numelts);
    cout<<"résultat d'ordre 2:"<<PP2<<endl;
    return 0;
}
