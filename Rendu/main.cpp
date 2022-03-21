#include <iostream>
#include "Project.hpp"

using namespace std;

int main()
{

    Maillage M(0,14.,0,14.,100,100);
    M.savecoord("Coorneu.txt");
    M.savenumtri("Numtri.txt");
    Point P(0,2);
    vecteur VV =V(P);
    vecteur PP = resolution_1(0.1,50,M.sommets,M.numelts);
    cout<<"résultat d'ordre 1:"<<PP<<endl;
    vecteur PP2 = resolution_2(0.1,50,M.sommets,M.numelts);
    cout<<"résultat d'ordre 2:"<<PP2<<endl;
    return 0;
}
