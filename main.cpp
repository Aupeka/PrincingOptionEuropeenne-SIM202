#include "Projet.hpp"
//#include "matrix.hpp"
#include <iostream>

//Declaration Variables Globales
double r;

//Déclaration de Xi (matrice défini par un tableau)
double Xi[4];

/*
//Définition de A(x) (matrice défini par un tableau)
double A(vecteur x){
	double A[4];
	A[0] = Xi[0]*x(1)*x(1); A[1] = Xi[1]*x(1)*x(2); A[2] = Xi[2]*x(1)*x(2); A[3] = Xi[3]*x(2)*x(2);
}

//Définition de V(x)
vecteur V(const vecteur& x){
	vecteur V(2, 0.);
	V[1] = (Xi[0] + Xi[2]/2 - r)*x(1);
	V[2] = (Xi[3] + Xi[1]/2 - r)*x(2);

	return V;

}
*/
int main(){

//Remplissage Xi

//Xi[0] = 1; Xi[1] = 0; Xi[2] = 0; Xi[3] = 1; //L'identité au début 

//double a_1 = ...;
//double a_2 = ...;

int T = 50;	

double sigma_2 = 3; //Variance de l'actif, volatilité
double r = 0.04; //Taux d'interet

//Discrétisation --> Methode des éléments finis

double delta_T = 0.01; //Pas de récurrence
double a = 2.; //Limite du domaine

//Maillage

int Nbpt_cote = 50; //Nb de coupes de [0,a]
int Nbpt = Nbpt_cote*Nbpt_cote;
Maillage Maillage_b_s(0,a,0,a,Nbpt_cote,Nbpt_cote);

//Renomage donnée importante --> [Corneu = Sommets] & [Numtri = Numels]

//vector<Point> Coorneu = Maillage_b_s.sommets; //Coordonnée de tous les sommets
//list<Numeros> Numtri = Maillage_b_s.numelts; //Pour chaque triangle (en ligne) donne le numéro des 3 sommets



/*
// /!\ Penser à mettre Numtri (+1)

//Matrice de Masse
matrice_sym M(Nbpt);

//Matrice K
matrice_sym K(Nbpt);

//Matrice B
matrice_nonsym B(Nbpt);
*/

	//Construction du maillage d'ANN202

/*
Maillage M(0,1,0,1,20,20);
M.savecoord("Coorneu.txt");
M.savenumtri("Numtri.txt");
*/

return 0;

}







