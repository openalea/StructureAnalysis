#include <iostream>

#include "Estimateur.hpp"


using namespace std;


// **************************************** class Estimateur ****************************************

// Constructeur
Estimateur::Estimateur (  )
{
    
}

// méthode virtuelle (pour les classes filles) qui lance l'estimation
void Estimateur::fit ()
{
    cout << "Type d'estimateur non sélectionné" << endl;
}







// **************************************** class MCEM ****************************************


// Constructeur
MCEM::MCEM (  )
{
    

}


// méthode qui lance l'estimation via MCEM
void MCEM::fit ()
{

    // On initialise les paramètres initiaux


    int m = 0;                   // tour de boucle numéro (m) de l'algo EM
    double critere_arret = 0;    // arrète l'algo EM si ce critère devient inférieur ou égal à seuil_arret_EM
    while ( (m < max_iteration_EM) & (critere_arret > seuil_arret_EM) )
    {


        // Etape E




        // Etape M




        m++;
    }

}