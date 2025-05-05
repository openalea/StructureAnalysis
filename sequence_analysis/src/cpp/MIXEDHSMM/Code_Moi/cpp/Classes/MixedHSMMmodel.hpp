#ifndef MIXEDHSMMMODEL_HPP_INCLUDED
#define MIXEDHSMMMODEL_HPP_INCLUDED


/*##################################################################################################################
Objet contenant les paramètres du modèle :
- les paramètres Beta
- les paramètres sigma^2
##################################################################################################################*/

class MixedHSMMmodel
{

    private :
    /*************************************** Attributs *****************************************/
    double***** beta;     // tableau de double tel que : beta[gamma][l][i][k][j] = (j+1)-ième composante du vecteur \Beta^{\gamma+1,l+1}_{i,k}
                          // avec 0 <= gamma <= 3, [si gamma==0 : 0 <= l <= L-1] et [si gamma!=0 : 0 <= l <= 0], 0 <= i <= Card(E)-1, 0 <= k <= ça dépend
                          // si l'indice n'existe pas (k par exemple), alors indice=0 fixé (k=0 par exemple)

    double**** sigma2;   // tableau de double tel que : sigma2[gamma][l][i][k] = \sigma^2_{\gamma,l,i,k}
                          // avec 0 <= gamma <= 3, [si gamma==0 : 0 <= l <= L-1] et [si gamma!=0 : 0 <= l <= 0], 0 <= i <= Card(E)-1, 0 <= k <= ça dépend
                          // si l'indice n'existe pas (k par exemple), alors indice=0 fixé (k=0 par exemple)



   /*******************************************************************************************/


    public :
    /*************************************** Méthodes ****************************************/

    // Constructeur
    MixedHSMMmodel ( );







    // acces aux atributs
    

    /*****************************************************************************************/




};


#endif
