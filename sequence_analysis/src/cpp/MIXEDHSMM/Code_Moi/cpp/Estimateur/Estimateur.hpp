#ifndef ESTIMATEUR_HPP_INCLUDED
#define ESTIMATEUR_HPP_INCLUDED


class Estimateur
{

    private :
    /*************************************** Attributs *****************************************/

    // le modèle
    

   /*******************************************************************************************/


    public :
    /*************************************** Méthodes ****************************************/

    // Constructeur
    Estimateur ( );

    // méthode virtuelle (pour les classes filles) qui lance l'estimation
    virtual void fit ();





    // acces aux atributs
    

    /*****************************************************************************************/

};




class MCEM : public Estimateur
{

    private :
    /*************************************** Attributs *****************************************/
    int max_iteration_EM;    // nombre maximum de tour de boucle de l'algo EM  (0 <= m < max_iteration_EM)
    double seuil_arret_EM;   // seuil qui arrète l'algo EM si un critère d'arrêt devient inférieur ou égal à ce seuil

    double***** beta;     // tableau de double tel que : beta[gamma][l][i][k][j] = (j+1)-ième composante du vecteur \Beta^{\gamma+1,l+1}_{i,k}
                          // avec 0 <= gamma <= 3, [si gamma==0 : 0 <= l <= L-1] et [si gamma!=0 : 0 <= l <= 0], 0 <= i <= Card(E)-1, 0 <= k <= ça dépend
                          // si l'indice n'existe pas (k par exemple), alors indice=0 fixé (k=0 par exemple)

   /*******************************************************************************************/


    public :
    /*************************************** Méthodes ****************************************/

    // Constructeur
    MCEM ( );

    // méthode qui lance l'estimation
    virtual void fit ();







    // acces aux atributs
    

    /*****************************************************************************************/

};





#endif