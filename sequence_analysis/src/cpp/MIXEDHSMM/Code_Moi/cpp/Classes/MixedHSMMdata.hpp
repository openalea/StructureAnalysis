#ifndef MIXEDHSMMDATA_HPP_INCLUDED
#define MIXEDHSMMDATA_HPP_INCLUDED


/*##################################################################################################################
Objet contenant les données :
- les séquences observées ((Y^q_t)_t)_q, avec Y^q_t = ( Y^q_{t,1}, Y^q_{t,2}, ..., Y^q_{t,L} )  (L variables observées)
- les covariables (X^q)_q [non dynamique], ou bien ((X^q_t)_t)_q [dynamique]
##################################################################################################################*/

class MixedHSMMdata
{

    private :
    /*************************************** Attributs *****************************************/
    int nb_sequences;                // (Q) nombre de séquences
    int* lenght_sequences;           // liste d'entiers telle que : length_sequences[q]= m^q = la longueur de la séquence numéro q+1 (avec 0 <= q <= Q-1)

    int nb_variables_observed;       // (L) nombre de variables observées
    bool* type_variables_observed;   // liste de booléen telle que : type_variables_observed[l] = type de la variable observée Y^q_{t,(l+1)} (avec 0 <= l <= L-1), avec type=true si variable continue (double), et type=false si variable discrète (int)

    // les covariables sont assumées être fournies sous forme d'un tableau disjonctif complet (pour l'instant elles ne sont pas centrées-réduites, mais si c'est le cas, alors le type de chaque colonne de X^q sera forcément un double) 
    bool covariates_are_dynamics;    // true = covariables dynamiques, et false = covariables non dynamiques
    int nb_covariates;               // (J) nombre de colonnes de X^q_t si dynamique, ou X^q si non dynamique
    bool* type_covariates;           // liste de booléens telle que : type_covariates[j] = type de X^q_{t,j+1} (avec 0 <= j <= J-1), c'est à dire la (j+1)-ième colonne de X^q_t (ou X^q si covariables non dynamique), avec type=true si continue (double), et type=false si indicatrice (int)


    void*** variables_observed;      // tableau contenant des int ou des doubles tel que : variables_observed[q][t][l] = Y^{q+1}_{t,(l+1)}  (avec 0 <= q <= Q-1, 0 <= t <= m^q, 0 <= l <= L-1)
    void*** covariates;              // tableau contenant des int ou des doubles tel que : covariates[q][t][j] = X^{q+1}_{t,j+1}  (avec 0 <= q <= Q-1, 0 <= j <= J-1, [si dynamique : 0 <= t <= m^q ] ou [si non dynamique : 0 <= t <= 0 ]


    
   /*******************************************************************************************/


    public :
    /*************************************** Méthodes ****************************************/

    // Constructeur
    MixedHSMMdata ( );







    // acces aux atributs
    int get_nb_sequences() const { return nb_sequences; }
    int* get_lenght_sequences() const { return lenght_sequences; }

    int get_nb_variables_observed() const { return nb_variables_observed; }
    bool* get_type_variables_observed() const { return type_variables_observed; }

    bool get_covariates_are_dynamics() const { return covariates_are_dynamics; }
    int get_nb_covariates() const { return nb_covariates; }
    bool* get_type_covariates() const { return type_covariates; }

    void*** get_variables_observed() const { return variables_observed; }
    void*** get_covariates() const { return covariates; }


    /*****************************************************************************************/




};


#endif