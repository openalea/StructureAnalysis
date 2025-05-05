#include <iostream>
#include "MixedHSMMdata.hpp"
#include "MixedHSMMmodel.hpp"

using namespace std;


// Prototype avec les valeurs par défaut
int* nombreDeSecondes(int heures, int minutes = 0, int secondes = 0);



// Main
int main()
{

    int T[10][10][10];

    for (int i=0; i<10; i++)
    {
        for (int j=0; j<10; j++)
        {
            for (int k=0; k<10; k++)
            {
                T[i][j][k] = i+j+k;
            }
        }
    }

    cout << T[2][5][3] << endl;
    cout << T << endl;

    ////////


    int ***t = new int**[10];

    for (int i=0; i<10; i++)
    {
        t[i] = new int*[10];
        for (int j=0; j<10; j++)
        {
            t[i][j] = new int[10];
        }
    }

    for (int i=0; i<10; i++)
    {
        for (int j=0; j<10; j++)
        {
            for (int k=0; k<10; k++)
            {
                t[i][j][k] = i+j+k;
            }
        }
    }

    cout << t[2][5][3] << endl;
    cout << t << endl;


    ///////////////

    
    int *W = new int[10];
    
    for (int i=0; i<10; i++)
    {
        W[i] = i+100;
    }

    int a = 0;

    cout << W[a] << endl;

    delete [] W;
    W = nullptr;

    cout << W[a] << endl;


    ////////

    /*
    int* b;
    int heures = 5;
    b = nombreDeSecondes(heures);

    cout << *b << endl;

    *b = 9;
    cout << *b << endl;
    */
    



    return 0;
}




// Définition de la fonction, SANS les valeurs par défaut
int* nombreDeSecondes(int heures, int minutes, int secondes)
{
    int total = 0;

    total = heures * 60 * 60;
    total += minutes * 60;
    total += secondes;

    int* ptr = new int;
    *ptr = 5;

    int a = 7;
    int* ptr2 = &a;


    return ptr;
}

