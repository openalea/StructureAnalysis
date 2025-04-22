#include <iostream>
#include "classe.hpp"

using namespace std;


int main()
{

    double tab[5];
    for (int i=0; i<5; i++)
    {
        tab[i] = 1;
    }

    CovariateData nom (false, 5, 3, tab);


    return 0;
}
