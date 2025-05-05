#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <iostream>

namespace py = pybind11;
using namespace std;



// *********************************** Ajouter une fonction *********************************

int add (int i=1, int j=2) {
    return i + j;
}


// pour effectuer des contrôle sur les indices / sécurié (je crois uniquement pour data)
// ne fonctionne pas pour le multidimensionel (en fait c'est plus compliqué)
int fois_deux (py::array_t<double> data) {    // data est un np.array (dtype=np.float64)
    py::buffer_info buf = data.request();
    double* ptr = static_cast<double*> (buf.ptr);
    int const x = buf.shape[0];

    /*
    int const y1 = buf.ndim;       // dimension de data
    int const y2 = buf.shape[k];   // longueur de la k-ième diemnsion  (avec 0 <= k <= buf.ndim -1)
    */

    for (int i=0; i<x; i++) {
        ptr[i] = ptr[i]*2;         // ptr[i][j][k] ne fonctionne pas pour multidimensionel
    }

    return x;
}


// accès direct à l'objet, sans contrôle. équivalent à pointeur brut en c++
double fois_deux_bis (py::array_t<double>& data) {    // data est un np.array (dtype=np.float64)
    auto mat = data.mutable_unchecked<3>();
    int const n1 = mat.shape(0);
    int const n2 = mat.shape(1);
    int const n3 = mat.shape(2);

    /*
    data.mutable_unchecked<N>()      // lecture et écriture    (N = nombre de dimensions connue à la compilation)
    data.unchecked<N>()              // lecture seule

    int const y1 = mat.ndim();       // dimension de data
    int const y2 = mat.shape(k);     // longueur de la k-ième diemnsion  (avec 0 <= k <= mat.ndim() -1)
    */

    double sum = 0;
    for (int i=0; i<n1; i++)
    {
        for (int j=0; j<n2; j++)
        {
            for (int k=0; k<n3; k++)
            {
                sum += mat(i, j, k);
                mat(i, j, k) = mat(i, j, k)*2;
            }
        }
    }

    return sum;
}





    
double*** numpy_to_tableau (py::array_t<double>& np_array) {    // np_array est un np.array 3D (dtype=np.float64)
    py::buffer_info buf = np_array.request();
    double* ptr = static_cast<double*> (buf.ptr);               // np_array est équivalent à un tableau 1D en mémoire, et est pointé par ptr

    int const n1 = buf.shape[0]
    int const n2 = buf.shape[1]
    int const n3 = buf.shape[2]

    tableau = new double***[n1];
    for (int i=0, i<n1, i++) {
        

    }




}




// ******************************************************************************************



// *********************************** Ajouter une Class ************************************









// ******************************************************************************************


PYBIND11_MODULE(example, m) {
    m.def("add", &add, "A function that adds two numbers", py::arg("i")=1, py::arg("j")=2);
    m.def("fois_deux", &fois_deux, "blabla");
    m.def("fois_deux_bis", &fois_deux_bis, "blabla");
}






