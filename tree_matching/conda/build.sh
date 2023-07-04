#!/bin/bash

if [[ -d build ]]; then
    rm -rf build
fi
mkdir build
cd build



echo "****** CMAKE CONFIG"

export GMPDIR=${PREFIX}

cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} \
      -DCMAKE_PREFIX_PATH=${PREFIX} \
      -DCMAKE_BUILD_TYPE=Release  \
      -DPython3_EXECUTABLE=${PYTHON} \
       ${SYSTEM_DEPENDENT_ARGS[@]} \
      -LAH .. 

echo
echo "****** COMPILE"
export VERBOSE=1
make -j${CPU_COUNT} 
echo "****** INSTALL CXX LIB"
make install

echo
echo "****** INSTALL PYTHON LIB"
cd ..

pip install --prefix=${PREFIX} 

echo "****** END OF BUILD PROCESS"
