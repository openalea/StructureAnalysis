rm src/py/wrapper/*.cpp
rm src/py/wrapper/*.h
source activate vplants && scons; source deactivate vplants
cd src/cpp && source activate statiskit-dev && python AutoWIG.py; source deactivate statiskit-dev && cd ../..
source activate vplants && scons; source deactivate vplants