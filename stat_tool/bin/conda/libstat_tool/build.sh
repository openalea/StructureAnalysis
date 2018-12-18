set -ve

scons cpp -f SConstructWIG --prefix=$PREFIX -j$CPU_COUNT --visibility=default

set +ve