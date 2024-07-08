set -ve

scons py -f SConstructWIG --prefix=$PREFIX -j$CPU_COUNT --visibility=default

set +ve