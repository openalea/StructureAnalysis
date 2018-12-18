echo ON

scons py -f SConstructWIG --prefix=%LIBRARY_PREFIX% -j%CPU_COUNT% --arch=%ARCH%
if errorlevel 1 exit 1

echo OFF