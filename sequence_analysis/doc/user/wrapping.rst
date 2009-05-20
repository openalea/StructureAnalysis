########
WRAPPING
########


.. todo:: check all wrapper with optional arguments and switch to BOOST_OVERLOAD ? 

100% means either fully done or check the comments in the wrapping files. Maybe some functions have not been exported because irrelevant or not enough information.

known issues
============

* ComputeSelfTransition: the prototype is protected. Only the empty argument case has been wrapped in export_markovian.cpp
