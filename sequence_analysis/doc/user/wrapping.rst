########
WRAPPING
########


.. todo:: check all wrapper with optional arguments and switch to BOOST_OVERLOAD ? 


wrapping to be done
===================


* Remains to be done or check carefully: Estimation, Simulate and Extract*
* Fully done : Compare and all other "simple" functions.

* time_events: estimation ? 
* renewal: function to extct forward, backward and so on.
* renewal_data: estimation ?

known issues
============

* ComputeSelfTransition: the prototype is protected. Only the empty argument case has been wrapped in export_markovian.cpp
