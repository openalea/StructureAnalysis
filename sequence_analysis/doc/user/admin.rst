
Current Developments
####################

.. contents::



documentation
=============

The sphinx documentation is now used to manage the documentation so as to 
include a reference guide as well as a User guide (and this administration
section).

docstrings need to be written.


package CPP structure
=====================

independant classes:

*        backward;
*        Correlation;
*        function;
*        length_bias;
*        nb_event;
*        non_homogenous_markov;
*        nonparametric_sequence_process;
*        renewal;
*        renewal_iterator;
*        self_transition;
*        semi_markov_iterator
*        sequence_characteritics
*        top_paramerters
*        variable_order_chain_data
*        variable_order_markov_iterator

   
.. graphviz::
 
    digraph G1{
        tops -> sequences;
        non_homogenous_markov_data->markovian_sequences -> sequences;
        semi_markov_data->markovian_sequences -> sequences;
        variable_order_markov_dataa->markovian_sequences -> sequences;
    }


.. graphviz::

    digraph G2{
        hidden_semi_markov->semi_markov
    }

.. graphviz::

    digraph G3{
        hidden_variable_semi_markov->variable_order_markov
    }

.. graphviz::

    digraph G4{
        renewal_data->time_events
    }


Wrapper and validation (tests)
==============================

work in progress
----------------

=============== =================== =================== =============== =========== ==========
CPP file        Boost wrapping      python module       docstring       test        sphinx doc
=============== =================== =================== =============== =========== ==========
Sequences       80%                 sequences.py(10%)   10%             50%         0%
Tops            10%
=============== =================== =================== =============== =========== ==========


Known issues
------------
* Sequence: reverse method does not seem to work
* Sequences:several moving average:find apporpriate name
* Sequences: see wrapping.rst for a list of methods that have not yet beem implemented

