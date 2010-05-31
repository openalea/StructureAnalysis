.. _stat_tool_syntax:

.. |leq|   unicode:: U+02264 
.. |geq|   unicode:: U+02265 

.. testsetup:: *
      
    from openalea.stat_tool import *


.. .. include:: alias.rst

File Syntax, STAT module
########################

.. contents::


An ASCII file format is defined for each of the following object type of the STAT module:


4.7.1 type COMPOUND
===================

A compound (or stopped-sum) distribution is defined as the distribution of the sum of n independent and identically distributed random variables :math:`X_i` where :math:`n` is the value taken by the random variable :math:`N`. The distribution of :math:`N` is referred to as the sum distribution while the distribution of the :math:`X_i` is referred to as the elementary distribution. Consider the following example:

.. literalinclude:: syntax_compound.dat

you can load this data by using


.. doctest::
    
    >>> compound = Compound('doc/user/syntax_compound.dat')

The first line gives the distribution type. The parametric sum distribution and the parametric elementary distribution are then defined in subsequent lines according to the syntactic form defined for the type DISTRIBUTION.

4.7.2 type CONVOLUTION
======================

The distribution of the sum of independent random variables is the convolution of the distributions of these elementary random variables. Consider the following example:

.. literalinclude:: syntax_convolution.dat

you can load this data by using


.. doctest::
    
    >>> convolution = Convolution('doc/user/syntax_convolution.dat')


The first line gives the distribution type and the number of elementary distributions (2 or 3). The elementary parametric distributions are then defined in subsequent lines according to the syntactic form defined for the type DISTRIBUTION.


4.7.3 type DISTRIBUTION, type RENEWAL
=====================================

The available parametric discrete distributions are the binomial distribution, the Poisson distribution, the negative binomial distribution and the uniform (rectangular) distribution with an additional shift parameter which defines the lower bound to the range of possible values. The name of the distribution is first given, then the name of each parameter followed by its actual value as shown in the following examples:

.. literalinclude:: syntax_distribution.dat

you can load this data by using


.. doctest::
       
    >>> binomial = Distribution('doc/user/syntax_distribution.dat')

INF_BOUND and SUP_BOUND are integer-valued parameters while PARAMETER and PROBABILITY are real-valued parameters.

For every parametric distributions, the following constraint applies to the shift parameter:
0 |leq| INF_BOUND |leq| 200

For a BINOMIAL or a UNIFORM distribution, the following constraint applies to the parameters INF_BOUND and SUP_BOUND which define the range of possible values:
0 < SUP_BOUND - INF_BOUND |leq| 500

For a BINOMIAL distribution, the following constraint applies to the probability of 'success':
0 |leq| PROBABILITY |leq| 1

For a POISSON distribution, the following constraint applies to the parameter (which is equal to the mean):
0 < PARAMETER _ 200

For a NEGATIVE_BINOMIAL distribution, the following constraints apply to the parameters:
0 < PARAMETER
0 < PROBABILITY < 1
PARAMETER (1 - PROBABILITY) / PROBABILITY |leq| 200.

Pour une loi de type UNIFORM, les contraintes suivantes sur les paramètres doivent être respectées :
0 _ SUP_BOUND - INF_BOUND _ 500

A renewal process is built from a discrete parametric distribution (BINOMIAL, POISSON or NEGATIVE_BINOMIAL) termed the inter-event distribution which represents the time interval between consecutive events. Hence, the types DISTRIBUTION and RENEWAL share the same ASCII file format.

4.7.4 type HIDDEN_MARKOV
========================
A hidden Markov chain is constructed from an underlying Markov chain and nonparametric observation (or state-dependent) distributions. Consider the following example::

    HIDDEN_MARKOV_CHAIN

    2 STATES
    ORDER 1

    INITIAL_PROBABILITIES
    0.8  0.2

    TRANSITION_PROBABILITIES
    0.6  0.4
    0.1  0.9

    OBSERVATION_PROBABILITIES

    2 SPACES

    SPACE 1

    STATE 0
    OBSERVATION 0 : 1.0

    STATE 1
    OBSERVATION 0 : 0.2
    OBSERVATION 1 : 0.8

    SPACE 2

    STATE 0
    OBSERVATION 0 : 0.2
    OBSERVATION 1 : 0.4
    OBSERVATION 2 : 0.4
    
    STATE 1
    OBSERVATION 0 : 0.8
    OBSERVATION 1 : 0.1
    OBSERVATION 2 : 0.1

The first line gives the object type. The underlying Markov chain is then defined on subsequent lines according to the syntactic form defined for the type MARKOV. The observation (or state-dependent) probabilities relating the output processes to the non-observable state process are then defined. Since the process is 'hidden', at least one possible output should be observable in more than one state.

4.7.5 type HIDDEN_SEMI-MARKOV
=============================

A hidden semi-Markov chain is constructed from an underlying semi-Markov chain (first-order Markov chain representing transition between distinct states and state occupancy distributions associated to the non-absorbing states) and nonparametric observation (or state-dependent) distributions. The state occupancy distributions are defined as objects of type DISTRIBUTION with the additional constraint that the minimum time spent in a given state is 1 (INF_BOUND |leq| 1). Consider the following example::

    HIDDEN_SEMI-MARKOV_CHAIN

    4 STATES

    INITIAL_PROBABILITIES
    0.8  0.2  0.0  0.0

    TRANSITION_PROBABILITIES
    0.0  0.6  0.4  0.0
    0.0  0.0  0.7  0.3
    0.0  0.2  0.0  0.8
    0.0  0.0  0.0  1.0

    STATE 0 OCCUPANCY_DISTRIBUTION
    NEGATIVE_BINOMIAL  INF_BOUND : 2  PARAMETER : 3.2  PROBABILITY : 0.4

    STATE 1 OCCUPANCY_DISTRIBUTION
    BINOMIAL  INF_BOUND : 1  SUP_BOUND : 12  PROBABILITY : 0.6

    STATE 2 OCCUPANCY_DISTRIBUTION
    POISSON  INF_BOUND : 1  PARAMETER : 5.4

    OBSERVATION_PROBABILITIES
    
    1 SPACE

    SPACE 1

    STATE 0
    OBSERVATION 0 : 1.0

    STATE 1
    OBSERVATION 0 : 0.3
    OBSERVATION 1 : 0.6
    OBSERVATION 2 : 0.1

    STATE 2
    OBSERVATION 0 : 0.2
    OBSERVATION 1 : 0.4
    OBSERVATION 2 : 0.4

    STATE 3
    OBSERVATION 2 : 1.0

Note that absorbing states such as state 3 :math:`(p_{33}=1)` are by nature Markovian. It is also possible to define nonabsorbing Markovian states such as state 2 :math:`(0 < p_{22} < 1)`. In this case, the resulting model is a hybrid hidden Markov/semi--Markov chain.

The first line gives the object type. The underlying semi-Markov chain (embedded first-order Markov chain and state occupancy distributions associated to the nonabsorbing states) is then defined on subsequent lines according to the syntactic form defined for the type SEMI-MARKOV. The observation (or state-dependent) probabilities relating the output processes to the non-observable state process are then defined. Since the process is 'hidden', at least one possible output should be observable in more than one state.

4.7.6 type HISTOGRAM
====================
The syntactic form of the type HISTOGRAM consists in giving, in a first column, the values in increasing order and, in a second column, the corresponding frequencies. If a value is not given, the corresponding frequency is assumed to be null. Consider the following example:

.. literalinclude:: syntax_histogram.dat

you can load this data by using


.. doctest::

    >>> histogram = Histogram('doc/user/syntax_histogram.dat')

4.7.7 type MARKOV
=================
Consider the following example of an homogeneous Markov chain::

    MARKOV_CHAIN

    2 STATES
    ORDER 2

    INITIAL_PROBABILITIES
    0.8  0.2

    TRANSITION_PROBABILITIES
    0.6  0.4
    0.1  0.9
    0.3  0.7
    0.2  0.8

The first line gives the object type. Then, the number of states (between 2 and 15) and the order (between 1 and 4) are defined on the two subsequent lines. On the next lines, the initial probabilities and the transition probabilities are given. Since, the initial probabilities and the transition probabilities for a given memory constitute distributions, the elements of a line should sum to one.

It is also possible to define observation (or state-dependent) probabilities if each possible output can be observed in a single state. With this restriction, the state space corresponds to a partition of the output space and the overall process is a lumped process::

    OBSERVATION_PROBABILITIES

    2 SPACES
    
    SPACE 1
    
    STATE 0
    OBSERVATION 0 : 1.0

    STATE 1
    OBSERVATION 1 : 0.2
    OBSERVATION 2 : 0.8

    SPACE 2

    STATE 0
    OBSERVATION 0 : 0.7
    OBSERVATION 1 : 0.3
    
    STATE 1
    OBSERVATION 2 : 0.6
    OBSERVATION 3 : 0.4

Consider the following example of a non-homogeneous Markov chain::

    NON-HOMOGENEOUS_MARKOV_CHAIN
    
    3 STATES
    ORDER 1
    
    INITIAL_PROBABILITIES
    0.5  0.3  0.2
    
    TRANSITION_PROBABILITIES
    0.6  0.2  0.2
    0.1  0.8  0.1
    0.2  0.1  0.7
    
    STATE 0 HOMOGENEOUS
    
    STATE 1 NON-HOMOGENEOUS
    MONOMOLECULAR FUNCTION  PARAMETER 1 : 0.99  PARAMETER 2 : -0.34  PARAMETER 3 : 0.3
    
    STATE 2 NON-HOMOGENEOUS
    LOGISTIC FUNCTION  PARAMETER 1 : 0. 99  PARAMETER 2 : 2.8  PARAMETER 3 : 0.2

The first line gives the object type. Then, the initial probabilities and the transition probabilities are given in the same way as for an homogeneous Markov chain. The non-homogeneous / homogeneous character is then defined state by state. In the case of a non-homogeneous transition distribution, the function :math:`p_{ii}(t)` represents the self-transition in state `i` as a function of the index parameter `t`. The corresponding transition distribution defined in the transition probability matrix gives the relative weights of the probabilities of leaving state `i`.

For a MONOMOLECULAR function :math:`\left(p_{ii}(t)=a+b \exp{(-ct)}\right)`, the following constraints apply:

* 0 |leq| PARAMETER 1 |leq| 1
* 0 |leq| PARAMETER 1 + PARAMETER 2 |leq| 1
* PARAMETER 3 > 0
    
For a MONOMOLECULAR function :math:`\left(p_{ii}(t)=a/ \{ 1+b \exp{(-ct)}\}\right)`, the following constraints apply:

* 0 |leq| PARAMETER 1 |leq| 1
* 0 |leq| PARAMETER 1 / (1. + PARAMETER 2) |leq| 1
* PARAMETER 3 > 0

4.7.8 type MIXTURE
==================
A mixture is a parametric model of classification where each elementary distribution or component represents a class with its associated weight. Consider the following example::

    MIXTURE 2 DISTRIBUTIONS

    DISTRIBUTION 1  WEIGHT : 0.3
    BINOMIAL  INF_BOUND : 2  SUP_BOUND : 5  PROBABILITY : 0.8

    DISTRIBUTION 2  WEIGHT : 0.7
    NEGATIVE_BINOMIAL  INF_BOUND : 5  PARAMETER : 3.2  PROBABILITY : 0.4

The first line gives the distribution type and the number of components of the mixture (between 2 and 4). The components are then defined on two lines, the first one giving the associated weight and the second one giving the definition of the elementary parametric distribution according to the syntactic form defined for the type DISTRIBUTION. The weights should sum to one.

4.7.9 type SEMI-MARKOV
======================
A semi-Markov chain is constructed from a first-order Markov chain representing transition between distinct states and state occupancy distributions associated to the nonabsorbing states. The state occupancy distributions are defined as objects of type DISTRIBUTION with the additional constraint that the minimum time spent in a given state is at least 1 (INF_BOUND |leq| 1). Consider the following example::

    SEMI-MARKOV_CHAIN

    4 STATES

    INITIAL_PROBABILITIES
    0.8  0.2  0.0  0.0

    TRANSITION_PROBABILITIES
    0.0  0.6  0.4  0.0
    0.0  0.0  0.7  0.3
    0.0  0.2  0.0  0.8
    0.0  0.0  0.0  1.0

    STATE 0 OCCUPANCY_DISTRIBUTION
    NEGATIVE_BINOMIAL  INF_BOUND : 2  PARAMETER : 3.2  PROBABILITY : 0.4

    STATE 1 OCCUPANCY_DISTRIBUTION
    BINOMIAL  INF_BOUND : 1  SUP_BOUND : 12  PROBABILITY : 0.6

    STATE 2 OCCUPANCY_DISTRIBUTION
    POISSON  INF_BOUND : 1  PARAMETER : 5.4

The first line gives the object type while the second line gives the number of states (between 2 and 15). The embedded first-order Markov chain is then defined on subsequent lines by its initial probabilities and its transition probabilities (note that, unlike for the type MARKOV, the order should not be specified). Since this embedded Markov chain represents only transitions between distinct states, the self-transitions (i.e. elements of the main diagonal) should be equal to zero except in the case of absorbing states where the self-transitions are equal to one (e.g. state 3 in the above example). The state occupancy distributions are then defined for each nonabsorbing state according to the syntactic form defined for the type DISTRIBUTION with the additional constraint that time spent in a given state is at least 1 (INF_BOUND |leq| 1). Like for the type MARKOV, observation (or state-dependent) probabilities can be defined in order to specify a lumped process (with the restriction that each possible output can be observed in a single state).

Note that absorbing states such as state 3 :math:`(p_{33}=1)` are by nature Markovian. It is also possible to define nonabsorbing Markovian states such as state 2 :math:`(0 < p_{22} < 1)`. In this case, the resulting model is a hybrid hidden Markov/semi--Markov chain.


4.7.10 type SEQUENCES
=====================
The syntactic form of the type SEQUENCES is constituted of a header giving the number and the type of variables and of the sequence. Consider the following example of univariate sequences::

    1 VARIABLE

    VARIABLE 1 : STATE

    1 0 0 0 1 1 2 0 2 2 2 1 1 0 1 0 1 1 1 1 0 1 1 1 \
    0 1 2 2 2 1

    0 0 0 1 1 0 2 0 2 2 2 1 1 1 1 0 1 0 0 0 0 0

The type STATE is the generic type. The character '\' enables to continue a sequence on the following line.

Consider the following example of multivariate sequences::

    2 VARIABLES

    VARIABLE 1 : STATE
    VARIABLE 2 : STATE
    
    1 0 | 0 0 | 1 0 | 2 0 | 2 1 | 2 1 | 1 0 | 1 0 | 1 0 | 0 1 | 0 1 | 1 1 \
    0 1 | 2 0 | 2 1
    
    0 0 | 0 0 | 1 0 | 2 0 | 2 1 | 1 1 | 1 0 | 1 0 | 0 0 | 0 0

The character '|' enables to separate successive vectors.

Consider the following example of sequences with an explicit index parameter of type POSITION::

    2 VARIABLES

    VARIABLE 1 : POSITION
    VARIABLE 2 : STATE

    10 1 | 12 0 | 13 1 | 14 2 | 15 2 | 20 2 | 22 1 | 23 1 | 27 1 | 30 0 | 31 0 | 32 1 \
    35 1 | 37 0 | 40 1 | 45

    5 0 | 7 0 | 10 0 | 11 0 | 15 1 | 18 1 | 20 0 | 21 0 | 22 0 | 25 0 | 25

This explicit index parameter is given as a first variable and the other variables (at least one) should be of type STATE. The index values should be increasing along sequences and the sequence ends with a final index value.

The explicit index parameter of type POSITION can be replaced by inter-position intervals::

    2 VARIABLES

    VARIABLE 1 : POSITION_INTERVAL
    VARIABLE 2 : STATE

    10 1 | 2 0 | 1 1 | 1 2 | 1 2 | 5 2 | 2 1 | 1 1 | 4 1 | 3 0 | 1 0 | 1 1 \
    3 1 | 2 0 | 3 1 | 5

    5 0 | 2 0 | 3 0 | 1 0 | 4 1 | 3 1 | 2 0 | 1 0 | 1 0 | 3 0 | 0

Consider the following example of sequences with an explicit index parameter of type TIME::

    2 VARIABLES

    VARIABLE 1 : TIME
    VARIABLE 2 : STATE
    
    3 1 | 7 4 | 10 8 | 14 10 | 18 15 | 21 16 | 25 18 | 28 19 | 31 20 | 35 22 | 39 23 | 42 24 \
    45 25 | 49 25
    
    3 1 | 7 2 | 10 6 | 14 9 | 18 13 | 21 14 | 25 15 | 28 16 | 31 17 | 35 17
    
The only difference with the explicit index parameter of type POSITION is that the index values should be strictly increasing along sequences and that no final index value is required.

The explicit index parameter of type TIME can be replaced by time intervals::

    2 VARIABLES
    
    VARIABLE 1 : TIME_INTERVAL
    VARIABLE 2 : STATE
    
    3 1 | 4 4 | 3 8 | 4 10 | 4 15 | 3 16 | 4 18 | 3 19 | 3 20 | 4 22 | 4 23 | 3 24 \
    3 25 | 4 25
    
    3 1 | 4 2 | 3 6 | 4 9 | 4 13 | 3 14 | 4 15 | 3 16 | 3 17 | 4 17

4.7.11 type TIME_EVENTS
=======================
The syntactic form of data of type {time interval between two observation dates, number of events occurring between these two observation dates} consists in giving, in a first column, the time interval between two observation dates (length of the observation period), in a second column, the number of events occurring between these two observation dates and, in a third column, the corresponding frequency. The time interval between two observation dates should be given in increasing order and then, for each possible time interval, the number of events should be given in increasing order. This is equivalent of giving successively the frequency distribution of the number of events for each possible time interval between two observation dates, ranked in increasing order.

::

    # frequency distribution of the number of events for an observation period of length 20
    20  2   1
    20  3   2
    20  4   4
    20  5   12
    20  6   14
    20  7   6
    20  8   2
    20  9   1

::

    #frequency distribution of the number of events for an observation period of length 30
    30  3   1
    30  5   2
    30  6   4
    30  7   12
    30  8   14
    30  9   6
    30  10  2
    30  12  1

4.7.12 type TOPS
================
Consider the following example::

    2 VARIABLES

    VARIABLE 1 : POSITION
    VARIABLE 2 : NB_INTERNODE
    
    10 5 | 12 5 | 13 6 | 13 8 | 15 7 | 20 10 | 22 11 | 23 11 | 27 15 | 30 16 | 31 15 | 32 17 \
    35 16 | 37 18 | 40 19 | 45
    
    5 2 | 7 4 | 10 5 | 11 6 | 15 7 | 18 8 | 20 9 | 21 11 | 22 11 | 25 12 | 25

The syntactic form of the type TOPS is a variant of the syntactic form of the type SEQUENCES. 'Tops' can be seen as sequences with an explicit index parameter of type POSITION. This index parameter represents the position of successive offspring shoots along the parent shoot and a final index value gives the number of internodes of the parent shoot. The second variable of type NB_INTERNODE gives the number of internodes of the offspring shoots.

The explicit index parameter of type POSITION can be replaced by inter-position intervals::

    2 VARIABLES
    
    VARIABLE 1 : POSITION_INTERVAL
    VARIABLE 2 : NB_INTERNODE
    
    10 5 | 2 5 | 1 6 | 0 8 | 2 7 | 5 10 | 2 11 | 1 11 | 4 15 | 3 16 | 1 15 | 1 17 \
    3 16 | 2 18 | 3 19 | 5
    
    5 2 | 2 4 | 3 5 | 1 6 | 4 7 | 3 8 | 2 9 | 1 11 | 1 11 | 3 12 | 0

4.7.13 type TOP_PARAMETERS
==========================
A model of 'tops' is defined by three parameters, namely the growth probability of the parent shoot, the growth probability of the offspring shoots (both in the sense of Bernoulli processes) and the growth rhythm ratio offspring shoots / parent shoot. Consider the following example::

    TOP_PARAMETERS
    
    PROBABILITY : 0.7
    AXILLARY_PROBABILITY : 0.6
    RHYTHM_RATIO : 0.8
    
The following constraints apply to the parameters:
    
* 0.05 |leq| PROBABILITY |leq| 1
* 0.05 |leq| AXILLARY_PROBABILITY |leq| 1
* 1/3 |leq| RHYTHM_RATIO |leq| 3

4.7.14 type VECTOR_DISTANCE
===========================
The parameters of definition of a distance between vectors are the number of variables, the distance type (ABSOLUTE_VALUE or QUADRATIC) if there is more than one variable, the variable types (NUMERIC, SYMBOLIC, ORDINAL or CIRCULAR), and eventually the weights of the variables (default behaviour: the variables have the same weight), and in the symbolic case, explicit distances between symbols (default behaviour: 0 / 1 for mismatch / match). Consider the following example:

.. literalinclude:: syntax_vector_distance.dat

you can load this data by using

.. doctest::

    >>> vector_distance = VectorDistance('doc/user/syntax_vector_distance.dat')

4.7.15 type VECTORS
===================

.. warning:: Check this example and data files which can do be loaded.

In the syntactic form of the type VECTORS, each row corresponds to an individual and each column corresponds to a variable. Consider the following example:

.. literalinclude:: syntax_vectors.dat

you can load this data by using

.. doctest::
    :options: +SKIP

    >>> vector = Vectors('doc/user/syntax_vectors.dat')

.. here below are aliases that won't appear in the output documents


