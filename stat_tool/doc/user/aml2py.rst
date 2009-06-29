AML syntax versus pythonic syntax
=================================

This is a simple example to illustrate that the old AML syntax can still be used, as well as the new object-python syntax

.. literalinclude:: aml2py.py

The two Merge calls return the same output::


    histogram - sample size: 424
    mean: 13.1179   variance: 34.6858   standard deviation: 5.88947
    coefficient of skewness: -0.142591   coefficient of kurtosis: -0.997803
    mean absolute deviation: 4.99081   coefficient of concentration: 0.255743
    information: -1284.16 (-3.02869)


