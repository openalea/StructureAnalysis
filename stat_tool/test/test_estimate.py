"""estimate tests

:Author: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""

__version__ = "$Id$"

from pathlib import Path

from openalea.stat_tool._stat_tool import CUMUL_THRESHOLD

from openalea.stat_tool import _stat_tool, distribution
from openalea.stat_tool.compound import Compound
from openalea.stat_tool.data_transform import ExtractDistribution, Shift
from openalea.stat_tool.distribution import Binomial
from openalea.stat_tool.enums import compound_type, distribution_identifier_type
from openalea.stat_tool.estimate import Estimate, likelihood_penalty_type
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.simulate import Simulate

from openalea.stat_tool.distribution import set_seed

from importlib.resources import files, as_file

import pytest

@pytest.fixture
def path():
    return Path(__file__).parent

@pytest.fixture
def data_path():
    # Data directory within module
    return files('openalea.stat_tool.data')

@pytest.fixture
def meri1(data_path):
    _meri1 = Histogram(str((data_path / "meri1.his")))
    return _meri1

@pytest.fixture
def meri5(data_path):
    _meri5 = Histogram(str((data_path / "meri5.his")))
    return _meri5

@pytest.fixture
def peup2(data_path):
    _peup2 = Histogram(str((data_path /  "peup2.his")))
    return _peup2

@pytest.fixture
def compound1(path):
    _compound1 = Compound(str((path / "data" / "compound1.cd")))
    return _compound1

def test_nonparametric(meri1):
    epsilon = 1.0 - CUMUL_THRESHOLD
    h = meri1
    e = h.estimate_nonparametric()
    # find likelihood value
    strh = str(h)
    word = "information:"
    i = strh.find(word)
    j = strh.find("(", i)
    val = float(strh[i + len(word) : j])
    assert e
    assert abs(e.likelihood() - val) < epsilon

def test_nb(peup2):
    """negativebinomial"""
    h = peup2
    assert h.estimate_parametric("NB")

def test_binomial(meri5):
    """binomial distribution"""
    h = meri5
    assert h.estimate_parametric("B")

def test_poisson():
    """poisson distribution
    >>> p = distribution.poisson(0, 10)
    >>> h = p.simulate(1000)
    """
    set_seed(0)
    p = distribution.Poisson(0, 10)
    h = p.simulate(1000)

    assert h.estimate_parametric("P")

def test_binomial_estimation():
    """estimate binomial distribution"""
    from openalea.stat_tool.enums import distribution_identifier_type as dist_type

    p = distribution.Binomial(2, 12, 0.7)
    set_seed(0)
    h = p.simulate(1000)
    d1 = h.estimate_parametric("B")
    d2 = h.default_parametric_estimation(dist_type["B"])
    assert d1
    assert d2

def test_uniform_estimation():
    """Estimate Uniform distribution"""
    from openalea.stat_tool.enums import distribution_identifier_type as dist_type

    p = distribution.Uniform(2, 12)
    set_seed(0)
    h = p.simulate(1000)
    d1 = h.estimate_parametric("U")
    d2 = h.default_parametric_estimation(dist_type["U"])
    assert d1
    assert d2

def test_Poisson_estimation():
    """Estimate Binomial distribution"""
    from openalea.stat_tool.enums import distribution_identifier_type as dist_type

    p = distribution.Poisson(0, 10)
    set_seed(0)
    h = p.simulate(1000)
    d1 = h.estimate_parametric("P")
    d2 = h.default_parametric_estimation(dist_type["P"])
    assert d1
    assert d2

def test_negative_binomial_estimation():
    """Estimate Negative Binomial distribution"""
    from openalea.stat_tool.enums import distribution_identifier_type as dist_type

    p = distribution.NegativeBinomial(2, 4.5, 0.6)
    set_seed(0)
    h = p.simulate(1000)
    d1 = h.estimate_parametric("NB")
    d2 = h.default_parametric_estimation(dist_type["NB"])
    assert d1
    assert d2

def test_mixture_1(peup2):
    distributions = ["B", "NB", "NB", "NB"]
    h = peup2
    m1 = h.estimate_DiscreteMixture(distributions, NbComponent="Estimated")
    assert m1

    types = []
    for d in distributions:
        temp = distribution_identifier_type[d]
        types.append(temp)

    c = h.discrete_mixture_estimation2(
        types, 0, True, True, likelihood_penalty_type["AIC"]
    )

    # slight variations in quantiles
    # CPL : Core dumped!!!!
    #assert str(c)[0:1336] == str(m1)[0:1336]

def test_mixture_2(peup2):
    h = peup2
    m2 = h.estimate_DiscreteMixture([Binomial(0, 10, 0.5), "NB"])
    assert m2

def test_convolution():
    elementary = Histogram("data/nothofagus_antarctica_bud_2.his")
    total = Histogram("data/nothofagus_antarctica_shoot_2.his")

    convol1 = Estimate(
        Shift(total, 1),
        "CONVOLUTION",
        Estimate(elementary, "NP"),
        NbIteration=100,
        Estimator="PenalizedLikelihood",
        Weight=0.5,
    )

    convol2 = total.shift(1).estimate_convolution(
        elementary.estimate_nonparametric(),
        NbIteration=100,
        Estimator="PenalizedLikelihood",
        Weight=0.5,
    )

    assert convol1 and convol2
    assert convol1 == convol2

def test_compound_two_distribution(compound1):
    """test to be checked"""
    cdist1 = compound1
    set_seed(0)
    chisto1 = Simulate(cdist1, 200)

    cdist2 = Estimate(
        chisto1,
        "COMPOUND",
        ExtractDistribution(cdist1, "Sum"),
        "Sum",
        InitialDistribution=ExtractDistribution(cdist1, "Elementary"),
    )

    # If we call the method directly, we need to provide
    # the default values and perform a conversion.
    # Default is LIKELIHOOD -1 -1.0 SECOND_DIFFERENCE ZERO, which
    #  corresponds to 0, -1,-1,1,0
    # In addition because the type is 's', the 2 distributions
    # must be reversed.

    cdist3 = chisto1.compound_estimation1(
        ExtractDistribution(cdist1, "Elementary"),
        ExtractDistribution(cdist1, "Sum"),
        compound_type["Sum"],
        _stat_tool.LIKELIHOOD,
        -1,
        -1.0,
        _stat_tool.SECOND_DIFFERENCE,
        _stat_tool.ZERO,
    )
    assert str(cdist2) == str(cdist3)

