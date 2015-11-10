# 

"""
#########################################################################
#
#  apricot tree follow-up
#
#  dates : 21-mars 24-mars 28-mars 31-mars 03-avr 07-avr 10-avr 14-avr 18-avr
#          09-mai 12-mai 16-mai 19-mai 23-mai 26-mai 30-mai 02-juin.
#          = 0  3  7 10 13 17 20 24 28 49 52 56 59 63 66 70 73
#
#  Data : Evelyne Costes.
#
#########################################################################
"""

from openalea.sequence_analysis import *
from openalea.sequence_analysis import get_shared_data


seq1 = Sequences(get_shared_data("abricotier_suivi_11.seq"))
ipe = IndexParameterExtract(seq1, 3, MaxIndex=24)

Plot(seq1, ViewPoint="Data")
Plot(ipe, ViewPoint="Data")
    
Display(seq1, ViewPoint="Data", Format="Line")

timev3_7 = TimeScaling(TimeEvents(seq1, 3, 7),  12)
hist1 = ExtractHistogram(timev3_7, "NbEvent", 48)
hist2 = ExtractHistogram(TimeEvents(get_shared_data("abri13.ren")), "NbEvent", 16)

Plot(hist1, hist2)
Plot(timev3_7)

# "Equilibrium" or "Ordinary"
mytype = "Equilibrium"
mytype= "Ordinary"
renew3_7 = Estimate(timev3_7, mytype)
renew3_7_s = Estimate(timev3_7, mytype, NbIteration=10000, Estimator="PenalizedLikelihood", Weight=0.6)
renew3_7_np = Estimate(timev3_7, mytype, NbIteration=10000)
renew3_7_p = Estimate(timev3_7, mytype, Estimator="Parametric")


Plot(renew3_7_p)
Plot(ExtractDistribution(renew3_7, "InterEvent"), 
        ExtractDistribution(renew3_7_s, "InterEvent"),
        ExtractDistribution(renew3_7_np, "InterEvent"),
        ExtractDistribution(renew3_7_p, "InterEvent"), 
        Distribution("NB", 1, 1, 1. / 11.0191))
Plot(ExtractDistribution(renew3_7, "NbEvent", 48),
        ExtractDistribution(renew3_7_s, "NbEvent", 48),
        ExtractDistribution(renew3_7_np, "NbEvent", 48),
        ExtractDistribution(renew3_7_p, "NbEvent", 48))

"""
#########################################################################
#
#  coffee tree follow-up
#
#  Data : Christian Cilas.
#
#########################################################################
"""
seq11 = Sequences(get_shared_data("cafe_ortho1.seq"))
seq12 = Sequences(get_shared_data("cafe_ortho2.seq"))
seq13 = Sequences(get_shared_data("cafe_ortho3.seq"))
seq14 = Sequences(get_shared_data("cafe_ortho4.seq"))
seq15 = Sequences(get_shared_data("cafe_ortho5.seq"))
seq16 = Sequences(get_shared_data("cafe_ortho6.seq"))

seq10 = Merge(seq11, seq12, seq13, seq14, seq15, seq16)


timev1_40 = RenewalData(seq10, 0, 39)
timev41_80 = RenewalData(seq10, 40, 79)
timev1_80 = RenewalData(seq10, 0, 79)
timev81_120 = RenewalData(seq10, 80, 119)
timev121_167 = RenewalData(seq10, 120, 166)
timev81_167 = RenewalData(seq10, 80, 166)
timev1_167 = RenewalData(seq10, 0, 166)


Plot(timev1_40)

renew1_40_i = Estimate(timev1_40, InitialInterEvent=Distribution("U", 1, 38), NbIteration=10000)
renew1_40_is = Estimate(timev1_40, InitialInterEvent=Distribution("U", 1, 50), NbIteration=10000, Estimator="PenalizedLikelihood", Weight=0.4)


Plot(renew1_40_is)

Display(renew1_40_is)
renew1_40_c = Estimate(timev1_40, "Equilibrium", NbIteration=10000)

Plot(ExtractDistribution(renew1_40_i, "InterEvent"), ExtractDistribution(renew1_40_is, "InterEvent"), ExtractDistribution(renew1_40_c, "InterEvent"))
Plot(Fit(Merge(Shift(ExtractHistogram(timev1_40, "Backward"), 1), ExtractHistogram(timev1_40, "Forward")), ExtractDistribution(renew1_40_is,  "Forward")))


timev10 = Simulate(renew1_40_c, "Equilibrium", 200, 40)
timev11 = Simulate(renew1_40_c, "Equilibrium", Simulate(Mixture(0.5, Distribution("U", 40, 40), 0.5, Distribution("U", 50, 50)), 200))

Plot(timev11)

Display(timev11)

timev14 = TimeSelect(timev11, 40)
timev15 = NbEventSelect(timev11, 0, 2)
timev16 = Merge(timev10, timev11)


