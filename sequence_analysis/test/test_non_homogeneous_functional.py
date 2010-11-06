#########################################################################
#
#  Floraison des tiges de vanille decrites du sommet vers la base.
#  deux echantillons : apex mort (m) ou decapite (d).
#
#  VARIABLE 1 : non-fleuri (0) / fleuri (1).
#
#########################################################################
from openalea.sequence_analysis import Sequences, Estimate, NonhomogeneousMarkov
from openalea.sequence_analysis import ComputeSelfTransition, Plot
from openalea.sequence_analysis import get_shared_data
seq_m = Sequences(get_shared_data("vanille_m.seq"))
ComputeSelfTransition(seq_m)
Plot(seq_m, "SelfTransition")
Plot(seq_m, "Intensity")
Plot(seq_m, "FirstOccurrence")

mc_m = Estimate(seq_m, "NONHOMOGENEOUS_MARKOV", "MONOMOLECULAR", "VOID")
Plot(mc_m, "SelfTransition")
Plot(mc_m, "Intensity")

