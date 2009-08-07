"""common enumerate

:Author: Thomas Cokelaer <Thomas.Cokelaer@inria.fr>

"""
__version__ = "$Id: enumerate.py 6708 2009-08-05 14:20:27Z cokelaer $"

import _sequence_analysis


type_dict = {"Pearson": 0,
             "Spearman": 1,
             "Kendall": 2,
             "Spearman2": 3}


norm_type = {"Approximated": 0,
                 "Exact": 1, }

type_dict2 = {
        "Length"            : _sequence_analysis.LENGTH,
        "NbRun"             : _sequence_analysis.NB_RUN,
        "NbOccurrence"      : _sequence_analysis.NB_OCCURRENCE,
        "FirstOccurrence"   : _sequence_analysis.FIRST_OCCURRENCE,
        "Mean"              : _sequence_analysis.SEQUENCE_MEAN,
        "Cumul"             : _sequence_analysis.SEQUENCE_CUMUL
        }    