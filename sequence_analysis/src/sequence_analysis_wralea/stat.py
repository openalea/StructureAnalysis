###############################################################################
# -*- python -*-
#
#       amlPy function implementation
#
#       Copyright or (C) or Copr. 2006 INRIA - CIRAD - INRA  
#
#       File author(s): Christophe Pradal <christophe.prada@cirad.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
# 
#       OpenAlea WebSite : http://openalea.gforge.inria.fr
#
###############################################################################

__doc__="""
amlPy functions
"""

__license__= "Cecill-C"
__revision__=" $Id: stat.py 9410 2010-08-11 16:04:02Z cokelaer $ "

#//////////////////////////////////////////////////////////////////////////////

from openalea.core import *
from openalea.stat_tool import *
from openalea.sequence_analysis import *

import types

# !!! Do not plot in nosetests !!!
import sys

DISABLE_PLOT = False
if("nosetests" in sys.argv[0]):
    DISABLE_PLOT = True




#//////////////////////////////////////////////////////////////////////////////
# Adapters
def adapt2list(arg):
    if arg is None:
        return []
    elif type(arg) == types.ListType:
        return arg
    else:
        return [arg]

#//////////////////////////////////////////////////////////////////////////////
# Input/output functions
#//////////////////////////////////////////////////////////////////////////////
@add_docstring(Fit)
def py_fit(histo, dist):
    return (Fit(histo, dist),)


@add_docstring(Compound)
def py_compound_ascii( filename ):
    if filename and filename != '':
        return (Compound(filename),)

@add_docstring(Convolution)
def py_convolution_ascii( filename ):
    if filename and filename != '':
        return (Convolution(filename),)

@add_docstring(Distribution)
def py_distribution_ascii( filename ):
    if filename and filename != '':
        return (Distribution(filename),)

@add_docstring(Histogram)
def py_histogram_ascii(filename):
    if filename and filename != '':
        return (Histogram(filename),)


#//////////////////////////////////////////////////////////////////////////////

@add_docstring(Compound)
def py_compound( sum_dist, dist ):
    if sum_dist and dist:
        return (Compound(sum_dist, dist),)

#//////////////////////////////////////////////////////////////////////////////

@add_docstring(Convolution)
def py_convolution( list_of_dist ):
    l=adapt2list(list_of_dist)
    if l:
        return Convolution(*l)

@add_docstring(Mixture)
def py_mixture( list_of_dist ):
    l=adapt2list(list_of_dist)
    if l:
        arg=[]
        for w,d in l:
            arg.append(w)
            arg.append(d)
    return Mixture(*arg)

#//////////////////////////////////////////////////////////////////////////////

@add_docstring(Binomial)
def py_dist_binomial( inf_bound= 0, sup_bound = 10, proba = 0.5 ):
    return (Distribution("BINOMIAL",inf_bound, sup_bound, proba),)

#//////////////////////////////////////////////////////////////////////////////

@add_docstring(Poisson)
def py_dist_poisson( inf_bound= 0, param = 10 ):
    """test"""
    return (Distribution("P",inf_bound, param),)

#//////////////////////////////////////////////////////////////////////////////

@add_docstring(Distribution)
def py_dist_negativebinomial( inf_bound= 0, param = 10., proba = 0.5 ):
    return (Distribution("NB",inf_bound, param, proba),)

#//////////////////////////////////////////////////////////////////////////////

def py_dist_uniform( inf_bound= 0, sup_bound = 10 ):
    return (Distribution("U",inf_bound, sup_bound),)

#//////////////////////////////////////////////////////////////////////////////

def py_hiddenmarkov( filename, length=20 ):

    if filename and filename != '':
        return (HiddenMarkov(filename,length),)

#//////////////////////////////////////////////////////////////////////////////

def py_hiddensemimarkov( filename, length=20, counting=True ):

    if filename and filename != '':
        return (HiddenSemiMarkov(filename,length,counting),)

#//////////////////////////////////////////////////////////////////////////////

def py_markov( filename, length=20 ):

    if filename and filename != '':
        return (Markov(filename,length),)

#//////////////////////////////////////////////////////////////////////////////

def py_mixture( filename ):

    if filename and filename != '':
        return (Mixture(filename),)

#//////////////////////////////////////////////////////////////////////////////

class PyRenewal(Node):
    def __init__(self):
        Node.__init__(self)
        self.Types=["Equilibriun","Ordinary"]
        self.add_input(name="inter_event")
        self.add_input(name="Type", interface = IEnumStr(self.Types), value = self.Types[0])
        self.add_input(name="ObservationTime", interface=IInt, value=20)

    def __call__(self, inputs):
        inter_event= self.get_input("inter_event")
        Type=self.get_input("Type")
        obs_time=self.get_input("ObservationTime")
        return (Renewal(inter_event,Type,obs_time),)

class PyRenewalAscii(Node):
    def __init__(self):
        Node.__init__(self)
        self.Types=["Equilibriun","Ordinary"]
        self.add_input(name="filename", interface=IFileStr)
        self.add_input(name="Type", interface = IEnumStr(self.Types), value = self.Types[0])
        self.add_input(name="ObservationTime", interface=IInt, value=20)

    def __call__(self, inputs):
        filename= self.get_input("filename")
        Type=self.get_input("Type")
        obs_time=self.get_input("ObservationTime")
        if filename:
            return (Renewal(filename,Type,obs_time),)

#//////////////////////////////////////////////////////////////////////////////

def py_semimarkov( filename, length=20, counting=True ):

    if filename and filename != '':
        return (SemiMarkov(filename,length,counting),)

#//////////////////////////////////////////////////////////////////////////////

@add_docstring(Histogram)
def py_histogram(seq):
    return (Histogram(seq),)

#//////////////////////////////////////////////////////////////////////////////

@add_docstring(Plot)
def py_plot( obj, fig_id=1, title='', viewpoint="v", **kargs):
    if DISABLE_PLOT:
        return
    if obj:
        if type(obj) == types.ListType:
            status = Plot(*obj, FigureId=fig_id, ViewPoint=viewpoint, **kargs)
        else:
            Plot(obj, FigureId=fig_id, ViewPoint=viewpoint, **kargs)

def py_plot_segprofile(seq, ind, nb_seg, model, output):
    if DISABLE_PLOT:
        return
    if seq:
        args= [seq, ind, nb_seg]
        model = model.split(',')
        model = map(lambda x:x.strip(),model)
        args.extend(model)
        kwds={'ViewPoint':'SegmentProfile', 'Output':output}
        Plot(*args, **kwds)



#//////////////////////////////////////////////////////////////////////////////

class PyObjectFromFile(Node):
    """ 
    Construct a model (Distribution, vectors, ...) from a file.
    """

    Types= {"BINOMIAL" : Distribution,
            "POISSON" : Distribution,
            "NEGATIVE_BINOMIAL" : Distribution,
            "UNIFORM" : Distribution,
            "COMPOUND" : Compound,
            "CONVOLUTION" : Convolution,
            "MIXTURE" : Mixture,
            "HISTOGRAM" : Histogram,
            "SEQUENCES" : Sequences,
            "TIME_EVENTS" : TimeEvents,
            "TOPS" : Tops,
            "VECTOR_DISTANCE" : VectorDistance,
            "VECTORS" : Vectors,
            
            }
    def __init__(self):
        Node.__init__(self)
        self.add_input(name="filename", interface=IFileStr)
        self.add_input(name="Type", 
                       interface = IEnumStr(self.Types.keys()), 
                       value = "BINOMIAL")
        self.add_output(name='Model')
        self.set_caption("import model")

    def __call__(self, inputs):
        filename= self.get_input("filename")
        name=self.get_input("Type")
        self.set_caption(name)

        klass = self.Types.get(name)
        if filename:
            return (klass(filename),)

#//////////////////////////////////////////////////////////////////////////////

@add_docstring(Sequences)
def py_sequences(seq=[], identifiers=[], indexParameter='Position'):
    if seq:
        if identifiers:
            return (Sequences(seq,identifier,indexParameter),)
        else:
            return(Sequences(seq),)

#//////////////////////////////////////////////////////////////////////////////
# Compare family : TODO

@add_docstring(Compare)
def py_compare_frequency( histos, Type ):
    _types = {"Numeric" : "N","Ordinal" : "O", "Symbolic" : "S"}
    _type=_types[Type]
    if histos is None:
        return None,
    args = list(adapt2list(histos))
    args.append(_type)
    
    return Compare(*args),
        


@add_docstring(Compare)
def py_compare_vectors( vectors, vector_distance ):
    if vectors is None or vector_distance is None:
        return None,
    return Compare(vectors, vector_distance),

#//////////////////////////////////////////////////////////////////////////////

@add_docstring(ComparisonTest)
def py_comparisontest(Type, histo1, histo2):
    if Type in ['F', 'T', 'W'] and histo1 and histo2:
        return ComparisonTest(Type, histo1, histo2)

#//////////////////////////////////////////////////////////////////////////////

@add_docstring(Estimate)
def py_estimate_dist( histo,
                      distribution,
                      MinInfBound,
                      InfBoundStatus,
                      ):
    if not histo: return
    if distribution == 'NON-PARAMETRIC':
        return Estimate(histo, distribution)
    else:
        return Estimate(histo, 
                        distribution, 
                        MinInfBound=MinInfBound, 
                        InfBoundStatus=InfBoundStatus)


@add_docstring(Estimate)
def py_estimate_mixture( histo,
                         components,
                         MinInfBound,
                         InfBoundStatus,
                         DistInfBoundStatus,
                         NbComponents,
                         Penalty,):

    if not histo: return
    args = [histo]
    args.append("MIXTURE")

    mix = components.split(',')
    mix = map(lambda x:x.strip(),mix)
    args.extend(mix)

    kwds = {}
    kwds['MinInfBound'] = MinInfBound
    kwds['InfBoundStatus'] = InfBoundStatus
    kwds['DistInfBoundStatus'] = DistInfBoundStatus
    kwds['NbComponent'] = NbComponents

    if NbComponents == "Estimated":
        kwds['Penalty'] = Penalty

    return Estimate(*args, **kwds)

@add_docstring(Estimate)
def py_estimate_conv( histo,
                      dist,
                      Estimator,
                      InitialDistribution,
                      MinInfBound,
                      NbIteration,
                      Penalty,
                      Weight,
                      Outside,
                    ):

    if not histo or not dist: return
    
    kwds = {}
    kwds['Estimator']=Estimator
    kwds['MinInfBound']=MinInfBound
    kwds['NbIteration']=NbIteration
    if InitialDistribution:
        kwds['InitialDistribution'] = InitialDistribution

    if Estimator == "PenalizedLikelihood":
        kwds['Penalty'] = Penalty
        kwds['Weight'] = Weight
        kwds['Outside'] = Outside
    
    return Estimate(histo,"CONVOLUTION", dist,*kwds)

@add_docstring(Estimate)
def py_estimate_compound( histo,
                          dist,
                          type,
                          Estimator,
                          InitialDistribution,
                          MinInfBound,
                          NbIteration,
                          Penalty,
                          Weight,
                          Outside,
                          ):

    if not histo or not dist: return
    
    kwds = {}
    kwds['Estimator']=Estimator
    kwds['MinInfBound']=MinInfBound
    kwds['NbIteration']=NbIteration
    if InitialDistribution:
        kwds['InitialDistribution'] = InitialDistribution

    if Estimator == "PenalizedLikelihood":
        kwds['Penalty'] = Penalty
        kwds['Weight'] = Weight
        kwds['Outside'] = Outside
    
    return Estimate(histo, "COMPOUND", dist, type, *kwds)



"""
def py_estimate_distrib( histo,
                      distribution,
                      mixtures,
                      unknow,
                      MinInfBound,
                      InfBoundStatus,
                      DistInfBoundStatus,
                      NbComponents,
                      Penalty,
                      Parametric,
                      InitialDistribution ):
    if not histo: return
    args = [histo]
    args.append(distribution)

    if distribution == 'MIXTURE':
        mix = mixtures.split(',')
        mix = map(lambda x:x.strip(),mix)
        args.extend(mix)

    if distribution in ['CONVOLUTION','COMPOUND'] and unknow != '': 
        args.append(unknow)

    kwds = {}
    if MinInfBound != 0:
        kwds['MinInfBound'] = MinInfBound

    if distribution != 'NON-PARAMETRIC':
        if InfBoundStatus != 'Free':
            kwds['InfBoundStatus'] = InfBoundStatus

    if distribution == 'MIXTURE':
        if DistInfBoundStatus != 'Fixed':
            kwds['DistInfBoundStatus'] = DistInfBoundStatus
        if NbComponents == "Estimated":
            kwds['NbComponent'] = NbComponents
            kwds['Penalty'] = Penalty

    if distribution in ['CONVOLUTION','COMPOUND']:
        kwds['Parametric'] = Parametric
        if InitialDistribution:
            kwds['InitialDistribution'] = InitialDistribution
            if 'MinInfBound' in kwds:
                del kwds['MinInfBound']

    return Estimate(*args, **kwds)
"""

@add_docstring(Merge)
def py_merge( data ):
    
    if data:
        if type(data) == types.ListType:
            return Merge(*data)
        else:
            return Merge(data)

@add_docstring(ExtractData)
def py_extractdata( model ):
    if model is not None:
        return ExtractData(model)

@add_docstring(ExtractHistogram)
def py_extract_histogram( model, choice):
    if model is not None:
        return ExtractHistogram(model, choice)

@add_docstring(ExtractHistogram)
def py_extract_histogram_sequences( model, input_data, choice, variable):
    if model is not None:
        if input_data=='Sequences':
            try:
                return ExtractHistogram(model, choice, variable)
            except:
                return ExtractHistogram(model, choice)
        elif input_data=='Vectors':
            return ExtractHistogram(model, variable)



@add_docstring(ExtractDistribution)
def py_extract_distribution( model, key, index):
    if model is not None:
        try:
            return ExtractDistribution(model, key, index)
        except:
            return ExtractDistribution(model, key)

@add_docstring(Shift)
def py_shift(obj, param=0):
    if obj: return Shift(obj,param)

@add_docstring(Shift)
def py_shiftn(obj, variable=0, param=0):
    if obj: return Shift(obj,variable, param)

@add_docstring(Cluster)
def py_cluster(obj, mode, variable_rank, step, information_ratio, limits):
    if obj is None: return

    if mode == "Step":
        if variable_rank:
            return Cluster(obj, mode, variable_rank, step)
        else:
            return Cluster(obj, mode, step)
    elif mode == 'Information':
        return Cluster(obj, mode, information_ratio)
    elif mode == "Limit":
        _limits = adapt2list(limits)
        if variable_rank:
            return Cluster(obj, mode, variable_rank, _limits)
        else:
            return Cluster(obj, mode, _limits)

@add_docstring(Simulate)
def py_simulate_dist(obj, size=100):
    if obj is None:
        return
    return Simulate(obj,size)

@add_docstring(SelectVariable)
def py_select_variable(obj, variables, mode):
    if not obj:
        return
    vars = adapt2list(variables)
    return SelectVariable(obj, vars, Mode=mode),

@add_docstring(ValueSelect)
def py_value_select(obj, variables, value):
    if not obj:
        return
    return ValueSelect(obj, variables, value),

@add_docstring(SelectIndividual)
def py_select_individual(obj, individuals, mode):
    if not obj:
        return
    vars = adapt2list(individuals)
    return SelectIndividual(obj, vars, Mode=mode),

@add_docstring(Segmentation)
def py_segmentation(seq, ind, nb_segment, change_points,  model_list, model, NbSegment, Output):
    if not seq:
        return
    
    kwds={}
    args= [seq, ind]
    if change_points:
        args.append(change_points)
    else:
        args.append(nb_segment)
    model = model.split(',')
    model = map(lambda x:x.strip(),model)
    args.extend(model)

    if change_points:
        kwds['Output']=Output
        return Segmentation(*args, **kwds),
    else:
        kwds['NbSegment']=NbSegment
        if NbSegment != 'Estimated':
            kwds['Output']=Output
        return Segmentation(*args, **kwds),


@add_docstring(Segmentation)
def py_segmentation_sample(seq, nb_segment, model_list, model, Output):
    if not seq:
        return
    kwds={}
    args= [seq, nb_segment]

    model = model.split(',')
    model = map(lambda x:x.strip(),model)
    args.extend(model)

    kwds['Output']=Output
    return Segmentation(*args, **kwds)

@add_docstring(ComputeCorrelation)
def py_compute_correlation(seq,  MaxLag,  Type,  Normalization):
    if not seq:
        return 
    return ComputeCorrelation(seq,  MaxLag=MaxLag,  Type=Type,  Normalization=Normalization)

@add_docstring(ComputeCorrelation)
def py_compute_correlation_mult(seq,  auto,  var1,  var2,  MaxLag,  Type,  Normalization):
    if not seq:
        return 
    if auto:
        return ComputeCorrelation(seq,  var1,  MaxLag=MaxLag,  Type=Type,  Normalization=Normalization)
    else:
        return ComputeCorrelation(seq,  var1, var2,  MaxLag=MaxLag,  Type=Type,  Normalization=Normalization)

@add_docstring(PointwiseAverage)
def py_pointwise_average(seq,  StandardDeviation,  Output,  dirname,  FileName,  Format):
    if not seq:
        return
    
    kwds={}
    kwds['StandardDeviation'] = StandardDeviation
    kwds['Output'] = Output
    if dirname and FileName:
        FileName = os.path.join(dirname, FileName)
        kwds[FileName] = FileName
        kwds['Format'] = Format
        
    return PointwiseAverage(seq,  **kwds), 

@add_docstring(Vectors)
def py_vectors(seq, IndexVariable = False):
    if not seq:
        return
    return Vectors(seq, IndexVariable=IndexVariable),


@add_docstring(Regression)
def py_regression(vec, regressionModel, explanatoryVariable, responseVariable, filter, frequencies, distribution, span, Algorithm, Weighting ):
    if vec is None or explanatoryVariable is None or responseVariable is None:
        return

    if regressionModel == 'Linear':
        reg = Regression(vec, regressionModel, explanatoryVariable, responseVariable)
    elif regressionModel == 'MovingAverage':
        if filter:
            reg = Regression(vec, regressionModel, explanatoryVariable, responseVariable, filter, Algorithm=Algorithm)
        elif frequencies:
            reg = Regression(vec, regressionModel, explanatoryVariable, responseVariable, frequencies, Algorithm=Algorithm)
        elif distribution:
            reg = Regression(vec, regressionModel, explanatoryVariable, responseVariable, distribution, Algorithm=Algorithm)
        else:
            raise "Missing filter parameters: Please, give us filter, frequencies or distribution to estimate a MovingAverage regression."
    else:
        reg = Regression(vec, regressionModel, explanatoryVariable, responseVariable, span, Weighting=Weighting)
    return reg,

@add_docstring(Cumulate)
def py_cumulate(seq):
    if seq is not None:
        return Cumulate(seq),

@add_docstring(ToHistogram)
def py_to_histogram(obj):
    if obj is not None:
        return ToHistogram(obj),

@add_docstring(ToDistribution)
def py_to_distribution(obj):
    if obj is not None:
        return ToDistribution(obj),

@add_docstring(ExtractData)
def py_extract_data(obj):
    if obj is not None:
        return ExtractData(obj),
