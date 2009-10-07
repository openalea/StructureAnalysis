"""Output functions

:Author: Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>
"""
__version__ = "$Id$"

import plot
import os
import glob
import sys
import error
# Output functions
from _stat_tool import FORWARD_DYNAMIC_PROGRAMMING, FORWARD_BACKWARD_SAMPLING,\
    GENERALIZED_VITERBI, FORWARD_BACKWARD_SAMPLING


# exported in sequence_analysis but hard-coded here to prevent dependency on 
# sequence analysis
NB_SEGMENTATION = 10
NB_STATE_SEQUENCE = 10

def add_doc(function):
    """a simple decorator to replace f's docstring by a new one
    
    The new one is the docstring of the function's name capitalized. 
    E.g: if function's name is display, then 
        display.__doc__ = Display.__doc__
    """
    name = function.__name__
     
    function.__doc__ = eval(name.capitalize()).__doc__
    return function
     
def Display(obj, *args, **kargs):
    """ASCII output of an object of the STAT module 

    ASCII output of sets of sequences or tops (ViewPoint="Data"): the format
    "Column" corresponds to the ASCII file syntax for objects of type sequences 
    or tops. For a given value of the index parameter, the different variables 
    are successively displayed. With the format "Line", the univariate sequence
    for each variable are displayed on consecutive lines. In the case of
    univariate sequences, the two formats give the same output.

    ASCII output of a (frequency) distribution and the associate hazard or
    survival rates (ViewPoint="Survival"): It is assumed that the (frequency)
    distribution represents lifetime and the hazard or survival rates are
    deduced from this lifetime distribution.

    ASCII output of the state profile given by the smoothed probabilities 
    :math:`P(S_t=j|X_0^{\tau-1}=x_0^\tau)` as a
    function of the index parameter `t` computed from the parameters of a hidden
    Markovian model for the sequence :math:`x_0^\tau` (ViewPoint="StateProfile").

    :Parameters:
      * `obj` - object to display,
      * `vec` (`_Vectors`), 
      * `seq` (`_Sequences`, `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData`, `_Tops`),
      * `dist` (`_Distribution`, `_MixtureDist`, `_Convolution`, `_Compound`),
      * `histo` (`_Histogram`, `_MixtureData`, `_ConvolutionData`, `_CompoundData`),
      * `hmc` (`_HiddenMarkov`),
      * `hsmc` (`_HiddenSemiMarkov`),
      * `identifier` (int) - identifier of a sequence.

    :Keywords:
      * ViewPoint (string): point of view on the object ("Survival" or "Data"
        or "StateProfile"). This optional argument can be set at
          * "Data" only if the first argument is of type `_Vectors`,
            `_Sequences`, `_DiscreteSequences`, `_MarkovData`,
            `_SemiMarkovData` or `_Tops`,
          * "Survival" only if the first argument is of type `_Distribution`,
            `_MixtureDist`, `_Convolution`, `_Compound`, `_Histogram`,
            `_MixtureData`, `_ConvolutionData` or `_CompoundData`
          * "StateProfile" only if the first argument is of type `_HiddenMarkov` 
            or `_HiddenSemiMarkov`.
      * Detail (int): level of detail: 1 (default value) or 2. 
        This optional argument cannot be used if the optional argument 
        ViewPoint is set at "Survival" or "StateProfile".
      * Format (string): format of sequences (only relevant for multivariate 
        sequences): "Column" (default value) or "Line". This optional argument
        can only be used if the  optional argument ViewPoint is set at "Data",
        and hence, if the first argument is of type `_Vectors`, `_Sequences`,
        `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData` or `_Tops`. 

    :Returns:
        A string

    :Examples:

    .. doctest::
        :options: +SKIP 
    
        >>> from openalea.stat_tool.output import Display
        >>> Display(obj, Detail=2)
        >>> Display(vec, ViewPoint="Data", Detail=2)
        >>> Display(seq, ViewPoint="Data", Format="Line", Detail=2)
        >>> Display(dist, ViewPoint="Survival")
        >>> Display(histo, ViewPoint="Survival")
        >>> Display(hmc, identifier, ViewPoint="StateProfile")
        >>> Display(hsmc, identifier, ViewPoint="StateProfile") 

    .. seealso::
        :func:`~openalea.stat_tool.output.Plot`,
        :func:`~openalea.stat_tool.output.Save`.

    """
    return obj.display(*args, **kargs)
    


def Plot(obj, *args, **kargs):
    """     
    Graphical output of an object of the STAT module using the GNUPLOT software.

    In the case of Markovian models or sequences, the graphical outputs are 
    grouped as follows:
      * "SelfTransition": add outgoing server thunderbirdself-transition probability as a function of the
        index parameter (non-homogeneous Markov chain),
      * "Observation": observation distributions attached to each state of the
        underlying (semi-)Markov chain (lumped processes or hidden Markovian
        processes),
      * "Intensity": (empirical) probabilities of states/outputs as a function
        of the index parameter,
      * "FirstOccurrence": (frequency) distributions of the time-up to the first
        occurrence of a state/output (or first-passage time in a state/output
        distributions),
      * "Recurrence" (frequency) distributions of the recurrence time in a
        state/output,
      * "Sojourn": (frequency) distributions of the sojourn time in a 
        state/output (or state/output occupancy distributions). For the 
        frequency distributions extracted from sequences, the sojourn times in
        the last visited states which are considered as censored are isolated.
      * "Counting": counting (frequency) distributions (either distributions of
        the number of runs (or clumps) of a state/output per sequence or
        distributions of the number of occurrences of a state/output per
        sequence).

    Graphical output of a (frequency) distribution and the associate hazard or
    survival rates (ViewPoint="Survival"): It is assumed that the (frequency)
    distribution represents lifetime and the hazard or survival rates are
    deduced from this lifetime distribution.
    Graphical output of the state profile given by the smoothed probabilities as
    a function of the index parameter t computed from the parameters of a hidden
    Markovian model for the sequence (ViewPoint="StateProfile"). 

    :Parameters:
      * obj1 (`_Distribution`, `_Mixture`, `_Convolution`, `_Compound`,
        `_DistributionData`, `_MixtureData`, `_ConvolutionData`,
        `_CompoundData`,`_Renewal`, `_TimeEvents`, `_RenewalData`,
        `_Sequences`, `_DistanceMatrix`, ` _TopParameters`, `_Tops`),
      * vec1 (`_Vectors`): values,
      * vecn (`_Vectors`): vectors,
      * variable (int): variable index,
      * obj2: (`_Markov`, `_SemiMarkov`, `_HiddenMarkov`, `_HiddenSemiMarkov`,
        `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData`): Markovian model
        for discrete univariate sequences  or discrete univariate sequences,
      * obj3: (`_Markov`, `_SemiMarkov`, `_HiddenMarkov`, `_HiddenSemiMarkov`,
        `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData`): Markovian model
        for discrete multivariate sequences or  discrete multivariate sequences,
      * type (string): type of graphical outputs in the case of Markovian models
        or sequences: "SelfTransition", "Observation", "Intensity",
        "FirstOccurrence", "Recurrence", "Sojourn" or "Counting",
      * dist1, dist2, ... (`_Distribution`, `_Mixture`, `_Convolution`, `_Compound`),
      * histo1, histo2, ... (`_DistributionData`, `_MixtureData`, `_ConvolutionData`,
        `_CompoundData`),
      * seq (`_Sequences`, `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData`,
        `_Tops`),
      * dist (`_Distribution`, `_Mixture`, `_Convolution`, `_Compound`),
      * histo (`_DistributionData`, `_MixtureData`, `_ConvolutionData`, 
        `_CompoundData`),
      * hmc (_HiddenMarkov),
      * hsmc (_HiddenSemiMarkov),
      * identifier (int): identifier of a sequence.

    :Keywords:
    
      * ViewPoint (string): point of view on the object ("Data" or "Survival" 
        or "StateProfile"). This optional argument can be set at :
          * "Data" only if the first mandatory argument is of type sequences, 
            discrete_sequences, markov_data, semi-markov_data or tops, 
          * "Survival" only if the first mandatory argument is of type distribution, 
            mixture, convolution, compound, histogram, mixture_data, 
            convolution_data or compound_data 
          * "StateProfile" only if the first mandatory argument is of type 
            hidden_markov or hidden_semi-markov.
      * Title (string): graphic title (the default: no title). 

    :Returns:
        Nothing.
        
    :Examples:

    .. doctest::
        :options: +SKIP

        >>> from openalea.stat_tool.output import Display
        >>> Plot(obj1, Title="Distribution")
        >>> Plot(vec1, Title="Values")
        >>> Plot(vecn, variable, Title="Vectors")
        >>> Plot(variable)
        >>> Plot(obj2, type, Title="Sequences")
        >>> Plot(type)
        >>> Plot(obj3, type, variable, Title="Multivariate sequences")
        >>> Plot(type, variable)
        >>> Plot(dist1, dist2,..., Title="Family of distributions")
        >>> Plot(histo1, histo2,..., Title="Family of frequency distributions")
        >>> Plot(seq, ViewPoint="Data")
        >>> Plot(dist, ViewPoint="Survival", Title="Survival rates")
        >>> Plot(histo, ViewPoint="Survival", Title="Survival rates")
        >>> Plot(hsmc, identifier, ViewPoint="StateProfile", Title="Smoothed probabilities") 

    .. seealso::
        :func:`~openalea.stat_tool.output.Display`,
        :func:`~openalea.stat_tool.output.Save`.
    """
    
    
    return obj.plot(*args, **kargs)

def Save(obj, *args, **kargs):
    """     
    Saving of an object of the STAT module in a file.

    Saving of sets of sequences or 'tops' (ViewPoint="Data"): the format "Column" 
    corresponds to the ASCII file syntax for objects of type _Sequences or _Tops. 
    For a given value of the index parameter, the different variables are
    successively written. With the format "Line", the univariate sequence for
    each variable are written on consecutive lines. In the case of univariate
    sequences, the two formats give the same file.

    Saving of a (frequency) distribution and the associate hazard or survival rates 
    (ViewPoint="Survival"): It is assumed that the (frequency) distribution represents 
    lifetime and the hazard or survival rates are deduced from this lifetime
    distribution.

    Saving of the state profile given by the smoothed probabilities as a function of 
    the index parameter t computed from the parameters of a hidden Markovian model 
    for the sequence (ViewPoint="StateProfile"). 

    .. note:: The persistence mechanism is implemented by the Save function. 
    
    
    :Parameters:
    
      * obj: object of the STAT module (except objects of type vector_distance),
      * file_name (string),
      * histo (_Histogram, _MixtureData, _ConvolutionData, _CompoundData),
      * vec (_Vectors),
      * timev (_TimeEvents, _RenewalData),
      * seq (_Sequences, _DiscreteSequences, _MarkovData, _SemiMarkovData, _Tops).
      * dist (_Distribution, _Mixture, _Convolution, _Compound),
      * hmc (_HiddenMarkov),
      * hsmc (_HiddenSemiMarkov). 

    :Keywords:
    
      * ViewPoint (string): point of view on the object ("Data" or "Survival" or "StateProfile").
        This optional argument can be set at :
          * "Data" only if the first argument is of type `_Vectors`, `_Sequences`, 
            `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData` or `_Tops`, 
          * "Survival" only if the first argument is of type `_Distribution`,
            `_Mixture`, `_Convolution`, `_Compound`, `_Histogram`, `_MixtureData`, 
            `_ConvolutionData` or `_CompoundData` 
          * "StateProfile" only if the first argument is of type `_HiddenMarkov or 
            `_HiddenSemiMarkov`.
      * Detail (int): level of detail: 1 (default value) or 2. 
        This optional argument can only be used if the optional argument ViewPoint 
        is not set, or if the optional argument ViewPoint is set at "Data" and 
        if the first mandatory argument is of type `_Vectors`, `_Sequences`, 
        `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData` or `_Tops`.
      * Format (string): file format: "ASCII" (default format), "Binary" or "SpreadSheet". 
        These file formats cannot be specified if the optional argument ViewPoint
        is set at "Data". The optional argument Format can only be set at "Binary"
        if the optional argument ViewPoint is not set.
      * Format (string): format of sequences (only relevant for multivariate sequences): 
        "Column" (default value) or "Line". This optional argument can only be used if the 
        optional argument ViewPoint is set at "Data", and hence, if the first argument is of 
        type `_Vectors`, `_Sequences`, `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData` 
        or `_Tops`. 
      * Sequence (int): identifier of a sequence. This optional argument can only be used 
        if the optional argument ViewPoint is set at "StateProfile", and hence, if the first 
        mandatory argument is of type `_HiddenMarkov` or `_HiddenSemiMarkov`. 
  
    :Returns:
        No object returned.
        
    :Examples:

    .. doctest::
        :options: +SKIP

        >>> Save(obj, file_name, Format="ASCII", Detail=2)
        >>> Save(histo, file_name, ViewPoint="Data")
        >>> Save(vec, file_name, ViewPoint="Data", Detail=2)
        >>> Save(timev, file_name, ViewPoint="Data")
        >>> Save(seq, file_name, ViewPoint="Data", Format="Line", Detail=2)
        >>> Save(dist, file_name, ViewPoint="Survival", Format="SpreadSheet")
        >>> Save(histo, file_name, ViewPoint="Survival", Format="SpreadSheet")
        >>> Save(hmc, ViewPoint="StateProfile", Sequence=1, Format="SpreadSheet")
        >>> Save(hsmc, ViewPoint="StateProfile", Sequence=1, Format="SpreadSheet") 

    .. seealso::
        :func:`~openalea.stat_tool.output.Display`,
        :func:`~openalea.stat_tool.output.Plot`.
        
    .. todo::
        In the statInterface, Format is used for ViewPoint=="Data" need to be
        clarified
    """

    return obj.save(*args, **kargs)



class StatInterface(object):
    """ Abstract base class for stat_tool objects """
    
    
   

    def old_plot(self, *args, **kargs):
        """ Old AML style plot """
        #todo: to be replace by correct enumerate but depends on sequence_analysis
        output_type = { 
               "ChangePoint" : 0,
               "Segment" : 1
               }
        title = kargs.get("Title", "")
        ViewPoint = kargs.get("ViewPoint", "")
        suffix = kargs.get("Suffix", "")
        params = kargs.get("Params", ())
        output = kargs.get("Output", 0)

        data = bool(ViewPoint.lower() == "data")
        survival = bool(ViewPoint.lower() == "survival")
        stateprofile = bool(ViewPoint.lower() == "stateprofile")
        segmentprofile = bool(ViewPoint.lower() == "segmentprofile")

        
        import tempfile
        prefix = tempfile.mktemp()
        
        if(data):
            try:
                self.plot_data_write(prefix, title)
            except AttributeError:
                raise AttributeError("%s has not 'data' viewpoint"
                                     % (str(type(self))))
        elif(survival):
            try:
                self.survival_plot_write(prefix, title)
            except AttributeError:
                raise AttributeError("%s has not 'survival' viewpoint"
                                     % (str(type(self))))

        elif(stateprofile):
            try:
                self.state_profile_plot_write(prefix, title, *params)
            except AttributeError:
                raise AttributeError("%s has not 'state_profile' viewpoint"
                                     % (str(type(self))))
        elif (segmentprofile):
            try:
                
                error.CheckType([args[0],args[1]], 
                                [int,int])
                if len(args)==2:
                    error.CheckType([args[2]], [[list, str]])
                    models = []
                    for model in args[2]:
                        models.append(model_type[args[2]])
                else:
                    models = [3] #Gaussian todo: check this is correct
                output = output_type[output]
                self.segment_profile_write(prefix, args[0], args[1], models,
                                           output, title)
            except AttributeError:
                raise AttributeError("%s has not 'segment_profile' viewpoint"
                                     % (str(type(self))))
        elif(args):
            self.plot_write(prefix, title, list(args))
        else:            
            self.plot_write(prefix, title)
            

        plot_file = prefix + suffix + ".plot"

        f = open(plot_file, "a")
        f.write("pause -1")
        f.close()
        if("win32" in sys.platform):
            # replace file separators
            f = open(plot_file, "r")
            ct = f.read()
            f.close()
            ctrp = ct.replace('\\', '\\\\')
            ctrp = ctrp.replace(',\\\\', ',\\')
            f = open(plot_file, "w")
            f.write(ctrp)
            f.close()
            print plot_file, "\n"
		
        try:
            import Gnuplot
            command = Gnuplot.GnuplotOpts.gnuplot_command
        except ImportError:
            if("win32" in sys.platform):
                command = "pgnuplot.exe"
            else:
                command = "gnuplot"

        if(not plot.DISABLE_PLOT):
            os.system("%s %s"%(command, plot_file))
        
        #  for f in glob.glob(prefix+"*"):
        #     os.remove(f)

    def plot_print(self, *args, **kargs):
        """ Old AML style print into .ps file """

        title = kargs.get("Title", "")
        ViewPoint = kargs.get("ViewPoint", "")
        suffix = kargs.get("Suffix", "")
        params = kargs.get("Params", ())

        survival = bool(ViewPoint.lower() == "survival")
        stateprofile = bool(ViewPoint.lower() == "stateprofile")
        
        import tempfile
        prefix = tempfile.mktemp()

        if(survival):
            try:
                self.survival_plot_write(prefix, title)
            except AttributeError:
                raise AttributeError("%s has not 'survival' viewpoint"
                                     % (str(type(self))))

        elif(stateprofile):
            try:
                self.state_profile_plot_write(prefix, title, *params)
            except AttributeError:
                raise AttributeError("%s has not 'state_profile' viewpoint"
                                     % (str(type(self))))


        elif(args):
            self.plot_write(prefix, title, list(args))
        else:
            self.plot_write(prefix, title)

        plot_file = prefix + suffix + ".print"
        print "Graph printed into file:",  prefix + suffix + ".ps"
        f = open(plot_file, "r")
        f.readline()
        contents = f.read()
        f.close()
        f = open(plot_file, "w+")
        f.write("set terminal postscript color \n")
        f.write(contents)
        f.close()
        f = open(plot_file, "a")
        f.write("pause -1")
        f.close()            
        if("win32" in sys.platform):
            # replace file separators
            f = open(plot_file, "r")
            ct = f.read()
            f.close()
            ctrp = ct.replace('\\', '\\\\')
            ctrp = ctrp.replace(',\\\\', ',\\')
            f = open(plot_file, "w")
            f.write(ctrp)
            f.close()
            print plot_file, "\n"

        try:
            import Gnuplot
            command = Gnuplot.GnuplotOpts.gnuplot_command
        except ImportError:
            if("win32" in sys.platform):
                command = "pgnuplot.exe"
            else:
                command = "gnuplot"

        if(not plot.DISABLE_PLOT):
            os.system("%s %s"%(command, plot_file))
        
        for f in glob.glob(prefix+"*"):
            if f != prefix + suffix + ".ps":
                os.remove(f)
    @add_doc
    def plot(self, *args, **kargs):
        
        title = kargs.get("Title", "")
        ViewPoint = kargs.get("ViewPoint", "")
        params = kargs.get("Params", ())
        groups = kargs.get("Groups", ())

        survival = bool(ViewPoint.lower() == "survival")
        stateprofile = bool(ViewPoint.lower() == "stateprofile")
        segmentprofile = bool(ViewPoint.lower() == "segmentprofile")
        data = bool(ViewPoint.lower() == "data")
        
        try:
            if (survival):
                plotable = self.survival_get_plotable(*params)

            elif (stateprofile):
                plotable = self.stateprofile_get_plotable(*params)

            elif (segmentprofile):
                plotable = self.segmentprofile_get_plotable(*params)
            else:
                if (args):
                    if len(args)==1 and type(args[0])==int:
                        plotable = self.get_plotable_list()
                    else:
                        plotable = self.get_plotable_list(list(args), *params)
                else:
                    plotable = self.get_plotable(*params)
            
            plotter = plot.get_plotter()
        except:
            import warnings
            warnings.warn("Cannot use new plotter. Use old style plot.")
            plotable = None
    
        if plot.DISABLE_PLOT:
            return

        if(plotable is not None):
            plotter.plot(plotable, title, groups, *args, **kargs)
        else:
            self.old_plot(*args, **kargs)
            
    @add_doc
    def display(self, Detail=1, ViewPoint='v', Format='c', **kargs):
        
        format_map = {'c':'c', 'l':'l', "Column": 'c', "Line":'l'}
        viewpoint_map = {'v':'v', "Data":"d", "Survival":'s',
                         "SegmentProfile":'q', "StateProfile":'p'}
          
        segmentations_map = {
                "DynamicProgramming": FORWARD_DYNAMIC_PROGRAMMING,
                "ForwardBackwardSampling": FORWARD_BACKWARD_SAMPLING
                }
        
        state_seq_map = { "GeneralizedViterbi": GENERALIZED_VITERBI,
                         "ForwardBackwardSampling": FORWARD_BACKWARD_SAMPLING
                         }
            
        
        # Detail level
        exhaustive = False
        if Detail == 2:
            exhaustive = True
        elif Detail != 1:
            raise ValueError("Detail can be either 1 or 2. You gave %s" 
                             % Detail)
        
        Format = error.CheckDictKeys(Format, format_map)
        ViewPoint = error.CheckDictKeys(ViewPoint, viewpoint_map)
        
        # those lines are not used yet. Depend on sequence_analysis
         
        #Segmentations = error.ParseKargs(kargs, "Segmentations", 
        #                                 "DynamicProgramming",
        #                                 segmentations_map ) 
        #NbSegmentation = kargs.get("NbSegmentation", NB_SEGMENTATION)
        
        #state_sequence = error.ParseKargs(kargs, "StateSequences", 
        #                                  "GeneralizedViterbi", state_seq_map)
        
        #nb_state_sequence = kargs.get("NbStateSequence", NB_STATE_SEQUENCE)
        # ViewPoint

        
        # Survival
        if ViewPoint == 's':
            try:
                output = self.survival_ascii_write()
            except AttributeError:
                raise AttributeError("""
                %s has not 'survival' viewpoint""" % (str(type(self))))

        # Data
        elif ViewPoint == "d":
            try:
                # try with format argument
                output = self.ascii_data_write(exhaustive, Format)

            except Exception, e:
                try:
                    output = self.ascii_data_write(exhaustive)

                except AttributeError:
                    raise AttributeError("""
                        %s has not 'data' viewpoint""" % (str(type(self))))

        # StatProfile
        elif ViewPoint == "p":
            try:
                output = self.state_profile_ascii_write() # A completer
            except AttributeError:
                raise AttributeError("""
                    %s has not 'stateprofile' viewpoint""" % (str(type(self))))
                
        else:
            output = self.ascii_write(exhaustive)

        return output 
    
    @add_doc
    def save(self, filename, Detail=2, ViewPoint="", Format="ASCII" ):
        
        # Detail level
        if(Detail>1):
            exhaustive = True
        else:
            exhaustive = False
            
        if(Format.lower() == "spreadsheet"):
            self.spreadsheet_write(filename)
            #f = open(filename, 'w')
            #f.write(outstr)
            #f.close()

        elif(Format.lower() == "data"):
            self.file_ascii_data_write(filename, exhaustive)        
        else:
            self.file_ascii_write(filename, exhaustive)
        
