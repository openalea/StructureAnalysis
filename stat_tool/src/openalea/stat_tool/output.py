__doc__ = """ Output functions """
__docformat__ = " restructuredtext "



import plot
import os
import glob

# Output functions


def Display(obj, *args, **kargs):
    """     
    ASCII output of an object of the STAT module 

    Usage
    -----
     * ``Display(obj, Detail=2)``
     * ``Display(vec, ViewPoint="Data", Detail=2)``
     * ``Display(seq, ViewPoint="Data", Format="Line", Detail=2)``

     * ``Display(dist, ViewPoint="Survival")``
     * ``Display(histo, ViewPoint="Survival")``

     * ``Display(hmc, identifier, ViewPoint="StateProfile")``
     * ``Display(hsmc, identifier, ViewPoint="StateProfile")`` 

    Parameters
    ----------
     * obj: object to display,
     * vec (`_Vectors`), 
     * seq (`_Sequences`, `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData`, `_Tops`),

     * dist (`_Distribution`, `_MixtureDist`, `_Convolution`, `_Compound`),
     * histo (`_Histogram`, `_MixtureData`, `_ConvolutionData`, `_CompoundData`),

     * hmc (`_HiddenMarkov`),
     * hsmc (`_HiddenSemiMarkov`),
     * identifier (int): identifier of a sequence. 

    Keywords
    --------

      * ViewPoint (string): point of view on the object ("Survival" or "Data" or "StateProfile"). \
        This optional argument can be set at : 
        * "Data" only if the first argument is of type `_Vectors`, `_Sequences`, \
              `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData` or `_Tops`, 
        * "Survival" only if the first argument is of type `_Distribution`, `_MixtureDist`, \
              `_Convolution`, `_Compound`, `_Histogram`, `_MixtureData`, \
              `_ConvolutionData` or `_CompoundData`
        * "StateProfile" only if the first argument is of type `_HiddenMarkov` or \
              `_HiddenSemiMarkov`.

      * Detail (int): level of detail: 1 (default value) or 2. \
      This optional argument cannot be used if the optional argument ViewPoint is set at \
      "Survival" or "StateProfile".
    
      * Format (string): format of sequences (only relevant for multivariate sequences): \
      "Column" (default value) or "Line". This optional argument can only be used if the \
      optional argument ViewPoint is set at "Data", and hence, if the first argument is of\ 
      type `_Vectors`, `_Sequences`, `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData`\ 
      or `_Tops`. 

    Description
    -----------
    ASCII output of sets of sequences or tops (ViewPoint="Data"): the format "Column" 
    corresponds to the ASCII file syntax for objects of type sequences or tops. For a given 
    value of the index parameter, the different variables are successively displayed. 
    With the format "Line", the univariate sequence for each variable are displayed on 
    consecutive lines. In the case of univariate sequences, the two formats give the same 
    output.

    ASCII output of a (frequency) distribution and the associate hazard or survival rates 
    (ViewPoint="Survival"): It is assumed that the (frequency) distribution represents 
    lifetime and the hazard or survival rates are deduced from this lifetime distribution.

    ASCII output of the state profile given by the smoothed probabilities as a function of 
    the index parameter t computed from the parameters of a hidden Markovian model for the 
    sequence (ViewPoint="StateProfile"). 

    Return
    ------
    A string

    See Also
    --------
    `Plot`, `Save`.
    """

    return obj.display(*args, **kargs)


def Plot(obj, *args, **kargs):
    """     
    Graphical output of an object of the STAT module using the GNUPLOT software.

    Usage
    -----
      * Plot(obj1, Title="Distribution")
      * Plot(vec1, Title="Values")
      * Plot(vecn, variable, Title="Vectors")
      * Plot(variable)

      * Plot(obj2, type, Title="Sequences")
      * Plot(type)
      * Plot(obj3, type, variable, Title="Multivariate sequences")
      * Plot(type, variable)

      * Plot(dist1, dist2,..., Title="Family of distributions")
      * Plot(histo1, histo2,..., Title="Family of frequency distributions")

      * Plot(seq, ViewPoint="Data")

      * Plot(dist, ViewPoint="Survival", Title="Survival rates")
      * Plot(histo, ViewPoint="Survival", Title="Survival rates")

      * Plot(hmc, identifier, ViewPoint="StateProfile",
      * Title="Smoothed probabilities")
      * Plot(hsmc, identifier, ViewPoint="StateProfile", Title="Smoothed probabilities") 

    Paramters
    ---------
      * obj1 (`_Distribution`, `_Mixture`, `_Convolution`, `_Compound`,\
      `_DistributionData`, `_MixtureData`, `_ConvolutionData`, `_CompoundData`, \
      `_Renewal`, `_TimeEvents`, `_RenewalData`, `_Sequences`, `_DistanceMatrix`,\
      ` _TopParameters`, `_Tops`),
      * vec1 (`_Vectors`): values,
      * vecn (`_Vectors`): vectors,
      * variable (int): variable index,

      * obj2: (`_Markov`, `_SemiMarkov`, `_HiddenMarkov`, `_HiddenSemiMarkov`, \
      `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData`): Markovian model \
      for discrete univariate sequences  or discrete univariate sequences,
      * obj3: (`_Markov`, `_SemiMarkov`, `_HiddenMarkov`, `_HiddenSemiMarkov`,\
      `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData`): Markovian model \
      for discrete multivariate sequences or  discrete multivariate sequences,
      * type (string): type of graphical outputs in the case of Markovian models or \
      sequences: "SelfTransition", "Observation", "Intensity", "FirstOccurrence", \
      "Recurrence", "Sojourn" or "Counting",

      * dist1, dist2, ... (`_Distribution`, `_Mixture`, `_Convolution`, `_Compound`),
      * histo1, histo2, ... (`_DistributionData`, `_MixtureData`, `_ConvolutionData`, \
      `_CompoundData`),

      * seq (`_Sequences`, `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData`, `_Tops`),

      * dist (`_Distribution`, `_Mixture`, `_Convolution`, `_Compound`),
      * histo (`_DistributionData`, `_MixtureData`, `_ConvolutionData`, `_CompoundData`),

      * hmc (_HiddenMarkov),
      * hsmc (_HiddenSemiMarkov),
      * identifier (int): identifier of a sequence. 

    Keywords
    --------
      * ViewPoint (string): point of view on the object ("Data" or "Survival" or "StateProfile"). \
      This optional argument can be set at :

        * "Data" only if the first mandatory argument is of type sequences, 
          discrete_sequences, markov_data, semi-markov_data or tops, 
        * "Survival" only if the first mandatory argument is of type distribution, 
          mixture, convolution, compound, histogram, mixture_data, convolution_data 
          or compound_data 
        * "StateProfile" only if the first mandatory argument is of type hidden_markov 
          or hidden_semi-markov.

      * Title (string): graphic title (the default: no title). 

    Return
    ------
      Nothing.

    Description
    -----------
      In the case of Markovian models or sequences, the graphical outputs are grouped as follows:
        * "SelfTransition": self-transition probability as a function of the index parameter (non-homogeneous Markov chain),
        * "Observation": observation distributions attached to each state of the underlying (semi-)Markov chain (lumped processes or hidden Markovian processes),
        * "Intensity": (empirical) probabilities of states/outputs as a function of the index parameter,
        * "FirstOccurrence": (frequency) distributions of the time-up to the first occurrence of a state/output (or first-passage time in a state/output distributions),
        * "Recurrence" (frequency) distributions of the recurrence time in a state/output,
        * "Sojourn": (frequency) distributions of the sojourn time in a state/output (or state/output occupancy distributions). For the frequency distributions extracted from sequences, the sojourn times in the last visited states which are considered as censored are isolated.
        * "Counting": counting (frequency) distributions (either distributions of the number of runs (or clumps) of a state/output per sequence or distributions of the number of occurrences of a state/output per sequence).

    Background
    ----------
    Graphical output of a (frequency) distribution and the associate hazard or survival 
    rates (ViewPoint="Survival"): It is assumed that the (frequency) distribution 
    represents lifetime and the hazard or survival rates are deduced from this 
    lifetime distribution.
    Graphical output of the state profile given by the smoothed probabilities as a function 
    of the index parameter t computed from the parameters of a hidden Markovian model for 
    the sequence (ViewPoint="StateProfile"). 

    See Also
    --------
    `Display`, `Save`
    """

    return obj.plot(*args, **kargs)


def Save(obj, *args, **kargs):
    """     
    Saving of an object of the STAT module in a file.

    Usage
    -----
      * ``Save(obj, file_name, Format="ASCII", Detail=2)``

      * ``Save(histo, file_name, ViewPoint="Data")``
      * ``Save(vec, file_name, ViewPoint="Data", Detail=2)``
      * ``Save(timev, file_name, ViewPoint="Data")``
      * ``Save(seq, file_name, ViewPoint="Data", Format="Line", Detail=2)``

      * ``Save(dist, file_name, ViewPoint="Survival", Format="SpreadSheet")``
      * ``Save(histo, file_name, ViewPoint="Survival", Format="SpreadSheet")``
    
      * ``Save(hmc, ViewPoint="StateProfile", Sequence=1, Format="SpreadSheet")``
      * ``Save(hsmc, ViewPoint="StateProfile", Sequence=1, Format="SpreadSheet")`` 

    Parameters
    ----------
      * obj: object of the STAT module (except objects of type vector_distance),
      * file_name (string),

      * histo (_Histogram, _MixtureData, _ConvolutionData, _CompoundData),
      * vec (_Vectors),
      * timev (_TimeEvents, _RenewalData),
      * seq (_Sequences, _DiscreteSequences, _MarkovData, _SemiMarkovData, _Tops).

      * dist (_Distribution, _Mixture, _Convolution, _Compound),

      * hmc (_HiddenMarkov),
      * hsmc (_HiddenSemiMarkov). 

    Keywords
    --------

      * ViewPoint (string): point of view on the object ("Data" or "Survival" or "StateProfile").
      This optional argument can be set at :
         * "Data" only if the first argument is of type `_Vectors`, `_Sequences`, \
              `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData` or `_Tops`, 
         * "Survival" only if the first argument is of type `_Distribution`, `_Mixture`, \
              `_Convolution`, `_Compound`, `_Histogram`, `_MixtureData`, \
              `_ConvolutionData or `_CompoundData 
         * "StateProfile" only if the first argument is of type `_HiddenMarkov or \
              `_HiddenSemiMarkov`.

      * Detail (int): level of detail: 1 (default value) or 2. 
      This optional argument can only be used if the optional argument ViewPoint 
      is not set, or if the optional argument ViewPoint is set at "Data" and 
      if the first mandatory argument is of type `_Vectors`, `_Sequences`, 
      `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData` or `_Tops`.

      * Format (string): file format: "ASCII" (default format), "Binary" or "SpreadSheet". 
      These file formats cannot be specified if the optional argument ViewPoint is set at 
      "Data". The optional argument Format can only be set at "Binary" if the optional 
      argument ViewPoint is not set.
    
      * Format (string): format of sequences (only relevant for multivariate sequences): 
      "Column" (default value) or "Line". This optional argument can only be used if the 
      optional argument ViewPoint is set at "Data", and hence, if the first argument is of 
      type `_Vectors`, `_Sequences`, `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData` 
      or `_Tops`. 

      * Sequence (int): identifier of a sequence. This optional argument can only be used 
      if the optional argument ViewPoint is set at "StateProfile", and hence, if the first 
      mandatory argument is of type `_HiddenMarkov` or `_HiddenSemiMarkov`. 

    
    Description
    -----------
    Saving of sets of sequences or 'tops' (ViewPoint="Data"): the format "Column" 
    corresponds to the ASCII file syntax for objects of type _Sequences or _Tops. 
    For a given value of the index parameter, the different variables are successively 
    written. With the format "Line", the univariate sequence for each variable 
    are written on consecutive lines. In the case of univariate sequences, 
    the two formats give the same file.

    Saving of a (frequency) distribution and the associate hazard or survival rates 
    (ViewPoint="Survival"): It is assumed that the (frequency) distribution represents 
    lifetime and the hazard or survival rates are deduced from this lifetime distribution.

    Saving of the state profile given by the smoothed probabilities as a function of 
    the index parameter t computed from the parameters of a hidden Markovian model 
    for the sequence (ViewPoint="StateProfile"). 

    Background
    ----------
    The persistence mechanism is implemented by the Save function. 

    Return
    ------
    No object returned.

    See Also
    --------

    `Display`, `Plot`
    """

    return obj.save(*args, **kargs)



class StatInterface(object):
    """ Abstract base class for stat_tool objects """

    def old_plot(self, title="", *args):
        """ Old AML style plot """

        import tempfile
        prefix = tempfile.mktemp()
        self.plot_write(prefix, title)
        
        plot_file = prefix + ".plot"

        # Add an infinite pause in the command file
        f = open(plot_file, 'a')
        f.write("pause -1")
        f.close()
        
        # call gnuplot
        os.system("gnuplot %s"%(plot_file))
        
        for f in glob.glob(prefix+"*"):
            os.remove(f)


    def plot(self, title="", *args):
        __doc__ = Plot.__doc__


        try:
            plotable = self.get_plotable()
            plotter = plot.get_plotter()

        except AttributeError:
            import warnings
            warnings.warn("Use old style plot.")
            
            plotable = None

        except ImportError:
            import warnings
            warnings.warn("No Plotter available. Use old style plot.")
            plotable = None

        
        if(plotable):
            plotter.plot(plotable, title, *args)
        else:
            self.old_plot(title, *args)
            

    def display(self, Detail=1, ViewPoint="", Format=""):
        __doc__ = Plot.__doc__

        # Detail level
        if(Detail>1):
            exaustive = True
        else:
            exaustive = False

        # ViewPoint

        # Survival
        if(ViewPoint.lower() == "survival"):
            try:
                return self.survival_ascii_write()
            except AttributeError:
                raise AttributeError("%s has not 'survival' viewpoint"%(str(type(self))))

        # Data
        elif(ViewPoint.lower() == "data"):
            try:
                # try with format argument
                return self.ascii_data_write(exaustive, format)

            except Exception, e:
                try:
                    return self.ascii_data_write(exaustive)

                except AttributeError:
                    raise AttributeError("%s has not 'data' viewpoint"%(str(type(self))))
                

        # StatProfile
        elif(ViewPoint.lower() == "stateprofile"):
            try:
                return self.state_profile_ascii_write() # A completer
            except AttributeError:
                raise AttributeError("%s has not 'stateprofile' viewpoint"%(str(type(self))))
                
        else:
            return self.ascii_write(exaustive)



    def save(self, filename, Detail=2, ViewPoint="", Format="ASCII" ):
        __doc__ = Save.__doc__

        if(Format.lower() == "spreadsheet"):
            outstr = self.spreadsheet_write(filename)

        else:
            outstr = self.display(Detail, ViewPoint, Format)

            f = open(filename, 'w')
            f.write(outstr)
            f.close()
        

################################################################################

class Test:

    
    def get_mixture(self):
        
        from distribution import Uniform
        from mixture import Mixture
 
        d1 = Uniform(0, 10)
        d2 = Uniform(10, 20)
        d3 = Uniform(20, 30)

        m = Mixture(0.1, d1, 0.2, d2, 0.7, d3)
        
        return m

    def get_mixture_2(self):

        from histogram import Histogram 

        h = Histogram("../../../test/meri2.his")
        m = h.estimate_mixture(["B", "NB"])
        return m


    def test_old_plot(self):
        m = self.get_mixture()
        #m.old_plot()


    def test_plot_mixture_1(self):
        m = self.get_mixture()
        m.plot()

    def test_plot_mixture_2(self):
        m = self.get_mixture_2()
        m.plot()

    def test_plot_mixture_data(self):
        from mixture import Mixture
        from simulate import Simulate
        from distribution import Distribution

        mixt1 = Mixture(0.6, Distribution("B", 2, 18, 0.5), 0.4, 
                        Distribution("NB", 10, 10, 0.5))
        mixt_histo1 = Simulate(mixt1, 200)

        mixt1.plot()
        mixt_histo1.plot()


    def test_plot_convolution(self):
        
        from convolution import Convolution
        from simulate import Simulate
        from estimate import Estimate
        from data_transform import Shift
        from histogram import Histogram
        
        convol1 = Convolution("../../../test/convolution1.conv")
        convol1.plot()

        histo_b2 = Histogram("../../../test/nothofagus_antarctica_bud_2.his")
        histo_s2 = Histogram("../../../test/nothofagus_antarctica_shoot_2.his")
        
        convol31 = Estimate(Shift(histo_s2, 1), "CONVOLUTION", Estimate(histo_b2, "NP"), 
                            NbIteration=100, Estimator="PenalizedLikelihood", Weight=0.5)
        
        convol31.plot()


    def test_plot_convolution_data(self):
        
        from convolution import Convolution
        from simulate import Simulate
        
        convol1 = Convolution("../../../test/convolution1.conv")
        convol_histo1 = Simulate(convol1, 200)
        convol_histo1.plot()

    
    def test_plot_distribution_set(self):
        
        from distribution import Distribution 
        d1 = Distribution("B", 2, 18, 0.5) 
        d2 = Distribution("NB", 10, 10, 0.5)
        d3 = Distribution("U", 10, 20)
        
        Plot(d1, d2, d3)


