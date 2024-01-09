#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Output functions

.. topic:: output.py summary

    A module dedicated to Output functions (plot, display, save)

    :Code status: mature
    :Documentation status: to be completed
    :Authors:
        * Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>
        * Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: output.py 15183 2013-11-06 10:35:50Z jbdurand $

"""
__version__ = "$Id: output.py 15183 2013-11-06 10:35:50Z jbdurand $"

import plot
import os
import glob
import sys
import error
# Output functions
from openalea.stat_tool._stat_tool import FORWARD_DYNAMIC_PROGRAMMING, \
    FORWARD_BACKWARD_SAMPLING, GENERALIZED_VITERBI, \
    FORWARD_BACKWARD_SAMPLING

try:
    from openalea.sequence_analysis.enums import NB_STATE_SEQUENCE, \
        NB_SEGMENTATION
except:
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
    :math:`P(S_t=j|X_0^{\\tau-1}=x_0^\\tau)` as a
    function of the index parameter `t` computed from the parameters of a hidden
    Markovian model for the sequence :math:`x_0^\\tau` (ViewPoint="StateProfile").

    :Parameters:

    * `obj` - object to display,
    * `vec` (`_Vectors`),
    * `seq` (`_Sequences`, `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData`, `_Tops`),
    * `dist` (`_Distribution`, `_MixtureDist`, `_Convolution`, `_Compound`),
    * `histo` (`_FrequencyDistribution`, `_DiscreteMixtureData`, `_ConvolutionData`, `_CompoundData`),
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
        `_MixtureDist`, `_Convolution`, `_Compound`, `_FrequencyDistribution`,
        `_DiscreteMixtureData`, `_ConvolutionData` or `_CompoundData`
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
        `_DiscreteDistributionData`, `_DiscreteMixtureData`, `_ConvolutionData`,
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
      * histo1, histo2, ... (`_DiscreteDistributionData`, `_DiscreteMixtureData`, `_ConvolutionData`,
        `_CompoundData`),
      * seq (`_Sequences`, `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData`,
        `_Tops`),
      * dist (`_Distribution`, `_Mixture`, `_Convolution`, `_Compound`),
      * histo (`_DiscreteDistributionData`, `_DiscreteMixtureData`, `_ConvolutionData`,
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
      * nbcol (int): number of columns in the output figure
      * Show:
      * legend_size: 10
      * legend_nbcol: 2
      * legend_loc: best
      * legend: True/False


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
      * histo (_FrequencyDistribution, _DiscreteMixtureData, _ConvolutionData, _CompoundData),
      * vec (_Vectors),
      * timev (_TimeEvents, _RenewalData),
      * seq (_Sequences, _DiscreteSequences, _MarkovData, _SemiMarkovData, _Tops).
      * dist (_Distribution, _Mixture, _Convolution, _Compound),
      * hmc (_HiddenMarkov),
      * hsmc (_HiddenSemiMarkov).

    :Keywords:

      * ViewPoint (string): point of view on the object ("Data" or "Survival" or "StateProfile").
        This optional argument can be set at :
          * "Data" only if the first argument is of type `_Sequences`,
            `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData` or `_Tops`,
          * "Survival" only if the first argument is of type `_Distribution`,
            `_Mixture`, `_Convolution`, `_Compound`, `_FrequencyDistribution`, `_DiscreteMixtureData`,
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
        type `_Sequences`, `_DiscreteSequences`, `_MarkovData`, `_SemiMarkovData`
        or `_Tops`. If the first argument is of type `_Vectors`, use Format="Data" to actually
        save the data rather than their summary.
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
        >>> Save(vec, file_name, Format="Data")
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
                        try:
                            from openalea.sequence_analysis.enums import model_type
                            models.append(model_type[args[2]])
                        except:
                            pass
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

        Title = kargs.get("Title", "")

        params = kargs.get("Params", ())
        groups = kargs.get("Groups", ())

        possible_modes = {'Blocking':False, 'NonBlocking':True}
        Mode = error.ParseKargs(kargs, 'Mode', 'Blocking', possible=possible_modes)

        viewpoint_map = {'v':'v',
                         "Data":"d",
                         "Survival":'s',
                         "SegmentProfile":'q',
                         "StateProfile":'p'}
        ViewPoint = error.ParseKargs(kargs, "ViewPoint", "v", possible=viewpoint_map)


        #todo: check the compatibilities between options
        """
        if ((output_option) && ((view_point != 'q') ||
       ((args[0].tag() != AMObjType::SEQUENCES) && (args[0].tag()
       != AMObjType::MARKOVIAN_SEQUENCES) &&
        (args[0].tag() != AMObjType::VARIABLE_ORDER_MARKOV_DATA) &&
        (args[0].tag() != AMObjType::SEMI_MARKOV_DATA) &&
        (args[0].tag() != AMObjType::NONHOMOGENEOUS_MARKOV_DATA))) &&
      ((view_point != 'p') || ((args[0].tag() != AMObjType::HIDDEN_SEMI_MARKOV)
      &&
        (args[0].tag() != AMObjType::MARKOVIAN_SEQUENCES) &&
        (args[0].tag() != AMObjType::VARIABLE_ORDER_MARKOV_DATA) &&
        (args[0].tag() != AMObjType::SEMI_MARKOV_DATA) &&
        (args[0].tag() != AMObjType::NONHOMOGENEOUS_MARKOV_DATA)))) {
    status = false;
    genAMLError(ERRORMSG(INCOMPATIBLE_OPTIONS_s) , "Plot");
  }


  if ((config) && (view_point != 'p') && ((args[0].tag() ==
      AMObjType::MARKOVIAN_SEQUENCES) ||
       (args[0].tag() == AMObjType::HIDDEN_VARIABLE_ORDER_MARKOV)
           || (args[0].tag() == AMObjType::HIDDEN_SEMI_MARKOV) ||
       (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA)
        || (args[0].tag() == AMObjType::SEMI_MARKOV_DATA))) {
    variable = args[1].val.i;

    switch (args[0].tag()) {

    case AMObjType::MARKOVIAN_SEQUENCES : {
      seq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
      if ((variable <= seq->get_nb_variable()) && (seq->get_characteristics(variable - 1))) {
        status = false;
        genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Plot");
      }
      break;
    }

    case AMObjType::HIDDEN_VARIABLE_ORDER_MARKOV : {
      hmarkov = (HiddenVariableOrderMarkov*)((STAT_model*)args[0].val.p)->pt;
      if ((variable <= hmarkov->get_nb_output_process()) &&
          (hmarkov->get_nonparametric_process(variable))) {
        status = false;
        genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Plot");
      }
      break;
    }

    case AMObjType::HIDDEN_SEMI_MARKOV : {
      hsmarkov = (HiddenSemiMarkov*)((STAT_model*)args[0].val.p)->pt;
      if ((variable <= hsmarkov->get_nb_output_process()) &&
          (hsmarkov->get_nonparametric_process(variable))) {
        status = false;
        genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Plot");
      }
      break;
    }

    case AMObjType::VARIABLE_ORDER_MARKOV_DATA : {
      seq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
      if ((variable < seq->get_nb_variable()) && (seq->get_characteristics(variable))) {
        status = false;
        genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Plot");
      }
      break;
    }

    case AMObjType::SEMI_MARKOV_DATA : {
      seq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
      if ((variable < seq->get_nb_variable()) && (seq->get_characteristics(variable))) {
        status = false;
        genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Plot");
      }
      break;
    }
    }
  }

        """
        try:
            from openalea.sequence_analysis.enums import output_display
        except:
            from openalea.stat_tool.enums import output_display



        if kargs.get('Output'):
            try:
                Output = None
                Output = error.ParseKargs(kargs, "Output", 'Segment', output_display)
            except:
                print 'warning could not import output_display from sequence_analysis'
        else:
            try:
                from openalea.sequence_analysis.enums import output_display
            except:
                from openalea.stat_tool.enums import output_display
            Output = None

        if Output is None:
            if ViewPoint == 'q':
                Output = output_display['Segment']
            elif ViewPoint == 'p':
                Output = output_display['State']
        elif (ViewPoint == 'q' and Output not in [output_display['ChangePoint'], output_display['Segment']]) \
            or (ViewPoint == 'p'  and Output not in [output_display['State'], output_display['InState'], output_display['OutState']]):
            raise ValueError(" INCOMPATIBLE_OPTIONS between ViewPoint and Output")


        #calling the plot functions from here
        try:
            if ViewPoint=='s':
                from openalea.stat_tool.enums import histogram_types
                from openalea.stat_tool.enums import model_distribution_types
                #todo is *params needed or not?
                if type(self) in model_distribution_types:
                    #equivalent to dist->suvival_plot_write(error, Plot_prefix, title)
                    plotable = self.survival_get_plotable(*params)
                elif type(self) in histogram_types:
                    #equivalent to histo->survival_plot_write(error , Plot_prefix , title)
                    output = self.survival_get_plotable(*params)
                else:
                    raise ValueError("""(%s) has no survival point. Use another
                        Viewpoint or use a first argument in DISTRIBUTION or MIXTURE or
                        CONVOLUTION or COMPOUND or FREQUENCY_DISTRIBUTION or
                        MIXTURE_DATA or CONVOLUTION_DATA or COMPOUND_DATA"""
                        % str(type(self)))

            elif ViewPoint=='p':
                #print 'viewpoint = state-profile'
                Plot_prefix=''
                plotable = None
                from openalea.sequence_analysis._sequence_analysis import \
                    _HiddenVariableOrderMarkov, _HiddenSemiMarkov
                if type(self) == _HiddenVariableOrderMarkov:
                    plotable = self.state_profile_plotable_write(args[0])
                elif type(self) == _HiddenSemiMarkov:
                    if len(args)==0:
                        raise SyntaxError("expect an identifier (Plot(hsmc25, 1, ViewPoint='StateProfile')")
                    elif len(args)==1:
                        identifier = args[0]
                    else:
                        #print 'iiiiiiiiiiiiiii'
                        raise SyntaxError("expect only one identifier Plot(hsmc25, 1, ViewPoint='StateProfile'")
                    plotable = self.state_profile_plotable_write(identifier, Output)
                else:
                    #todo 3 args required
                    from openalea.sequence_analysis._sequence_analysis import _MarkovianSequences, _VariableOrderMarkovData, _SemiMarkovData
                    # , _NonhomogeneousMarkovData
                    assert type(self) in [_MarkovianSequences, _VariableOrderMarkovData,
                                          _SemiMarkovData] #, _NonhomogeneousMarkovData]
                    if type(args[1])==_HiddenVariableOrderMarkov:
                        plotable = args[1].state_profile_plotable_write2(self , args[0]);
                    elif type(args[1])==_HiddenSemiMarkov:
                        plotable = args[1].state_profile_plotable_write2(self , args[0], Output)
                    else:
                        raise TypeError("expect HiddenVariableOrderMarkov or HiddenSemiMarkov")

                if plotable == None:
                    try:
                        plotable = self.stateprofile_get_plotable(*params)
                    except:
                        pass

            elif ViewPoint=='q':
                from openalea.sequence_analysis._sequence_analysis import _Sequences, _MarkovianSequences, _VariableOrderMarkovData, _SemiMarkovData
                if type(self) not in [_Sequences, _MarkovianSequences, _VariableOrderMarkovData, _SemiMarkovData]:
                    raise TypeError('object must be in SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA')

                try:
                    self.nb_variable
                except:
                    raise ValueError("object has no nb_variable. check that it is a sequence")
                nb_variable = self.nb_variable
                assert len(args)>=2
                error.CheckType([args[0], args[1]], [[int], [int]])
                #construct model_type
                from openalea.sequence_analysis.enums import model_type
                types = []
                for i in range(0, nb_variable):
                    error.CheckType([args[i+2]], [str])
                    if i == 0:
                        types.append(model_type[args[i+2]])
                        #Multinomial or Poisson or Ordinal or Gaussian or
                        # Mean or Variance or MeanVariance
                        if args[i+2] in ["Mean", "MeanVariance"]:
                            for j in range(1, nb_variable):
                                types.append(types[i])
                            break
                    else:
                        # Multinomial or Poisson or Ordinal or Gaussian
                        # or Variance
                        types.append(model_type[args[i+2]])
                #seq->segment_profile_plot_write(
                #         error , Plot_prefix , args[1].val.i ,
                #           args[2].val.i , model_type , output , title);

                plotable = self.segment_profile_plotable_write(args[0], args[1],
                                       types, Output)

            #data viewPoint
            elif ViewPoint == 'd':
                from openalea.sequence_analysis._sequence_analysis import _SemiMarkovData, _MarkovianSequences, _Sequences, _Tops
                #, _NonHomogeneousMarkovData, 
                if type(self) in [_SemiMarkovData, _MarkovianSequences, _Sequences,
                                  _Tops]: #, _NonHomogeneousMarkovData, 
                    #status = seq->plot_data_write(error , Plot_prefix , title);
                    plotable = self.get_plotable_data(*params)
            elif ViewPoint == 'v':
                # plot_write(error , Plot_prefix , title);

                if args:
                    #sequence case:
                    #todo: make it looser: observation, intensity INTENSITY?
                    choices = ["SelfTransition", "Observation", "Intensity",
                         "FirstOccurrence", "Recurrence", "Sojourn", "Counting"]
                    if args[0] in choices:
                        multiplotset = self.get_plotable()
                        viewpoints = [x for x in multiplotset.viewpoint]
                        plotable = []
                        try:
                            from openalea.sequence_analysis import enums
                        except:
                            raise ImportError("sequence analysis not installed !!")

                        if len(args)==1:
                            variable = 0
                        elif len(args)==2:
                            variable = args[1]
                        for index, xx in enumerate(viewpoints):
                            if xx==enums.markovian_sequence_type[args[0]]:
                                if multiplotset.variable[index]==variable:
                                    plotable.append(multiplotset[index])
                    elif len(args) == 1 and type(args[0]) == str:
                        raise SyntaxError("first argument must be in %s and second arg (int) may be provided." % choices)
                    elif len(args)==1 and type(args[0])==int:
                        from openalea.stat_tool._stat_tool import _Vectors
                        if type(self)==_Vectors:
                            #Plot(vector, 1)
                            multiplotset = self.get_plotable()
                            viewpoints = [x for x in multiplotset.viewpoint]
                            plotable = []
                            try:
                                from openalea.sequence_analysis import enums
                            except:
                                raise ImportError("sequence analysis not installed !!")
                            plotable = [multiplotset[args[0]]]
                        else:
                            #Plot(hist1, hist2, hist3)
                            plotable = self.get_plotable_list()
                    elif len(args)==1:
                        #e.g., list of histograms
                        plotable = self.get_plotable_list(list(args), *params)
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
            plotter.plot(plotable, Title, groups, *args, **kargs)
        else:
            self.old_plot(*args, **kargs)

    @add_doc
    def display(self, *args, **kargs):
        format_map = {
                      'c':'c',
                      'l':'l',
                      'Column': 'c',
                      'Line': 'l'}
        viewpoint_map = {'v':'v',
                         "Data":"d",
                         "Survival":'s',
                         "SegmentProfile":'q',
                         "StateProfile":'p'}

        segmentations_map = {
                "DynamicProgramming": FORWARD_DYNAMIC_PROGRAMMING,
                "ForwardBackwardSampling": FORWARD_BACKWARD_SAMPLING
                }

        state_seq_map = { "GeneralizedViterbi": GENERALIZED_VITERBI,
                         "ForwardBackwardSampling": FORWARD_BACKWARD_SAMPLING
                         }

        # Detail level
        Detail = error.ParseKargs(kargs, "Detail", 1, [1, 2])
        if Detail == 2:
            exhaustive = True
        else:
            exhaustive = False
        Format = error.ParseKargs(kargs, "Format", "c", format_map)
        ViewPoint = error.ParseKargs(kargs, "ViewPoint", "v", viewpoint_map)
        NbStateSequence = error.ParseKargs(kargs, "NbStateSequence",
                                           NB_STATE_SEQUENCE)

        error.CheckType([NbStateSequence], [[int, float]])
        NbSegmentation = error.ParseKargs(kargs, "NbSegmentation",
                                          NB_SEGMENTATION)
        error.CheckType([NbSegmentation], [[int, float]])


        StateSequence = error.ParseKargs(kargs, "StateSequence",
                                          "GeneralizedViterbi", state_seq_map)
        Segmentation = error.ParseKargs(kargs, "Segmentation",
                                        "DynamicProgramming",
                                        segmentations_map)
        #todo it seems that by default, segmentation = FORWARD_DYNAMIC_PROGRAMMING ,


        # !! in AML, Output is not set y default, i.e. equivalent to None
        # the ParseKargs does not accept None sinc we provide the list of
        # possible keys in output_display (which do not contain None)
        #, so we first need to check the presence of Output in the kargs
        # then, to give a default value!=None. But be aware that tis default
        # value is a dummy variable that is not used.
        try:
            from openalea.sequence_analysis.enums import output_display
        except:
            from openalea.stat_tool.enums import output_display

        if kargs.get('Output'):
            try:
                Output = None
                Output = error.ParseKargs(kargs, "Output", 'Segment', output_display)
            except:
                print 'warning could not import output_display from sequence_analysis'
        else:
            try:
                from openalea.sequence_analysis.enums import output_display
            except:
                from openalea.stat_tool.enums import output_display
            Output = None

        if Output is None:
            if ViewPoint == 'q':
                Output = output_display['Segment']
            elif ViewPoint == 'p':
                Output = output_display['State']
        elif (ViewPoint == 'q' and Output not in [output_display['ChangePoint'], output_display['Segment']]) \
            or (ViewPoint == 'p'  and Output not in [output_display['State'], output_display['InState'], output_display['OutState']]):
            raise ValueError(" INCOMPATIBLE_OPTIONS between ViewPoint and Output")




        #check arguments compatibilities

        if Detail == 2 and ViewPoint not in ['v', 'd']:
            raise ValueError("incompatible options")
        if Format == 'l' and ViewPoint != 'd':
            raise ValueError("incompatible options")
        """if segmentations_option or nb_segmentation_option)  and \
           (view_point!='q' or args[0] not in
            (
                (args[0].tag() != AMObjType::SEQUENCES)
                && (args[0].tag() != AMObjType::MARKOVIAN_SEQUENCES)
                && (args[0].tag() != AMObjType::VARIABLE_ORDER_MARKOV_DATA)
                && (args[0].tag() != AMObjType::SEMI_MARKOV_DATA)
                && (args[0].tag() != AMObjType::NONHOMOGENEOUS_MARKOV_DATA)
            )
        )
  if (((state_sequences_option) || (nb_state_sequence_option)) && ((view_point != 'p') ||
       ((args[0].tag() != AMObjType::HIDDEN_VARIABLE_ORDER_MARKOV) &&
        (args[0].tag() != AMObjType::HIDDEN_SEMI_MARKOV) &&
        (args[0].tag() != AMObjType::MARKOVIAN_SEQUENCES) &&
        (args[0].tag() != AMObjType::VARIABLE_ORDER_MARKOV_DATA) &&
        (args[0].tag() != AMObjType::SEMI_MARKOV_DATA) &&
        (args[0].tag() != AMObjType::NONHOMOGENEOUS_MARKOV_DATA)))) {
    status = false;
    genAMLError(ERRORMSG(INCOMPATIBLE_OPTIONS_s) , "Display");
  }
  if ((output_option) && ((view_point != 'q') ||
       ((args[0].tag() != AMObjType::SEQUENCES) && (args[0].tag() != AMObjType::MARKOVIAN_SEQUENCES) &&
        (args[0].tag() != AMObjType::VARIABLE_ORDER_MARKOV_DATA) &&
        (args[0].tag() != AMObjType::SEMI_MARKOV_DATA) &&
        (args[0].tag() != AMObjType::NONHOMOGENEOUS_MARKOV_DATA))) &&
      ((view_point != 'p') || ((args[0].tag() != AMObjType::HIDDEN_SEMI_MARKOV) &&
        (args[0].tag() != AMObjType::MARKOVIAN_SEQUENCES) &&
        (args[0].tag() != AMObjType::VARIABLE_ORDER_MARKOV_DATA) &&
        (args[0].tag() != AMObjType::SEMI_MARKOV_DATA) &&
        (args[0].tag() != AMObjType::NONHOMOGENEOUS_MARKOV_DATA)))) {
    status = false;
    genAMLError(ERRORMSG(INCOMPATIBLE_OPTIONS_s) , "Display");
  }
        """



        # ---------------- ViewPoint
        # 1-Survival
        if ViewPoint == 's':
            from openalea.stat_tool.enums import histogram_types
            from openalea.stat_tool.enums import model_distribution_types

            if type(self) in model_distribution_types:
                output = self.survival_ascii_write()
            elif type(self) in histogram_types:
                output = self.survival_ascii_write()
            else:
                raise ValueError("""(%s) has no survival point. Use another
                Viewpoint or use a first argument in DISTRIBUTION or MIXTURE or
                CONVOLUTION or COMPOUND or FREQUENCY_DISTRIBUTION or
                MIXTURE_DATA or CONVOLUTION_DATA or COMPOUND_DATA"""
                % str(type(self)))

        # Data
        elif ViewPoint == "d":
            try:
                #todo checkType
                # Markovian_Sequences, VOMData, SMData,
                # or Nonhomogenous_Markov_data
                output = self.ascii_data_write(exhaustive, Format)

            except Exception, e:
                #for vectors only
                #todo checkType
                try:
                    output = self.ascii_data_write(exhaustive)

                except AttributeError:
                    raise AttributeError("""
                        %s has not 'data' viewpoint""" % (str(type(self))))

        # StatProfile
        elif ViewPoint == 'p':
            try:
                from openalea.sequence_analysis._sequence_analysis import \
                    _HiddenVariableOrderMarkov, _HiddenSemiMarkov, \
                    _MarkovianSequences, _VariableOrderMarkovData, \
                    _SemiMarkovData# , _NonhomogenousMarkovData
            except:
                raise ImportError("openalea.sequence_analysis not found")
            assert len(args)>=1
            assert type(args[0]) is int

            if type(self) == _HiddenVariableOrderMarkov:
                assert len(args)==1
                self.state_profile_ascii_write(args[0], StateSequence,
                                               NbStateSequence)
            elif type(self) == _HiddenSemiMarkov:
                assert len(args)==1
                ## HERE HERE
                #args[0] is identifier
                output = self.state_profile_ascii_write(
                    args[0], Output, StateSequence, NbStateSequence)
            elif type(self) in [_MarkovianSequences, VariableOrderMarkovData,
                                _SemiMarkovData]:#, _NonhomogenousMarkovData]:
                assert len(args)==2
                if type(args[1]) == _HiddenVariableOrderMarkov:
                    args[1].state_profile_write(self , args[0], 'a',
                                     StateSequence , NbStateSequence)
                elif type(args[1]) == _HiddenSemiMarkov:
                    args[1].state_profile_write(self , args[0] , output ,
                                      'a' , StateSequence , NbStateSequence)
                else:
                    raise ValueError(
                    """Display with state profile requires 3 arugments and
                    second one must be HIDDEN_VARIABLE_ORDER_MARKOV or
                    HIDDEN_SEMI_MARKOV""")
            else:
                raise ValueError("Wrong arguments combinaison. Check them")

        #segment profile
        elif ViewPoint == 'q':
            try:
                self.nb_variable
            except:
                raise ValueError("object has no nb_variable. check that it is a sequence")
            nb_variable = self.nb_variable
            assert len(args)>=2
            error.CheckType([args[0], args[1]], [[int],[int]])
            #construct model_type
            try:
                from openalea.sequence_analysis.enums import model_type
            except:
                pass
            types = []
            for i in range(0, nb_variable):
                error.CheckType([args[i+2]], [str])
                if i==0:
                    types.append(model_type[args[i+2]])
                    #Multinomial or Poisson or Ordinal or Gaussian or
                    # Mean or Variance or MeanVariance
                    if args[i+2] in ["Mean","MeanVariance"]:
                        for j in range(1, nb_variable):
                             types.append(types[i])
                        break
                else:
                    # Multinomial or Poisson or Ordinal or Gaussian
                    # or Variance
                    types.append(model_type[args[i+2]])
            output = self.segment_profile_write(args[0], args[1], types, Output,
                                       'a', Segmentation, NbSegmentation)
        elif ViewPoint == 'v':
            from openalea.stat_tool.enums import all_stat_tool_types
            try:
                from openalea.sequence_analysis.enums import all_sequences_types
            except:
                pass
            if type(self) in all_stat_tool_types:
                output = self.ascii_write(exhaustive)
            elif type(self) in all_sequences_types:
                output = self.ascii_write(exhaustive)
            else:
                raise TypeError("wrong input type.")

        return output

    # TODO
    # Save(hsmc, ViewPoint="StateProfile", Sequence=1, Format="SpreadSheet")
    @add_doc
    def save(self, filename, Detail=2, ViewPoint="", Sequence="", Format="ASCII" ):

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
