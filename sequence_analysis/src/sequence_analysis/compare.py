#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""

.. topic:: compare.py summary

    A module dedicated to Compare functions

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: compare.py 9035 2010-06-01 19:17:27Z cokelaer $
"""
__version__ = "$Id: compare.py 9035 2010-06-01 19:17:27Z cokelaer $"


from openalea.stat_tool.comparison import compare_histo, compare_vectors
from openalea.stat_tool._stat_tool import I_DEFAULT
from openalea.stat_tool._stat_tool import _VectorDistance, _Vectors
from openalea.sequence_analysis._sequence_analysis import _SemiMarkov, \
     _HiddenSemiMarkov, _VariableOrderMarkov, _HiddenVariableOrderMarkov

from openalea.stat_tool.enums import histogram_types, bool_type, \
    format_type, algorithm_type

from enums import begin_aligned_map

from enums import sequence_alignment_first_arg, output_sequence, \
    markov_model_comparison_first_arg, indel_cost_map, ms_vomd_smd_nhmd,\
    markov_model_for_sequences_first_arg, markovian_algorithms, \
    markov_model_for_sequences_second_arg

from openalea.stat_tool import error


def _compare_markovian_sequences(obj, *args, **kargs):
    """
    :Example:

    .. doctest::
        :options: +SKIP
        
        >>> Compare(mc1, length_histo1, mc2, length_histo2,...,  FileName="result")
        >>> Compare(mc1, mc2,..., nb_seq, length, FileName="result")
        >>> Compare(mc1, seqm1, mc2, seqm2,..., nb_seq, FileName="result")
        >>> Compare(smc1, length_histo1, smc2, length_histo2,...,   FileName="result")
        >>> Compare(smc1, smc2,..., nb_seq, length, FileName="result")
        >>> Compare(smc1, seqm1, smc2, seqm2,..., nb_seq, FileName="result")
        >>> Compare(hmc1, length_histo1, hmc2, length_histo2,...,    FileName="result")
        >>> Compare(hmc1, hmc2,..., nb_seq, length, FileName="result")
        >>> Compare(hmc1, seqm1, hmc2, seqm2,..., nb_seq, FileName="result")
        >>> Compare(hsmc1, length_histo1, hsmc2, length_histo2,...,    FileName="result")
        >>> Compare(hsmc1, hsmc2,..., nb_seq, length, FileName="result")
        >>> Compare(hsmc1, seqm1, hsmc2, seqm2,..., nb_seq,  FileName="result")
    """

    error.CheckType([obj], [markov_model_comparison_first_arg])

    first_list = []
    second_list = []
    nb_seq = None
    length = None

    filename = kargs.get("Filename", None)

    from enums import ms_vomd_smd_list
    # Type of arg0 is same as type of obj., so we have the following case
    # >>> Compare(hsmc1, hsmc2, ,..., nb_seq,  FileName="result")
    if type(args[0]) == type(obj):

        #first_list.append(obj)
        for arg in args:
            if (isinstance(arg, int)) and nb_seq is None:
                nb_seq = arg
            elif (isinstance(arg, int)) and length is None:
                length = arg
            else:
                first_list.append(arg)

        return obj.divergence_computation_length(first_list, nb_seq,
                                                 length, filename)
    # Case where second arguments is Markovian and alternates with obj's type
    elif type(args[0]) in ms_vomd_smd_list:
        #first_list.append(obj)

        for arg in args:
            if (isinstance(arg, int)) and nb_seq is None:
                nb_seq = arg
            elif type(arg) == type(obj):
                first_list.append(arg)
            else:
                second_list.append(arg)
        return obj.divergence_computation_sequences(first_list, second_list,
                                                     nb_seq, filename)

    # Case where second arguments is histogram and then
    # alternates with obj's type
    elif (isinstance(arg[0], histogram_types)):
        hlength = []
        for arg in args:
            if type(arg) == type(obj):
                first_list.append(arg)
            else:
                hlength.append(arg)
        return obj.divergence_computation_histogram(first_list, hlength,
                                                     filename)

    else:
        raise Exception("case not handled. ")




def _compare_sequences(seq, *args, **kargs):
    """compare function related to sequences"""
    #int indel_cost = ADAPTATIVE , algorithm = AGGLOMERATIVE;
    #double indel_factor , transposition_factor = TRANSPOSITION_FACTOR;

    from openalea.sequence_analysis._sequence_analysis import \
        INDEL_FACTOR_1, INDEL_FACTOR_N, TRANSPOSITION_FACTOR

    error.CheckType([seq], [sequence_alignment_first_arg])

    ref_identifier = kargs.get("RefSequence", I_DEFAULT)
    test_identifier = kargs.get("TestSequence", I_DEFAULT)
    Begin = error.ParseKargs(kargs, "Begin", "Aligned", begin_aligned_map)
    End = error.ParseKargs(kargs, "End", "Aligned", begin_aligned_map)
    FileName = kargs.get("FileName", None)
    Format = error.ParseKargs(kargs, "Format", 'ASCII', format_type)
    AlignmentFormat = error.ParseKargs(kargs, "AlignmentFormat", 'ASCII',
                                       format_type)
    AlignmentFileName = kargs.get("AlignmentFileName", None)
    IndelCost = error.ParseKargs(kargs, "IndelCost", "Adaptative",
                                 indel_cost_map)
    IndelFactor = kargs.get("IndelFactor", INDEL_FACTOR_1)
    Transposition = error.ParseKargs(kargs, "Transposition", False, bool_type)
    TranspositionFactor = error.ParseKargs(kargs, "TranspositionFactor",
                                           TRANSPOSITION_FACTOR)
    Algorithm = error.ParseKargs(kargs, "Algorithm", "Agglomerative",
                                 algorithm_type)
    # check all int and float cases
    error.CheckType([ref_identifier, test_identifier, IndelFactor,
                     TranspositionFactor, TranspositionFactor],
                     [int, int, [int, float], [int, float], [int, float]])



    # case 2 of AML
    if len(args) == 1:

        if isinstance(args[0], _VectorDistance):

            Output = error.ParseKargs(kargs, "Output", "DistanceMatrix",
                                      output_sequence)

            if Output == 'm':
                if not error.ParseKargs(kargs, IndelFactor):
                    IndelFactor = INDEL_FACTOR_1


                if kargs.get("Algorithm", None):
                    raise ValueError("Algorithm cannot be used in this context")
                #todo: othet error case ?
                #if ((!transposition_option) && (transposition_factor_option)
                #if ((!file_name_option) && (format_option))
                #if ((!alignment_file_name_option) && (alignment_format_option)

                dist_matrix = seq.alignment_vector_distance(args[0],
                                    ref_identifier, test_identifier, Begin, End,
                                    IndelCost, IndelFactor, Transposition,
                                    TranspositionFactor, FileName, Format,
                                    AlignmentFileName, AlignmentFormat)
                return dist_matrix

            elif Output == 's':
                #check errors
                # if (ref_sequence_option):
                #if (test_sequence_option):
                # if (transposition_option)
                # if (transposition_factor_option)
                # if (format_option)
                # if (alignment_file_name_option)
                # if (alignment_format_option)

                # overwrite default values if Output == 's'
                if not error.ParseKargs(kargs, IndelFactor):
                    IndelFactor = INDEL_FACTOR_N


                sequence = seq.multiple_alignment(args[0], Begin, End,
                                    IndelCost, IndelFactor, Algorithm,
                                     FileName)

                if hasattr(seq, 'markovian_sequences'):
                    return sequence.markovian_sequences()
                else:
                    return sequence

    else: #case 1 of AML
        dist_matrix = seq.alignment(ref_identifier, test_identifier, Begin ,
                                  End , FileName , Format , AlignmentFileName,
                                  AlignmentFormat)
        return dist_matrix


def _compare_markovian_models_for_sequences(obj, *args, **kargs):
    """
    require test for the first, second and fourth if statement

    """

    error.CheckType([obj], [ms_vomd_smd_nhmd])

    Filename = kargs.get("Filename", None)
    Algorithm = error.ParseKargs(kargs, "Algorithm", 'Forward',
                                 markovian_algorithms)
    markov_list = []

    for arg in args:
        if type(arg) in markov_model_for_sequences_second_arg:
            markov_list.append(arg)

    if isinstance(args[0], _HiddenVariableOrderMarkov):
        return obj.comparison_hidden_variable_order_markov(markov_list,
                                                            Algorithm, Filename)
    if isinstance(args[0], _HiddenSemiMarkov):
        return obj.comparison_hidden_semi_markov(markov_list,
                                                  Algorithm, Filename)
    elif isinstance(args[0], _VariableOrderMarkov):
        return obj.comparison_variable_order_markov(markov_list, Filename)
    # obj should be a sequence
    elif isinstance(args[0], _SemiMarkov):
        return obj.comparison_semi_markov(markov_list, Filename)


def Compare(arg1, *args, **kargs):
    """Comparison functions factory

    * Comparison of frequency distributions.
    * Comparison of vectors.
    * Comparison of sequences.
    * Comparison of Markovian models for sequences.
    * Comparison of Markovian models.

    :Usage:

    .. doctest::
        :options: +SKIP
    
        >>> Compare(histo1, histo2,..., type, FileName="result", Format="ASCII")

        >>> Compare(vec, vector_distance)

        >>> Compare(seq, vector_distance, RefSequence=3, TestSequence=8,
                Begin="Free", End="Free", Transposition=True,
                FileName="result", Format="SpreadSheet",
                AlignmentFileName="alignment",
                AlignmentFormat="ASCII")
        >>> Compare(seq, RefSequence=3, TestSequence=8,
                Begin="Free", End="Free",  FileName="result",
                Format="SpreadSheet", AlignmentFileName="alignment",
                AlignmentFormat="ASCII")

        >>> Compare(discrete_seq, mc1, mc2,..., FileName="result")
        >>> Compare(discrete_seq, smc1, smc2,..., FileName="result")
        >>> Compare(discrete_seq, hmc1, hmc2,..., Algorithm="Viterbi",  FileName="result")
        >>> Compare(discrete_seq, hsmc1, hsmc2,..., Algorithm="Viterbi", FileName="result")

        >>> Compare(mc1, length_histo1, mc2, length_histo2,...,  FileName="result")
        >>> Compare(mc1, mc2,..., nb_seq, length, FileName="result")
        >>> Compare(mc1, seqm1, mc2, seqm2,..., nb_seq, FileName="result")
        >>> Compare(smc1, length_histo1, smc2, length_histo2,...,   FileName="result")
        >>> Compare(smc1, smc2,..., nb_seq, length, FileName="result")
        >>> Compare(smc1, seqm1, smc2, seqm2,..., nb_seq, FileName="result")
        >>> Compare(hmc1, length_histo1, hmc2, length_histo2,...,    FileName="result")
        >>> Compare(hmc1, hmc2,..., nb_seq, length, FileName="result")
        >>> Compare(hmc1, seqm1, hmc2, seqm2,..., nb_seq, FileName="result")
        >>> Compare(hsmc1, length_histo1, hsmc2, length_histo2,...,    FileName="result")
        >>> Compare(hsmc1, hsmc2,..., nb_seq, length, FileName="result")
        >>> Compare(hsmc1, seqm1, hsmc2, seqm2,..., nb_seq,  FileName="result")

    :Arguments:

    * histo1, histo2, ... (histogram, mixture_data, convolution_data, compound_data),
    * type (string): variable type ("NUMERIC" ("N"), "ORDINAL" ("O") or "SYMBOLIC" ("S")).

    * vec (vectors),
    * vector_distance (vector_distance).

    * seq (sequences, discrete_sequences, markov_data, semi-markov_data),

    * discrete_seq (discrete_sequences, markov_data, semi-markov_data),
    * mc1, mc2, ... (markov),
    * smc1, smc2, ... (semi-markov),
    * hmc1, hmc2, ... (hidden_markov),
    * hsmc1, hsmc2, ... (hidden_semi-markov).

    * mc1, mc2, ... (markov),
    * smc1, smc2, ... (semi-markov),
    * hmc1, hmc2, ... (hidden_markov),
    * hsmc1, hsmc2, ... (hidden_semi-markov),
    * length_histo1, length_histo2, ... (histogram, mixture_data,
      convolution_data, compound_data): frequency distribution of lengths
      of generated sequences,
    * nb_seq (int): number of generated sequences,
    * length (int): length of generated sequences,
    * seqm1, seqm2, ... (discrete_sequences, markov_data, semi-markov_data).


    :Optional Arguments:

    **vector case**

    * FileName (string): name of the result file,
    * Format (string): format of the result file: "ASCII" (default format) or
      "SpreadSheet". This optional argument can only be used in conjunction with
      the optional argument FileName.

    **sequence case**

    * RefSequence (int): identifier of the reference sequence,
    * TestSequence (int): identifier of the test sequence,
    * Begin(STRING: "Free" or "Fixed" (the default). If this optional argument
      is set at "Free", any space at the beginning of the alignment contrid=bute
      a weight of 0 (beginning-space free alignment).
    * End(STRING: "Free" or "Fixed" (the default). If this optional argument is
      set at "Free", any space at the beginning of the alignment contrid=bute a
      weight of 0 (end-space free alignment).
    * Transposition (bool): use of the transposition operation (default value:
      False). This optional argument requires the second mandatory argument
      being of type vector_distance.
    * FileName (string): name of result file,
    * Format (string): format of result file: "ASCII" (default format) or
      "SpreadSheet". This optional argument can only be used in conjunction
      with the optional argument FileName.
    * AlignmentFileName (string): file name of the sequences of edit operations
      (deletion/insertion/exact matching and eventually substitution and
      transposition) resulting from sequence alignments,
    * AlignmentFormat (string): format of the file of sequences of edit
      operations: "ASCII" (default format) or "Binary". This optional argument
      can only be used in conjunction with the optional argument
      AlignmentFileName.

    **markovian for sequences case**

    * Algorithm (string): type of algorithm: "Forward" (the default) or
      "Viterbi". This optional argument applies only with models of type
      hidden_markov or hidden_semi-markov,
    * FileName (string): name of result file.

    **markovian case**

    * FileName (string): name of result file. If this optional argument is used,
      some complementary results, with respect to the returned object of type
      distance_matrix, are saved on a file.



    :Returned Object:

    **distribution and markovian for sequences cases**

    No object returned.

    **vectors, sequences, markovian cases**

    An object of type distance_matrix is returned.


    :Description:

    If number of alignments <= 30, the alignments are displayed in the shell window.

    The result of comparisons is displayed in the shell window.


    :Background:

    **sequences case**

    The result of the comparison of two sequences takes the form of the sequence
    of edit operations required for transforming the reference sequence into the
    test sequence with associated costs:

    * deletion (`d`): deletion of an element of the reference sequence,
    * insertion (`i`): insertion of an element of the test sequence,
    * exact matching (`m`): matching of an element of the test sequence
      with the same element in the reference sequence: null cost,
    * substitution (`s`) (mismatching): replacement of an element of the
      reference sequence by another element of the test sequence: The associated
      cost is the standardized distance between the two elements.
    * transposition (two successive `t`) or swap: interchange of adjacent
      elements in the sequence with the additional constraint that each element
      can participate in no more than one swap. This edit operation applies only
      in the case where two successive element in the reference sequence exactly
      match two successive elements in the test sequence taken in reverse order.

    A purely structural alignment consists in allowing only exact matching,
    insertion and deletion. In this case, the argument *vector_distance* specifying
    the local distance between elements is not required.

    **markovian for sequences case**

    In the case of Markov chains (type markov) or semi-Markov chains (type
    semi-markov), the comparison relies on the likelihood of a sequence for the
    different models being compared. In the case of hidden Markov chains (type
    hidden_markov) or hidden semi-Markov chains (hidden_semi-markov), the
    comparison relies either on the likelihood of a sequence (of every possible
    state sequences that can generate the observed sequence: Algorithm="Forward"),
    or on the likelihood of the state sequence that best explains the observed
    sequence (Algorithm="Viterbi").

    **markovian case**

    The comparison of Markovian models relies on the Kullback-Leibler directed
    divergence (or cross-entropy or discriminant information). For each model
    being compared, a sample of sequences is generated and the log-likelihoods
    of these sequences for the different models are computed (including the
    reference model used for simulation). The dissimilarity measure is the
    `divergence` from the reference model to the target model on the
    basis of log-likelihoods of the sequences normalized by their cumulative
    length. This procedure is repeated for each model.


    .. seealso::

        :func:`~openalea.stat_tool.vectors.VectorDistance`,
        :func:`~openalea.stat_tool.cluster.Clustering`.
        :func:`~openalea.sequence_analysis.data_transform.ComputeStateSequences`.
        :func:`~openalea.stat_tool.estimate.Estimate`


    """

    p1 = arg1


    # COMPARE 1
    if type(p1) in histogram_types:
        return compare_histo(arg1, *args, **kargs)
    # COMPARE 2
    elif isinstance(p1, _Vectors):
        return compare_vectors(arg1, *args, **kargs)
    #COMPARE 3
    elif type(p1) in sequence_alignment_first_arg and len(args)==0:
        return _compare_sequences(arg1, *args, **kargs)
    #COMPARE3 bis
    elif type(p1) in sequence_alignment_first_arg and \
            isinstance(args[0], _VectorDistance):
        return _compare_sequences(arg1, *args, **kargs)
    #Compare 4
    elif type(p1) in markov_model_for_sequences_first_arg and \
            type(args[0]) in markov_model_for_sequences_second_arg:
        return _compare_markovian_models_for_sequences(arg1, *args, **kargs)
    #COMPARE 5
    elif type(p1) in markov_model_comparison_first_arg:
        return _compare_markovian_sequences(arg1, *args, **kargs)


    raise Exception("Error in Compare. No case corresponding to your command."
                    "Check your arguments.")


