/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2015 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: sequence_label.cpp 18070 2015-04-23 10:52:51Z guedon $
 *
 *       Forum for V-Plants developers:
 *
 *  ----------------------------------------------------------------------------
 *
 *                      GNU General Public Licence
 *
 *       This program is free software; you can redistribute it and/or
 *       modify it under the terms of the GNU General Public License as
 *       published by the Free Software Foundation; either version 2 of
 *       the License, or (at your option) any later version.
 *
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS For A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */



namespace sequence_analysis {


/****************************************************************
 *
 *  Mots cles (format des fichiers) :
 */


const char *SEQ_word[] = {
  "MARKOV_CHAIN" ,
  "EQUILIBRIUM_MARKOV_CHAIN" ,
  "HIDDEN_MARKOV_CHAIN" ,
  "EQUILIBRIUM_HIDDEN_MARKOV_CHAIN" ,

  "SEMI-MARKOV_CHAIN" ,
  "EQUILIBRIUM_SEMI-MARKOV_CHAIN" ,
  "HIDDEN_SEMI-MARKOV_CHAIN" ,
  "EQUILIBRIUM_HIDDEN_SEMI-MARKOV_CHAIN" ,

  "NONHOMOGENEOUS_MARKOV_CHAIN" ,
  "HOMOGENEOUS" ,
  "NONHOMOGENEOUS" ,

  "OCCUPANCY_DISTRIBUTION" ,

  "INDEX_PARAMETER"

//  "TOP_PARAMETERS" ,
//  "AXILLARY_PROBABILITY" ,
//  "RHYTHM_RATIO"
};


const char *SEQ_index_parameter_word[] = {
  " " ,
  "TIME" ,
  "TIME_INTERVAL" ,
  "POSITION" ,
  "POSITION_INTERVAL"
};



/****************************************************************
 *
 *  Labels :
 */


const char *SEQ_label[] = {
  "log-likelihood of the state sequence" ,
  "log-likelihood of the state sequences" ,
  "log-likelihood of the observed sequences" ,
  "information of the sequences in the iid case" ,

  "smoothed" ,
  "observed" ,
  "theoretical" ,
  "smoothed observed probabilities" ,

  "ordinary renewal process" ,
  "equilibrium renewal process" ,
  "time between 2 observation" ,
  "inter-event" ,
  "backward" ,
  "forward" ,
  "recurrence time" ,
  "length-biased" ,
  "inter-event time within the observation period" ,
  "inter-event time censored on both ends" ,
  "inter-event time censored on one end" ,
  "complete inter-event time" ,
  "time up to event" ,
  "number of event" ,
  "during" ,
  "time unit" ,
  "mixture of number of event distributions" ,
  "no-event probability" ,
  "event probability" ,

  "Markov chain" ,
  "hidden Markov chain" ,
  "semi-Markov chain" ,
  "hidden semi-Markov chain" ,

  "maximum order" ,
  "memory tree" ,
  "transition tree" ,
  "memory transition matrix" ,
  "non-terminal" ,
  "terminal" ,
  "completion" ,
  "completed" ,
  "confidence intervals for transition probabilities" ,
  "free transient parameter" ,
  "free transient parameters" ,
  "recommended maximum order" ,
  "pruning threshold" ,
  "initial counts" ,
  "transition counts" ,
  "maximum transition count difference" ,
  "log-likelihoods" ,
  "count" ,
  "delta" ,
  "Krichevsky-Trofimov" ,
  "likelihood ratio test" ,

  "self-transition" ,
  "asymptote" ,

  "occupancy distribution" ,

  "probability of no-occurrence of " ,
  "time up to the first occurrence" ,
  "time up to the first occurrence of " ,
  "probability of leaving " ,
  "absorption probability of " ,
  "biased" ,
  "occupancy" ,
  "complete/censored state occupancy weights",
  "sojourn time" ,
  "initial run" ,
  "final run" ,
  "mixture of " ,
  "number of runs" ,
  "number of runs of " ,
  "number of occurrences" ,
  "number of occurrences of " ,
  "per sequence" ,
  "per length" ,
  "missing value" ,
  "words" ,

  "state probabilities" ,
  "posterior state sequence probability" ,
  "posterior state sequence probability log ratio" ,
  "state begin" ,
  "posterior state probabilities" ,
  "posterior in state probabilities" ,
  "posterior out state probabilities" ,
  "conditional entropy" ,
  "marginal entropy" ,
  "sum of marginal entropies" ,
  "partial state sequence entropy" ,
  "state sequence entropy" ,
  "state sequence divergence" ,
  "upper bound" ,
  "number of state sequences" ,
  "maximum posterior state probabilities" ,
  "maximum posterior in state probabilities" ,
  "maximum posterior out state probabilities" ,
  "likelihood ratio" ,

  "correlation function" ,
  "partial" ,
  "auto" ,
  "cross-" ,
  "Pearson" ,
  "Spearman" ,
  "Kendall" ,
  "rank" ,
  "lag" ,
  "maximum lag" ,
  "white noise" ,
  "randomness 95% confidence limit" ,
  "pair frequency" ,

  "index" ,

  "simulated" ,
  "sequence" ,
  "sequences" ,
  "vertex identifier" ,
  "index parameter" ,
  "minimum index parameter" ,
  "maximum index parameter" ,
  "time" ,
  "time interval" ,
  "position" ,
  "position interval" ,
  "length" ,
  "sequence length" ,
  "cumulative length" ,
  "shift" ,

  "alignment length" ,
  "aligned on" ,
  "maximum gap length" ,
  "alignment coding" ,
  "consensus" ,

  "optimal" ,
  "change point" ,
  "change points" ,
  "mean change-point amplitude" ,
  "segment" ,
  "segments" ,
  "segment sample size" ,
  "global standard deviation" ,
  "piecewise linear function" ,
  "number of segments" ,
  "posterior probability" ,
  "penalty" ,
  "number of segmentations" ,
  "segmentations" ,
  "segmentation log-likelihood" ,
  "log-likelihood of all the possible segmentations" ,
  "posterior change-point probabilities" ,
  "posterior segment probabilities" ,
  "segmentation entropy" ,
  "first-order dependency entropy" ,
  "change-point entropy" ,
  "uniform entropy" ,
  "segmentation divergence" ,
  "begin conditional entropy" ,
  "end conditional entropy" ,
  "maximum change-point likelihood" ,
  "maximum segment likelihood" ,
  "maximum posterior change-point probabilities" ,
  "maximum posterior segment probabilities" ,
  "ambiguity"

//  "top" ,
//  "tops" ,
//  "number of internode"
};



/****************************************************************
 *
 *  Messages d'erreur pour l'analyse des fichiers :
 */


const char *SEQ_parsing[] = {
  "time data not ordered" ,
  "time data too large" ,
  "number of event data not ordered" ,

  "bad state" ,
  "bad number of memories" ,

  "time index not ordered" ,
  "position not ordered" ,
  "position not allowed" ,
  "bad maximum sequence length: should be greater than 1"
};



/****************************************************************
 *
 *  Messages d'erreur de traitement :
 */


const char *SEQ_error[] = {
  "only time interval censored on both ends: choose a longer observation period" ,
  "incompatible renewal data" ,
  "maximum number of events too small: choose a longer observation period" ,
  "average number of events too small: choose a longer observation period" ,
  "time unit too large" ,
  "number of complete intervals too small" ,
  "complete time interval: bad minimum value" ,
  "forward recurrence time interval: bad minimum value" ,
  "no event observation period: bad minimum value" ,
  "bad inter-event mean computation method" ,
  "initial inter-event distribution support incompatible with data" ,

  "bad time between two observation dates" ,
  "bad minimum time" ,
  "bad maximum time" ,
  "bad minimum number of events" ,
  "bad maximum number of events" ,
  "empty renewal data structure" ,
  "time between two observation dates too short" ,
  "time between two observation dates too long" ,

  "bad model structure" ,
  "single state component" ,
  "bad number of states" ,
  "missing state" ,
  "bad order" ,
  "bad minimum order" ,
  "bad maximum order" ,
  "too many parameters" ,
  "overlap of values observed in the different states" ,
  "bad model type" ,
  "bad output process type" ,
  "no parametric output process" ,
  "bad minimum number of state sequences" ,
  "bad number of state sequences" ,
  "average state occupancy too short" ,
  "bad number of sequences" ,
  "bad sequence identifier" ,
  "bad sequence identifiers" ,
  "bad reference sequence identifier" ,
  "bad test sequence identifier" ,
  "bad sequence length" ,
  "sequence length too short" ,
  "sequence length too long" ,
  "cumulative sequence length too long" ,
  "states not represented" ,
  "failure in the computation of the optimal state sequences" ,
  "reference model" ,
  "target model" ,
  "number of failures for the Kullback-Leibler divergence estimation" ,

  "vertex identifier not allowed: change the sample order" ,
  "bad vertex identifier" ,
  "bad index parameter type" ,
  "bad index parameter" ,
  "bad variable indices" ,
  "bad variable lag" ,
  "bad date order" ,
  "bad begin index parameter" ,
  "bad end index parameter" ,
  "bad minimum sequence length" ,
  "bad maximum sequence length" ,
  "bad maximum run length" ,
  "bad minimum index parameter" ,
  "bad maximum index parameter" ,
  "bad number of selected values" ,
  "unequal index intervals: should be equal" ,
  "bad number of sequences: should be > 1" ,
  "bad position transform step" ,
  "bad length" ,
  "bad value" ,
  "bad correlation coefficient type" ,
  "bad frequency" ,
  "bad maximum lag" ,
  "incompatible with other correlation functions" ,
  "bad differencing order" ,
  "initial run histograms already built" ,
  "bad run length" ,
  "too high number of possible words" ,
  "bad minimum frequency: should be positive" ,

  "state sequences not in the data" ,
  "characteristics not computed" ,
  "consecutive values from 0" ,
  "non-existing characteristic distribution" ,
  "non-existing forward sojourn time distribution" ,
  "incompatible with model" ,
  "sequence incompatible with model" ,

  "too many alignment" ,
  "bad insertion/deletion factor: should be greater than 0.5" ,
  "bad transposition factor: should be between 0 and 2" ,

  "forbidden output" ,
  "bad number of segments" ,
  "bad minimum number of segments" ,
  "bad maximum number of segments" ,
  "bad change point" ,
  "segmentation failure" ,
  "bad number of segmentations" ,
  "bad change-point model"

//  "bad position" ,
//  "bad number of internodes" ,
//  "bad top identifier" ,
//  "bad main axe number of internodes: should be greater than the last position" ,
//  "bad minimum position" ,
//  "bad maximum position" ,
//  "bad neighborhood" ,
//  "not enough neighbors" ,
//  "equality of growth probabilities not possible" ,
//  "bad number of tops" ,
//  "bad number of trials" ,
//  "bad number of axillary shoots per node"
};


};  // namespace sequence_analysis
