# -*- python -
#
#       OpenAlea.Container
#
#       Copyright 2012 INRIA - CIRAD - INRA
#
#       File author(s):  Jonathan Legrand <jonathan.legrand@ens-lyon.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite: http://openalea.gforge.inria.fr
#
################################################################################
"""This module helps to use clustering and standardization methods on graphs."""

__license__ = "Cecill-C"
__revision__ = ""

import warnings, numpy as np, copy, math, time
from numpy import ndarray
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import gzip, cPickle as pickle

from sklearn.cluster import SpectralClustering, Ward, DBSCAN
from sklearn import metrics

def distance_matrix_from_vector(data, var_type, no_dist_index = []):
    """
    Function creating a distance matrix based on a vector (list) of values.
    Each values are attached to an individual.

    :Parameters:
     - `data` (list) - vector (list) of values
     - `var_type` (str) - type of variable

    :Returns:
     - `dist_mat` (np.array) - distance matrix
    """
    implemented_var_types = ["numeric", "ordinal"]
    N = len(data)
    dist_mat = np.zeros( shape = [N,N], dtype=float )

    if var_type.lower() == "numeric":
        for i in xrange(N):
            for j in xrange(i+1,N):# we skip when i=j because in that case the distance is 0.
                if data[i] is None or data[j] is None or i in no_dist_index or j in no_dist_index:
                    dist_mat[i,j] = dist_mat[j,i] = None
                else:
                    dist_mat[i,j] = dist_mat[j,i] = abs(data[i]-data[j])

    elif var_type.lower() == "ordinal":
        rank = data.argsort() # In case of ordinal variables, observed values are replaced by ranked values.
        for i in xrange(N):
            for j in xrange(i+1,N): # we skip when i=j because in that case the distance is 0.
                if data[i] is None or data[j] is None or i in no_dist_index or j in no_dist_index:
                    dist_mat[i,j] = dist_mat[j,i] = None
                else:
                    dist_mat[i,j] = dist_mat[j,i] = abs(rank[i]-rank[j])

    else:
        raise ValueError("Unknown variable type '{}'... should either be {}".format(var_type, implemented_var_types))
    return dist_mat


def mad_based_outlier(points, thresh=3.5):
    """
    The median absolute deviation (MAD) is a robust measure of the variability of a univariate sample of quantitative data.
    For a univariate data set X1, X2, ..., Xn, the MAD is defined as the median of the absolute deviations from the data's median:
    $$ MAD(X) = median_i(|X_i - median_j(X_j)|) $$
    In order to use the MAD as a consistent estimator for the estimation of the standard deviation $\sigma$, one takes
    $$ \hat{\sigma}=K . MAD, $$
    where $K$ is a constant scale factor, which depends on the distribution.
    """
    isdict = False
    if isinstance(points, dict):
        keys, points = points.keys(), points.values()
        isdict = True
    if not isinstance(points, ndarray):
        points = np.array(points)
    # - Detect the presence of NaNs unsupported by np.sum
    nans_index = np.isnan(points)
    points_no_nans = np.array([p for n,p in enumerate(points) if not nans_index[n]])
    if len(points_no_nans.shape) == 1:
        points_no_nans = points_no_nans[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points_no_nans - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)
    # - Compute a modified Z-score:
    modified_z_score = 0.6745 * diff / med_abs_deviation
    # - If NaNs were present in the begining, we introduce 'np.inf' so they are detected as outliers when returned !
    if sum(nans_index) != 0:
        for n in range(len(points)):
            if nans_index[n]:
                modified_z_score = list(modified_z_score[:n]) + [np.inf] + list(modified_z_score[n:])

    if isdict:
        return dict([ (keys[n], mzs > thresh) for n,mzs in enumerate(modified_z_score) ])
    else:
        return np.array(modified_z_score) > thresh


def vectorORarray(data, return_data=False):
    """
    Function checking if the provided `data` is a vector of value or a square matrix (array)
    """
    # -- Identifying case where numpy.array are 1D vectors:
    if isinstance(data,ndarray) and (len(data.shape)==1) or ((data.shape[0] == 1) and (data.shape[1] > 1)) or ((data.shape[1] == 1) and (data.shape[0] > 1)):
        data = data.tolist()

    if isinstance(data,list) or isinstance(data,tuple):
        if return_data:
            return 'vector', data
        else:
            return 'vector'

    if isinstance(data,ndarray) and (data.shape[0]==data.shape[1]) and ((data.shape[0]!=1)and(data.shape[1]!=1)):
        if return_data:
            return 'array', np.array(data)
        else:
            return 'array'


def missing_values(dist_mat):
    nan_index = np.isnan(dist_mat)
    if True in nan_index:
        nb_missing_values = len(np.where(nan_index is True)[0])
    else:
        nb_missing_values = 0.
    return nb_missing_values

def matrix_outliers2nan(dist_mat, outliers_index):
    dist_mat[outliers_index,:] = dist_mat[:,outliers_index] = np.nan
    return dist_mat

def vector_outliers2nan(vect, outliers_index):
    vect[outliers_index] = np.nan
    return vect

def standardisation_L1(dist_mat):
    N = dist_mat.shape[0]
    nb_missing_values = missing_values(dist_mat)
    return np.nansum(np.nansum(abs(dist_mat))) / (N*(N-1)-nb_missing_values)

def standardisation_L2(dist_mat):
    N = dist_mat.shape[0]
    nb_missing_values = missing_values(dist_mat)
    return np.nansum(np.nansum(dist_mat**2)) / (N*(N-1)-nb_missing_values)

def create_weight_matrix(N,weight,standard_distance_matrix):
    tmp_w_mat = np.zeros( shape = [N,N], dtype=float )
    tmp_w_mat.fill(weight)
    tmp_w_mat[np.where(np.isnan(standard_distance_matrix))]==np.nan
    return tmp_w_mat


def _within_cluster_distances(distance_matrix, clustering):
    """
    Function computing within cluster distance.
    $$ D(q) = \dfrac{ \sum_{i,j \in q; i \neq j} D(i,j) }{(N_{q}-1)N_{q}} ,$$
    where $D(i,j)$ is the distance matrix, $self._N$ is the total number of elements and $N_{q}$ is the number of elements found in clusters $q$.

    :Parameters:
     - `distance_matrix` (np.array) - distance matrix used to create the clustering.
     - `clustering` (list) - list giving the resulting clutering.

    :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
    """
    clusters_ids = list(set(clustering))
    nb_clusters = len(clusters_ids)
    nb_ids_by_clusters = [sum([i==q for i in clustering]) for q in clusters_ids]
    D_within = {}
    for n,q in enumerate(clusters_ids):
        index_q = [j for j,i in enumerate(clustering) if i==q]
        D_within[q] = 2. * sum( [distance_matrix[i,j] for i in index_q for j in index_q if j>i] ) / ( (nb_ids_by_clusters[n]-1) * nb_ids_by_clusters[n])

    if nb_clusters == 1:
        return D_within.values()[0]
    else:
        return D_within


def _between_cluster_distances(distance_matrix, clustering):
    """
    Function computing within cluster distance.
    $$ D(q) = \dfrac{ \sum_{i \in q} \sum_{j \not\in q} D(i,j) }{ (self._N - N_q) N_q }, $$
    where $D(i,j)$ is the distance matrix, $self._N$ is the total number of elements and $N_{q}$ is the number of elements found in clusters $q$.

    :Parameters:
     - `distance_matrix` (np.array) - distance matrix used to create the clustering.
     - `clustering` (list) - list giving the resulting clutering.

    :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
    """
    N = len(clustering)
    clusters_ids = list(set(clustering))
    nb_clusters = len(clusters_ids)
    nb_ids_by_clusters = [sum([i==q for i in clustering]) for q in clusters_ids]

    if 1 in nb_ids_by_clusters:
        raise ValueError("A cluster contain only one element!")

    D_between = {}
    for n,q in enumerate(clusters_ids):
        index_q = [j for j,i in enumerate(clustering) if i==q]
        index_not_q = [j for j,i in enumerate(clustering) if i!=q]
        D_between[q] = sum( [distance_matrix[i,j] for i in index_q for j in index_not_q] ) / ( (N-nb_ids_by_clusters[n]) * nb_ids_by_clusters[n])

    if nb_clusters == 1:
        return D_between.values()[0]
    else:
        return D_between


def _global_cluster_distances(distance_matrix, clustering):
    """
    Function computing global cluster distances, i.e. return the sum of within_cluster_distance and between_cluster_distance.

    :Parameters:
     - `distance_matrix` (np.array) - distance matrix used to create the clustering.
     - `clustering` (list) - list giving the resulting clutering.

    :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
    """
    w = _within_cluster_distances(distance_matrix, clustering)
    b = _between_cluster_distances(distance_matrix, clustering)
    N = len(clustering)
    clusters_ids = list(set(clustering))
    nb_ids_by_clusters = [len(np.where(clustering == q)[0]) for q in clusters_ids]

    gcd_w = sum( [ (nb_ids_by_clusters[q]*(nb_ids_by_clusters[q]-1))/float(sum([nb_ids_by_clusters[l]*(nb_ids_by_clusters[l]-1) for l in clusters_ids if l != q])) * w[q] for q in clusters_ids] )
    gcd_b = sum( [(N-nb_ids_by_clusters[q])*nb_ids_by_clusters[q]/float(sum([(N-nb_ids_by_clusters[l])*nb_ids_by_clusters[l] for l in clusters_ids if l != q])) * b[q] for q in clusters_ids] )
    return gcd_w, gcd_b


def _fn_distance2affinity(dist):
    return np.exp(-1 * dist / dist.std())


def CH_estimator(k,w,b,N):
    """
    Index of Calinski and Harabasz (1974).
    $$ \text{CH}(k) = \dfrac{B(k) / (k-1)}{W(k) / (n-k)} $$
    where $B(k)$ and $W(k)$ are the between- and within-cluster sums of squares, with $k$ clusters.
    The idea is to maximize $\text{CH}(k)$ over the number of clusters $k$. $\text{CH}(1)$ is not defined.

    :Parameters:
     - k: number of clusters;
     - w: global WITHIN cluster distances;
     - b: global BETWEEN cluster distances;
     - N: population size.
    """
    return (b/(k-1))/(w/(N-k))


def Hartigan_estimator(k,w_k,N):
    """
    Index of Hartigan (1975).
    Hartigan (1975) proposed the statistic:
    $$ \text{H}(k) = \left\{ \dfrac{W(k)}{W(k+1)} - 1 \right\} / (n-k-1)$$
    The idea is to start with $k=1$ and to add a cluster as long as $\text{H}(k)$ is sufficiently large.\\
    One can use an approximate $F$-distribution cut-off; instead Hartigan suggested that a cluster be added if $H(k) > 10$. Hence the estimated number of clusters is the smallest $k \geqslant 1$ such that $\text{H}(k) \leqslant 10$. This estimate is defined for $k=1$ and can potentially discriminate between one \textit{versus} more than one cluster.

    :Parameters:
     - k: number of clusters;
     - w_k: DICTIONARY of global WITHIN cluster distances (must contain k and k+1);
     - N: population size.
    """
    return (w_k[k]/w_k[k+1]-1)/(N-k-1)


def contiguous_ids(dic, starting_id = 0):
    """
    Check that the `dic` dictionary values are contiguous starting from `starting_id`, otherwise make them contiguous.
    """
    uniq, mini, maxi = list(set(dic.values())), min(dic.values()), max(dic.values())
    N_ids = len(uniq)
    if (N_ids == len(range(mini, maxi+1))) and (mini==starting_id):
        return dic
    elif (N_ids == len(range(mini, maxi+1))) and (mini!=starting_id):
        diff = mini - starting_id
        return dict([(k,v-diff) for k,v in dic.iteritems()])
    else:
        diffs = np.array(uniq) - range(starting_id, N_ids)
        return dict([(k,v-diffs[uniq.index(v)]) for k,v in dic.iteritems()])


def cluster2labels(clusters_dict):
    """
    Create a dict *key=clusters; *labels=ids from a dict *key=ids; *labels=clusters.
    """
    cluster2labels = {}
    for k,v in clusters_dict.iteritems():
        try:
            cluster2labels[v].append(k)
        except:
            cluster2labels[v] = []

    return cluster2labels


def cost_function(ids_1, ids_2, similarity = False):
    """
    The matching elements cost function (dissimilarity):
    $$ D(a,b) = 1 - (2* intersect(a,b)) / (N_a + N_b), $$
     where N_a and N_b are the number of elements in ensemble a and b.

    :Parameters:
     - `ids_1` (list|set) - list or set of id-like elements belonging to a cluster
     - `ids_2` (list|set) - list or set of id-like elements belonging to another cluster (from another clustering!)
     - `similarity` (bool) - if True, return similarity function instead of dissimilarity
    """
    if isinstance(ids_1, list):
        ids_1 = set(ids_1)
    if isinstance(ids_2, list):
        ids_2 = set(ids_2)
    if similarity:
        return float(2*len(ids_1 & ids_2))/(len(ids_1)+len(ids_2))
    else:
        return 1-float(2*len(ids_1 & ids_2))/(len(ids_1)+len(ids_2))


def ensemble_cost_function( cluster2labels_1, cluster2labels_2, similarity = False ):
    """
    Cost function between two clusters relating to the number of common ids in them.

    :Parameters:
     - `cluster2labels_1` (dict) - a dict *key=clusters; *labels=ids
     - `cluster2labels_2` (dict) - a dict *key=clusters; *labels=ids
     - `similarity` (bool) - if True, return similarity function instead of dissimilarity
    """
    cost_triplets = []
    for clusters_1, ids_1 in cluster2labels_1.iteritems():
        for clusters_2, ids_2 in cluster2labels_2.iteritems():
            if set(ids_1) & set(ids_2) != set([]):
                cost_triplets.append([clusters_1, clusters_2, cost_function(ids_1, ids_2, similarity)])

    return cost_triplets


def cluster_matching(clusters_dict_1, clusters_dict_2, return_all = False):
    """
    Function calling the BipartiteMatching C++ code to match clusters based on a Minimum cost flow algorithm over their ids composition.

    :Parameters:
     - `cluster2labels_1` (dict) - a dict *key=clusters; *labels=ids
     - `cluster2labels_2` (dict) - a dict *key=clusters; *labels=ids
     - `return_all` (bool) - if True, return score and unassociated groups
    """
    from openalea.tree_matching.bipartitematching import BipartiteMatching
    c2l1 = cluster2labels(clusters_dict_1)
    c2l2 = cluster2labels(clusters_dict_2)
    res = BipartiteMatching(c2l1.keys(), c2l2.keys(), ensemble_cost_function(c2l1,c2l2), [1.0 for i in xrange(len(c2l1))], [1.0 for i in xrange(len(c2l2))])
    match = res.match()

    if return_all:
        return match
    else:
        return match[1]


def clustering_naming(clustering_method, nb_clusters, global_distance_weights, global_distance_variables, outliers=None):
    if outliers is None:
        out_name = ""
    elif outliers == 'ignored':
        out_name = '-ignored_outliers'
    else:
        out_name = '-deleted_outliers'

    return str(clustering_method)+"_"+str(nb_clusters)+"_"+str([str(global_distance_weights[n])+"*"+str(global_distance_variables[n]) for n in xrange(len(global_distance_weights))])+out_name


def load_mvpd(fname):
    """
    Save the object on disk under `fname`.
    """
    print "Trying to open the mvpd_matrix file {}...".format(fname)
    t_start = time.time()
    with open(fname, 'rb') as input:
        loaded_mvpd = pickle.load(input)

    print "Time to load the mvpd_matrix: {}s".format(round(time.time()-t_start,3))

    tmp_obj = mvpd_matrix()
    unknown_att = set(loaded_mvpd.__dict__.keys())-set(tmp_obj.__dict__.keys())
    if unknown_att != set([]):
        print "Unknown attribute '{}' added to the object... please check versions!".format(unknown_att)

    unfound_att = set(tmp_obj.__dict__.keys())-set(loaded_mvpd.__dict__.keys())
    if unfound_att != set([]):
        print "Some attributes could not be found in the loaded object: '{}'".format(unfound_att)
        print "... please check versions!"

    return loaded_mvpd


class mvpd_matrix:
    """
    Class allowing to create a Multi-Variates Pairwise Distances Matrix by combining several variables given as vectors of values or (uni-variate) pairwise distance matrices.
    This class is usefull to create an MVPD matrix by combining different types of variables (informations) in a weigthed manner.
    The sum of weight must equal to one!
    We handle missing data by automatically re-weighting according to the presence or absence of data for each cmobined information.
    """
    def __init__(self, variable_list=[], variable_names=[], variable_types=[], variable_weights=[], variable_units=[], standardisation_method='L1', standardize_over=None):
        """
        Initialisation of the multivariate pairwise distance matrix object. The default standardisation method is to use the L1 metric.
        All '__init__' parameters can remain empty, but in such case you will have to add variables using 'self.add_variables'.
        If you provide 'variable_list', we expect you to provide 'variable_names' and 'variable_types' as well. All those lists should be of the same lenght !
        If you provide 'variable_weights', we will compute the multivariate pairwise distance matrix.

        :Parameters:
         - `variable_list` (list(list|array)) - list of vectors or arrays of values representing the variable observations to add;
         - `variable_names` (list(str)) - list of names to give to the variables;
         - `variable_types` (list(str)) - list of strings declaring the type of each used variable; takes values in ["Numeric", "Ordinal", "Interval"];
         - `variable_weights` (list(str)) - (Optional) list of strings giving the unit (dimension) of each variables;
         - `variable_units` (list(str)) - (Optional) list of strings giving the unit (dimension) of each variables;
         - `standardisation_method` (str) - (Optional) define which standarisation metric to apply to the data; takes value in ["L1", "L2"]; default is "L1";
        """
        # -- Initialisation:
        self._nb_var = len(variable_list)
        self.standardisation_method = standardisation_method

        # -- Checking provided information coherence:
        assert self._nb_var == len(variable_names)
        assert self._nb_var == len(variable_types)
        if len(variable_weights) != 0:
            assert self._nb_var == len(variable_weights)
        if len(variable_units) != 0:
            assert self._nb_var == len(variable_units)
        else:
            variable_units = [None for i in xrange(self._nb_var)]

        # -- Attributes saving 'initial' informations about variables:
        self._nb_values = None
        self._distance_matrix_dict = {} # save NON-standardized pairwise distance matrix with the varible name as key.
        self._distance_matrix_info = {}
        self._outliers_index_dict = None
        self._outliers_management = None
        # -- Attributes for caching information:
        self._var_data_dict = {}
        self._global_distance_matrix = None
        self._global_distance_ids = None
        self._global_distance_weights = None
        self._global_distance_variables = None
        # -- Attributes declaring the function used to transform the distance matrix into an affinity matrix (Spectral Clustering)
        #self._fn_distance2affinity = lambda dist: np.exp(-1 * dist / dist.std())
        # -- Attributes related to the clustering:
        self._method = None
        self._nb_clusters = None
        self._clustering = None # dict *keys=observation_ids, *values=clustering_labels
        self._full_tree = None
        self._clusterings_dict = None # dict *keys=`_nb_clusters`, *values=`_clustering`

        # -- If a list of variable observations, names and types are provided, we compute the NON-standardized pairwise distance matrix
        if self._nb_var != 0:
            for i in xrange(self._nb_var):
                self.add_variable(variable_list[i], variable_names[i], variable_types[i], variable_units[i])
            print "Done integrating variables {} as separate NON-standardised pairwise distance matrices.".format(variable_names)

        # -- If the weigths are also provided, we can compute the multi-variate pairwise distance matrix rigth now:
        if len(variable_weights) > 0:
            self.create_mvpd_matrix(variable_names, variable_weights, standardize_over)

        # -------------- DONE!


    def save_mvpd(self, fname):
        """
        Save the object on disk under `fname`.
        """
        import cPickle as pickle, gzip, time

        t_start = time.time()
        with open(fname, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

        return "{}s to save the mvpd_matrix under filename {}".format(round(time.time()-t_start,3), fname)


    def add_variable(self, var_data, var_name, var_type, var_unit = None, MAD_outliers_thres = None):
        """
        Add a distance matrix related to vertices properties form the graph.

        :Parameters:
         - `var_data` (list|array) -  a vector or array of values representing the variable to add;
         - `var_name` (str) - (Optional) name to give to the variable
         - `var_type` (str) - string or list of strings declaring the type of properties used, should be 'Numeric','Ordinal'
         - `var_unit` (str) - (Optional) unit (dimension) to give to the variable
         - `MAD_outliers_thres` (int) - (Optional) threshold used for the detection of outliers using MeanAbsDev.
        """
        # NEED to check if we already have a pairwise distance matrix and if the data to add have the adequate number of observations!!
        if len(self._distance_matrix_dict.keys())!=0:
            if vectorORarray(var_data) == "vector":
                nb_values = len(var_data)
            else:
                nb_values = var_data.shape[0]
            assert self._nb_values == nb_values

        if (var_name is None or var_name == ""):
            raise ValueError("You have to give a name to your variable if you want to use it!")

        # - We make sure self._distance_matrix_dict can receive the pairwise distance matrix with the name `var_name`:
        if self._distance_matrix_dict.has_key(var_name):
            raise KeyError("You already have a property named '{}'".format(var_name))

        print("Computing the pairwise distance matrix for the variable '{}'...".format(var_name))
        if vectorORarray(var_data) == 'vector':
            self._var_data_dict[var_name] = var_data # kept to be able to represent data dispersion in each cluster;
            self._distance_matrix_dict[var_name] = distance_matrix_from_vector(var_data, var_type) # used for clustering purposes;
        else:
            self._distance_matrix_dict[var_name] = var_data

        # - Adding its infos tho the list
        self._distance_matrix_info[var_name] = (var_type.lower(), var_unit)
        self._nb_values = self._distance_matrix_dict.values()[0].shape[0]

        if MAD_outliers_thres is not None:
            if vectorORarray(var_data) != "vector":
                raise TypeError("To detect outliers, on the basis of their MeanAbsDev, you need to provide the data as a vector of values!")
            else:
                self._outliers_index_dict[var_name] = mad_based_outlier(var_data, MAD_outliers_thres)

        return "Done adding and creating the pairwise distance matrix for the variable '{}'.".format(var_name)


    def remove_pairwise_distance_matrix(self, var_id):
        """
        Remove a pairwise distance matrix form the dictionary `self._distance_matrix_dict` and it's attached information in `self._distance_matrix_info`
        """
        self._distance_matrix_dict.pop(var_id)
        self._distance_matrix_info.pop(var_id)
        return "Done."

    def _standardisation(self, data, norm = 'L1', variable_types = None, outliers_index = []):
        """
        :Parameters:
         - `data` (np.array|list) - pairwise distance matrix (if np.array) or vector of values (if list or 1D array);
         - `norm` (str) - define which standarisation metric to apply to the data, takes value in ["L1", "L2"];
        :Optional:
         - `variable_types` (str) - used only if a vector of values is provided for `data`, takes value in ["Numeric", "Ordinal", "Interval"];
         - `outliers_index` (list) - list gathering the (previously detected) outliers position (index) in `data`;

        :Returns:
         - `standard_mat` (np.array) - standardized distance matrix
        """
        vOa, data = vectorORarray(data, True)
        if vOa == 'vector':
            distance_matrix = distance_matrix_from_vector(data, variable_types)
        elif vOa == 'array':
            distance_matrix = data
        else:
            raise ValueError("Can not convert the provided data.")

        # -- We create a copy of the `distance_matrix``as it should also contain the outliers values in the final pairwise distance matrix
        # (the ouliers are indeed removed only for the dispersion/standardisation measure computation).
        dist_mat = copy.copy(distance_matrix)
        # -- Handling outliers by setting their value to np.nan (will exclude them of standardisation value computation).
        if outliers_index != []:
            dist_mat[outliers_index,:] = dist_mat[:,outliers_index] = np.nan

        # -- Now we can start the standardisation:
        if norm.upper() == "L1":
            absd = standardisation_L1(dist_mat)
            return distance_matrix / absd
        elif norm.upper() == "L2":
            sd = standardisation_L2(dist_mat)
            return distance_matrix / sd
        else:
            raise ValueError("Wrong type of standardisation metric, choose 'L1' or 'L2'.")

    def create_mvpd_matrix(self, variable_names, variable_weights, standardize_over=None, ignore_outliers = False, delete_outliers = False, return_data = False):
        """
        Funtion creating the global weighted distance matrix.
        Provided `variable_names` should exist in `self._distance_matrix_dict`.
        Outliers can be computed when adding a vertex using MAD estimator.
        When using `standardize_over`, remember that the list index start at O !!

        :Parameters:
         - `variable_names` (list) - list of variables names to combine; the names should be keys in `self._distance_matrix_dict`;
         - `variable_weights` (list) - list of weightss used to create the global weighted distance matrix (MUST sum to 1!);
         - `standardize_over` (list) - variables index (as given in variable_names) to be used simultaneously for standardisation value computation;
         - `ignore_outliers` (bool) - if True ignore outliers when computing standardisation value;
         - `delete_outliers` (bool) - if True 'delete outliers' (values set to np.nan) of pairwise distance matrix;
         - `return_data` (bool) - if true the function return the sorted ids list `ids` & a dictionary of pairwise distance matrix `global_matrix`.
        """
        print "# -- Computing the multi-variate pairwise distance matrix..."
        if isinstance(variable_names,str):
            variable_names = [variable_names]
        if isinstance(variable_weights,int) or isinstance(variable_weights,float):
            variable_weights = [float(variable_weights)]
        assert len(variable_names) == len(variable_weights)
        assert math.fsum(variable_weights)==1.

        # -- Checking the presence of all requested information (i.e. variables names) in the dictionary `self._distance_matrix_dict`:
        for var_name in variable_names:
            if not self._distance_matrix_dict.has_key(var_name):
                raise KeyError("'{}' is not in the dictionary `self._distance_matrix_dict`".format(var_name))

        if ((not ignore_outliers) and (not delete_outliers)):
            self._outliers_management = None
        elif ignore_outliers:
            self._outliers_management = "ignored"
        else:
            self._outliers_management = "deleted"

        # -- Shortcut when asking for the same result:
        if variable_weights == self._global_distance_weights and variable_names == self._global_distance_variables and ids == self._global_distance_ids and self._outliers == outliers_management:
            if return_data:
                return self._global_distance_ids, self._global_distance_matrix
            else:
                print "- The global pairwise distance matrix was already computed!"
                return None
        else: # Otherwise, clean any possible clustering:
            self._method = None
            self._nb_clusters = None
            self._clustering = None
            self._full_tree = None
            self._clusterings_dict = None

        # By default we use all the given observations of the variables:
        ids = range(self._nb_values) ## SHOULD BE MODIFIED according to the way outliers should be handdled !!

        # -- Standardization step:
        standard_distance_matrix = {}
        mat_topo_dist_standard = []
        nb_var = 0
        if standardize_over is not None:
            self._standardize_over_variables(standardize_over, variable_names)
            print "- Found the following standardisation values: {}".format(self._standardisation_values)

        nb_std_var = 0
        for n, var_name in enumerate(variable_names):
            # reduce the matrix to the list of selected ids `ids`
            distance_matrix = self._distance_matrix_dict[var_name][ids,:][:,ids]
            # handle case where outliers need to be deleted of the distance matrix:
            if self._outliers_management == "deleted":
                distance_matrix = self._outliers_handler(distance_matrix, var_name)
            elif self._outliers_management == "ignored":
                outliers_index = outliers_index = self._outliers_index_dict[var_name]
            else:
                outliers_index = []
            # now we can compute the standardised version of the pairwise distance matrix:
            if ((standardize_over is not None) and (n in standardize_over)):
                standard_distance_matrix[var_name] = distance_matrix/self._standardisation_values[n]
            else:
                standard_distance_matrix[var_name] = self._standardisation(distance_matrix, self.standardisation_method, outliers_index)
            nb_var += 1

        # -- Checking for a simple case: to 'combine' only ONE pairwise distance matrix (no re-weighting to do !)
        if nb_var == 1:
            global_matrix = standard_distance_matrix.values()[0]
        else:
            # - Creating weight matrix and replacing nan by zeros for computation in standardized matrix.
            N = len(ids)
            weight_matrix, standardized_matrix = [], []
            standardized_matrix.extend( [np.nan_to_num(standard_distance_matrix[var_name]) for var_name in variable_names] )
            for n, var_name in enumerate(variable_names):
                weight_matrix.append(create_weight_matrix(N, variable_weights[n], standard_distance_matrix[var_name]))

            # Finally making the global weighted pairwise standard distance matrix:
            for w, w_m in zip(variable_weights, weight_matrix):
                binary = ~np.isnan(w_m)
                weight_matrix = [weight_matrix[r]*(binary+~binary*w) for r in xrange(nb_var)]

            global_matrix = np.zeros( shape = [N,N], dtype=float )
            for wei_mat, standard_mat in zip(weight_matrix, standardized_matrix):
                global_matrix += np.nan_to_num(wei_mat) * standard_mat

        # -- We update caching variables :
        self._global_distance_matrix = global_matrix
        self._global_distance_ids = ids
        self._global_distance_weights = variable_weights
        self._global_distance_variables = variable_names

        if return_data:
            return ids, global_matrix
        else:
            print "...done !"

    def _standardize_over_variables(self, standardize_over, variable_names):
        self._standardisation_values = {}
        if isinstance(standardize_over, list) and isinstance(standardize_over[0], int):
            standardize_over = [standardize_over]
        for var_list in standardize_over:
            print "Computing common standarisation value for variables {}...".format(var_list)
            big_vector = []
            for var_id in var_list:
                var_name = variable_names[var_id]
                big_vector.extend(self._outliers_handler(self._var_data_dict[var_name], var_name))
            # Now that we have created the 'big_vector' we can compute the common standardisation value:
            dist_mat = distance_matrix_from_vector(big_vector, self._distance_matrix_info[var_name][0])
            if self.standardisation_method.upper() == "L1":
                std_val = standardisation_L1(dist_mat)
            elif self.standardisation_method.upper() == "L2":
                std_val = standardisation_L2(dist_mat)
            else:
                raise ValueError("Wrong type of standardisation metric, choose 'L1' or 'L2'.")

            self._standardisation_values.update(dict([(var_id, std_val) for var_id in var_list]))
        return "Done."


    def _outliers_handler(self, data, var_name):
        if self._outliers_management is None: 
            outliers_index = []
        else:
            try:
                outliers_index = self._outliers_index_dict[var_name]
            except:
                print "No outliers defined for property '{}'... skipping that part.".format(var_name)

        if (vectorORarray(data)=='array'):
            return matrix_outliers2nan(data, outliers_index)
        elif (vectorORarray(data)=='vector'):
            return vector_outliers2nan(data, outliers_index)
        else:
            return "`data` type not understood!"


    def cluster(self, n_clusters, method = "ward"):
        """
        Actually run the clustering method depending on the selected method and the number of clusters.

        :Parameters:
         - `n_clusters` (int) - number of cluster to create;
         - `method` (str) - clustering method to use, "ward", "spectral" and "DBSCAN";
        """
        # -- Usual paranoia:
        if n_clusters is None:
            raise ValueError("You have to provide the number of clusters you want for the Ward method.")
        if self._global_distance_matrix is not None:
            global_matrix = copy.copy(self._global_distance_matrix)
        else:
            raise ValueError("No standardised pairwise distance matrix saved, please give one!")

        if method.lower() == "ward":
            clustering = Ward(n_clusters = n_clusters, compute_full_tree=False, connectivity=None).fit(global_matrix)
            clustering_labels = clustering.labels_

        if method.lower() == "spectral":
            affinity_matrix = _fn_distance2affinity(global_matrix)
            clustering = SpectralClustering(n_clusters = n_clusters, affinity='precomputed').fit(affinity_matrix)
            clustering_labels = list(clustering.labels_)

        self._clustering = dict([ (label,clustering_labels[n]) for n,label in enumerate(self._global_distance_ids) ])
        self._method = method.lower()
        self._nb_clusters = n_clusters

        return "Done computing {} clustering for {} clusters.".format(self._method, self._nb_clusters)


    def full_tree(self, range_clusters=None, method = "ward", global_matrix = None, connectivity = None):
        """
        Actually run the clustering method.
        :Parameters:
         - `range_clusters` (list) (optional) - if provided return a dict of clustering for each cluster number in it
         - `method` (str) (optional) - clustering method to use, "ward", "spectral" and "DBSCAN"
         - `global_matrix` (np.array) (optional) - distance matrix to cluster
         - `connectivity` (np.array) (optional) - connectivity matrix to use while clustering (if possible regarding the clustering method)
        """
        # -- Usual paranoia:
        if global_matrix is None:
            if self._global_distance_matrix is not None:
                global_matrix = copy.copy(self._global_distance_matrix)
            else:
                raise ValueError("No distance matrix saved, please give one!")

        if self._full_tree is None:
            if method.lower() == "ward":
                full_tree = Ward(compute_full_tree=True, connectivity=connectivity).fit(global_matrix)
            if method.lower() == "spectral":
                print('Not ready yet!')
                return None
            self._full_tree = full_tree
            self._method = method.lower()
        else:
            full_tree = self._full_tree
        # Stop HERE if no `range_clusters` provided.
        if range_clusters is None:
            return full_tree

        n_clusters = None
        if isinstance(range_clusters,int):
            n_clusters = copy.copy(range_clusters)
            assert n_clusters < full_tree.n_leaves_
            self._nb_clusters = n_clusters
            # If the clustering already exist, return the desired clustering as 'cluster' function would:
            if (self._clusterings_dict is not None) and self._clusterings_dict.has_key(n_clusters):
                print "Returning clustering previously computed by Cluster.full_tree()."
                clustering = self._clusterings_dict[n_clusters]
                self._clustering = clustering
                return clustering
            # If not we define a range of clusters to compute later:
            if n_clusters < 6 :
                m = 5
                while m*n_clusters > full_tree.n_leaves_:
                    m-=1
                range_clusters = range(2, m*n_clusters)
            elif n_clusters+5 < full_tree.n_leaves_:
                range_clusters = range(2, n_clusters+5)
            else:
                range_clusters = range(2, n_clusters)

        # Check: if we already computed that and if yes if it contain what we need (the clustering associated to each `n_clusters` in `range_clusters`) otherwise we compute it:
        if (self._clusterings_dict is None) or (sum([self._clusterings_dict.has_key(k) for k in range_clusters])!=len(range_clusters)):
            from sklearn.cluster.hierarchical import _hc_cut
            clusterings_dict = {}
            for n, c in enumerate(range_clusters):
                clust_labels = _hc_cut(c, full_tree.children_, full_tree.n_leaves_)
                clusterings_dict[c] = dict(zip(self._global_distance_ids, clust_labels))
                if n > 0 and clusterings_dict.has_key(range_clusters[n-1]):
                    clust_comp = ClustererComparison(clusterings_dict[range_clusters[n-1]], clusterings_dict[range_clusters[n]])
                    clusterings_dict[c] = clust_comp.relabel_clustering_2()

            self._clusterings_dict = clusterings_dict
        if n_clusters is not None:
            return self._clusterings_dict[n_clusters]
        else:
            return self._clusterings_dict


    def silhouette_estimators(self, clustering_method, k_min=4, k_max=15, beta = 1, plot_estimator = True):
        """
        Compute various estimators based on clustering results.

        :Parameters:
         - distance_matrix (np.array): distance matrix to be used for clustering;
         - clustering_method (str): clustering methods to be applyed, must be "Ward" or "Spectral"
         - clustering_range (list): range of clusters to consider.
        """
        clustering_labels={}
        sil = {}
        assert k_min>1

        for k in xrange(k_min, k_max+1):
            if clustering_method.lower() == "ward":
                clustering = Ward(n_clusters=k, compute_full_tree=True).fit(self._global_distance_matrix)
                clustering_labels[k] = clustering.labels_
            if clustering_method.lower() == "spectral":
                similarity = np.exp(-beta * self._global_distance_matrix / self._global_distance_matrix.std())
                clustering = SpectralClustering(n_clusters = k, precomputed=True).fit(self._global_distance_matrix)
                clustering_labels[k] = clustering.labels_

            sil[k] = metrics.silhouette_score(self._global_distance_matrix, clustering_labels[k], metric='euclidean')

        if plot_estimator:
            fig = plt.figure(figsize=(4,4),dpi=100)
            plt.plot(xrange(k_min, k_max+1), sil.values(), color='red')
            plt.title("Silhouette estimator")

        return sil


    def clustering_estimators(self, clustering_method, k_min=4, k_max=15, beta = 1, plot_estimator = True):
        """
        Compute various estimators based on clustering results.

        :Parameters:
         - distance_matrix (np.array): distance matrix to be used for clustering;
         - clustering_method (str): clustering methods to be applyed, must be "Ward" or "Spectral"
         - clustering_range (list): range of clusters to consider.
        """
        clustering_labels, w, N = {}, {}, {}
        CH, sil = {}, {}
        assert k_min>=1

        for k in xrange(k_min, k_max+1):
            if clustering_method.lower() == "ward":
                clustering = Ward(n_clusters=k, compute_full_tree=True).fit(self._global_distance_matrix)
                clustering_labels[k] = clustering.labels_
            if clustering_method.lower() == "spectral":
                similarity = np.exp(-beta * self._global_distance_matrix / self._global_distance_matrix.std())
                clustering = SpectralClustering(n_clusters = k, precomputed=True).fit(self._global_distance_matrix)
                clustering_labels[k] = clustering.labels_

            N[k] = len(clustering_labels[k])
            if k!=1:
                w[k], b = _global_cluster_distances(self._global_distance_matrix, clustering_labels[k])
                CH[k] = CH_estimator(k,w[k],b,N[k])
                sil[k] = metrics.silhouette_score(self._global_distance_matrix, clustering_labels[k], metric='euclidean')
            else:
                w[k] = _within_cluster_distances(self._global_distance_matrix, clustering_labels[k])

        Hartigan = dict( [(k, Hartigan_estimator(k,w,N[k])) for k in xrange(k_min, k_max)] )

        if plot_estimator:
            fig = plt.figure(figsize=(12,4),dpi=100)
            fig.add_subplot(131)
            plt.plot(xrange(k_min if k_min>1 else 2, k_max+1),CH.values())
            plt.title("Calinski and Harabasz estimator")
            fig.add_subplot(132)
            plt.plot(xrange(k_min, k_max), Hartigan.values(), color='green')
            plt.title("Hartigan estimator")
            fig.add_subplot(133)
            plt.plot(xrange(k_min if k_min>1 else 2, k_max+1), sil.values(), color='red')
            plt.title("Silhouette estimator")

        return CH, Hartigan, sil


    def group_clusters(self, regions2group, make_contiguous_ids = True):
        """
        """
        return group_regions(self, regions2group, region_names, make_contiguous_ids)


    def group_regions(self, regions2group, region_names=None, make_contiguous_ids = True):
        """
        Function modifying clusterer._clustering by grouping ids.
         :Parameters:
         - `regions2group` (list) - list of strings or ids naming regions (at least 2) to group, should be referenced in self.graph.graph_property()
        """
        assert len(regions2group)>=2
        if isinstance(regions2group[0], str):
            assert region_names is not None
            regions2group_index = [region_names.index(r_name) for r_name in regions2group]
        else:
            regions2group_index = regions2group

        # - Now we replace the regions2group ids with the min of their ids
        cid = min(regions2group_index)
        regions_left = list(set(regions2group_index)-set([cid]))
        clust = {}
        for k, v in self._clustering.iteritems():
            clust[k] = cid if v in regions_left else v

        if make_contiguous_ids:
            clust = contiguous_ids(clust, 0)

        self._clustering = clust
        self._nb_clusters = len(np.unique(self._clustering.values()))



class ClustererChecker:
    """
    Class allowing to analyse a clustered `mvpd_matrix` object.
    """
    def __init__(self, mvpd_matrix):
        # - Paranoia :
        assert mvpd_matrix.__class__.__name__ == 'mvpd_matrix'
        # - Initialisation
        self._clustering = mvpd_matrix._clustering.values()
        self._distance_matrix = mvpd_matrix._global_distance_matrix
        self._var_data_dict = dict([(var_name, mvpd_matrix._var_data_dict[var_name]) for var_name in mvpd_matrix._global_distance_variables])
        self._var_units_dict = dict([(var_name, mvpd_matrix._distance_matrix_info[var_name][1]) for var_name in mvpd_matrix._global_distance_variables])
        # - Precompute usefull values :
        self._N = len(self._clustering)
        self._clusters_ids = list(set(self._clustering))
        self._nb_clusters = len(self._clusters_ids)
        self._nb_ids_by_clusters = dict( [ (q, self._clustering.count(q)) for q in self._clusters_ids] )
        self._ids_by_clusters = dict( [(q, [mvpd_matrix._global_distance_ids[k] for k in [i for i,j in enumerate(self._clustering) if j==q]]) for q in self._clusters_ids] )
        self._sorted_ids_by_clusters = self.sorted_ids_by_clusters()

        # - Reformat inherited (usefull) info :
        self.info_clustering = {'method': mvpd_matrix._method,
                                'variables': mvpd_matrix._global_distance_variables,
                                'weights': mvpd_matrix._global_distance_weights}

        self.clustering_name = clustering_naming(mvpd_matrix._method, mvpd_matrix._nb_clusters, mvpd_matrix._global_distance_weights, mvpd_matrix._global_distance_variables)

        # - Check for undesirable situtations:
        if min(self._nb_ids_by_clusters.values()) <=1:
            raise ValueError("In the provided clustering a cluster has only one element, this is wrong!")
            print "Provided clustering: {}".format(self.clustering_name)

        # -- DONE!
        print "ClustererChecker object initialisation done!"


    def cluster_distance_matrix(self, round_digits = 3):
        """
        Function computing distance between clusters.
        For $\ell \eq q$  :
        \[ D(q,\ell) = \dfrac{ \sum_{i,j \in q; i \neq j} D(i,j) }{(N_{q}-1)N_{q}} , \]
        For $\ell \neq q$  :
        \[ D(q,\ell) = \dfrac{ \sum_{i \in q} \sum_{j \in \ell} D(i,j) }{N_{q} N_{\ell} } , \]
        where $D(i,j)$ is the distance matrix, $N_{q}$ and $N_{\ell}$ are the number of elements found in clusters $q$ and $\ell$.

        :Hidden parameters:
         - `self._distance_matrix` (np.array) - distance matrix used to create the clustering.
         - `self._clustering` (list) - list giving the resulting clutering.

        :WARNING: `distance_matrix` and `clustering` should obviously be ordered the same way!
        """
        D = np.zeros( shape = [self._nb_clusters,self._nb_clusters], dtype = float )
        for n,q in enumerate(self._clusters_ids):
            for m,l in enumerate(self._clusters_ids):
                if n==m:
                    index_q = [i for i,j in enumerate(self._clustering) if j==q]
                    D[n,m] = sum( [self._distance_matrix[i,j] for i in index_q for j in index_q if i!=j] ) / ( (self._nb_ids_by_clusters[n]-1) * self._nb_ids_by_clusters[n])
                if n>m:
                    index_q = [i for i,j in enumerate(self._clustering) if j==q]
                    index_l = [i for i,j in enumerate(self._clustering) if j==l]
                    D[n,m] = D[m,n]= sum( [self._distance_matrix[i,j] for i in index_q for j in index_l] ) / (self._nb_ids_by_clusters[n] * self._nb_ids_by_clusters[m])

        if round_digits is None:
            return D
        else:
            return np.round(D,int(round_digits))


    def plot_cluster_distances(self, print_clustering_name=True, savefig=None, print_values=False, **kwargs):
        """
        Display a heat-map of cluster distances with matplotlib.
        """
        cluster_distances = self.cluster_distance_matrix()
        numrows, numcols = cluster_distances.shape

        fig = plt.figure(figsize=(6,5))
        ax1 = fig.add_subplot(111)
        plt.imshow(cluster_distances, cmap=cm.winter, interpolation='nearest')
        plt.xticks(range(numcols), range(numcols), fontsize=14, fontweight='bold')
        plt.yticks(range(numrows), range(numrows), fontsize=14, fontweight='bold')

        if print_values:
            font = {'family': 'monospace', 'color': 'white', 'weight': 'bold', 'size': 16}
            alignment = {'horizontalalignment':'center', 'verticalalignment':'center'}
            for nrow in range(numrows):
                for ncol in range(numcols):
                    plt.text(nrow, ncol, cluster_distances[nrow,ncol], fontdict=font, **alignment)

        if kwargs.has_key('title'):
            if isinstance(kwargs['title'],str):
                plt.title(kwargs['title'])
            elif kwargs['title']:
                plt.title("Cluster distances heat-map")
        if print_clustering_name:
            plt.suptitle(self.clustering_name)

        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize=14)
        ax1.axes.tick_params(length=0)
        plt.tight_layout()
        if isinstance(savefig,str):
            plt.savefig(savefig)
            plt.close()
        else:
            plt.show()


    def cluster_diameters(self, round_digits = 3):
        """
        Function computing within cluster diameter, i.e. the max distance between two vertex from the same cluster.
        $$ \max_{i,j \in q} D(j,i) ,$$
        where $D(i,j)$ is the distance matrix and $q$ a cluster, .

        :Hidden parameters:
         - `self._distance_matrix` (np.array) - distance matrix used to create the clustering.
         - `self._clustering` (list) - list giving the resulting clutering.

        :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
        """
        if round_digits is None:
            round_digits = 13
        diameters = {}
        for q in self._clusters_ids:
            index_q = [i for i,j in enumerate(self._clustering) if j==q]
            diameters[q] = np.round(max([self._distance_matrix[i,j] for i in index_q for j in index_q]),round_digits)

        if self._nb_clusters == 1:
            return diameters.values()
        else:
            return diameters


    def cluster_separation(self, round_digits = 3):
        """
        Function computing within cluster diameter, i.e. the min distance between two vertex from two diferent clusters.
        $$ \min_{i \in q, j \not\in q} D(j,i) ,$$
        where $D(i,j)$ is the distance matrix and $q$ a cluster.

        :Hidden parameters:
         - `self._distance_matrix` (np.array) - distance matrix used to create the clustering.
         - `self._clustering` (list) - list giving the resulting clutering.

        :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
        """
        if round_digits is None:
            round_digits = 13
        separation = {}
        for q in self._clusters_ids:
            index_q = [i for i,j in enumerate(self._clustering) if j==q]
            index_not_q = [i for i,j in enumerate(self._clustering) if j!=q]
            separation[q] = np.round(min([self._distance_matrix[i,j] for i in index_q for j in index_not_q ]),round_digits)

        if self._nb_clusters == 1:
            return separation.values()
        else:
            return separation

    def within_cluster_distances(self):
        """
        Function computing within cluster distance.
        """
        return _within_cluster_distances(self._distance_matrix, self._clustering)


    def between_cluster_distances(self):
        """
        Function computing within cluster distance.
        """
        return _between_cluster_distances(self._distance_matrix, self._clustering)


    def global_cluster_distances(self):
        """
        Function computing global cluster distances, i.e. return the sum of within_cluster_distance and between_cluster_distance.
        """
        return _global_cluster_distances(self._distance_matrix, self._clustering)


    def __score_param(func):
        def wrapped_function(self, dict_labels_expert, groups2compare = []):
            """
            Wrapped function for clustering score computation according to the knowledge of the ground truth class assignments.

            :Parameters:
             - dict_labels_expert (dict) - expert defined regions / clusters in wich keys are labels
             - groups2compare (list) - pair(s) of groups id to compare, with first the expert id then the predicted id ex. [0,6] or [[0,6],[4,3]]
            """
            if groups2compare != []:
                compare_groups = True
                groups2compare = np.array(groups2compare, ndmin=2)
            else:
                compare_groups = False

            clustering_dict = dict(zip(self._vtx_list, self._clustering))
            not_found = []
            labels_true, labels_pred = [], []
            max1 = max(dict_labels_expert.values())+1
            max2 = max(clustering_dict.values())+1
            for k,v in dict_labels_expert.iteritems():
                if clustering_dict.has_key(k):
                    v2 = clustering_dict[k]
                    if compare_groups:
                        labels_true.append(v if (v in groups2compare[:,0]) else max1)
                        labels_pred.append(v2 if (v2 in groups2compare[:,1]) else max2)
                        #~ labels_pred.append(v2)
                    else:
                        labels_true.append(v)
                        labels_pred.append(v2)
                else:
                    not_found.append(k)

            if not_found != []:
                warnings.warn("These labels were not found in the clustering result: {}".format(not_found))

            return func(labels_true, labels_pred)
        return wrapped_function


    @__score_param
    def adjusted_rand_score(labels_true, labels_pred):
        """
        The Adjusted Rand Index (ARI) is a function that measures the similarity of the two assignments, ignoring permutations and with chance normalization.

        :Parameters:
         - `labels_true` (list) - knowledge of the ground truth class assignments
         - `labels_pred` (list) - clustering algorithm assignments of the same samples

        :Notes:
         - Random (uniform) label assignments have a ARI score close to 0.0.
         - Bounded range [-1, 1]. Negative values are bad (independent labelings), similar clusterings have a positive ARI, 1.0 is the perfect match score.
         - No assumption is made on the cluster structure: can be used to compare clustering algorithms such as k-means which assumes isotropic blob shapes with results of spectral clustering algorithms which can find cluster with "folded" shapes.
        """
        return metrics.adjusted_rand_score(labels_true, labels_pred)

    @__score_param
    def adjusted_mutual_info_score(labels_true, labels_pred):
        """
        The Mutual Information (NMI and AMI) is a function that measures the agreement of the two assignments, ignoring permutations.
        Adjusted Mutual Information (AMI) was proposed more recently than NMI and is normalized against chance.

        :Parameters:
         - `labels_true` (list) - knowledge of the ground truth class assignments
         - `labels_pred` (list) - clustering algorithm assignments of the same samples

        :Note:
         - Random (uniform) label assignments have a AMI score close to 0.0.
        """
        return metrics.adjusted_mutual_info_score(labels_true, labels_pred)

    @__score_param
    def normalized_mutual_info_score(labels_true, labels_pred):
        """
        The Mutual Information (NMI and AMI) is a function that measures the agreement of the two assignments, ignoring permutations.
        Normalized Mutual Information (NMI) is often used in the literature, but it is NOT normalized against chance.
        """
        return metrics.normalized_mutual_info_score(labels_true, labels_pred)

    @__score_param
    def homogeneity_score(labels_true, labels_pred):
        """
        Homogeneity: each cluster contains only members of a single class.
        Bounded below by 0.0 and above by 1.0 (higher is better).

        :Parameters:
         - `labels_true` (list) - knowledge of the ground truth class assignments
         - `labels_pred` (list) - clustering algorithm assignments of the same samples

        :Note:
         - homogeneity_score(a, b) == completeness_score(b, a)
        """
        return metrics.homogeneity_score(labels_true, labels_pred)

    @__score_param
    def completeness_score(labels_true, labels_pred):
        """
        Completeness: all members of a given class are assigned to the same cluster.
        Bounded below by 0.0 and above by 1.0 (higher is better).

        :Parameters:
         - `labels_true` (list) - knowledge of the ground truth class assignments
         - `labels_pred` (list) - clustering algorithm assignments of the same samples

        :Note:
         - homogeneity_score(a, b) == completeness_score(b, a)
        """
        return metrics.completeness_score(labels_true, labels_pred)

    @__score_param
    def v_measure_score(labels_true, labels_pred):
        """
        Harmonic mean of homogeneity and completeness_score is called V-measure.

        :Parameters:
         - `labels_true` (list) - knowledge of the ground truth class assignments
         - `labels_pred` (list) - clustering algorithm assignments of the same samples

        :Note:
         - `v_measure_score` is symmetric, it can be used to evaluate the agreement of two independent assignments on the same dataset.
        """
        return metrics.v_measure_score(labels_true, labels_pred)

    @__score_param
    def homogeneity_completeness_v_measure(labels_true, labels_pred):
        """
        Homogeneity, completensess and V-measure can be computed at once using homogeneity_completeness_v_measure

        :Parameters:
         - `labels_true` (list) - knowledge of the ground truth class assignments
         - `labels_pred` (list) - clustering algorithm assignments of the same samples
        """
        return metrics.homogeneity_completeness_v_measure(labels_true, labels_pred)


    def vertex2clusters_distance(self):
        """
        Compute the mean distance between a vertex and those from each group.
        $$ D(i,q) = \dfrac{ \sum_{i \neq j} D(i,j) }{ N_q } \: , \quad \forall i \in [1,self._N] \:, \: q \in [1,Q],$$
        where $D(i,j)$ is the distance matrix and $q$ a cluster.

        :Hidden parameters:
         - `self._distance_matrix` (np.array) - distance matrix used to create the clustering.
         - `self._clustering` (list) - list giving the resulting clutering.

        :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
        """
        index_q = {}
        for n,q in enumerate(self._clusters_ids):
            index_q[q] = [i for i,j in enumerate(self._clustering) if j==q]

        vertex_cluster_distance = {}
        for i in xrange(self._N):
            tmp = np.zeros( shape=[self._nb_clusters], dtype=float )
            for n,q in enumerate(self._clusters_ids):
                tmp[n] = sum([self._distance_matrix[i,j] for j in index_q[q] if i!=j])/(self._nb_ids_by_clusters[n]-1)

            vertex_cluster_distance[i] = tmp

        return vertex_cluster_distance


    def vertex_distance2cluster_center(self):
        """
        Compute the distance between a vertex and the center of its group.

        :Hidden parameters:
         - `self._distance_matrix` (np.array) - distance matrix used to create the clustering.
         - `self._clustering` (list) - list giving the resulting clutering.

        :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
        """
        D_iq = self.vertex2clusters_distance()
        vertex_distance2center = {}
        for i in xrange(self._N):
            vertex_distance2center[i] = D_iq[i][self._clustering[i]]

        return vertex_distance2center


    def ids_in_cluster(self,cluster_id):
        """
        Return a dictionary which keys are cluster ids and values are dictance-to-center sorted individuals.
        """
        return self._ids_by_clusters[cluster_id]

    def ids_by_clusters(self):
        """
        Return a dictionary which keys are cluster ids and values are dictance-to-center sorted individuals.
        """
        return self._ids_by_clusters

    def sorted_ids_by_clusters(self, return_distance=False):
        """
        Return a dictionary which keys are cluster ids and values are dictance-to-center sorted individuals.
        """
        vtx2center = self.vertex2clusters_distance()
        sorted_ids_by_clusters, sorted_dist_by_clusters = {}, {}
        for clust, ids in self._ids_by_clusters.iteritems():
            dist2center = [vtx2center[i][clust] for i in ids]
            sorted_ids_by_clusters[clust] = [ids[i] for i in np.argsort(dist2center)]
            sorted_dist_by_clusters[clust] = sorted(dist2center)

        self._sorted_ids_by_clusters = sorted_ids_by_clusters
        if return_distance:
            return self._sorted_ids_by_clusters, sorted_dist_by_clusters
        else:
            return self._sorted_ids_by_clusters


    def plot_vertex_distance2cluster_center(self, cluster_names=None, print_clustering_name=True, n_colors=None, savefig=None):
        """
        Plot the distance between a vertex and the center of its group.

        :Hidden parameters:
         - `self._distance_matrix` (np.array) - distance matrix used to create the clustering.
         - `self._clustering` (list) - list giving the resulting clutering.

        :WARNING: `distance_matrix` and `clustering` should obviously ordered the same way!
        """
        sorted_ids_by_clusters, sorted_dist_by_clusters = self.sorted_ids_by_clusters(True)

        # - Tricks to extend the colormap manually:
        if n_colors is None:
            n_colors = len(self._clusters_ids)
        else:
            assert n_colors >=len(self._clusters_ids)

        cmap = plt.get_cmap('jet')
        colors_i = np.concatenate((np.linspace(0, 1., n_colors), (0.,0.,0.,0.)))
        colors_rgba = cmap(colors_i)

        fig = plt.figure(figsize=(10,5))
        for clust_id, distances in sorted_dist_by_clusters.iteritems():
            if cluster_names is None:
                plt.plot( distances, '.-', label = "Cluster "+str(self._clusters_ids[clust_id]), figure=fig, color=tuple(colors_rgba[clust_id]) )
            else:
                plt.plot( distances, '.-', label = str(cluster_names[clust_id]), figure=fig, color=tuple(colors_rgba[clust_id]) )
            if print_clustering_name:
                plt.title(self.clustering_name)
            else:
                plt.title("Vertex distance to their cluster center")

            plt.xlabel("Ranked elements")
            plt.ylabel("Distance to cluster center")
            plt.axis([0,max(self._nb_ids_by_clusters.values()), 0, max(max(sorted_dist_by_clusters.values()))])

        plt.legend(ncol=3, framealpha=0.7, fontsize='small')
        plt.tight_layout()
        if isinstance(savefig,str):
            plt.savefig(savefig)
            plt.close()
        else:
            plt.show()


    def update_cluster_labels(self, old2new_labels, contruct_clustered_graph = True):
        """
        Change the cluster labels according to given info in `old2new_labels`.

        :Parameter:
         - `old2new_labels` (dict) - translation dictionary
        """
        assert isinstance(old2new_labels, dict)
        self.mvpd_matrix._clustering = [old2new_labels[k] for k in self._clustering]

        self.__init__( self.mvpd_matrix, contruct_clustered_graph )
        print "Clustering labels have been udpated !"


    def properties_boxplot_by_cluster(self, cluster_names=None, print_clustering_name=True, savefig=None):
        """
        Display boxplots of properties (used for clustering) by clusters.
        """
        import matplotlib.ticker as mticker
        ppts = [d for n,d in enumerate(self.info_clustering['variables']) if self.info_clustering['weights'][n]!=0]
        if 'topological' in ppts:
            ppts.remove('topological')

        N_ppts = len(ppts)
        fig = plt.figure(figsize=(3.5*N_ppts,5))
        for n,ppt in enumerate(ppts):
            ax = plt.subplot(1,N_ppts,n+1)
            data = [[v for k, v in enumerate(self._var_data_dict[ppt]) if k in self._ids_by_clusters[c]] for c in self._clusters_ids]
            plt.boxplot(data)
            try:
                plt.ylabel(ppt +' ('+ self._var_units_dict[ppt]+')',family='freesans', fontsize=13)
            except:
                plt.ylabel(ppt, fontsize=13)
            ax.set_yticklabels(ax.get_yticks(), fontsize=14)
            y_formatter = mticker.ScalarFormatter(useOffset=False)
            ax.yaxis.set_major_formatter(y_formatter)
            ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,3))
            if cluster_names is not None:
                xtickNames = plt.setp(ax, xticklabels=cluster_names)
                plt.setp(xtickNames, rotation=90, fontsize=9)
                fig.subplots_adjust(bottom=0.2)
            else:
                xtickNames = plt.setp(ax, xticklabels=["{}".format(n) for n in xrange(self._nb_clusters)])
                fig.subplots_adjust(bottom=0.05)
            plt.yticks(fontsize=9)
            if print_clustering_name:
                plt.title(" ")
                plt.suptitle(self.clustering_name, fontsize=14)

        plt.tight_layout()
        if isinstance(savefig,str):
            plt.savefig(savefig)
            plt.close()
        else:
            plt.show()

