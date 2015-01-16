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
__revision__ = " $Id: graph_clusterer.py 17387 2014-08-31 18:28:14Z jlegra02 $ "

import warnings, numpy as np, copy, math
from numpy import ndarray
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import gzip, cPickle as pickle

from openalea.container.temporal_graph_analysis import exist_relative_at_rank
from sklearn.cluster import SpectralClustering, Ward, DBSCAN
from sklearn import metrics
from scipy.sparse import csr_matrix
from openalea.container.temporal_graph_analysis import translate_keys_Image2Graph, add_graph_vertex_property_from_dictionary


def distance_matrix_from_vector(data, variable_types, no_dist_index = []):
    """
    Function creating a distance matrix based on a vector (list) of values.
    Each values are attached to an individual.

    :Parameters:
     - `data` (list) - vector (list) of values
     - `variable_types` (str) - type of variable

    :Returns:
     - `dist_mat` (np.array) - distance matrix
    """
    N = len(data)
    dist_mat = np.zeros( shape = [N,N], dtype=float )

    if variable_types == "Numeric":
        for i in xrange(N):
            for j in xrange(i+1,N):# we skip when i=j because in that case the distance is 0.
                if data[i] is None or data[j] is None or i in no_dist_index or j in no_dist_index:
                    dist_mat[i,j] = dist_mat[j,i] = None
                else:
                    dist_mat[i,j] = dist_mat[j,i] = abs(data[i]-data[j])

    if variable_types == "Ordinal":
        rank = data.argsort() # In case of ordinal variables, observed values are replaced by ranked values.
        for i in xrange(N):
            for j in xrange(i+1,N): # we skip when i=j because in that case the distance is 0.
                if data[i] is None or data[j] is None or i in no_dist_index or j in no_dist_index:
                    dist_mat[i,j] = dist_mat[j,i] = None
                else:
                    dist_mat[i,j] = dist_mat[j,i] = abs(rank[i]-rank[j])

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
    from numpy import ndarray
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
    if isinstance(data,list) or isinstance(data,tuple) or
     (isinstance(data,ndarray) and (((data.shape[0] == 1) and (data.shape[1] > 1)) or ((data.shape[1] == 1) and (data.shape[0] > 1)))):
        if return_data:
            return 'vector', data.tolist()
        else:
            return 'vector'

    if isinstance(data,ndarray) and (data.shape[0]==data.shape[1]) and ((data.shape[0]!=1)and(data.shape[1]!=1)):
        if return_data:
            return 'array', np.array(data)
        else:
            return 'array'


def standardisation(data, norm = 'L1', variable_types = None, outliers_index = []):
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

    if (norm.upper() != 'L1') and (norm.upper() != 'L2'):
        raise ValueError("Undefined standardisation metric")

    # -- Identifying case where numpy.array are 1D vectors:
    if isinstance(data,ndarray) and (((data.shape[0] == 1) and (data.shape[1] > 1)) or ((data.shape[1] == 1) and (data.shape[0] > 1))):
        data.tolist()

    # -- Creating the pairwise distance matrix if not provided in `data`:
    if isinstance(data,list):
        if isinstance(data[0],list) and len(data)==len(data[0]):
            data = np.array(data)
        elif not isinstance(data[0],list):
            distance_matrix = distance_matrix_from_vector(data, variable_types)
        else:
            raise ValueError("Can not convert the provided data.")

    if isinstance(data,ndarray) and (data.shape[0]==data.shape[1]) and ((data.shape[0]!=1)and(data.shape[1]!=1)):
        distance_matrix = data

    # -- We create a copy of the `distance_matrix``as it should also contain the outliers values in the final pairwise distance matrix
    # (the ouliers are indeed removed only for the dispersion/standardisation measure computation).
    dist_mat = copy.copy(distance_matrix)
    
    # -- Handling outliers by setting their value to np.nan (will exclude them of standardisation value computation).
    if outliers_index != []:
        dist_mat[outliers_index,:] = dist_mat[:,outliers_index] = np.nan

    # -- Now we can start the standardisation:
    nan_index = np.isnan(dist_mat)
    if True in nan_index:
        nb_missing_values = len(np.where(nan_index is True)[0])
    else:
        nb_missing_values = 0.

    N = distance_matrix.shape[0]
    if norm.upper() == "L1":
        absd = np.nansum(np.nansum(abs(dist_mat))) / (N*(N-1)-nb_missing_values)
        return distance_matrix / absd

    if norm.upper() == "L2":
        sd = np.nansum(np.nansum(dist_mat**2)) / (N*(N-1)-nb_missing_values)
        return distance_matrix / sd


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


class mvpd_matrix:
    """
    Class allowing to create a Multi-Variates Pairwise Distances Matrix by combining several variables given as vectors of values or (uni-variate) pairwise distance matrices.
    This class is usefull to create an MVPD matrix by combining different types of variables (informations) in a weigthed manner.
    The sum of weight must equal to one!
    We handle missing data by automatically re-weighting according to the presence or absence of data for each cmobined information.
    """
    def __init__(self, variable_list=[], variable_names=[], variable_types=[], variable_weights=[], variable_units=[], standardisation_method='L1'):
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

        # -- Variables saving 'initial' informations about variables:
        self._distance_matrix_dict = {} # save NON-standardized pairwise distance matrix with the varible name as key.
        self._distance_matrix_info = {} 

        # -- Variables for caching information:
        self._global_distance_matrix = None 
        self._global_distance_ids = None
        self._global_distance_weights = None
        self._global_distance_variables = None

        # -- If a list of variable observations, names and type are provided, we compute the NON-standardized pairwise distance matrix
        if self._nb_var != 0:
            for i in xrange(self._nb_var):
                self.add_variable(variable_list[i], variable_names[i], variable_types[i], variable_units[i])
            print 'Done integrating variables {} as separate NON-standardised pairwise distance matrices.'.format(variable_names)

        # -- If the weigths are also provided, we can compute the multi-variate pairwise distance matrix rigth now:
        if len(variable_weights) > 0:
            self.create_mvpd_matrix(variable_names, variable_weights)

        # -- DONE!
        print "mvpd_matrix object initialisation done!"


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
            if vectorORarray(var_data) == 'vector':
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
            self._distance_matrix_dict[var_name] = distance_matrix_from_vector(var_data, var_type)
        else:
            self._distance_matrix_dict[var_name] = var_data

        # - Adding its infos tho the list
        self._distance_matrix_info[var_name] = (var_type.lower(), var_unit)
        self._nb_values = self._distance_matrix_dict.values()[0].shape[0]
        
        if MAD_outliers_thres is not None:
            if vectorORarray(var_data) != 'vector':
                raise TypeError('To detect outliers, on the basis of their MeanAbsDev, you need to provide the data as a vector of values!')
            else:
                self._outliers_index_dict[var_name] = mad_based_outlier(var_data, MAD_outliers_thres)

        return 'Done adding and creating the pairwise distance matrix for variable '{}'.'.format(var_name)


    def remove_pairwise_distance_matrix(self, var_id):
        """
        Remove a pairwise distance matrix form the dictionary `self._distance_matrix_dict` and it's attached information in `self._distance_matrix_info`
        """
        self._distance_matrix_dict.pop(var_id)
        self._distance_matrix_info.pop(var_id)
        return 'Done.'


    def create_mvpd_matrix(self, variable_names, variable_weights, ignore_outliers = False, delete_outliers = False, return_data = False):
        """
        Funtion creating the global weighted distance matrix.
        Provided `variable_names` should exist in `self._distance_matrix_dict`

        :Parameters:
         - `variable_names` (list) - list of variables names to combine; the names should be keys in `self._distance_matrix_dict`;
         - `variable_weights` (list) - list of weightss used to create the global weighted distance matrix (MUST sum to 1!);
         - `ignore_outliers` (bool) - if True ignore outliers when computing standardisation value. Outliers are computed when adding a vertex
         - `delete_outliers` (bool) - if True 'delete outliers' (values set to np.nan) of pairwise distance matrix. Outliers are computed when adding a vertex
         - `return_data` (bool) - if true the function return the sorted ids list `ids` & a dictionary of pairwise distance matrix `global_matrix`
        """
        print "Computing the multi-variate pairwise distance matrix..."
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
            outliers_management = None
        elif ignore_outliers:
            outliers_management = 'ignored'
        else:
            outliers_management = 'deleted'

        # -- Shortcut when asking for the same result:
        if variable_weights == self._global_distance_weights and variable_names == self._global_distance_variables and ids == self._global_distance_ids and self._outliers == outliers_management:
            if return_data:
                return self._global_distance_ids, self._global_distance_matrix
            else:
                print "The global pairwise distance matrix was already computed!"
                return None
        else: # Otherwise, clean any possible clustering:
            self._method = None
            self._nb_clusters = None
            self._clustering = None
            self._full_tree = None
            self._clusterings_dict = None

        # By default we use all the given observations of the variables:
        ids = range(self._nb_values)
        ##SHOULD BE MODIFIED according to the way outliers should be handdled !!

        # -- Standardization step:
        standard_distance_matrix = {}
        mat_topo_dist_standard = []
        nb_var = 0
        for n, var_name in enumerate(variable_names):
            # reduce the matrix to the list of selected ids `ids`
            distance_matrix = self._distance_matrix_dict[var_name][ids,:][:,ids]
            # managing outliers if asked for:
            outliers_index = []
            if (ignore_outliers or delete_outliers):
                try:
                    outliers_index = self._outliers_index_dict[var_name]
                except:
                    print "No outliers defined for property '{}' !".format(var_name)
            if delete_outliers:
                distance_matrix[outliers_index,:] = distance_matrix[:,outliers_index] = np.nan
            # now we can compute the standardised version of the pairwise distance matrix:
            standard_distance_matrix[var_name] = standardisation(distance_matrix, self.standardisation_method, outliers_index)
            nb_var += 1

        # -- Checking for a simple case: to 'combine' only ONE pairwise distance matrix (no re-weighting to do !)
        if nb_var == 1:
            global_matrix = standard_distance_matrix.values()[0]
        else:
            # - Creating weight matrix and replacing nan by zeros for computation in standardized matrix.
            N = len(ids)
            weight_matrix, standardized_matrix = [], [], []
            standardized_matrix.extend( [np.nan_to_num(standard_distance_matrix[var_name]) for var_name in variable_names] )
            for n, var_name in enumerate(variable_names):
                weight_matrix.append(create_weight_matrix(N, variable_weights[n], standard_distance_matrix[var_name]))

            print("Creating the global pairwise weighted standard distance matrix...")
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
        self._outliers = outliers_management

        if return_data:
            return ids, global_matrix
        else:
            print "...done !"


    def cluster(self, n_clusters, method = "ward", ids = None, global_matrix = None, connectivity = None):
        """
        Actually run the clustering method.
        :Parameters:
         - `n_clusters` (int) - number of cluster to create
         - `method` (str) - clustering method to use, "ward", "spectral" and "DBSCAN"
         - `ids` (list) - list of ids
         - `global_matrix` (np.array) - distance matrix to cluster (ordered the same way than `ids`)
        """
        if global_matrix is None:
            if self._global_distance_matrix is not None:
                global_matrix = copy.copy(self._global_distance_matrix)
            else:
                raise ValueError("No distance matrix saved, please give one!")
        # -- Creating the list of vertices:
        if ids is None:
            ids = self._global_distance_ids
        elif ids == 'lineaged':
            ids = [k for k in self._global_distance_ids if k in self.graph.lineaged_vertex(False)]
        elif ids == 'fully lineaged':
            ids = [k for k in self._global_distance_ids if k in self.graph.lineaged_vertex(True)]
        else:
            assert isinstance(ids, list)
        # - Checking for unwanted ids:
        if ids is not None:
            id_not_in_labels = list(set(ids)-set(self.vtx_labels))
            if len(id_not_in_labels) != 0:
                warnings.warn("Some of the ids you provided has not been found in the graph : {}".format(id_not_in_labels))
                print ("Removing them...")
                ids = list(set(ids)-set(id_not_in_labels))
        # - Need to check if there is any changes in the `ids` list compared to the initial list used to create the pairwise distance matrix:
        if (set(self._global_distance_ids)-set(ids)) != set([]):
            ids_index = [self._global_distance_ids.index(v) for v in ids]
            global_matrix = global_matrix[ids_index,:][:,ids_index]

        if n_clusters is None:
            raise ValueError("You have to provide the number of clusters you want for the Ward method.")
        if method.lower() == "ward":
            clustering = Ward(n_clusters = n_clusters, compute_full_tree=True, connectivity=connectivity).fit(global_matrix)
            clustering_labels = clustering.labels_
        if method.lower() == "spectral":
            if connectivity is not None:
                print('Can not use connectivity matrix with Spectral clustering.')
            clustering = SpectralClustering(n_clusters = n_clusters, affinity='precomputed').fit(global_matrix)
            clustering_labels = list(clustering.labels_)

        self._clustering = dict([ (label,clustering_labels[n]) for n,label in enumerate(ids) ])
        self._method = method.lower()
        self._nb_clusters = n_clusters

        return "Done computing {} clustering for {} clusters.".format(self._method, self._nb_clusters)


    def full_tree(self, range_clusters=None, method = "ward", ids = None, global_matrix = None, connectivity = None):
        """
        Actually run the clustering method.
        :Parameters:
         - `range_clusters` (list) (optional) - if provided return a dict of clustering for each cluster number in it
         - `method` (str) (optional) - clustering method to use, "ward", "spectral" and "DBSCAN"
         - `ids` (list) (optional) - list of ids
         - `global_matrix` (np.array) (optional) - distance matrix to cluster (ordered the same way than `ids`)
         - `connectivity` (np.array) (optional) - connectivity matrix to use while clustering (if possible regarding the clustering method) (ordered the same way than `ids`)
        """
        if global_matrix is None:
            if self._global_distance_matrix is not None:
                global_matrix = copy.copy(self._global_distance_matrix)
            else:
                raise ValueError("No distance matrix saved, please give one!")
        # -- Creating the list of vertices:
        if ids is None:
            ids = self._global_distance_ids
        elif ids == 'lineaged':
            ids = [k for k in self._global_distance_ids if k in self.graph.lineaged_vertex(False)]
        elif ids == 'fully lineaged':
            ids = [k for k in self._global_distance_ids if k in self.graph.lineaged_vertex(True)]
        else:
            assert isinstance(ids, list)
        # - Checking for unwanted ids:
        if ids is not None:
            id_not_in_labels = list(set(ids)-set(self.vtx_labels))
            if len(id_not_in_labels) != 0:
                warnings.warn("Some of the ids you provided has not been found in the graph : {}".format(id_not_in_labels))
                print ("Removing them...")
                ids = list(set(ids)-set(id_not_in_labels))
        # - Need to check if there is any changes in the `ids` list compared to the initial list used to create the pairwise distance matrix:
        if (set(self._global_distance_ids)-set(ids)) != set([]):
            ids_index = [self._global_distance_ids.index(v) for v in ids]
            global_matrix = global_matrix[ids_index,:][:,ids_index]

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
                clusterings_dict[c] = dict(zip(ids, clust_labels))
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
        from sklearn import metrics
        from sklearn.cluster import spectral_clustering, Ward
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
        from sklearn import metrics
        from sklearn.cluster import spectral_clustering, Ward
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

