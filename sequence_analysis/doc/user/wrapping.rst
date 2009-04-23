########
WRAPPING
########


Tops
====

    friend Tops* tops_ascii_read(Format_error &error , const char *path , bool old_format);

    Tops();
    Tops(int nb_top , int *iidentifier , int *nb_position , bool init_flag = true);
    Tops(const Sequences &seq);
    Tops(const Tops &tops , int inb_sequence , int *index);
    Tops(const Tops &tops , bool model_flag = true , bool reverse_flag = false);
    Tops(int nb_sample , const Tops **ptops);
    ~Tops();
    Tops& operator=(const Tops &tops);

    Distribution_data* extract(Format_error &error , int position) const;

    Tops* shift(Format_error &error , int inb_internode) const;
    Tops* select_individual(Format_error &error , int inb_sequence ,
                            int *iidentifier , bool keep = true) const;
    Tops* reverse(Format_error &error) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_data_write(std::ostream &os , char format = 'c' , bool exhaustive = false) const;
    bool ascii_data_write(Format_error &error , const char *path ,
                          char format = 'c' , bool exhaustive = false) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;


    void build_nb_internode_histogram();

    Top_parameters* estimation(Format_error &error , int imin_position , int imax_position ,
                               int neighborhood , bool equal_probability = false) const;
    Top_parameters* estimation(Format_error &error , int neighborhood , bool equal_probability = false) const
    { return estimation(error , 1 , max_position , neighborhood , equal_probability); }

    // acces membres de la classe

    Top_parameters* get_top_parameters() const { return top_parameters; }
    Histogram* get_nb_internode() const { return nb_internode; }
    int get_max_position() const { return max_position; }
    Histogram* get_axillary_nb_internode(int position) const
    { return axillary_nb_internode[position]; }

Sequences wrapping
##################


constructors NOT done
=====================


*    Sequences(const Histogram &ihlength , int inb_variable , bool init_flag = false);
*    Sequences(const Renewal_data &timev);
*    Sequences(const Sequences &seq , int inb_sequence , int *index);
*    Sequences(const Sequences &seq , bool *segment_mean);
*    Sequences(const Sequences &seq , char transform = 'c' , int param = DEFAULT);


Public method not DONE
======================

*    Vectors* build_vectors(bool index_variable) const;
*    Vectors* extract_vectors(Format_error &error , int feature_type , int variable = I_DEFAULT , int value = I_DEFAULT) const;
*    Markovian_sequences* markovian_sequences(Format_error &error) const;
*    Tops* tops(Format_error &error) const;
*    bool check(Format_error &error , const char *pattern_label);
*    Time_events* extract_time_events(Format_error &error , int variable ,int begin_date , int end_date ,  int previous_date = I_DEFAULT , int next_date = I_DEFAULT) const;
*    Renewal_data* extract_renewal_data(Format_error &error , int variable ,  int begin_index , int end_index) const;
*    std::ostream& line_write(std::ostream &os) const;


*int min_index_parameter_computation() const;
*    int max_index_parameter_computation(bool last_position = false) const;

*   void marginal_histogram_computation(int variable);
*    double mean_computation(int variable) const;
*    double variance_computation(int variable , double mean) const;
*    double mean_absolute_deviation_computation(int variable , double mean) const;
*    double mean_absolute_difference_computation(int variable) const;
*    double skewness_computation(int variable , double mean , double variance) const;
*    double kurtosis_computation(int variable , double mean , double variance) const;
*    Histogram* value_index_interval_computation(Format_error &error , int variable , int value) const;
*    Correlation* correlation_computation(Format_error &error , int variable1 , int variable2 , int itype = PEARSON , int max_lag = I_DEFAULT , int normalization = EXACT) const;
*    Correlation* partial_autocorrelation_computation(Format_error &error , int variable ,  int itype = PEARSON , int max_lag = I_DEFAULT) const;

*    Distance_matrix* alignment(Format_error &error , std::ostream *os , const Vector_distance &ivector_dist ,
                               int ref_identifier = I_DEFAULT , int test_identifier = I_DEFAULT ,
                               bool begin_free = false , bool end_free = false , int indel_cost = ADAPTATIVE ,
                               double indel_factor = INDEL_FACTOR_1 , bool transposition_flag = false ,
                               double transposition_factor = TRANSPOSITION_FACTOR ,
                               const char *result_path = 0 , char result_format = 'a' ,
                               const char *alignment_path = 0 , char alignment_format = 'a') const;

*    Distance_matrix* alignment(Format_error &error , std::ostream *os , int ref_identifier = I_DEFAULT ,
                               int test_identifier = I_DEFAULT , bool begin_free = false , bool end_free = false ,
                               const char *result_path = 0 , char result_format = 'a' ,
                               const char *alignment_path = 0 , char alignment_format = 'a') const;

*    Sequences* multiple_alignment(Format_error &error , std::ostream &os , const Vector_distance &ivector_dist ,
                                  bool begin_free = false , bool end_free = false , int indel_cost = ADAPTATIVE ,
                                  double indel_factor = INDEL_FACTOR_N , int algorithm = AGGLOMERATIVE ,
                                  const char *path = 0) const;

*    Sequences* segmentation(Format_error &error , std::ostream &os , int iidentifier ,
                            int nb_segment , int *ichange_point , int *model_type ,
                            int output = SEQUENCE) const;
*    Sequences* segmentation(Format_error &error , std::ostream &os , int *nb_segment ,int *model_type , int iidentifier = I_DEFAULT ,  int output = SEQUENCE) const;
*    Sequences* segmentation(Format_error &error , std::ostream &os , int iidentifier , int max_nb_segment , int *model_type) const;

*    Sequences* hierarchical_segmentation(Format_error &error , std::ostream &os , int iidentifier , int max_nb_segment , int *model_type) const;

*    Sequences* segmentation(Format_error &error , int iidentifier , int nb_segment , const Vector_distance &ivector_dist , std::ostream &os , int output = SEGMENT) const;

*    bool segment_profile_write(Format_error &error , std::ostream &os , int iidentifier ,
                               int nb_segment , int *model_type , int output = SEGMENT ,
                               char format = 'a' , int segmentation = FORWARD_DYNAMIC_PROGRAMMING ,
                               int nb_segmentation = NB_SEGMENTATION) const;
*    bool segment_profile_write(Format_error &error , const char *path , int iidentifier ,
                               int nb_segment , int *model_type , int output = SEGMENT ,
                               char format = 'a' , int segmentation = FORWARD_DYNAMIC_PROGRAMMING ,
                               int nb_segmentation = NB_SEGMENTATION) const;
*    bool segment_profile_plot_write(Format_error &error , const char *prefix ,
                                    int iidentifier , int nb_segment , int *model_type ,
                                    int output = SEGMENT , const char *title = 0) const;

    // acces membres de la classe

*    Histogram* get_hlength() const { return hlength; }
*    Histogram* get_hindex_parameter() const { return hindex_parameter; }
*    Histogram* get_index_interval() const { return index_interval; }
*    int get_index_parameter(int iseq , int index) const    { return index_parameter[iseq][index]; }
*    Histogram* get_marginal(int variable) const { return marginal[variable]; }


