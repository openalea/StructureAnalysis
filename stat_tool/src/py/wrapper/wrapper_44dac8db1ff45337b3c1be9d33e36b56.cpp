#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::DistanceMatrix const volatile * get_pointer<class ::stat_tool::DistanceMatrix const volatile >(class ::stat_tool::DistanceMatrix const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_44dac8db1ff45337b3c1be9d33e36b56()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::stat_tool::DistanceMatrix * (::stat_tool::DistanceMatrix::*method_pointer_94512dd122a5535a9bd5e0926db1ee2d)(class ::stat_tool::StatError &, int , class ::std::vector< int, class ::std::allocator< int > > &, bool ) const = &::stat_tool::DistanceMatrix::select_individual;
    class ::stat_tool::DistanceMatrix * (::stat_tool::DistanceMatrix::*method_pointer_1d604885d9985299bd1f9b51dd023a91)(class ::stat_tool::StatError &) const = &::stat_tool::DistanceMatrix::symmetrize;
    class ::stat_tool::DistanceMatrix * (::stat_tool::DistanceMatrix::*method_pointer_e35d0838984b58eaa9ae519e30b67285)(class ::stat_tool::StatError &) const = &::stat_tool::DistanceMatrix::unnormalize;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::DistanceMatrix::*method_pointer_a92baa767dcd543298b28203a271abfa)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::DistanceMatrix::spreadsheet_write;
    bool  (::stat_tool::DistanceMatrix::*method_pointer_7b0a0590fbf553b7a09c3fbad02a8977)() const = &::stat_tool::DistanceMatrix::test_symmetry;
    void  (::stat_tool::DistanceMatrix::*method_pointer_3ffddd5defe953c8a5e3739abd7fd36a)(int , int , double , int , double , int , double , int , int , double , int , double , int ) = &::stat_tool::DistanceMatrix::update;
    void  (::stat_tool::DistanceMatrix::*method_pointer_09477275018c56d9ba6e229798e23b4e)(int , int , double , int ) = &::stat_tool::DistanceMatrix::update;
    class ::stat_tool::Dendrogram * (::stat_tool::DistanceMatrix::*method_pointer_66cb399f4c14583285a1271f9d559d7a)(enum ::stat_tool::hierarchical_strategy , enum ::stat_tool::linkage ) const = &::stat_tool::DistanceMatrix::agglomerative_hierarchical_clustering;
    class ::stat_tool::Dendrogram * (::stat_tool::DistanceMatrix::*method_pointer_c67f9204769859059a3b4645a1e5e9fc)() const = &::stat_tool::DistanceMatrix::divisive_hierarchical_clustering;
    bool  (::stat_tool::DistanceMatrix::*method_pointer_4622fdb2f29c5d65bb6728496f12a0fe)(class ::stat_tool::StatError &, bool , enum ::stat_tool::hierarchical_strategy , enum ::stat_tool::linkage , class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const, enum ::stat_tool::output_format ) const = &::stat_tool::DistanceMatrix::hierarchical_clustering;
    int  (::stat_tool::DistanceMatrix::*method_pointer_d161f244602a5b298026fb73f5b551df)() const = &::stat_tool::DistanceMatrix::get_nb_row;
    int  (::stat_tool::DistanceMatrix::*method_pointer_041820e3e60b5e0fbc7752e9bb7e6df0)() const = &::stat_tool::DistanceMatrix::get_nb_column;
    int  (::stat_tool::DistanceMatrix::*method_pointer_e5791fadaf98521793555cc5daec34c8)(int ) const = &::stat_tool::DistanceMatrix::get_row_identifier;
    int  (::stat_tool::DistanceMatrix::*method_pointer_23b018bd666f5a7aa18b4d1c5d7664e5)(int ) const = &::stat_tool::DistanceMatrix::get_column_identifier;
    double  (::stat_tool::DistanceMatrix::*method_pointer_a7a2fe5ef8f45a7e8390c15ba9beb01a)(int , int ) const = &::stat_tool::DistanceMatrix::get_distance;
    int  (::stat_tool::DistanceMatrix::*method_pointer_fa6b06e35d8c54d68b3adf4639f03059)(int , int ) const = &::stat_tool::DistanceMatrix::get_length;
    double  (::stat_tool::DistanceMatrix::*method_pointer_5349aedc446b5d518f24236b7d9249d9)(int , int ) const = &::stat_tool::DistanceMatrix::get_deletion_distance;
    int  (::stat_tool::DistanceMatrix::*method_pointer_e8f2840ac1865120853f42a83005fd76)(int , int ) const = &::stat_tool::DistanceMatrix::get_nb_deletion;
    double  (::stat_tool::DistanceMatrix::*method_pointer_f42c6ce805f95d9f8910e1e206c4bb1b)(int , int ) const = &::stat_tool::DistanceMatrix::get_insertion_distance;
    int  (::stat_tool::DistanceMatrix::*method_pointer_a1900738e606518c80126b5c8ea6db1c)(int , int ) const = &::stat_tool::DistanceMatrix::get_nb_insertion;
    int  (::stat_tool::DistanceMatrix::*method_pointer_5daf6ed78a085b47b0aa017698bef83b)(int , int ) const = &::stat_tool::DistanceMatrix::get_nb_match;
    double  (::stat_tool::DistanceMatrix::*method_pointer_235553c9f49a59bab65bb5c30642a358)(int , int ) const = &::stat_tool::DistanceMatrix::get_substitution_distance;
    int  (::stat_tool::DistanceMatrix::*method_pointer_9d4d705c8f5a518c860e5381dd7e6e1a)(int , int ) const = &::stat_tool::DistanceMatrix::get_nb_substitution;
    double  (::stat_tool::DistanceMatrix::*method_pointer_19395ee5a86456a4997ad475e9fa9a75)(int , int ) const = &::stat_tool::DistanceMatrix::get_transposition_distance;
    int  (::stat_tool::DistanceMatrix::*method_pointer_c1b986330e9753c99e851439c9c9e36d)(int , int ) const = &::stat_tool::DistanceMatrix::get_nb_transposition;
    int  (::stat_tool::DistanceMatrix::*method_pointer_599a5c6d226b513a8bff477874d05e29)() const = &::stat_tool::DistanceMatrix::get_label_size;
    bool  (::stat_tool::DistanceMatrix::*method_pointer_a9eacd07104557f79e066bb01dfd1632)() = &::stat_tool::DistanceMatrix::is_deletion;
    bool  (::stat_tool::DistanceMatrix::*method_pointer_e09ee9e598c25864b63bf9bd0600e075)() = &::stat_tool::DistanceMatrix::is_insertion;
    bool  (::stat_tool::DistanceMatrix::*method_pointer_d660715d707e599ebdcbc5fd61e33ab7)() = &::stat_tool::DistanceMatrix::is_match;
    bool  (::stat_tool::DistanceMatrix::*method_pointer_6771d0fd57555f84afba4c238fe29bcb)() = &::stat_tool::DistanceMatrix::is_substitution;
    bool  (::stat_tool::DistanceMatrix::*method_pointer_c771a762e34954c989927235bcef6981)() = &::stat_tool::DistanceMatrix::is_transposition;
    boost::python::class_< class ::stat_tool::DistanceMatrix, autowig::Held< class ::stat_tool::DistanceMatrix >::Type, boost::python::bases< class ::stat_tool::StatInterface > > class_44dac8db1ff45337b3c1be9d33e36b56("DistanceMatrix", "Distance matrix\n\n", boost::python::no_init);
    class_44dac8db1ff45337b3c1be9d33e36b56.def(boost::python::init<  >(""));
    class_44dac8db1ff45337b3c1be9d33e36b56.def(boost::python::init< class ::stat_tool::DistanceMatrix const &, enum ::stat_tool::matrix_transform  >(""));
    class_44dac8db1ff45337b3c1be9d33e36b56.def("select_individual", method_pointer_94512dd122a5535a9bd5e0926db1ee2d, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("symmetrize", method_pointer_1d604885d9985299bd1f9b51dd023a91, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("unnormalize", method_pointer_e35d0838984b58eaa9ae519e30b67285, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("spreadsheet_write", method_pointer_a92baa767dcd543298b28203a271abfa, boost::python::return_internal_reference<>(), "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("test_symmetry", method_pointer_7b0a0590fbf553b7a09c3fbad02a8977, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("update", method_pointer_3ffddd5defe953c8a5e3739abd7fd36a, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("update", method_pointer_09477275018c56d9ba6e229798e23b4e, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("agglomerative_hierarchical_clustering", method_pointer_66cb399f4c14583285a1271f9d559d7a, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("divisive_hierarchical_clustering", method_pointer_c67f9204769859059a3b4645a1e5e9fc, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("hierarchical_clustering", method_pointer_4622fdb2f29c5d65bb6728496f12a0fe, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("get_nb_row", method_pointer_d161f244602a5b298026fb73f5b551df, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("get_nb_column", method_pointer_041820e3e60b5e0fbc7752e9bb7e6df0, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("get_row_identifier", method_pointer_e5791fadaf98521793555cc5daec34c8, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("get_column_identifier", method_pointer_23b018bd666f5a7aa18b4d1c5d7664e5, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("get_distance", method_pointer_a7a2fe5ef8f45a7e8390c15ba9beb01a, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("get_length", method_pointer_fa6b06e35d8c54d68b3adf4639f03059, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("get_deletion_distance", method_pointer_5349aedc446b5d518f24236b7d9249d9, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("get_nb_deletion", method_pointer_e8f2840ac1865120853f42a83005fd76, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("get_insertion_distance", method_pointer_f42c6ce805f95d9f8910e1e206c4bb1b, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("get_nb_insertion", method_pointer_a1900738e606518c80126b5c8ea6db1c, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("get_nb_match", method_pointer_5daf6ed78a085b47b0aa017698bef83b, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("get_substitution_distance", method_pointer_235553c9f49a59bab65bb5c30642a358, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("get_nb_substitution", method_pointer_9d4d705c8f5a518c860e5381dd7e6e1a, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("get_transposition_distance", method_pointer_19395ee5a86456a4997ad475e9fa9a75, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("get_nb_transposition", method_pointer_c1b986330e9753c99e851439c9c9e36d, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("get_label_size", method_pointer_599a5c6d226b513a8bff477874d05e29, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("is_deletion", method_pointer_a9eacd07104557f79e066bb01dfd1632, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("is_insertion", method_pointer_e09ee9e598c25864b63bf9bd0600e075, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("is_match", method_pointer_d660715d707e599ebdcbc5fd61e33ab7, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("is_substitution", method_pointer_6771d0fd57555f84afba4c238fe29bcb, "");
    class_44dac8db1ff45337b3c1be9d33e36b56.def("is_transposition", method_pointer_c771a762e34954c989927235bcef6981, "");

    if(autowig::Held< class ::stat_tool::DistanceMatrix >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::DistanceMatrix >::Type, autowig::Held< class ::stat_tool::StatInterface >::Type >();
    }

}