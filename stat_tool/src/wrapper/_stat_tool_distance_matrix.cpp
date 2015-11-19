#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>
#include <stat_tool/distance_matrix.h>

void _stat_tool_distance_matrix()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::stat_tool::DistanceMatrix * (::stat_tool::DistanceMatrix::*method_pointer_1d604885d9985299bd1f9b51dd023a91)(class ::stat_tool::StatError &) const = &::stat_tool::DistanceMatrix::symmetrize;
        class ::stat_tool::DistanceMatrix * (::stat_tool::DistanceMatrix::*method_pointer_e35d0838984b58eaa9ae519e30b67285)(class ::stat_tool::StatError &) const = &::stat_tool::DistanceMatrix::unnormalize;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::DistanceMatrix::*method_pointer_e23319161b02558eb4a68ba86cf918de)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::DistanceMatrix::spreadsheet_write;
        bool (::stat_tool::DistanceMatrix::*method_pointer_7b0a0590fbf553b7a09c3fbad02a8977)() const = &::stat_tool::DistanceMatrix::test_symmetry;
        void (::stat_tool::DistanceMatrix::*method_pointer_3ffddd5defe953c8a5e3739abd7fd36a)(int, int, double, int, double, int, double, int, int, double, int, double, int) = &::stat_tool::DistanceMatrix::update;
        void (::stat_tool::DistanceMatrix::*method_pointer_09477275018c56d9ba6e229798e23b4e)(int, int, double, int) = &::stat_tool::DistanceMatrix::update;
        class ::stat_tool::Dendrogram * (::stat_tool::DistanceMatrix::*method_pointer_66cb399f4c14583285a1271f9d559d7a)(enum ::stat_tool::hierarchical_strategy, enum ::stat_tool::linkage) const = &::stat_tool::DistanceMatrix::agglomerative_hierarchical_clustering;
        class ::stat_tool::Dendrogram * (::stat_tool::DistanceMatrix::*method_pointer_c67f9204769859059a3b4645a1e5e9fc)() const = &::stat_tool::DistanceMatrix::divisive_hierarchical_clustering;
        int (::stat_tool::DistanceMatrix::*method_pointer_d161f244602a5b298026fb73f5b551df)() const = &::stat_tool::DistanceMatrix::get_nb_row;
        int (::stat_tool::DistanceMatrix::*method_pointer_041820e3e60b5e0fbc7752e9bb7e6df0)() const = &::stat_tool::DistanceMatrix::get_nb_column;
        int (::stat_tool::DistanceMatrix::*method_pointer_e5791fadaf98521793555cc5daec34c8)(int) const = &::stat_tool::DistanceMatrix::get_row_identifier;
        int (::stat_tool::DistanceMatrix::*method_pointer_23b018bd666f5a7aa18b4d1c5d7664e5)(int) const = &::stat_tool::DistanceMatrix::get_column_identifier;
        double (::stat_tool::DistanceMatrix::*method_pointer_a7a2fe5ef8f45a7e8390c15ba9beb01a)(int, int) const = &::stat_tool::DistanceMatrix::get_distance;
        int (::stat_tool::DistanceMatrix::*method_pointer_fa6b06e35d8c54d68b3adf4639f03059)(int, int) const = &::stat_tool::DistanceMatrix::get_length;
        double (::stat_tool::DistanceMatrix::*method_pointer_5349aedc446b5d518f24236b7d9249d9)(int, int) const = &::stat_tool::DistanceMatrix::get_deletion_distance;
        int (::stat_tool::DistanceMatrix::*method_pointer_e8f2840ac1865120853f42a83005fd76)(int, int) const = &::stat_tool::DistanceMatrix::get_nb_deletion;
        double (::stat_tool::DistanceMatrix::*method_pointer_f42c6ce805f95d9f8910e1e206c4bb1b)(int, int) const = &::stat_tool::DistanceMatrix::get_insertion_distance;
        int (::stat_tool::DistanceMatrix::*method_pointer_a1900738e606518c80126b5c8ea6db1c)(int, int) const = &::stat_tool::DistanceMatrix::get_nb_insertion;
        int (::stat_tool::DistanceMatrix::*method_pointer_5daf6ed78a085b47b0aa017698bef83b)(int, int) const = &::stat_tool::DistanceMatrix::get_nb_match;
        double (::stat_tool::DistanceMatrix::*method_pointer_235553c9f49a59bab65bb5c30642a358)(int, int) const = &::stat_tool::DistanceMatrix::get_substitution_distance;
        int (::stat_tool::DistanceMatrix::*method_pointer_9d4d705c8f5a518c860e5381dd7e6e1a)(int, int) const = &::stat_tool::DistanceMatrix::get_nb_substitution;
        double (::stat_tool::DistanceMatrix::*method_pointer_19395ee5a86456a4997ad475e9fa9a75)(int, int) const = &::stat_tool::DistanceMatrix::get_transposition_distance;
        int (::stat_tool::DistanceMatrix::*method_pointer_c1b986330e9753c99e851439c9c9e36d)(int, int) const = &::stat_tool::DistanceMatrix::get_nb_transposition;
        int (::stat_tool::DistanceMatrix::*method_pointer_599a5c6d226b513a8bff477874d05e29)() const = &::stat_tool::DistanceMatrix::get_label_size;
        bool (::stat_tool::DistanceMatrix::*method_pointer_a9eacd07104557f79e066bb01dfd1632)() = &::stat_tool::DistanceMatrix::is_deletion;
        bool (::stat_tool::DistanceMatrix::*method_pointer_e09ee9e598c25864b63bf9bd0600e075)() = &::stat_tool::DistanceMatrix::is_insertion;
        bool (::stat_tool::DistanceMatrix::*method_pointer_d660715d707e599ebdcbc5fd61e33ab7)() = &::stat_tool::DistanceMatrix::is_match;
        bool (::stat_tool::DistanceMatrix::*method_pointer_6771d0fd57555f84afba4c238fe29bcb)() = &::stat_tool::DistanceMatrix::is_substitution;
        bool (::stat_tool::DistanceMatrix::*method_pointer_c771a762e34954c989927235bcef6981)() = &::stat_tool::DistanceMatrix::is_transposition;
        boost::python::class_< class ::stat_tool::DistanceMatrix, std::shared_ptr< class ::stat_tool::DistanceMatrix >, boost::python::bases< class ::stat_tool::StatInterface > >("DistanceMatrix", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::stat_tool::DistanceMatrix const &, char >())
            .def("symmetrize", method_pointer_1d604885d9985299bd1f9b51dd023a91, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("unnormalize", method_pointer_e35d0838984b58eaa9ae519e30b67285, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("spreadsheet_write", method_pointer_e23319161b02558eb4a68ba86cf918de, boost::python::return_internal_reference<>())
            .def("test_symmetry", method_pointer_7b0a0590fbf553b7a09c3fbad02a8977)
            .def("update", method_pointer_3ffddd5defe953c8a5e3739abd7fd36a)
            .def("update", method_pointer_09477275018c56d9ba6e229798e23b4e)
            .def("agglomerative_hierarchical_clustering", method_pointer_66cb399f4c14583285a1271f9d559d7a, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("divisive_hierarchical_clustering", method_pointer_c67f9204769859059a3b4645a1e5e9fc, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_nb_row", method_pointer_d161f244602a5b298026fb73f5b551df)
            .def("get_nb_column", method_pointer_041820e3e60b5e0fbc7752e9bb7e6df0)
            .def("get_row_identifier", method_pointer_e5791fadaf98521793555cc5daec34c8)
            .def("get_column_identifier", method_pointer_23b018bd666f5a7aa18b4d1c5d7664e5)
            .def("get_distance", method_pointer_a7a2fe5ef8f45a7e8390c15ba9beb01a)
            .def("get_length", method_pointer_fa6b06e35d8c54d68b3adf4639f03059)
            .def("get_deletion_distance", method_pointer_5349aedc446b5d518f24236b7d9249d9)
            .def("get_nb_deletion", method_pointer_e8f2840ac1865120853f42a83005fd76)
            .def("get_insertion_distance", method_pointer_f42c6ce805f95d9f8910e1e206c4bb1b)
            .def("get_nb_insertion", method_pointer_a1900738e606518c80126b5c8ea6db1c)
            .def("get_nb_match", method_pointer_5daf6ed78a085b47b0aa017698bef83b)
            .def("get_substitution_distance", method_pointer_235553c9f49a59bab65bb5c30642a358)
            .def("get_nb_substitution", method_pointer_9d4d705c8f5a518c860e5381dd7e6e1a)
            .def("get_transposition_distance", method_pointer_19395ee5a86456a4997ad475e9fa9a75)
            .def("get_nb_transposition", method_pointer_c1b986330e9753c99e851439c9c9e36d)
            .def("get_label_size", method_pointer_599a5c6d226b513a8bff477874d05e29)
            .def("is_deletion", method_pointer_a9eacd07104557f79e066bb01dfd1632)
            .def("is_insertion", method_pointer_e09ee9e598c25864b63bf9bd0600e075)
            .def("is_match", method_pointer_d660715d707e599ebdcbc5fd61e33ab7)
            .def("is_substitution", method_pointer_6771d0fd57555f84afba4c238fe29bcb)
            .def("is_transposition", method_pointer_c771a762e34954c989927235bcef6981);
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::DistanceMatrix >, std::shared_ptr< class ::stat_tool::StatInterface > >();
}