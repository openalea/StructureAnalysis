#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>
#include <stat_tool/distance_matrix.h>

void _stat_tool_clusters()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::Clusters::*method_pointer_87add017690852d8bd2e555b142d5a09)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::Clusters::line_write;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::Clusters::*method_pointer_0a4671320bed551993a0769295d337a2)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool) const = &::stat_tool::Clusters::ascii_write;
        class ::stat_tool::MultiPlotSet * (::stat_tool::Clusters::*method_pointer_08d7f1c713815e52bb438c2daa011025)() const = &::stat_tool::Clusters::get_plotable;
        void (::stat_tool::Clusters::*method_pointer_bd5dc2ffdea2582897b9a29b13576dcf)() = &::stat_tool::Clusters::cluster_nb_pattern_computation;
        void (::stat_tool::Clusters::*method_pointer_f505de8e1c6a54ba9ede9c30f2f924b9)() = &::stat_tool::Clusters::pattern_distance_computation;
        void (::stat_tool::Clusters::*method_pointer_2d73e8b3a8ee536aab910acfc6be5f04)() = &::stat_tool::Clusters::cluster_distance_computation_1;
        void (::stat_tool::Clusters::*method_pointer_376fb47d457e5d4992064731cea86f9e)() = &::stat_tool::Clusters::cluster_distance_computation_2;
        class ::stat_tool::DistanceMatrix * (::stat_tool::Clusters::*method_pointer_5eb1a0b4af13514590f13802a9c6315e)() = &::stat_tool::Clusters::get_distance_matrix;
        int (::stat_tool::Clusters::*method_pointer_e8891968c97153348b2c67e5bb360b05)() const = &::stat_tool::Clusters::get_nb_pattern;
        int (::stat_tool::Clusters::*method_pointer_b7454a9815015c828f3d580714a2e873)() const = &::stat_tool::Clusters::get_nb_cluster;
        int (::stat_tool::Clusters::*method_pointer_074256310ab7527bbebbd928bf843e2e)(int) const = &::stat_tool::Clusters::get_cluster_nb_pattern;
        int (::stat_tool::Clusters::*method_pointer_595f8519a4355964993ca33b41bf9d5f)(int) const = &::stat_tool::Clusters::get_assignment;
        double (::stat_tool::Clusters::*method_pointer_26d1380771fb525e9d731ebe2b533dfb)(int, int) const = &::stat_tool::Clusters::get_pattern_distance;
        int (::stat_tool::Clusters::*method_pointer_31d24ecbdd745c40bc8084cb406ca6be)(int, int) const = &::stat_tool::Clusters::get_pattern_length;
        boost::python::class_< class ::stat_tool::Clusters, std::shared_ptr< class ::stat_tool::Clusters >, boost::python::bases< class ::stat_tool::DistanceMatrix > >("Clusters", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::stat_tool::DistanceMatrix const &, int >())
            .def(boost::python::init< class ::stat_tool::Clusters const & >())
            .def("line_write", method_pointer_87add017690852d8bd2e555b142d5a09, boost::python::return_internal_reference<>())
            .def("ascii_write", method_pointer_0a4671320bed551993a0769295d337a2, boost::python::return_internal_reference<>())
            .def("get_plotable", method_pointer_08d7f1c713815e52bb438c2daa011025, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("cluster_nb_pattern_computation", method_pointer_bd5dc2ffdea2582897b9a29b13576dcf)
            .def("pattern_distance_computation", method_pointer_f505de8e1c6a54ba9ede9c30f2f924b9)
            .def("cluster_distance_computation_1", method_pointer_2d73e8b3a8ee536aab910acfc6be5f04)
            .def("cluster_distance_computation_2", method_pointer_376fb47d457e5d4992064731cea86f9e)
            .def("get_distance_matrix", method_pointer_5eb1a0b4af13514590f13802a9c6315e, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_nb_pattern", method_pointer_e8891968c97153348b2c67e5bb360b05)
            .def("get_nb_cluster", method_pointer_b7454a9815015c828f3d580714a2e873)
            .def("get_cluster_nb_pattern", method_pointer_074256310ab7527bbebbd928bf843e2e)
            .def("get_assignment", method_pointer_595f8519a4355964993ca33b41bf9d5f)
            .def("get_pattern_distance", method_pointer_26d1380771fb525e9d731ebe2b533dfb)
            .def("get_pattern_length", method_pointer_31d24ecbdd745c40bc8084cb406ca6be);
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::Clusters >, std::shared_ptr< class ::stat_tool::DistanceMatrix > >();
}