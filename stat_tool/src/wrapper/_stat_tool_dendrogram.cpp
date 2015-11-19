#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>
#include <stat_tool/distance_matrix.h>

void _stat_tool_dendrogram()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::stat_tool::DistanceMatrix * (::stat_tool::Dendrogram::*method_pointer_05a4690bdeb350c696a0f8adad6e2758)() = &::stat_tool::Dendrogram::get_distance_matrix;
        enum ::stat_tool::cluster_scale (::stat_tool::Dendrogram::*method_pointer_c37315bb315357baa9ef3613b899d13d)() const = &::stat_tool::Dendrogram::get_scale;
        int (::stat_tool::Dendrogram::*method_pointer_b63ab2c9cca65c0eb3933c295e2683ae)() const = &::stat_tool::Dendrogram::get_nb_cluster;
        int (::stat_tool::Dendrogram::*method_pointer_5ee4829725c3552798934890b73511c1)(int) const = &::stat_tool::Dendrogram::get_cluster_nb_pattern;
        int (::stat_tool::Dendrogram::*method_pointer_25a8ada1965a580c82029a41414a73ea)(int, int) const = &::stat_tool::Dendrogram::get_cluster_pattern;
        int (::stat_tool::Dendrogram::*method_pointer_3700ed52828659f0a9e20a1fa6ec0195)(int) const = &::stat_tool::Dendrogram::get_parent;
        int (::stat_tool::Dendrogram::*method_pointer_4abc421a81aa59e7b1606226a101c4de)(int, int) const = &::stat_tool::Dendrogram::get_child;
        double (::stat_tool::Dendrogram::*method_pointer_8f969e3f51bb5dcc983b7455bdcbad1f)(int) const = &::stat_tool::Dendrogram::get_child_distance;
        double (::stat_tool::Dendrogram::*method_pointer_9a45befeeac854ef954f6c33ff6df67c)(int) const = &::stat_tool::Dendrogram::get_within_cluster_distance;
        double (::stat_tool::Dendrogram::*method_pointer_591298fbdda054f59d711ed83c7c1c17)(int) const = &::stat_tool::Dendrogram::get_between_cluster_distance;
        double (::stat_tool::Dendrogram::*method_pointer_20d69dbdb7e75621a70afde8ad1ea82f)(int) const = &::stat_tool::Dendrogram::get_max_within_cluster_distance;
        double (::stat_tool::Dendrogram::*method_pointer_4f6f95a33b3050ecb5f0e1df9ed0383a)(int) const = &::stat_tool::Dendrogram::get_min_between_cluster_distance;
        boost::python::class_< class ::stat_tool::Dendrogram, std::shared_ptr< class ::stat_tool::Dendrogram >, boost::python::bases< class ::stat_tool::StatInterface > >("Dendrogram", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::stat_tool::DistanceMatrix const &, enum ::stat_tool::cluster_scale >())
            .def(boost::python::init< class ::stat_tool::Dendrogram const & >())
            .def("get_distance_matrix", method_pointer_05a4690bdeb350c696a0f8adad6e2758, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_scale", method_pointer_c37315bb315357baa9ef3613b899d13d)
            .def("get_nb_cluster", method_pointer_b63ab2c9cca65c0eb3933c295e2683ae)
            .def("get_cluster_nb_pattern", method_pointer_5ee4829725c3552798934890b73511c1)
            .def("get_cluster_pattern", method_pointer_25a8ada1965a580c82029a41414a73ea)
            .def("get_parent", method_pointer_3700ed52828659f0a9e20a1fa6ec0195)
            .def("get_child", method_pointer_4abc421a81aa59e7b1606226a101c4de)
            .def("get_child_distance", method_pointer_8f969e3f51bb5dcc983b7455bdcbad1f)
            .def("get_within_cluster_distance", method_pointer_9a45befeeac854ef954f6c33ff6df67c)
            .def("get_between_cluster_distance", method_pointer_591298fbdda054f59d711ed83c7c1c17)
            .def("get_max_within_cluster_distance", method_pointer_20d69dbdb7e75621a70afde8ad1ea82f)
            .def("get_min_between_cluster_distance", method_pointer_4f6f95a33b3050ecb5f0e1df9ed0383a);
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::Dendrogram >, std::shared_ptr< class ::stat_tool::StatInterface > >();
}