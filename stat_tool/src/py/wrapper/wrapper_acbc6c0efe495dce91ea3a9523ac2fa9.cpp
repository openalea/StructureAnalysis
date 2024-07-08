#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::Dendrogram const volatile * get_pointer<class ::stat_tool::Dendrogram const volatile >(class ::stat_tool::Dendrogram const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_acbc6c0efe495dce91ea3a9523ac2fa9()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::stat_tool::DistanceMatrix * (::stat_tool::Dendrogram::*method_pointer_05a4690bdeb350c696a0f8adad6e2758)() = &::stat_tool::Dendrogram::get_distance_matrix;
    enum ::stat_tool::cluster_scale  (::stat_tool::Dendrogram::*method_pointer_c37315bb315357baa9ef3613b899d13d)() const = &::stat_tool::Dendrogram::get_scale;
    int  (::stat_tool::Dendrogram::*method_pointer_b63ab2c9cca65c0eb3933c295e2683ae)() const = &::stat_tool::Dendrogram::get_nb_cluster;
    int  (::stat_tool::Dendrogram::*method_pointer_5ee4829725c3552798934890b73511c1)(int ) const = &::stat_tool::Dendrogram::get_cluster_nb_pattern;
    int  (::stat_tool::Dendrogram::*method_pointer_25a8ada1965a580c82029a41414a73ea)(int , int ) const = &::stat_tool::Dendrogram::get_cluster_pattern;
    int  (::stat_tool::Dendrogram::*method_pointer_3700ed52828659f0a9e20a1fa6ec0195)(int ) const = &::stat_tool::Dendrogram::get_parent;
    int  (::stat_tool::Dendrogram::*method_pointer_4abc421a81aa59e7b1606226a101c4de)(int , int ) const = &::stat_tool::Dendrogram::get_child;
    double  (::stat_tool::Dendrogram::*method_pointer_8f969e3f51bb5dcc983b7455bdcbad1f)(int ) const = &::stat_tool::Dendrogram::get_child_distance;
    double  (::stat_tool::Dendrogram::*method_pointer_9a45befeeac854ef954f6c33ff6df67c)(int ) const = &::stat_tool::Dendrogram::get_within_cluster_distance;
    double  (::stat_tool::Dendrogram::*method_pointer_591298fbdda054f59d711ed83c7c1c17)(int ) const = &::stat_tool::Dendrogram::get_between_cluster_distance;
    double  (::stat_tool::Dendrogram::*method_pointer_20d69dbdb7e75621a70afde8ad1ea82f)(int ) const = &::stat_tool::Dendrogram::get_max_within_cluster_distance;
    double  (::stat_tool::Dendrogram::*method_pointer_4f6f95a33b3050ecb5f0e1df9ed0383a)(int ) const = &::stat_tool::Dendrogram::get_min_between_cluster_distance;
    boost::python::class_< class ::stat_tool::Dendrogram, autowig::Held< class ::stat_tool::Dendrogram >::Type, boost::python::bases< class ::stat_tool::StatInterface > > class_acbc6c0efe495dce91ea3a9523ac2fa9("Dendrogram", "Hierarchical clustering results\n\n", boost::python::no_init);
    class_acbc6c0efe495dce91ea3a9523ac2fa9.def(boost::python::init<  >(""));
    class_acbc6c0efe495dce91ea3a9523ac2fa9.def(boost::python::init< class ::stat_tool::DistanceMatrix const &, enum ::stat_tool::cluster_scale  >(""));
    class_acbc6c0efe495dce91ea3a9523ac2fa9.def(boost::python::init< class ::stat_tool::Dendrogram const & >(""));
    class_acbc6c0efe495dce91ea3a9523ac2fa9.def("get_distance_matrix", method_pointer_05a4690bdeb350c696a0f8adad6e2758, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_acbc6c0efe495dce91ea3a9523ac2fa9.def("get_scale", method_pointer_c37315bb315357baa9ef3613b899d13d, "");
    class_acbc6c0efe495dce91ea3a9523ac2fa9.def("get_nb_cluster", method_pointer_b63ab2c9cca65c0eb3933c295e2683ae, "");
    class_acbc6c0efe495dce91ea3a9523ac2fa9.def("get_cluster_nb_pattern", method_pointer_5ee4829725c3552798934890b73511c1, "");
    class_acbc6c0efe495dce91ea3a9523ac2fa9.def("get_cluster_pattern", method_pointer_25a8ada1965a580c82029a41414a73ea, "");
    class_acbc6c0efe495dce91ea3a9523ac2fa9.def("get_parent", method_pointer_3700ed52828659f0a9e20a1fa6ec0195, "");
    class_acbc6c0efe495dce91ea3a9523ac2fa9.def("get_child", method_pointer_4abc421a81aa59e7b1606226a101c4de, "");
    class_acbc6c0efe495dce91ea3a9523ac2fa9.def("get_child_distance", method_pointer_8f969e3f51bb5dcc983b7455bdcbad1f, "");
    class_acbc6c0efe495dce91ea3a9523ac2fa9.def("get_within_cluster_distance", method_pointer_9a45befeeac854ef954f6c33ff6df67c, "");
    class_acbc6c0efe495dce91ea3a9523ac2fa9.def("get_between_cluster_distance", method_pointer_591298fbdda054f59d711ed83c7c1c17, "");
    class_acbc6c0efe495dce91ea3a9523ac2fa9.def("get_max_within_cluster_distance", method_pointer_20d69dbdb7e75621a70afde8ad1ea82f, "");
    class_acbc6c0efe495dce91ea3a9523ac2fa9.def("get_min_between_cluster_distance", method_pointer_4f6f95a33b3050ecb5f0e1df9ed0383a, "");

    if(autowig::Held< class ::stat_tool::Dendrogram >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::Dendrogram >::Type, autowig::Held< class ::stat_tool::StatInterface >::Type >();
    }

}