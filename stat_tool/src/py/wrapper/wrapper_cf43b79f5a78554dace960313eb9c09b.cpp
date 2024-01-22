#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::Clusters const volatile * get_pointer<class ::stat_tool::Clusters const volatile >(class ::stat_tool::Clusters const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_cf43b79f5a78554dace960313eb9c09b()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::Clusters::*method_pointer_bd5dc2ffdea2582897b9a29b13576dcf)() = &::stat_tool::Clusters::cluster_nb_pattern_computation;
    void  (::stat_tool::Clusters::*method_pointer_f505de8e1c6a54ba9ede9c30f2f924b9)() = &::stat_tool::Clusters::pattern_distance_computation;
    void  (::stat_tool::Clusters::*method_pointer_2d73e8b3a8ee536aab910acfc6be5f04)() = &::stat_tool::Clusters::cluster_distance_computation_1;
    void  (::stat_tool::Clusters::*method_pointer_376fb47d457e5d4992064731cea86f9e)() = &::stat_tool::Clusters::cluster_distance_computation_2;
    class ::stat_tool::DistanceMatrix * (::stat_tool::Clusters::*method_pointer_5eb1a0b4af13514590f13802a9c6315e)() = &::stat_tool::Clusters::get_distance_matrix;
    int  (::stat_tool::Clusters::*method_pointer_e8891968c97153348b2c67e5bb360b05)() const = &::stat_tool::Clusters::get_nb_pattern;
    int  (::stat_tool::Clusters::*method_pointer_b7454a9815015c828f3d580714a2e873)() const = &::stat_tool::Clusters::get_nb_cluster;
    int  (::stat_tool::Clusters::*method_pointer_074256310ab7527bbebbd928bf843e2e)(int ) const = &::stat_tool::Clusters::get_cluster_nb_pattern;
    int  (::stat_tool::Clusters::*method_pointer_595f8519a4355964993ca33b41bf9d5f)(int ) const = &::stat_tool::Clusters::get_assignment;
    double  (::stat_tool::Clusters::*method_pointer_26d1380771fb525e9d731ebe2b533dfb)(int , int ) const = &::stat_tool::Clusters::get_pattern_distance;
    int  (::stat_tool::Clusters::*method_pointer_31d24ecbdd745c40bc8084cb406ca6be)(int , int ) const = &::stat_tool::Clusters::get_pattern_length;
    boost::python::class_< class ::stat_tool::Clusters, autowig::Held< class ::stat_tool::Clusters >::Type, boost::python::bases< class ::stat_tool::DistanceMatrix > > class_cf43b79f5a78554dace960313eb9c09b("Clusters", "Partitioning clustering results\n\n", boost::python::no_init);
    class_cf43b79f5a78554dace960313eb9c09b.def(boost::python::init<  >(""));
    class_cf43b79f5a78554dace960313eb9c09b.def(boost::python::init< class ::stat_tool::DistanceMatrix const &, int  >(""));
    class_cf43b79f5a78554dace960313eb9c09b.def(boost::python::init< class ::stat_tool::Clusters const & >(""));
    class_cf43b79f5a78554dace960313eb9c09b.def("cluster_nb_pattern_computation", method_pointer_bd5dc2ffdea2582897b9a29b13576dcf, "");
    class_cf43b79f5a78554dace960313eb9c09b.def("pattern_distance_computation", method_pointer_f505de8e1c6a54ba9ede9c30f2f924b9, "");
    class_cf43b79f5a78554dace960313eb9c09b.def("cluster_distance_computation_1", method_pointer_2d73e8b3a8ee536aab910acfc6be5f04, "");
    class_cf43b79f5a78554dace960313eb9c09b.def("cluster_distance_computation_2", method_pointer_376fb47d457e5d4992064731cea86f9e, "");
    class_cf43b79f5a78554dace960313eb9c09b.def("get_distance_matrix", method_pointer_5eb1a0b4af13514590f13802a9c6315e, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_cf43b79f5a78554dace960313eb9c09b.def("get_nb_pattern", method_pointer_e8891968c97153348b2c67e5bb360b05, "");
    class_cf43b79f5a78554dace960313eb9c09b.def("get_nb_cluster", method_pointer_b7454a9815015c828f3d580714a2e873, "");
    class_cf43b79f5a78554dace960313eb9c09b.def("get_cluster_nb_pattern", method_pointer_074256310ab7527bbebbd928bf843e2e, "");
    class_cf43b79f5a78554dace960313eb9c09b.def("get_assignment", method_pointer_595f8519a4355964993ca33b41bf9d5f, "");
    class_cf43b79f5a78554dace960313eb9c09b.def("get_pattern_distance", method_pointer_26d1380771fb525e9d731ebe2b533dfb, "");
    class_cf43b79f5a78554dace960313eb9c09b.def("get_pattern_length", method_pointer_31d24ecbdd745c40bc8084cb406ca6be, "");

    if(autowig::Held< class ::stat_tool::Clusters >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::Clusters >::Type, autowig::Held< class ::stat_tool::DistanceMatrix >::Type >();
    }

}