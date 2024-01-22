#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::ChainData const volatile * get_pointer<class ::stat_tool::ChainData const volatile >(class ::stat_tool::ChainData const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_8470eec064f35f569821a3873bd503cd()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    int  (::stat_tool::ChainData::*method_pointer_4f13feefcea2566cb8fdaa172ed5ab53)() const = &::stat_tool::ChainData::nb_parameter_computation;
    void  (::stat_tool::ChainData::*method_pointer_576a1605ed7f5cf090c80db6f04828bf)(class ::stat_tool::Chain &) const = &::stat_tool::ChainData::estimation;
    boost::python::class_< class ::stat_tool::ChainData, autowig::Held< class ::stat_tool::ChainData >::Type, boost::python::bases< class ::stat_tool::ChainReestimation< int > > > class_8470eec064f35f569821a3873bd503cd("ChainData", "Data structure corresponding to a Markov chain\n\n", boost::python::no_init);
    class_8470eec064f35f569821a3873bd503cd.def(boost::python::init< enum ::stat_tool::process_type , int , int , bool  >(""));
    class_8470eec064f35f569821a3873bd503cd.def(boost::python::init< class ::stat_tool::ChainData const & >(""));
    class_8470eec064f35f569821a3873bd503cd.def("nb_parameter_computation", method_pointer_4f13feefcea2566cb8fdaa172ed5ab53, "");
    class_8470eec064f35f569821a3873bd503cd.def("estimation", method_pointer_576a1605ed7f5cf090c80db6f04828bf, "");

    if(autowig::Held< class ::stat_tool::ChainData >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::ChainData >::Type, autowig::Held< class ::stat_tool::ChainReestimation< int > >::Type >();
    }

}