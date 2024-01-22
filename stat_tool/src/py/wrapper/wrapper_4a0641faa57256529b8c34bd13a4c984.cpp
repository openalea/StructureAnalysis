#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::ChainReestimation< int > const volatile * get_pointer<class ::stat_tool::ChainReestimation< int > const volatile >(class ::stat_tool::ChainReestimation< int > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_4a0641faa57256529b8c34bd13a4c984()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::ChainReestimation< int >::*method_pointer_3dbcd9b08f305ec1b045caa7f09f5049)() = &::stat_tool::ChainReestimation< int >::init;
    void  (::stat_tool::ChainReestimation< int >::*method_pointer_1d3bc5d3ddf85900be092a22c8302c66)() = &::stat_tool::ChainReestimation< int >::remove;
    boost::python::class_< class ::stat_tool::ChainReestimation< int >, autowig::Held< class ::stat_tool::ChainReestimation< int > >::Type > class_4a0641faa57256529b8c34bd13a4c984("_ChainReestimation_4a0641faa57256529b8c34bd13a4c984", "", boost::python::no_init);
    class_4a0641faa57256529b8c34bd13a4c984.def(boost::python::init< enum ::stat_tool::process_type , int , int , bool  >(""));
    class_4a0641faa57256529b8c34bd13a4c984.def("init", method_pointer_3dbcd9b08f305ec1b045caa7f09f5049, "");
    class_4a0641faa57256529b8c34bd13a4c984.def("remove", method_pointer_1d3bc5d3ddf85900be092a22c8302c66, "");

}