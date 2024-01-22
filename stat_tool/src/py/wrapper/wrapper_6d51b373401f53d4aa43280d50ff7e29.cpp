#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::Regression const volatile * get_pointer<class ::stat_tool::Regression const volatile >(class ::stat_tool::Regression const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_6d51b373401f53d4aa43280d50ff7e29()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::stat_tool::Vectors * (::stat_tool::Regression::*method_pointer_f0f8e6da2e0c5362a32a10014a11cb45)() const = &::stat_tool::Regression::get_vectors;
    int  (::stat_tool::Regression::*method_pointer_ceb101baacb35a3aad7ed11d70736b08)() const = &::stat_tool::Regression::get_nb_vector;
    double  (::stat_tool::Regression::*method_pointer_dded4e359e5850748ba155cde8aa6856)(int ) const = &::stat_tool::Regression::get_residual;
    boost::python::class_< class ::stat_tool::Regression, autowig::Held< class ::stat_tool::Regression >::Type, boost::python::bases< class ::stat_tool::StatInterface, class ::stat_tool::RegressionKernel > > class_6d51b373401f53d4aa43280d50ff7e29("Regression", "Regression function\n\n", boost::python::no_init);
    class_6d51b373401f53d4aa43280d50ff7e29.def(boost::python::init<  >(""));
    class_6d51b373401f53d4aa43280d50ff7e29.def(boost::python::init< enum ::stat_tool::parametric_function , int , int , class ::stat_tool::Vectors const & >(""));
    class_6d51b373401f53d4aa43280d50ff7e29.def(boost::python::init< class ::stat_tool::Regression const & >(""));
    class_6d51b373401f53d4aa43280d50ff7e29.def("get_vectors", method_pointer_f0f8e6da2e0c5362a32a10014a11cb45, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_6d51b373401f53d4aa43280d50ff7e29.def("get_nb_vector", method_pointer_ceb101baacb35a3aad7ed11d70736b08, "");
    class_6d51b373401f53d4aa43280d50ff7e29.def("get_residual", method_pointer_dded4e359e5850748ba155cde8aa6856, "");

    if(autowig::Held< class ::stat_tool::Regression >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::Regression >::Type, autowig::Held< class ::stat_tool::StatInterface >::Type >();
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::Regression >::Type, autowig::Held< class ::stat_tool::RegressionKernel >::Type >();
    }

}