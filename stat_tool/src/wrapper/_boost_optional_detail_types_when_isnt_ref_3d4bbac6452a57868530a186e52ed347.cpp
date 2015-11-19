#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _boost_optional_detail_types_when_isnt_ref_3d4bbac6452a57868530a186e52ed347()
{
        std::string boost_21ee8db290f35815a57c7bf74dca851d_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".boost");
        boost::python::object boost_21ee8db290f35815a57c7bf74dca851d_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(boost_21ee8db290f35815a57c7bf74dca851d_name.c_str()))));
        boost::python::scope().attr("boost") = boost_21ee8db290f35815a57c7bf74dca851d_module;
        boost::python::scope boost_21ee8db290f35815a57c7bf74dca851d_scope = boost_21ee8db290f35815a57c7bf74dca851d_module;        std::string optional_detail_8b9e69d15fe354e087a1a4561a039920_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".optional_detail");
        boost::python::object optional_detail_8b9e69d15fe354e087a1a4561a039920_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(optional_detail_8b9e69d15fe354e087a1a4561a039920_name.c_str()))));
        boost::python::scope().attr("optional_detail") = optional_detail_8b9e69d15fe354e087a1a4561a039920_module;
        boost::python::scope optional_detail_8b9e69d15fe354e087a1a4561a039920_scope = optional_detail_8b9e69d15fe354e087a1a4561a039920_module;
        boost::python::class_< struct ::boost::optional_detail::types_when_isnt_ref<std::locale>, std::shared_ptr< struct ::boost::optional_detail::types_when_isnt_ref<std::locale> > >("_TypesWhenIsntRef_3d4bbac6452a57868530a186e52ed347", boost::python::no_init);
}