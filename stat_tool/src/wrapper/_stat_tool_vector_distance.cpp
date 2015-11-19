#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>
#include <stat_tool/discrete_mixture.h>
#include <stat_tool/vectors.h>

void _stat_tool_vector_distance()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        void (::stat_tool::VectorDistance::*method_pointer_c48cb58c82df5fd0845f4da056f382ed)(int, double) const = &::stat_tool::VectorDistance::dispersion_update;
        int (::stat_tool::VectorDistance::*method_pointer_f522981ef750530a886bbf3e448c387b)() const = &::stat_tool::VectorDistance::get_nb_variable;
        enum ::stat_tool::metric (::stat_tool::VectorDistance::*method_pointer_8366a2e9ef805d4e98c43bcf3593ddfa)() const = &::stat_tool::VectorDistance::get_distance_type;
        enum ::stat_tool::variable_type (::stat_tool::VectorDistance::*method_pointer_74ca394b8dd25ed5a5e109512c3edd88)(int) const = &::stat_tool::VectorDistance::get_var_type;
        double (::stat_tool::VectorDistance::*method_pointer_4a8076132c285fa6800f866ed97a0fd3)(int) const = &::stat_tool::VectorDistance::get_weight;
        double (::stat_tool::VectorDistance::*method_pointer_99c788dc725954bbba868c6f3b38ad3a)(int) const = &::stat_tool::VectorDistance::get_dispersion;
        int (::stat_tool::VectorDistance::*method_pointer_a1b01637456f5084b89fff4122e32f10)(int) const = &::stat_tool::VectorDistance::get_nb_value;
        double (::stat_tool::VectorDistance::*method_pointer_610596ad20295a5cac0e85a6352d7c9e)(int, int, int) const = &::stat_tool::VectorDistance::get_category_distance;
        int (::stat_tool::VectorDistance::*method_pointer_97cf13043c2f5e40babb6f83fc923071)(int) const = &::stat_tool::VectorDistance::get_period;
        boost::python::class_< class ::stat_tool::VectorDistance, std::shared_ptr< class ::stat_tool::VectorDistance >, boost::python::bases< class ::stat_tool::StatInterface > >("VectorDistance", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::stat_tool::VectorDistance const & >())
            .def("dispersion_update", method_pointer_c48cb58c82df5fd0845f4da056f382ed)
            .def("get_nb_variable", method_pointer_f522981ef750530a886bbf3e448c387b)
            .def("get_distance_type", method_pointer_8366a2e9ef805d4e98c43bcf3593ddfa)
            .def("get_var_type", method_pointer_74ca394b8dd25ed5a5e109512c3edd88)
            .def("get_weight", method_pointer_4a8076132c285fa6800f866ed97a0fd3)
            .def("get_dispersion", method_pointer_99c788dc725954bbba868c6f3b38ad3a)
            .def("get_nb_value", method_pointer_a1b01637456f5084b89fff4122e32f10)
            .def("get_category_distance", method_pointer_610596ad20295a5cac0e85a6352d7c9e)
            .def("get_period", method_pointer_97cf13043c2f5e40babb6f83fc923071);
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::VectorDistance >, std::shared_ptr< class ::stat_tool::StatInterface > >();
}