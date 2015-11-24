#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _mbstate_t()
{

        boost::python::class_< ::__mbstate_t, std::shared_ptr< ::__mbstate_t > >("MbstateT", boost::python::no_init);
}