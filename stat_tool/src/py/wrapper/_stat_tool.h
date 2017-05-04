#ifndef AUTOWIG__STAT_TOOL
#define AUTOWIG__STAT_TOOL

#include <boost/python.hpp>
#include <type_traits>
#include <stat_tool/distribution_reestimation.hpp>
#include <stat_tool/stat_label.h>
#include <stat_tool/markovian.h>
#include <stat_tool/chain_reestimation.h>
#include <stat_tool/discrete_mixture.h>
#include <stat_tool/vectors.h>
#include <stat_tool/reestimation.h>
#include <stat_tool/convolution.h>
#include <stat_tool/quantile_computation.hpp>
#include <stat_tool/reestimation.hpp>
#include <stat_tool/continuous_parametric_estimation.hpp>
#include <stat_tool/mixture.h>
#include <stat_tool/distance_matrix.h>
#include <stat_tool/plotable.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/distribution.h>
#include <stat_tool/plotable.h>
#include <stat_tool/chain_reestimation.hpp>
#include <stat_tool/stat_tools.h>
#include <stat_tool/curves.h>
#include <stat_tool/regression.h>

namespace autowig
{
     template<class T> struct Held {
        typedef T* Type;
        static bool const is_class = false;
    };
}

/*namespace autowig
{
    class Wrap_57e8416c0462569e8916250db94d4d5a : public ::stat_tool::StatInterface, public boost::python::wrapper< class ::stat_tool::StatInterface >
    {
        public:
            
            virtual bool  plot_write(class ::stat_tool::StatError & param_0, char const * param_1, char const * param_2) const
            { return this->get_override("plot_write")(param_0, param_1, param_2); }
            virtual bool  spreadsheet_write(class ::stat_tool::StatError & param_0, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const param_1) const
            { return this->get_override("spreadsheet_write")(param_0, param_1); }
            virtual bool  ascii_write(class ::stat_tool::StatError & param_0, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const param_1, bool  param_2) const
            { return this->get_override("ascii_write")(param_0, param_1, param_2); }
            virtual class ::std::basic_ostream< char, struct ::std::char_traits< char > > & ascii_write(class ::std::basic_ostream< char, struct ::std::char_traits< char > > & param_0, bool  param_1) const
            { return this->get_override("ascii_write")(param_0, param_1); }
            virtual class ::std::basic_ostream< char, struct ::std::char_traits< char > > & line_write(class ::std::basic_ostream< char, struct ::std::char_traits< char > > & param_0) const
            { return this->get_override("line_write")(param_0); }

        protected:
            

        private:
            

    };
}*/
#endif