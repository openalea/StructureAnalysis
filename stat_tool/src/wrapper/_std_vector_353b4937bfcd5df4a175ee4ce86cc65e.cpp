#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_vector_353b4937bfcd5df4a175ee4ce86cc65e()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_1e960f79d6585bc38acb49374308d0c7)(unsigned long, class ::stat_tool::MultiPlot const &) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::assign;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_76875c9b3be0579c995d5b0ee76b1c08)(class ::std::initializer_list<stat_tool::MultiPlot>) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::assign;
        class ::__gnu_cxx::__normal_iterator<stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_8059c159071652fb9cdccefc0ba6d391)() = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::begin;
        class ::__gnu_cxx::__normal_iterator<const stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_bc4eac8133b15068b668327d84a951a6)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::begin;
        class ::__gnu_cxx::__normal_iterator<stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_9dc69e78cfdb5f779dd853ac92dd75f5)() = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::end;
        class ::__gnu_cxx::__normal_iterator<const stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_5f0f7a3db95652469f5cc55f4686083f)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::end;
        class ::std::reverse_iterator<__gnu_cxx::__normal_iterator<stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > > (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_d93258add7da541fbc00d162701508d8)() = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::rbegin;
        class ::std::reverse_iterator<__gnu_cxx::__normal_iterator<const stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > > (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_660d403f896555e1ad6354d186162488)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::rbegin;
        class ::std::reverse_iterator<__gnu_cxx::__normal_iterator<stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > > (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_015221517c585664bfc7487ed2d2ee2e)() = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::rend;
        class ::std::reverse_iterator<__gnu_cxx::__normal_iterator<const stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > > (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_d0323ce6cefb50ec8d053d396b4b3083)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::rend;
        class ::__gnu_cxx::__normal_iterator<const stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_2dc64eb8c9715ae6bafcd5b74fc23f02)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::cbegin;
        class ::__gnu_cxx::__normal_iterator<const stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_5da50d465a915a1bb2331a30314b0108)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::cend;
        class ::std::reverse_iterator<__gnu_cxx::__normal_iterator<const stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > > (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_9bcdb356d0155f7ab4588cf3f067b664)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::crbegin;
        class ::std::reverse_iterator<__gnu_cxx::__normal_iterator<const stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > > (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_bbedac2192ea58c2a06a2b7862368874)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::crend;
        unsigned long (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_b27a21b3a1b359ad9e1525ba95889172)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::size;
        unsigned long (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_bfbaf1f6ef5f5de0894d99e9cb08d1e8)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::max_size;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_d3d1d6bcd36552b3a9d00a12c7fb5914)(unsigned long) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::resize;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_2716e006579e5128b03fb3359850bb6f)(unsigned long, class ::stat_tool::MultiPlot const &) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::resize;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_787a19fa0c825310b116e3af10c0d1e1)() = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::shrink_to_fit;
        unsigned long (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_578e2742f84d5f72ab2df2174f633fca)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::capacity;
        bool (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_a608e3fca0bb505fb88bb9ecf9f98e83)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::empty;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_9cb2aa9a2b0a551a81a050b90ba90aba)(unsigned long) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::reserve;
        class ::stat_tool::MultiPlot & (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_b33e92976e59532b99563f931c194509)(unsigned long) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::operator[];
        class ::stat_tool::MultiPlot const & (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_fe27993b65a854c987b983ceff076a8a)(unsigned long) const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::operator[];
        class ::stat_tool::MultiPlot & (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_948799718e2b5295a98fef05cc84c6b9)(unsigned long) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::at;
        class ::stat_tool::MultiPlot const & (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_a2a8a118e903577ba230286a13e4a516)(unsigned long) const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::at;
        class ::stat_tool::MultiPlot & (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_f59ba0fe4c39568c8e96bc9199f317ca)() = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::front;
        class ::stat_tool::MultiPlot const & (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_312d887aa45e5a9186391111122428ca)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::front;
        class ::stat_tool::MultiPlot & (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_be9d3b9c32c85cf8b39f812f835ad2b4)() = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::back;
        class ::stat_tool::MultiPlot const & (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_42114a1d3c9d5702857cb80e08e89ade)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::back;
        class ::stat_tool::MultiPlot * (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_2c3088b9e89d521da686542b25d2069d)() = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::data;
        class ::stat_tool::MultiPlot const * (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_8f520dbcc16658198b3da9bbae7c71ad)() const = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::data;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_42da52736f6b52c3a8874835f3062487)(class ::stat_tool::MultiPlot const &) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::push_back;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_7022ca3e7c5557f089da3a139313c76a)() = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::pop_back;
        class ::__gnu_cxx::__normal_iterator<stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_9edccd27c5bb5ce9b90e9dc9947619e2)(class ::__gnu_cxx::__normal_iterator<const stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > >, class ::stat_tool::MultiPlot const &) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::insert;
        class ::__gnu_cxx::__normal_iterator<stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_e1fb6f6b8b96535798bbb7df31a1f874)(class ::__gnu_cxx::__normal_iterator<const stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > >, class ::std::initializer_list<stat_tool::MultiPlot>) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::insert;
        class ::__gnu_cxx::__normal_iterator<stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_46454bdf2ac550eeacf22ef913eb94fb)(class ::__gnu_cxx::__normal_iterator<const stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > >, unsigned long, class ::stat_tool::MultiPlot const &) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::insert;
        class ::__gnu_cxx::__normal_iterator<stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_9112e2f9e15c5a19990c92d1ccf96c0c)(class ::__gnu_cxx::__normal_iterator<const stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > >) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::erase;
        class ::__gnu_cxx::__normal_iterator<stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_b3d2931b21e55770a2494f64fd93a90e)(class ::__gnu_cxx::__normal_iterator<const stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > >, class ::__gnu_cxx::__normal_iterator<const stat_tool::MultiPlot *, std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > >) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::erase;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_3cd5f3b75e0f57bebbecbe69f20456a4)(class ::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > &) = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::swap;
        void (::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::*method_pointer_3468843386c55add8cf2668640bb3df0)() = &::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >::clear;
        boost::python::class_< class ::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> >, std::shared_ptr< class ::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > > >("_Vector_353b4937bfcd5df4a175ee4ce86cc65e", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::std::allocator<stat_tool::MultiPlot> const & >())
            .def(boost::python::init< unsigned long, class ::std::allocator<stat_tool::MultiPlot> const & >())
            .def(boost::python::init< unsigned long, class ::stat_tool::MultiPlot const &, class ::std::allocator<stat_tool::MultiPlot> const & >())
            .def(boost::python::init< class ::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > const & >())
            .def(boost::python::init< class ::std::vector<stat_tool::MultiPlot, std::allocator<stat_tool::MultiPlot> > const &, class ::std::allocator<stat_tool::MultiPlot> const & >())
            .def(boost::python::init< class ::std::initializer_list<stat_tool::MultiPlot>, class ::std::allocator<stat_tool::MultiPlot> const & >())
            .def("assign", method_pointer_1e960f79d6585bc38acb49374308d0c7)
            .def("assign", method_pointer_76875c9b3be0579c995d5b0ee76b1c08)
            .def("begin", method_pointer_8059c159071652fb9cdccefc0ba6d391)
            .def("begin", method_pointer_bc4eac8133b15068b668327d84a951a6)
            .def("end", method_pointer_9dc69e78cfdb5f779dd853ac92dd75f5)
            .def("end", method_pointer_5f0f7a3db95652469f5cc55f4686083f)
            .def("rbegin", method_pointer_d93258add7da541fbc00d162701508d8)
            .def("rbegin", method_pointer_660d403f896555e1ad6354d186162488)
            .def("rend", method_pointer_015221517c585664bfc7487ed2d2ee2e)
            .def("rend", method_pointer_d0323ce6cefb50ec8d053d396b4b3083)
            .def("cbegin", method_pointer_2dc64eb8c9715ae6bafcd5b74fc23f02)
            .def("cend", method_pointer_5da50d465a915a1bb2331a30314b0108)
            .def("crbegin", method_pointer_9bcdb356d0155f7ab4588cf3f067b664)
            .def("crend", method_pointer_bbedac2192ea58c2a06a2b7862368874)
            .def("size", method_pointer_b27a21b3a1b359ad9e1525ba95889172)
            .def("max_size", method_pointer_bfbaf1f6ef5f5de0894d99e9cb08d1e8)
            .def("resize", method_pointer_d3d1d6bcd36552b3a9d00a12c7fb5914)
            .def("resize", method_pointer_2716e006579e5128b03fb3359850bb6f)
            .def("shrink_to_fit", method_pointer_787a19fa0c825310b116e3af10c0d1e1)
            .def("capacity", method_pointer_578e2742f84d5f72ab2df2174f633fca)
            .def("empty", method_pointer_a608e3fca0bb505fb88bb9ecf9f98e83)
            .def("reserve", method_pointer_9cb2aa9a2b0a551a81a050b90ba90aba)
            .def("__setitem__", method_pointer_b33e92976e59532b99563f931c194509, boost::python::return_internal_reference<>())
            .def("__getitem__", method_pointer_fe27993b65a854c987b983ceff076a8a, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("at", method_pointer_948799718e2b5295a98fef05cc84c6b9, boost::python::return_internal_reference<>())
            .def("at", method_pointer_a2a8a118e903577ba230286a13e4a516, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("front", method_pointer_f59ba0fe4c39568c8e96bc9199f317ca, boost::python::return_internal_reference<>())
            .def("front", method_pointer_312d887aa45e5a9186391111122428ca, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("back", method_pointer_be9d3b9c32c85cf8b39f812f835ad2b4, boost::python::return_internal_reference<>())
            .def("back", method_pointer_42114a1d3c9d5702857cb80e08e89ade, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("data", method_pointer_2c3088b9e89d521da686542b25d2069d, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("data", method_pointer_8f520dbcc16658198b3da9bbae7c71ad, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("push_back", method_pointer_42da52736f6b52c3a8874835f3062487)
            .def("pop_back", method_pointer_7022ca3e7c5557f089da3a139313c76a)
            .def("insert", method_pointer_9edccd27c5bb5ce9b90e9dc9947619e2)
            .def("insert", method_pointer_e1fb6f6b8b96535798bbb7df31a1f874)
            .def("insert", method_pointer_46454bdf2ac550eeacf22ef913eb94fb)
            .def("erase", method_pointer_9112e2f9e15c5a19990c92d1ccf96c0c)
            .def("erase", method_pointer_b3d2931b21e55770a2494f64fd93a90e)
            .def("swap", method_pointer_3cd5f3b75e0f57bebbecbe69f20456a4)
            .def("clear", method_pointer_3468843386c55add8cf2668640bb3df0);
}