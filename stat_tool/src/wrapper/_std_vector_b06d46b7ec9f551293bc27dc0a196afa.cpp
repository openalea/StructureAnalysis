#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_vector_b06d46b7ec9f551293bc27dc0a196afa()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_208fbea201e25b8884bf6c5e64906edf)(unsigned long, enum ::stat_tool::process_distribution const &) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::assign;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_1322ab8f67355d2a828de626b055bf81)(class ::std::initializer_list<stat_tool::process_distribution>) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::assign;
        class ::__gnu_cxx::__normal_iterator<stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_7dfbb96ba70b572fb120198ac05314e3)() = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::begin;
        class ::__gnu_cxx::__normal_iterator<const stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_c264f8bde7dd58698339ce6e4145ab15)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::begin;
        class ::__gnu_cxx::__normal_iterator<stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_cfb19c6829735503b9746a91ae297ead)() = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::end;
        class ::__gnu_cxx::__normal_iterator<const stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_9aeef1267d3651d08c0e83559b8139b1)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::end;
        class ::std::reverse_iterator<__gnu_cxx::__normal_iterator<stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > > (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_5c782e385dec5db8886f47dd7f0b93bb)() = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::rbegin;
        class ::std::reverse_iterator<__gnu_cxx::__normal_iterator<const stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > > (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_cd7e286989085631bb4350334fdc19e4)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::rbegin;
        class ::std::reverse_iterator<__gnu_cxx::__normal_iterator<stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > > (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_c9cc74ed7fe153799a6a7cb43161c309)() = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::rend;
        class ::std::reverse_iterator<__gnu_cxx::__normal_iterator<const stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > > (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_d05389da916f5cb0b2ee2c27de13357d)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::rend;
        class ::__gnu_cxx::__normal_iterator<const stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_d0b72163dd185c51aaab3d15e80c67a7)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::cbegin;
        class ::__gnu_cxx::__normal_iterator<const stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_d34205575084545c87346e6022bb46f9)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::cend;
        class ::std::reverse_iterator<__gnu_cxx::__normal_iterator<const stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > > (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_cb139ac60d245079a0410998c8695bbb)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::crbegin;
        class ::std::reverse_iterator<__gnu_cxx::__normal_iterator<const stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > > (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_121616a191925aa0a5ac7a0aded6c260)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::crend;
        unsigned long (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_2976530ccb5657b6b6c0a54e0f5c9bf8)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::size;
        unsigned long (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_1aacab760834510ea8d4153fc713ba92)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::max_size;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_70ee84578153504299edef13598e052d)(unsigned long) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::resize;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_bd273cd8cf6052bab101f81c2e8e2af8)(unsigned long, enum ::stat_tool::process_distribution const &) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::resize;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_7e1a87abab0e55769d56b92f87de9364)() = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::shrink_to_fit;
        unsigned long (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_c4a00f13846e5f839d297a023ef32936)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::capacity;
        bool (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_ca8611fc69aa580da7d1319be7a1b2af)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::empty;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_a51bd8eeb45656289869d33118b75250)(unsigned long) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::reserve;
        enum ::stat_tool::process_distribution const & (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_fdeddab37b135a37b33c6270959703f7)(unsigned long) const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::operator[];
        enum ::stat_tool::process_distribution const & (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_563d298b97c95b388324b19f0f74ab68)(unsigned long) const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::at;
        enum ::stat_tool::process_distribution const & (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_73ae975684b25694bbeee2d13496bdba)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::front;
        enum ::stat_tool::process_distribution const & (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_6acc6f10ad2e54568819075f729606fa)() const = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::back;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_dd51d9a7dfbd547ea2f1392263c449cb)(enum ::stat_tool::process_distribution const &) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::push_back;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_2d363ed2fc5950f68a03d934d57851cc)() = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::pop_back;
        class ::__gnu_cxx::__normal_iterator<stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_162d49c09a54599f91a9a46de35d3d14)(class ::__gnu_cxx::__normal_iterator<const stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > >, enum ::stat_tool::process_distribution const &) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::insert;
        class ::__gnu_cxx::__normal_iterator<stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_abdb6a9576965e2db5d88d95340257bd)(class ::__gnu_cxx::__normal_iterator<const stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > >, class ::std::initializer_list<stat_tool::process_distribution>) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::insert;
        class ::__gnu_cxx::__normal_iterator<stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_3416c40ef9925177aafabc4732e41eb1)(class ::__gnu_cxx::__normal_iterator<const stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > >, unsigned long, enum ::stat_tool::process_distribution const &) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::insert;
        class ::__gnu_cxx::__normal_iterator<stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_761bc5712055510b8c735d5f0d3cef9f)(class ::__gnu_cxx::__normal_iterator<const stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > >) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::erase;
        class ::__gnu_cxx::__normal_iterator<stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_f7f15843e71f5e08ac9de17d0717e8f4)(class ::__gnu_cxx::__normal_iterator<const stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > >, class ::__gnu_cxx::__normal_iterator<const stat_tool::process_distribution *, std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > >) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::erase;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_f617e6b977e651a6ae8985dd17f99bda)(class ::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > &) = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::swap;
        void (::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::*method_pointer_823a942d3fe95fb9bee873042664f0b9)() = &::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >::clear;
        boost::python::class_< class ::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> >, std::shared_ptr< class ::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > > >("_Vector_b06d46b7ec9f551293bc27dc0a196afa", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::std::allocator<stat_tool::process_distribution> const & >())
            .def(boost::python::init< unsigned long, class ::std::allocator<stat_tool::process_distribution> const & >())
            .def(boost::python::init< unsigned long, enum ::stat_tool::process_distribution const &, class ::std::allocator<stat_tool::process_distribution> const & >())
            .def(boost::python::init< class ::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > const & >())
            .def(boost::python::init< class ::std::vector<stat_tool::process_distribution, std::allocator<stat_tool::process_distribution> > const &, class ::std::allocator<stat_tool::process_distribution> const & >())
            .def(boost::python::init< class ::std::initializer_list<stat_tool::process_distribution>, class ::std::allocator<stat_tool::process_distribution> const & >())
            .def("assign", method_pointer_208fbea201e25b8884bf6c5e64906edf)
            .def("assign", method_pointer_1322ab8f67355d2a828de626b055bf81)
            .def("begin", method_pointer_7dfbb96ba70b572fb120198ac05314e3)
            .def("begin", method_pointer_c264f8bde7dd58698339ce6e4145ab15)
            .def("end", method_pointer_cfb19c6829735503b9746a91ae297ead)
            .def("end", method_pointer_9aeef1267d3651d08c0e83559b8139b1)
            .def("rbegin", method_pointer_5c782e385dec5db8886f47dd7f0b93bb)
            .def("rbegin", method_pointer_cd7e286989085631bb4350334fdc19e4)
            .def("rend", method_pointer_c9cc74ed7fe153799a6a7cb43161c309)
            .def("rend", method_pointer_d05389da916f5cb0b2ee2c27de13357d)
            .def("cbegin", method_pointer_d0b72163dd185c51aaab3d15e80c67a7)
            .def("cend", method_pointer_d34205575084545c87346e6022bb46f9)
            .def("crbegin", method_pointer_cb139ac60d245079a0410998c8695bbb)
            .def("crend", method_pointer_121616a191925aa0a5ac7a0aded6c260)
            .def("size", method_pointer_2976530ccb5657b6b6c0a54e0f5c9bf8)
            .def("max_size", method_pointer_1aacab760834510ea8d4153fc713ba92)
            .def("resize", method_pointer_70ee84578153504299edef13598e052d)
            .def("resize", method_pointer_bd273cd8cf6052bab101f81c2e8e2af8)
            .def("shrink_to_fit", method_pointer_7e1a87abab0e55769d56b92f87de9364)
            .def("capacity", method_pointer_c4a00f13846e5f839d297a023ef32936)
            .def("empty", method_pointer_ca8611fc69aa580da7d1319be7a1b2af)
            .def("reserve", method_pointer_a51bd8eeb45656289869d33118b75250)
            .def("__getitem__", method_pointer_fdeddab37b135a37b33c6270959703f7, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("at", method_pointer_563d298b97c95b388324b19f0f74ab68, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("front", method_pointer_73ae975684b25694bbeee2d13496bdba, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("back", method_pointer_6acc6f10ad2e54568819075f729606fa, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("push_back", method_pointer_dd51d9a7dfbd547ea2f1392263c449cb)
            .def("pop_back", method_pointer_2d363ed2fc5950f68a03d934d57851cc)
            .def("insert", method_pointer_162d49c09a54599f91a9a46de35d3d14)
            .def("insert", method_pointer_abdb6a9576965e2db5d88d95340257bd)
            .def("insert", method_pointer_3416c40ef9925177aafabc4732e41eb1)
            .def("erase", method_pointer_761bc5712055510b8c735d5f0d3cef9f)
            .def("erase", method_pointer_f7f15843e71f5e08ac9de17d0717e8f4)
            .def("swap", method_pointer_f617e6b977e651a6ae8985dd17f99bda)
            .def("clear", method_pointer_823a942d3fe95fb9bee873042664f0b9);
}