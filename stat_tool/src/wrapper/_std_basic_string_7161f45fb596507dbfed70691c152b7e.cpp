#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_basic_string_7161f45fb596507dbfed70691c152b7e()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        unsigned long (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_4d9f8cb0a6ca50e1b12aa434620a0df9)() const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::size;
        unsigned long (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_7da97470fe6e50afbd7d156ccb5e185c)() const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::length;
        unsigned long (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_f96880c40cd5581c9b43b4e6f16dd6f7)() const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::max_size;
        void (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_e85c72d22ed853e5a26c83f61939f0a5)(unsigned long, char) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::resize;
        void (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_c0a134e0c521551fa3bc026e11341ab3)(unsigned long) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::resize;
        void (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_c05a1012af4f583791eed75518272dd9)() = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::shrink_to_fit;
        unsigned long (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_d2e28dce8283515487542d3414a7198f)() const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::capacity;
        void (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_a389ccce937b57179e772a77ef41b59d)(unsigned long) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::reserve;
        void (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_8f5643ab2c3c504abc0124fe73bb974a)() = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::clear;
        bool (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_3b1638775ee45c42ace83d53561b439e)() const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::empty;
        char const & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_ec0b90f799fa5980af1270b29e6639a9)(unsigned long) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::operator[];
        char const & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_9e976617ada95eebac521a71277f857d)(unsigned long) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::at;
        char const & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_ddb4d5308ba558ffae352ba8df9f4cd9)() const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::front;
        char const & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_1ad749a54e055794bf8043fd85bec066)() const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::back;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_97eda6673d5250d6a77fee856230becd)(class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::operator+=;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_5baf1ee718835552b22bababb6df8b1c)(char) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::operator+=;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_544f8339f4ae546eb6eee6a646a54f85)(class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::append;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_6d1bf569c74859d78cd105f65f255c19)(class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, unsigned long, unsigned long) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::append;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_03c5040018505d949dcfb28af3b542b2)(unsigned long, char) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::append;
        void (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_da708dd35c0758cb9c0b77d4fd09b30d)(char) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::push_back;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_e5836a6fd43e51c3947c2385b8bec184)(class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::assign;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_c39f221dc38f5d8fb5e9b5352fc65be3)(class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, unsigned long, unsigned long) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::assign;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_551afc8548a2529e8a1709d95f963351)(unsigned long, char) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::assign;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_c766678d25a256ba80381939b7e23f17)(unsigned long, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::insert;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_980f32d36f265639ae31b5459cd04db7)(unsigned long, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, unsigned long, unsigned long) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::insert;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_7219968409505d3b86186534121e5926)(unsigned long, unsigned long, char) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::insert;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_eb9cbe3a084556c49e2c0fcc05b52707)(unsigned long, unsigned long) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::erase;
        void (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_4b15a8695d4d5ec48c1dc38d8064951b)() = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::pop_back;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_b47c862f8d2c5ec29306464701a8eb44)(unsigned long, unsigned long, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::replace;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_60f59d5bb5265d79a4930ea3eb4f6756)(unsigned long, unsigned long, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, unsigned long, unsigned long) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::replace;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > & (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_cce999fa9ec45975a93ce3bcc2df98f0)(unsigned long, unsigned long, unsigned long, char) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::replace;
        void (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_2803a1d2a6c15604b2bc0449a50ad6b0)(class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > &) = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::swap;
        unsigned long (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_85acc06d2d5457b4b6f1573d70885cc7)(class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, unsigned long) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::find;
        unsigned long (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_03c6e22de70955d5a2e71d844f3f80a7)(char, unsigned long) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::find;
        unsigned long (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_29fa605e6a875ca08a7de75977ae02ea)(class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, unsigned long) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::rfind;
        unsigned long (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_acd38ec160b351708d11cb4732825209)(char, unsigned long) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::rfind;
        unsigned long (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_3ee8fed5cda75fe3870a2900b1e20350)(class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, unsigned long) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::find_first_of;
        unsigned long (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_dcecca6be8395a1d8b062a559e4b76a3)(char, unsigned long) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::find_first_of;
        unsigned long (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_8be153d847635c9c89f427aeabecb1d0)(class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, unsigned long) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::find_last_of;
        unsigned long (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_2f2ed413410a5e2d87842baee08d3006)(char, unsigned long) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::find_last_of;
        unsigned long (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_189542fc8b1355b09de7b0a53eb42b49)(class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, unsigned long) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::find_first_not_of;
        unsigned long (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_1f771f012afc526aa23bc5227ce7478c)(char, unsigned long) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::find_first_not_of;
        unsigned long (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_7b0f5c3243e15c789d83ac53b4299e61)(class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, unsigned long) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::find_last_not_of;
        unsigned long (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_ce9ad1c08b4b56938c140a0956c3c4ee)(char, unsigned long) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::find_last_not_of;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_3f7192510311542d90639b28d300da58)(unsigned long, unsigned long) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::substr;
        int (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_b9ab99358759507793b52fb4462b1b92)(class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::compare;
        int (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_9ca46c97061951f59b012c6faf285629)(unsigned long, unsigned long, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::compare;
        int (::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::*method_pointer_c36f7b6379a250c9b76bfeee452e90cc)(unsigned long, unsigned long, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, unsigned long, unsigned long) const = &::std::basic_string<char, std::char_traits<char>, std::allocator<char> >::compare;
        boost::python::class_< class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::shared_ptr< class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > > >("_BasicString_7161f45fb596507dbfed70691c152b7e", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::std::allocator<char> const & >())
            .def(boost::python::init< class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const & >())
            .def(boost::python::init< class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, unsigned long, unsigned long >())
            .def(boost::python::init< class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, unsigned long, unsigned long, class ::std::allocator<char> const & >())
            .def(boost::python::init< class ::std::initializer_list<char>, class ::std::allocator<char> const & >())
            .def("size", method_pointer_4d9f8cb0a6ca50e1b12aa434620a0df9)
            .def("length", method_pointer_7da97470fe6e50afbd7d156ccb5e185c)
            .def("max_size", method_pointer_f96880c40cd5581c9b43b4e6f16dd6f7)
            .def("resize", method_pointer_e85c72d22ed853e5a26c83f61939f0a5)
            .def("resize", method_pointer_c0a134e0c521551fa3bc026e11341ab3)
            .def("shrink_to_fit", method_pointer_c05a1012af4f583791eed75518272dd9)
            .def("capacity", method_pointer_d2e28dce8283515487542d3414a7198f)
            .def("reserve", method_pointer_a389ccce937b57179e772a77ef41b59d)
            .def("clear", method_pointer_8f5643ab2c3c504abc0124fe73bb974a)
            .def("empty", method_pointer_3b1638775ee45c42ace83d53561b439e)
            .def("__getitem__", method_pointer_ec0b90f799fa5980af1270b29e6639a9, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("at", method_pointer_9e976617ada95eebac521a71277f857d, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("front", method_pointer_ddb4d5308ba558ffae352ba8df9f4cd9, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("back", method_pointer_1ad749a54e055794bf8043fd85bec066, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("__iadd__", method_pointer_97eda6673d5250d6a77fee856230becd, boost::python::return_internal_reference<>())
            .def("__iadd__", method_pointer_5baf1ee718835552b22bababb6df8b1c, boost::python::return_internal_reference<>())
            .def("append", method_pointer_544f8339f4ae546eb6eee6a646a54f85, boost::python::return_internal_reference<>())
            .def("append", method_pointer_6d1bf569c74859d78cd105f65f255c19, boost::python::return_internal_reference<>())
            .def("append", method_pointer_03c5040018505d949dcfb28af3b542b2, boost::python::return_internal_reference<>())
            .def("push_back", method_pointer_da708dd35c0758cb9c0b77d4fd09b30d)
            .def("assign", method_pointer_e5836a6fd43e51c3947c2385b8bec184, boost::python::return_internal_reference<>())
            .def("assign", method_pointer_c39f221dc38f5d8fb5e9b5352fc65be3, boost::python::return_internal_reference<>())
            .def("assign", method_pointer_551afc8548a2529e8a1709d95f963351, boost::python::return_internal_reference<>())
            .def("insert", method_pointer_c766678d25a256ba80381939b7e23f17, boost::python::return_internal_reference<>())
            .def("insert", method_pointer_980f32d36f265639ae31b5459cd04db7, boost::python::return_internal_reference<>())
            .def("insert", method_pointer_7219968409505d3b86186534121e5926, boost::python::return_internal_reference<>())
            .def("erase", method_pointer_eb9cbe3a084556c49e2c0fcc05b52707, boost::python::return_internal_reference<>())
            .def("pop_back", method_pointer_4b15a8695d4d5ec48c1dc38d8064951b)
            .def("replace", method_pointer_b47c862f8d2c5ec29306464701a8eb44, boost::python::return_internal_reference<>())
            .def("replace", method_pointer_60f59d5bb5265d79a4930ea3eb4f6756, boost::python::return_internal_reference<>())
            .def("replace", method_pointer_cce999fa9ec45975a93ce3bcc2df98f0, boost::python::return_internal_reference<>())
            .def("swap", method_pointer_2803a1d2a6c15604b2bc0449a50ad6b0)
            .def("find", method_pointer_85acc06d2d5457b4b6f1573d70885cc7)
            .def("find", method_pointer_03c6e22de70955d5a2e71d844f3f80a7)
            .def("rfind", method_pointer_29fa605e6a875ca08a7de75977ae02ea)
            .def("rfind", method_pointer_acd38ec160b351708d11cb4732825209)
            .def("find_first_of", method_pointer_3ee8fed5cda75fe3870a2900b1e20350)
            .def("find_first_of", method_pointer_dcecca6be8395a1d8b062a559e4b76a3)
            .def("find_last_of", method_pointer_8be153d847635c9c89f427aeabecb1d0)
            .def("find_last_of", method_pointer_2f2ed413410a5e2d87842baee08d3006)
            .def("find_first_not_of", method_pointer_189542fc8b1355b09de7b0a53eb42b49)
            .def("find_first_not_of", method_pointer_1f771f012afc526aa23bc5227ce7478c)
            .def("find_last_not_of", method_pointer_7b0f5c3243e15c789d83ac53b4299e61)
            .def("find_last_not_of", method_pointer_ce9ad1c08b4b56938c140a0956c3c4ee)
            .def("substr", method_pointer_3f7192510311542d90639b28d300da58)
            .def("compare", method_pointer_b9ab99358759507793b52fb4462b1b92)
            .def("compare", method_pointer_9ca46c97061951f59b012c6faf285629)
            .def("compare", method_pointer_c36f7b6379a250c9b76bfeee452e90cc);
}