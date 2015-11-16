#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _std_vector_a0d6e5da7dec5b598e07f046e4fb632f()
{
        std::string std_a5e4e9231d6351ccb0e06756b389f0af_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
        boost::python::object std_a5e4e9231d6351ccb0e06756b389f0af_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(std_a5e4e9231d6351ccb0e06756b389f0af_name.c_str()))));
        boost::python::scope().attr("std") = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        boost::python::scope std_a5e4e9231d6351ccb0e06756b389f0af_scope = std_a5e4e9231d6351ccb0e06756b389f0af_module;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_4ccc81500a895368a9e8d22d68286743)(unsigned long, int const &) = &::std::vector<int, std::allocator<int> >::assign;
        class ::__gnu_cxx::__normal_iterator<int *, std::vector<int, std::allocator<int> > > (::std::vector<int, std::allocator<int> >::*method_pointer_5106c97036bd5ad2bad92be071857d00)() = &::std::vector<int, std::allocator<int> >::end;
        class ::std::reverse_iterator<__gnu_cxx::__normal_iterator<int *, std::vector<int, std::allocator<int> > > > (::std::vector<int, std::allocator<int> >::*method_pointer_86e5e0d0812c5b54ac297b513ffe6167)() = &::std::vector<int, std::allocator<int> >::rbegin;
        class ::std::reverse_iterator<__gnu_cxx::__normal_iterator<int *, std::vector<int, std::allocator<int> > > > (::std::vector<int, std::allocator<int> >::*method_pointer_1b88a2e6f9815faab55b4ef0f8230a71)() = &::std::vector<int, std::allocator<int> >::rend;
        class ::std::reverse_iterator<__gnu_cxx::__normal_iterator<const int *, std::vector<int, std::allocator<int> > > > (::std::vector<int, std::allocator<int> >::*method_pointer_c20f755c5b5f52389a163f3d4e40c117)() const = &::std::vector<int, std::allocator<int> >::rend;
        class ::__gnu_cxx::__normal_iterator<const int *, std::vector<int, std::allocator<int> > > (::std::vector<int, std::allocator<int> >::*method_pointer_3d6c2e66e8b95ccea89f1915bc5f1940)() const = &::std::vector<int, std::allocator<int> >::cbegin;
        class ::std::reverse_iterator<__gnu_cxx::__normal_iterator<const int *, std::vector<int, std::allocator<int> > > > (::std::vector<int, std::allocator<int> >::*method_pointer_96a98e3b7b6f568598fd9238b4270030)() const = &::std::vector<int, std::allocator<int> >::crbegin;
        class ::std::reverse_iterator<__gnu_cxx::__normal_iterator<const int *, std::vector<int, std::allocator<int> > > > (::std::vector<int, std::allocator<int> >::*method_pointer_e0a06df8d2185638a7cd78334ec0db17)() const = &::std::vector<int, std::allocator<int> >::crend;
        unsigned long (::std::vector<int, std::allocator<int> >::*method_pointer_e07ec54b5ecf51aba9ea64933865d66f)() const = &::std::vector<int, std::allocator<int> >::size;
        unsigned long (::std::vector<int, std::allocator<int> >::*method_pointer_71bf07dd40255055b8a324e995732a46)() const = &::std::vector<int, std::allocator<int> >::max_size;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_33512a4a00db51a98e4d13123a3bad29)(unsigned long) = &::std::vector<int, std::allocator<int> >::resize;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_178ab65394295db6a364a19c28967f04)(unsigned long, int const &) = &::std::vector<int, std::allocator<int> >::resize;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_a7b649cf3de95b3ab2863b316b860dfc)() = &::std::vector<int, std::allocator<int> >::shrink_to_fit;
        unsigned long (::std::vector<int, std::allocator<int> >::*method_pointer_2449152b3f41528eb23d2589039de7c2)() const = &::std::vector<int, std::allocator<int> >::capacity;
        bool (::std::vector<int, std::allocator<int> >::*method_pointer_6763fbcbb2be52318e0f476bfa741210)() const = &::std::vector<int, std::allocator<int> >::empty;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_b4515a904899566390be0b4ff0194a91)(unsigned long) = &::std::vector<int, std::allocator<int> >::reserve;
        int const & (::std::vector<int, std::allocator<int> >::*method_pointer_96e6e91e2dc15b239c742513ae954e9d)(unsigned long) const = &::std::vector<int, std::allocator<int> >::operator[];
        int const & (::std::vector<int, std::allocator<int> >::*method_pointer_206368d540105e6b838fbd180d49ba5c)(unsigned long) const = &::std::vector<int, std::allocator<int> >::at;
        int const & (::std::vector<int, std::allocator<int> >::*method_pointer_638ebf46a23d50798c1d1afc6006bcc9)() const = &::std::vector<int, std::allocator<int> >::front;
        int const & (::std::vector<int, std::allocator<int> >::*method_pointer_37735e7cc42857a58f3b52a6b27205b7)() const = &::std::vector<int, std::allocator<int> >::back;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_3b9d531e6aa554fdaa7229f51472adc2)(int const &) = &::std::vector<int, std::allocator<int> >::push_back;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_8d12d434085a5f54a8c38754cb8ce092)() = &::std::vector<int, std::allocator<int> >::pop_back;
        class ::__gnu_cxx::__normal_iterator<int *, std::vector<int, std::allocator<int> > > (::std::vector<int, std::allocator<int> >::*method_pointer_fe3e52f4314e5fbda2618a82fe86cf95)(class ::__gnu_cxx::__normal_iterator<const int *, std::vector<int, std::allocator<int> > >, int const &) = &::std::vector<int, std::allocator<int> >::insert;
        class ::__gnu_cxx::__normal_iterator<int *, std::vector<int, std::allocator<int> > > (::std::vector<int, std::allocator<int> >::*method_pointer_695ff395f758583d887f199400f6b196)(class ::__gnu_cxx::__normal_iterator<const int *, std::vector<int, std::allocator<int> > >, class ::std::initializer_list<int>) = &::std::vector<int, std::allocator<int> >::insert;
        class ::__gnu_cxx::__normal_iterator<int *, std::vector<int, std::allocator<int> > > (::std::vector<int, std::allocator<int> >::*method_pointer_897934f2c090570dab07d88d1777798c)(class ::__gnu_cxx::__normal_iterator<const int *, std::vector<int, std::allocator<int> > >, unsigned long, int const &) = &::std::vector<int, std::allocator<int> >::insert;
        class ::__gnu_cxx::__normal_iterator<int *, std::vector<int, std::allocator<int> > > (::std::vector<int, std::allocator<int> >::*method_pointer_02f4ed4cb7f45720a79289600e969e2e)(class ::__gnu_cxx::__normal_iterator<const int *, std::vector<int, std::allocator<int> > >) = &::std::vector<int, std::allocator<int> >::erase;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_730a17c027895110a75f98866282c7c5)(class ::std::vector<int, std::allocator<int> > &) = &::std::vector<int, std::allocator<int> >::swap;
        void (::std::vector<int, std::allocator<int> >::*method_pointer_9d23ebaeda6d5d49b49b9958e15b5e51)() = &::std::vector<int, std::allocator<int> >::clear;
        boost::python::class_< class ::std::vector<int, std::allocator<int> >, std::shared_ptr< class ::std::vector<int, std::allocator<int> > > >("_Vector_a0d6e5da7dec5b598e07f046e4fb632f", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::std::allocator<int> const & >())
            .def(boost::python::init< unsigned long, class ::std::allocator<int> const & >())
            .def(boost::python::init< unsigned long, int const &, class ::std::allocator<int> const & >())
            .def(boost::python::init< class ::std::vector<int, std::allocator<int> > const & >())
            .def("assign", method_pointer_4ccc81500a895368a9e8d22d68286743)
            .def("end", method_pointer_5106c97036bd5ad2bad92be071857d00)
            .def("rbegin", method_pointer_86e5e0d0812c5b54ac297b513ffe6167)
            .def("rend", method_pointer_1b88a2e6f9815faab55b4ef0f8230a71)
            .def("rend", method_pointer_c20f755c5b5f52389a163f3d4e40c117)
            .def("cbegin", method_pointer_3d6c2e66e8b95ccea89f1915bc5f1940)
            .def("crbegin", method_pointer_96a98e3b7b6f568598fd9238b4270030)
            .def("crend", method_pointer_e0a06df8d2185638a7cd78334ec0db17)
            .def("size", method_pointer_e07ec54b5ecf51aba9ea64933865d66f)
            .def("max_size", method_pointer_71bf07dd40255055b8a324e995732a46)
            .def("resize", method_pointer_33512a4a00db51a98e4d13123a3bad29)
            .def("resize", method_pointer_178ab65394295db6a364a19c28967f04)
            .def("shrink_to_fit", method_pointer_a7b649cf3de95b3ab2863b316b860dfc)
            .def("capacity", method_pointer_2449152b3f41528eb23d2589039de7c2)
            .def("empty", method_pointer_6763fbcbb2be52318e0f476bfa741210)
            .def("reserve", method_pointer_b4515a904899566390be0b4ff0194a91)
            .def("__getitem__", method_pointer_96e6e91e2dc15b239c742513ae954e9d, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("at", method_pointer_206368d540105e6b838fbd180d49ba5c, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("front", method_pointer_638ebf46a23d50798c1d1afc6006bcc9, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("back", method_pointer_37735e7cc42857a58f3b52a6b27205b7, boost::python::return_value_policy< boost::python::return_by_value >())
            .def("push_back", method_pointer_3b9d531e6aa554fdaa7229f51472adc2)
            .def("pop_back", method_pointer_8d12d434085a5f54a8c38754cb8ce092)
            .def("insert", method_pointer_fe3e52f4314e5fbda2618a82fe86cf95)
            .def("insert", method_pointer_695ff395f758583d887f199400f6b196)
            .def("insert", method_pointer_897934f2c090570dab07d88d1777798c)
            .def("erase", method_pointer_02f4ed4cb7f45720a79289600e969e2e)
            .def("swap", method_pointer_730a17c027895110a75f98866282c7c5)
            .def("clear", method_pointer_9d23ebaeda6d5d49b49b9958e15b5e51);
}