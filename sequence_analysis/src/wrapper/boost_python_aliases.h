


// synopsis: converts a boost python list into a CPP array
// INPUT is the input variable name
// TYPE is the type of the CPP array to be returned into the variable data
// VARIABLE_NAME allows to use this macros several time inside the same function
#define CREATE_ARRAY(INPUT, TYPE, VARIABLE_NAME)\
  int VARIABLE_NAME##_size = len(INPUT);\
  sequence_analysis::wrap_util::auto_ptr_array<TYPE> \
      VARIABLE_NAME(new TYPE[VARIABLE_NAME##_size]);\
  for (int i = 0; i < VARIABLE_NAME##_size; i++)\
      VARIABLE_NAME[i] = boost::python::extract<TYPE>(INPUT[i]);


// prototype not included since the name cannot be overloaded
// INPUT the input argument
// METHOD_NAME: name of the method that will be called on INPUT
// OUTPUT_TYPE: type returned by METHOD_NAME
// optional list of arguments to be used by METHOD_NAME
#define SIMPLE_METHOD_TEMPLATE_1(INPUT, METHOD_NAME, OUTPUT_TYPE, ...)\
    StatError error; \
    OUTPUT_TYPE* ret;\
    ret = INPUT.METHOD_NAME(error, __VA_ARGS__ );\
    if (!ret)\
      sequence_analysis::wrap_util::throw_error(error);\
    return ret;\

#define SIMPLE_METHOD_TEMPLATE_0(INPUT, METHOD_NAME, OUTPUT_TYPE)\
    StatError error; \
    OUTPUT_TYPE* ret;\
    ret = INPUT.METHOD_NAME(error);\
    if (!ret)\
      sequence_analysis::wrap_util::throw_error(error);\
    return ret;\

#define FOOTER_OS \
    if (!ret) \
      sequence_analysis::wrap_util::throw_error(error);\
    cout << os.str() << endl;\
    return ret;\

#define FOOTER \
    if (!ret) \
      sequence_analysis::wrap_util::throw_error(error);\
    return ret;\

// !! don't use if TYPE=bool
#define HEADER(TYPE) \
    StatError error; \
    TYPE* ret;\


#define HEADER_OS(TYPE) \
    StatError error; \
    TYPE* ret;\
    std::stringstream os;\



// INPUT PARAMETER VALIDATION -------------



#define CHECK(VALUE, MIN, MAX) \
	if (VALUE< MIN || VALUE>=MAX)\
     {\
    	 PyErr_SetString(PyExc_TypeError,(error_message.str()).c_str());\
    	 throw_error_already_set();\
     }\


// quick alias to avoid writting the return_value_policy...
#define DEF_RETURN_VALUE(NAME, REFERENCE, ARGS, DOCSTRING) \
    .def(NAME, REFERENCE, return_value_policy< manage_new_object >(), ARGS, DOCSTRING)

#define DEF_RETURN_VALUE_NO_ARGS(NAME, REFERENCE, DOCSTRING) \
    .def(NAME, REFERENCE, return_value_policy< manage_new_object >(), DOCSTRING)


