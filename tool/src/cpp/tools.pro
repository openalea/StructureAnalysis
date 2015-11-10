HEADERS	+= \
        gparser.h \
	gsmbtable.h \
	util_tuple.h \
	config.h \
	std.h \
	bfstream.h \
	rcobject.h \
	util_polymath.h \
        gscanner.h \
	template2_parser.h \
	timer.h \
	errormsg.h \
	readline.h \
	dirnames.h \
	util_matrix.h \
	util_vector.h \
	util_math.h \

SOURCES	+= \
	timer.cpp \
	errormsg.cpp \
	readline.cpp \
	dirnames.cpp \
	util_matrix.cpp \
	util_vector.cpp \
	util_math.cpp 

TEMPLATE	= lib
CONFIG	+= debug

TARGET = tools

unix:OBJECT_DIR = .libs

