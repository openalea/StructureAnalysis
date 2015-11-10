# Microsoft Developer Studio Project File - Name="tools" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=tools - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "tools.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "tools.mak" CFG="tools - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "tools - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "tools - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "tools - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MDd /W3 /O1 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /FD /c
# ADD CPP /nologo /MDd /GR /O2 /I "$(QTDIR)\include" /I ".." /D "NDEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "UNICODE" /D "QT_THREAD_SUPPORT" /D "BISON_HPP" /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "tools - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MDd /W3 /Gm /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /FD /GZ /c
# ADD CPP /nologo /MDd /GR /ZI /Od /I ".." /I "$(QTDIR)\include" /D "_DEBUG" /D "WITHOUT_GLUT" /D "RWOUT" /D "WIN32" /D "_MBCS" /D "_LIB" /D "UNICODE" /D "QT_THREAD_SUPPORT" /D "BISON_HPP" /FR /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "tools - Win32 Release"
# Name "tools - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\dirnames.cpp
# End Source File
# Begin Source File

SOURCE=.\errormsg.cpp
# End Source File
# Begin Source File

SOURCE=.\readline.cpp
# End Source File
# Begin Source File

SOURCE=.\timer.cpp
# End Source File
# Begin Source File

SOURCE=.\util_enviro.cpp
# End Source File
# Begin Source File

SOURCE=.\util_glut.cpp
# End Source File
# Begin Source File

SOURCE=.\util_math.cpp
# End Source File
# Begin Source File

SOURCE=.\util_matrix.cpp
# End Source File
# Begin Source File

SOURCE=.\util_vector.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\all_tools
# End Source File
# Begin Source File

SOURCE=.\all_tools.h
# End Source File
# Begin Source File

SOURCE=.\bfstream.h
# End Source File
# Begin Source File

SOURCE=.\config.h
# End Source File
# Begin Source File

SOURCE=.\dirnames.h
# End Source File
# Begin Source File

SOURCE=.\errormsg.h
# End Source File
# Begin Source File

SOURCE=.\gparser.h
# End Source File
# Begin Source File

SOURCE=.\gscanner.h
# End Source File
# Begin Source File

SOURCE=.\gsmbtable.h
# End Source File
# Begin Source File

SOURCE=.\rcobject.h
# End Source File
# Begin Source File

SOURCE=.\readline.h
# End Source File
# Begin Source File

SOURCE=.\std.h
# End Source File
# Begin Source File

SOURCE=.\template2_parser.h
# End Source File
# Begin Source File

SOURCE=.\template2_scanner.h
# End Source File
# Begin Source File

SOURCE=.\timer.h
# End Source File
# Begin Source File

SOURCE=.\tools_namespace.h
# End Source File
# Begin Source File

SOURCE=.\util_assert.h
# End Source File
# Begin Source File

SOURCE=.\util_enviro.h
# End Source File
# Begin Source File

SOURCE=.\util_gl.h
# End Source File
# Begin Source File

SOURCE=.\util_glut.h
# End Source File
# Begin Source File

SOURCE=.\util_math.h
# End Source File
# Begin Source File

SOURCE=.\util_matrix.h
# End Source File
# Begin Source File

SOURCE=.\util_polymath.h
# End Source File
# Begin Source File

SOURCE=.\util_string.h
# End Source File
# Begin Source File

SOURCE=.\util_tuple.h
# End Source File
# Begin Source File

SOURCE=.\util_types.h
# End Source File
# Begin Source File

SOURCE=.\util_vector.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# Begin Source File

SOURCE=.\template2_parser.y
# End Source File
# Begin Source File

SOURCE=.\template2_scanner.l
# End Source File
# End Group
# Begin Group "Interfaces"

# PROP Default_Filter "ui"
# End Group
# Begin Group "RogueWave"

# PROP Default_Filter ""
# Begin Group "RW Source Files"

# PROP Default_Filter "*.cpp"
# Begin Source File

SOURCE=.\rw_cstring.cpp
# End Source File
# Begin Source File

SOURCE=.\rw_date.cpp
# End Source File
# Begin Source File

SOURCE=.\rw_time.cpp
# End Source File
# End Group
# Begin Group "RW Header Files"

# PROP Default_Filter "*.h"
# Begin Source File

SOURCE=.\rw_bitvec.h
# End Source File
# Begin Source File

SOURCE=.\rw_comp.h
# End Source File
# Begin Source File

SOURCE=.\rw_cstring.h
# End Source File
# Begin Source File

SOURCE=.\rw_date.h
# End Source File
# Begin Source File

SOURCE=.\rw_defs.h
# End Source File
# Begin Source File

SOURCE=.\rw_hash.h
# End Source File
# Begin Source File

SOURCE=.\rw_hdict.h
# End Source File
# Begin Source File

SOURCE=.\rw_hset.h
# End Source File
# Begin Source File

SOURCE=.\rw_list.h
# End Source File
# Begin Source File

SOURCE=.\rw_locale.h
# End Source File
# Begin Source File

SOURCE=.\rw_macro.h
# End Source File
# Begin Source File

SOURCE=.\rw_queue.h
# End Source File
# Begin Source File

SOURCE=.\rw_slist.h
# End Source File
# Begin Source File

SOURCE=.\rw_stack.h
# End Source File
# Begin Source File

SOURCE=.\rw_time.h
# End Source File
# Begin Source File

SOURCE=.\rw_tokenizer.h
# End Source File
# Begin Source File

SOURCE=.\rw_vector.h
# End Source File
# Begin Source File

SOURCE=.\rw_zone.h
# End Source File
# End Group
# End Group
# End Target
# End Project
