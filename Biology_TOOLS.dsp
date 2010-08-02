# Microsoft Developer Studio Project File - Name="Biology_TOOLS" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=Biology_TOOLS - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Biology_TOOLS.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Biology_TOOLS.mak" CFG="Biology_TOOLS - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Biology_TOOLS - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "Biology_TOOLS - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "Biology_TOOLS - Win32 Release"

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
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /Yu"stdafx.h" /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /Yu"stdafx.h" /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "Biology_TOOLS - Win32 Debug"

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
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /Yu"stdafx.h" /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /FR /FD /GZ /c
# SUBTRACT CPP /YX
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "Biology_TOOLS - Win32 Release"
# Name "Biology_TOOLS - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Group "Matrix_Class"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\Matrix_Class\bandmat.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\cholesky.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\evalue.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\fft.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\hholder.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\jacobi.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\myexcept.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newfft.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmat1.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmat2.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmat3.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmat4.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmat5.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmat6.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmat7.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmat8.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmat9.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmatex.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmatnl.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmatrm.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\solution.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\sort.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\submat.cpp
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\svd.cpp
# End Source File
# End Group
# Begin Group "Random_Class"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\Random_Class\extreal.cpp
# End Source File
# Begin Source File

SOURCE=.\Random_Class\newran.cpp
# End Source File
# End Group
# Begin Source File

SOURCE=.\Biology_TOOLS.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Group "Matrix_Class_H"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\Matrix_Class\boolean.h
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\controlw.h
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\include.h
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\myexcept.h
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmat.h
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmatap.h
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmatio.h
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmatnl.h
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmatrc.h
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\newmatrm.h
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\precisio.h
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\solution.h
# End Source File
# Begin Source File

SOURCE=.\Matrix_Class\tmt.h
# End Source File
# End Group
# Begin Group "Tree_Class_H"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\Tree_Class\tree_msvc.h
# End Source File
# End Group
# Begin Group "Random_Class_H"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\Random_Class\extreal.h
# End Source File
# Begin Source File

SOURCE=.\Random_Class\newran.h
# End Source File
# End Group
# Begin Group "tool_box"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\TOOLS_Class\tools.h
# End Source File
# End Group
# Begin Source File

SOURCE=.\comparativeAnalysis.h
# End Source File
# Begin Source File

SOURCE=.\metaAnalysis.h
# End Source File
# Begin Source File

SOURCE=.\phylogeny.h
# End Source File
# Begin Source File

SOURCE=.\statistics.h
# End Source File
# Begin Source File

SOURCE=.\StdAfx.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# Begin Source File

SOURCE=.\meta_data.txt
# End Source File
# Begin Source File

SOURCE=.\old_source_lambda.txt
# End Source File
# Begin Source File

SOURCE=.\ReadMe.txt
# End Source File
# End Target
# End Project
