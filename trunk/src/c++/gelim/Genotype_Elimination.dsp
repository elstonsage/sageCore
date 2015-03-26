# Microsoft Developer Studio Project File - Name="Genotype Elimination" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=Genotype Elimination - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Genotype_Elimination.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Genotype_Elimination.mak" CFG="Genotype Elimination - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Genotype Elimination - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "Genotype Elimination - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=xicl6.exe
RSC=rc.exe

!IF  "$(CFG)" == "Genotype Elimination - Win32 Release"

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
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /GX /O1 /I "..\..\include" /I "..\include" /I "d:\program files\xtl" /D "_WINDOWS" /D "HAS_SSTREAM" /D "MSDOS" /D "WIN32" /D "NDEBUG" /FD /c
# ADD BASE RSC /l 0x409
# ADD RSC /l 0x409
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=xilink6.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"Release\geno_elim.lib"
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Desc=Moving library to its final resting place
PostBuild_Cmds=copy .\Release\*.lib ..\lib
# End Special Build Tool

!ELSEIF  "$(CFG)" == "Genotype Elimination - Win32 Debug"

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
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /Gi /GX /ZI /Od /I "." /I "..\include" /D "_WINDOWS" /D "MSDOS" /D "HAS_SSTREAM" /D "WIN32" /D "_DEBUG" /FR /YX /FD /c
# ADD BASE RSC /l 0x409
# ADD RSC /l 0x409
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=xilink6.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"Debug\geno_elim.lib"
# Begin Special Build Tool
SOURCE="$(InputPath)"
PostBuild_Desc=Moving library to its final resting place
PostBuild_Cmds=copy .\Debug\*.lib ..\lib
# End Special Build Tool

!ENDIF 

# Begin Target

# Name "Genotype Elimination - Win32 Release"
# Name "Genotype Elimination - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=.\geno_eliminate.cpp
# End Source File
# Begin Source File

SOURCE=.\ML_Remap.cpp
# End Source File
# Begin Source File

SOURCE=.\ped_marker.cpp
# End Source File
# Begin Source File

SOURCE=.\ped_phenotypes.cpp
# End Source File
# Begin Source File

SOURCE=.\pedigree_region.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h"
# Begin Source File

SOURCE=..\include\geno_elim\allele_set.h
# End Source File
# Begin Source File

SOURCE=..\include\geno_elim\bitfield.h
# End Source File
# Begin Source File

SOURCE=..\include\geno_elim\geno_eliminate.h
# End Source File
# Begin Source File

SOURCE=..\include\geno_elim\genotype.h
# End Source File
# Begin Source File

SOURCE=..\include\geno_elim\ML_Remap.h
# End Source File
# Begin Source File

SOURCE=..\include\geno_elim\ped_marker.h
# End Source File
# Begin Source File

SOURCE=..\include\geno_elim\ped_phenotypes.h
# End Source File
# Begin Source File

SOURCE=..\include\geno_elim\pedigree_region.h
# End Source File
# End Group
# End Target
# End Project
