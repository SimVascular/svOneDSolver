# Hey emacs, this is a -*- makefile -*-

# Copyright (c) Stanford University, The Regents of the University of
#               California, and others.
#
# All Rights Reserved.
#
# See Copyright-SimVascular.txt for additional details.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# This is where we want to define things that will be useful for
# multiple packages.  Vtk, for example, is used by both the level set
# package for front extraction and by the utils package to implement
# general vtk helper functions.

# Important notes:
#   - set TARGETDIR
#   - specify inclusion of { Discrete, Parasolid, Irit } via SV_USE_*
#     vars

# ----------------------------------------------------
# by default, no check for dependancies when compiling
# ----------------------------------------------------

NO_DEPEND = 1

# -----------------------------------------------------------
# CLUSTER = { x64_cygwin, x64_linux, x64_macosx }
# -----------------------------------------------------------

CLUSTER = x64_cygwin
#CLUSTER = x64_linux
#CLUSTER = x64_macosx

# ---------------------------------------------------------------------
# CXX_COMPILER_VERSION = { icpc, msvc-19.0, msvc-19.16, clang,
#                          mingw-gcc, gcc}
# FORTRAN_COMPILER_VERSION = { ifort, mingw-gfortran, gfortran }
# ---------------------------------------------------------------------

# NOTE: CXX_COMPILER_VERSION and FORTRAN_COMPILER_VERSION
#       should be replaced with additional variables

# Windows 10
SV_COMPILER = msvc
#SV_COMPILER_VERSION = 19.0
#CXX_COMPILER_VERSION = msvc-19.0
SV_COMPILER_VERSION = 19.16
CXX_COMPILER_VERSION = msvc-19.16
FORTRAN_COMPILER_VERSION = ifort

# Ubuntu 18
# CLUSTER=x64_linux
# CXX_COMPILER_VERSION=gcc
# FORTRAN_COMPILER_VERSION=gfortran
# SV_COMPILER=gnu
# SV_COMPILER_VERSION=7.3

# optionally override with cluster options
# -----------------------------------------------------------------------

ifeq ($(LOCAL_DIR_CLUSTER_OVERRIDES),1)
-include cluster_overrides.mk
else
-include $(TOP)/cluster_overrides.mk
endif

# -------
# globals
# -------

SV_USE_SHARED = 1

# ----------------------------------------------
# Control inclusion of leslib
# {binary and dummy} are mutually exclusive opts
# ----------------------------------------------

SV_USE_BINARY_LESLIB = 0
SV_USE_DUMMY_LESLIB = 1

# --------------------------------------------------------
# Control inclusion of svLS
# {binary, dummy, source code} are mutually exclusive opts
# --------------------------------------------------------

SV_USE_DUMMY_SVLS = 0
SV_USE_SOURCE_CODE_SVLS = 1

# -----------------------------------------------------
# Compile with zlib
# -----------------------------------------------------

SV_USE_ZLIB = 1

# -----------------------------------------------------
# Compile with 3-D Solver and Related Programs
# -----------------------------------------------------

SV_USE_SOLVERIO = 1
SV_USE_THREEDSOLVER = 1
SV_USE_PRESOLVER = 1
SV_USE_POSTSOLVER = 1

# -----------------------------------------------------
# Compile Flowsolver Modules
# -----------------------------------------------------

SV_THREEDSOLVER_USE_CORONARY = 1
SV_THREEDSOLVER_USE_CLOSEDLOOP = 1
SV_THREEDSOLVER_USE_VARWALL = 1
SV_THREEDSOLVER_USE_VTK = 1

# -----------------------------------------------------
# Compile with MPI
# -----------------------------------------------------

SV_USE_MPI = 1
SV_USE_OPENMPI = 0
SV_USE_MPICH = 0

# by default, only build with real mpi on windows
ifeq ($(CLUSTER), x64_cygwin)
  SV_USE_DUMMY_MPI = 0
endif
ifeq ($(CLUSTER), x64_linux)
  SV_USE_DUMMY_MPI = 1
endif
ifeq ($(CLUSTER), x64_macosx)
  SV_USE_DUMMY_MPI = 1
endif

# -----------------------------------------------------
# Compile with VTK
# -----------------------------------------------------

SV_USE_VTK = 1
SV_USE_VTK_SHARED = 0

# -----------------------------------------------------
# Compile with sparse, metis, nspcg
# -----------------------------------------------------

SV_USE_SPARSE = 1
SV_USE_METIS = 1
SV_USE_NSPCG = 1

# -----------------------------------------------------
# Compile with Optimization
# -----------------------------------------------------

MAKE_OPTIMIZED = 1
LINK_WITH_DEBUG = 1

# -----------------------------------------------------
# Static link
# -----------------------------------------------------

SV_STATIC_BUILD = 1

# if you need to override anything above for a given site, do it here
# -----------------------------------------------------------------------

ifeq ($(LOCAL_DIR_SITE_OVERRIDES),1)
-include site_overrides.mk
else
-include $(TOP)/site_overrides.mk
endif

# ----------------
# Target directory
# ----------------

TARGETDIR = .

#SV_EXTERNALS_VERSION_NUMBER = 2019.02
SV_EXTERNALS_VERSION_NUMBER = 2019.06
SV_VTK_OPENGL_VERSION=gl2

ifeq ($(CLUSTER), x64_cygwin)
    SV_LOWERCASE_CMAKE_BUILD_TYPE=release
    SV_CMAKE_BUILD_TYPE=Release
    OPEN_SOFTWARE_BINARIES_TOPLEVEL = C:/cygwin64/usr/local/svsolver/ext/$(SV_EXTERNALS_VERSION_NUMBER)/$(SV_LOWERCASE_CMAKE_BUILD_TYPE)/$(SV_VTK_OPENGL_VERSION)/bin/$(SV_COMPILER)/$(SV_COMPILER_VERSION)/x64
endif

ifeq ($(CLUSTER), x64_linux)
    SV_LOWERCASE_CMAKE_BUILD_TYPE=release
    SV_CMAKE_BUILD_TYPE=Release
    OPEN_SOFTWARE_BINARIES_TOPLEVEL = /usr/local/svsolver/ext/$(SV_EXTERNALS_VERSION_NUMBER)/$(SV_LOWERCASE_CMAKE_BUILD_TYPE)/$(SV_VTK_OPENGL_VERSION)/bin/$(SV_COMPILER)/$(SV_COMPILER_VERSION)/x64
endif

ifeq ($(CLUSTER), x64_macosx)
    SV_LOWERCASE_CMAKE_BUILD_TYPE=release
    SV_CMAKE_BUILD_TYPE=Release
    OPEN_SOFTWARE_BINARIES_TOPLEVEL = /usr/local/svsolver/ext/$(SV_EXTERNALS_VERSION_NUMBER)/$(SV_LOWERCASE_CMAKE_BUILD_TYPE)/$(SV_VTK_OPENGL_VERSION)/bin/$(SV_COMPILER)/$(SV_COMPILER_VERSION)/x64
endif

# -------------------------------------------
#   Release version numbers for SimVascular
# -------------------------------------------

SV_USE_WIN32_REGISTRY=1
SV_REGISTRY_TOPLEVEL=SIMVASCULAR

# if you need to override anything above, stuff it in global_overrides.mk
# -----------------------------------------------------------------------

ifeq ($(LOCAL_DIR_GLOBAL_OVERRIDES),1)
-include global_overrides.mk
else
-include $(TOP)/global_overrides.mk
endif

SV_MAJOR_VERSION ?= $(shell date +"%Y")
SV_MINOR_VERSION ?= $(shell date +"%m")
SV_PATCH_VERSION ?= $(shell date +"%d")
SV_MAJOR_VERSION_TWO_DIGIT ?= $(shell date +"%y")
SV_MAJOR_VER_NO = "$(SV_MAJOR_VERSION_TWO_DIGIT).$(SV_MINOR_VERSION)"
SV_FULL_VER_NO = "$(SV_MAJOR_VERSION_TWO_DIGIT).$(SV_MINOR_VERSION).$(SV_PATCH_VERSION)"

ifeq ($(CLUSTER),x64_cygwin)
  SV_VERSION  = SimVascular
  SV_PLATFORM = x64
  SV_POSTFIX=
  SV_OS=windows
endif
ifeq ($(CLUSTER),x64_linux)
  SV_VERSION  = simvascular
  SV_PLATFORM = x64
  SV_POSTFIX=
  SV_OS=linux
endif
ifeq ($(CLUSTER),x64_macosx)
  SV_VERSION  = simvascular
  SV_PLATFORM = x64
  SV_POSTFIX=
  SV_OS=macosx
endif

# --------------
# Global defines
# --------------

GLOBAL_DEFINES = -DSV_VERSION=\"$(SV_VERSION)\" \
                 -DSV_MAJOR_VERSION=\"$(SV_MAJOR_VERSION)\" \
                 -DSV_MINOR_VERSION=\"$(SV_MINOR_VERSION)\" \
                 -DSV_PATCH_VERSION=\"$(SV_PATCH_VERSION)\" \
                 -DSV_MAJOR_VER_NO=\"$(SV_MAJOR_VER_NO)\" \
                 -DSV_FULL_VER_NO=\"$(SV_FULL_VER_NO)\" \
                 -DSV_REGISTRY_TOPLEVEL=\"$(SV_REGISTRY_TOPLEVEL)\"

ifeq ($(SV_USE_WIN32_REGISTRY), 1)
  GLOBAL_DEFINES += -DSV_USE_WIN32_REGISTRY
endif

ifeq ($(SV_STATIC_BUILD),1)
  GLOBAL_DEFINES += -DSV_STATIC_LINK -DSV_STATIC_BUILD
endif

ifeq ($(SV_GLOBALS_SHARED),1)
  GLOBAL_DEFINES += -DSV_GLOBALS_SHARED
endif

ifeq ($(CLUSTER), x64_cygwin)
   GLOBAL_DEFINES += -DSV_USE_NOTIMER -DWINDOWS -DWIN32
endif

ifeq ($(CLUSTER), x64_linux)
   GLOBAL_DEFINES += -DSV_USE_NOTIMER -DUNIX
endif

ifeq ($(CLUSTER), x64_macosx)
   GLOBAL_DEFINES += -DSV_USE_NOTIMER -DUNIX
endif

ifeq ($(SV_USE_ZLIB),1)
  GLOBAL_DEFINES += -DSV_USE_ZLIB
endif

ifeq ($(SV_USE_VTK),1)
  GLOBAL_DEFINES += -DSV_USE_VTK
endif

# ----------------------------------
# Platform-specific compiler options
# ----------------------------------

ifeq ($(CLUSTER), x64_cygwin)
  ifeq ($(CXX_COMPILER_VERSION), msvc-19.0)
	include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/compiler.vs19.10.x64_cygwin.mk
  endif
  ifeq ($(CXX_COMPILER_VERSION), msvc-19.16)
	include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/compiler.vs19.16.x64_cygwin.mk
  endif
  ifeq ($(FORTRAN_COMPILER_VERSION), ifort)
	include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/compiler.ifort.x64_cygwin.mk
        GLOBAL_DEFINES += -DSV_WRAP_FORTRAN_IN_CAPS_NO_UNDERSCORE
  endif
  ifeq ($(CXX_COMPILER_VERSION), mingw-gcc)
	include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/compiler.mingw-gcc.x64_cygwin.mk
  endif
  ifeq ($(FORTRAN_COMPILER_VERSION), mingw-gfortran)
	include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/compiler.mingw-gfortran.x64_cygwin.mk
        GLOBAL_DEFINES += -DSV_WRAP_FORTRAN_IN_LOWERCASE_WITH_UNDERSCORE
  endif
endif

ifeq ($(CLUSTER), x64_linux)
  ifeq ($(CXX_COMPILER_VERSION), icpc)
	include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/compiler.icpc.x64_linux.mk
  endif
  ifeq ($(FORTRAN_COMPILER_VERSION), ifort)
	include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/compiler.ifort.x64_linux.mk
        GLOBAL_DEFINES += -DSV_WRAP_FORTRAN_IN_LOWERCASE_WITH_UNDERSCORE
  endif
  ifeq ($(CXX_COMPILER_VERSION), gcc)
	include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/compiler.gcc.x64_linux.mk
  endif
  ifeq ($(FORTRAN_COMPILER_VERSION), gfortran)
	include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/compiler.gfortran.x64_linux.mk
        GLOBAL_DEFINES += -DSV_WRAP_FORTRAN_IN_LOWERCASE_WITH_UNDERSCORE
  endif
endif

ifeq ($(CLUSTER), x64_macosx)
  ifeq ($(CXX_COMPILER_VERSION), clang)
	include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/compiler.clang.x64_macosx.mk
  endif
  ifeq ($(FORTRAN_COMPILER_VERSION), ifort)
	include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/compiler.ifort.x64_macosx.mk
        GLOBAL_DEFINES += -DSV_WRAP_FORTRAN_IN_LOWERCASE_WITH_UNDERSCORE
  endif
  ifeq ($(CXX_COMPILER_VERSION), gcc)
	include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/compiler.gcc.x64_macosx.mk
  endif
  ifeq ($(FORTRAN_COMPILER_VERSION), gfortran)
	include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/compiler.gfortran.x64_macosx.mk
        GLOBAL_DEFINES += -DSV_WRAP_FORTRAN_IN_LOWERCASE_WITH_UNDERSCORE
  endif
endif

# --------------------------------
# build directory for object files
# --------------------------------

BUILD_DIR = obj/$(CLUSTER)/$(CXX_COMPILER_VERSION)-$(FORTRAN_COMPILER_VERSION)
LIB_BUILD_DIR = $(CLUSTER)/$(CXX_COMPILER_VERSION)-$(FORTRAN_COMPILER_VERSION)
BUILD_MPI_DIR = obj/$(CLUSTER)/$(CXX_COMPILER_VERSION)-$(FORTRAN_COMPILER_VERSION)/$(MPI_NAME)
LIB_MPI_BUILD_DIR = $(CLUSTER)/$(CXX_COMPILER_VERSION)-$(FORTRAN_COMPILER_VERSION)

# ---------------------
# Local lib directories
# ---------------------

LIBDIRS =
SHARED_LIBDIRS =
EXECDIRS =
LOCAL_INCDIR = $(TOP)/../Code/Source/Include/Make

ifeq ($(SV_USE_THREEDSOLVER),1)
     LIBDIRS += ../Code/FlowSolvers/ThreeDSolver
     SOLVERIO_INCDIR = -I $(TOP)/../Code/FlowSolvers/ThreeDSolver/SolverIO -I $(TOP)/../Code/FlowSolvers/Include/Make
     THREEDSOLVER_INCDIR = -I $(TOP)/../Code/FlowSolvers/ThreeDSolver
     EXECDIRS += ../Code/FlowSolvers/ThreeDSolver
endif

ifeq ($(SV_USE_SOLVERIO),1)
     LIBDIRS += ../Code/FlowSolvers/ThreeDSolver/SolverIO
     THREEDSOLVER_INCDIR = -I $(TOP)/../Code/FlowSolvers/ThreeDSolver
endif

SUBDIRS         = $(LIBDIRS) $(EXECDIRS)

# -------------------------
# Local include directories
# -------------------------

LOCAL_SUBDIRS   = $(LIBDIRS) $(SHARED_LIBDIRS) ../Code/Source/Include ../Code/Source/Include/Make
LOCAL_INCDIR    := $(foreach i, ${LOCAL_SUBDIRS}, -I$(TOP)/$(i))

# include path to find libs when linking
ifeq ($(CLUSTER), x64_cygwin)
  TOP_WINDOWS_STYLE=$(shell cygpath -m $(TOP))
  GLOBAL_LFLAGS 	 += $(LIBPATH_COMPILER_FLAG)$(TOP_WINDOWS_STYLE)/Lib/$(LIB_BUILD_DIR)
else
  GLOBAL_LFLAGS 	 += $(LIBPATH_COMPILER_FLAG)$(TOP)/Lib/$(LIB_BUILD_DIR)
endif

LFLAGS 	 = $(GLOBAL_LFLAGS)

#
# ThirdParty software that must be built
# from source if used.
#

# -----------------------------------------
# ***  Optional Open Source Packages    ***
# ***   (less restrictive licenses)     ***
# *** (e.g. MIT or BSD or Apache 2.0)   ***
# -----------------------------------------

# ------
# Sparse
# ------

ifeq ($(SV_USE_SPARSE),1)
  THIRD_PARTY_LIBDIRS += ../Code/ThirdParty/sparse
  SPARSE_TOP = $(TOP)/../Code/ThirdParty/sparse
  SPARSE_INCDIR  = -I $(SPARSE_TOP)
  SPARSE_LIBS    = $(SVLIBFLAG)_simvascular_thirdparty_sparse$(LIBLINKEXT)
endif

# ----
# zlib
# ----

ifeq ($(SV_USE_ZLIB),1)
  SV_LIB_THIRDPARTY_ZLIB_NAME=_simvascular_thirdparty_zlib
  THIRD_PARTY_LIBDIRS += ../Code/ThirdParty/zlib
  ZLIB_TOP = $(TOP)/../Code/ThirdParty/zlib
  ZLIB_INCDIR  = -I $(ZLIB_TOP)
  ZLIB_LIBS    = $(SVLIBFLAG)$(SV_LIB_THIRDPARTY_ZLIB_NAME)$(LIBLINKEXT)
endif

# -----------------------------------------
# ***  Optional Open Source Packages    ***
# ***           (LGPL code)             ***
# -----------------------------------------

# ------
# NSPCG
# ------

ifeq ($(SV_USE_NSPCG),1)
  THIRD_PARTY_LIBDIRS += ../Code/ThirdParty/nspcg
  NSPCG_TOP = $(TOP)/../Code/ThirdParty/nspcg
  NSPCG_INCDIR  = -I $(NSPCG_TOP)
  NSPCG_LIBS    = $(SVLIBFLAG)_simvascular_thirdparty_nspcg$(LIBLINKEXT)
endif

# -----------------------------------------
# ***  Optional Open Source Packages    ***
# ***  (not free for commercial use)    ***
# -----------------------------------------

# ----
# svLS
# ----

ifeq ($(SV_USE_DUMMY_SVLS),1)
    SVLS_DEFS   = 
    SVLS_INCDIR = -I ../svLS
    SVLS_LIBS   = $(SVLIBFLAG)_simvascular_dummy_svLS$(LIBLINKEXT)
endif

ifeq ($(SV_USE_SOURCE_CODE_SVLS),1)
    SVLS_DEFS   = 
    SVLS_INCDIR = -I ../svLS
    SVLS_LIBS   = $(SVLIBFLAG)_simvascular_svLS_$(MPI_NAME)$(LIBLINKEXT)
endif

# -----
# Metis
# -----

ifeq ($(SV_USE_METIS),1)
  THIRD_PARTY_LIBDIRS += ../Code/ThirdParty/metis
  METIS_TOP = $(TOP)/../Code/ThirdParty/metis
  METIS_INCDIR  = -I $(METIS_TOP)
  METIS_LIBS    = $(SVLIBFLAG)_simvascular_thirdparty_metis$(LIBLINKEXT)
endif


#
# ThirdParty software included from externals
#

# ---------------------------------------
# ***  Required Open Source Packages  ***
# ***  (no commercial restrictions)   ***
# ---------------------------------------

# ---------------------
# Visualization toolkit
# ---------------------

ifeq ($(SV_USE_VTK),1)

  ifeq ($(CLUSTER), x64_cygwin)
	include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/vtk.$(SV_VTK_OPENGL_VERSION).x64_cygwin.mk
  endif

  ifeq ($(CLUSTER), x64_linux)
	include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/vtk.$(SV_VTK_OPENGL_VERSION).x64_linux.mk
  endif

  ifeq ($(CLUSTER), x64_macosx)
	include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/vtk.$(SV_VTK_OPENGL_VERSION).x64_macosx.mk
  endif

endif

# -----------------------------------------
# ***  Optional Open Source Packages    ***
# ***   (less restrictive licenses)     ***
# *** (e.g. MIT or BSD or Apache 2.0)   ***
# -----------------------------------------

# -----
# MPI
# -----

ifeq ($(SV_USE_MPI),1)

ifeq ($(SV_USE_DUMMY_MPI),1)

  MPI_NAME      = nompi
  MPI_TOP       = ../dummyMPI
  MPI_INCDIR    = -I $(MPI_TOP)
  MPI_LIBS      = $(SVLIBFLAG)_simvascular_dummy_mpi$(LIBLINKEXT)
  MPI_SO_PATH   = 
  MPIEXEC_PATH  = 
  MPIEXEC       =

else

  MPI_NAME ?= mpi

  ifeq ($(CLUSTER), x64_cygwin)
	include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/msmpi.x64_cygwin.mk
  endif

  # on linux, use the OS installed version of mpich2
  ifeq ($(CLUSTER), x64_linux)
    ifeq ($(SV_USE_OPENMPI),1)
      include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/openmpi.x64_linux.mk
    endif
    ifeq ($(SV_USE_MPICH),1)
      include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/mpich.x64_linux.mk
    endif
  endif

  ifeq ($(CLUSTER), x64_macosx)
    ifeq ($(SV_USE_OPENMPI),1)
      include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/openmpi.x64_macosx.mk
    endif
    ifeq ($(SV_USE_MPICH),1)
      include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/mpich.x64_macosx.mk
    endif
  endif

endif

endif

# -----------------------------------------
# ***  Optional Open Source Packages    ***
# ***           (GPL code)              ***
# -----------------------------------------

# --------------------------------------
# ***  Optional Commercial Packages  ***
# --------------------------------------

# ------
# LesLib
# ------

ifeq ($(SV_USE_BINARY_LESLIB),1)

  ifeq ($(CLUSTER), x64_cygwin)
	  include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/leslib-1.5.x64_cygwin.mk
  endif

  ifeq ($(CLUSTER), x64_linux)
	  include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/leslib-1.5.x64_linux.mk
  endif

  #No leslib for mac osx

endif

ifeq ($(SV_USE_DUMMY_LESLIB),1)

  LESLIB_INCDIR = 
  LESLIB_LIBS   = 

  ifeq ($(CLUSTER), x64_cygwin)
    LESLIB_DEFS   = -DACUSIM_NT -DACUSIM_WIN -DACUSIM_WIN64
    LESLIB_LIBS   = 
  endif

  ifeq ($(CLUSTER), x64_linux)
    LESLIB_DEFS   = -DACUSIM_LINUX
  endif

  ifeq ($(CLUSTER), x64_macosx)
    LESLIB_DEFS   = -DACUSIM_LINUX
  endif

endif

# here's your chance to override package locations
# ------------------------------------------------

ifeq ($(LOCAL_DIR_PKG_OVERRIDES),1)
-include pkg_overrides.mk
else
-include $(TOP)/pkg_overrides.mk
endif

# --------------------------------
# define rules for file extensions
# --------------------------------

ifeq ($(CLUSTER), x64_cygwin)
	  include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/rules.x64_cygwin.mk
endif

ifeq ($(CLUSTER), x64_linux)
	  include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/rules.x64_linux.mk
endif

ifeq ($(CLUSTER), x64_macosx)
	  include $(TOP)/MakeHelpers/$(SV_EXTERNALS_VERSION_NUMBER)/rules.x64_macosx.mk
endif
