# Hey emacs, this is a -*- makefile -*-

# Copyright (c) 2014-2015 The Regents of the University of California.
# All Rights Reserved. 
#
# Portions of the code Copyright (c) 2009-2011 Open Source Medical 
# Software Corporation, University of California, San Diego.
# All rights reserved.
#
# Portions copyright (c) 1999-2007 Stanford University,
#   Nathan Wilson, Ken Wang, Charles Taylor.
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
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

TOP = ..

include $(TOP)/Makefile.cygwin

#
#  The following are helper scripts to create PC
#  distributions.  Note that some of the paths are hardcoded.
#

SV_CLUSTER = $(CLUSTER)

ifeq ($(CLUSTER),x64_cygwin)
  include $(TOP)/wix-3.11.x64_cygwin.mk
endif

SV_TIMESTAMP := $(shell date +"%Y-%m-%d")
SV_TIMESTAMP_PERIOD := $(shell date +"%y.%m.%d")
DIST_DIR_TOP = package/svonedsolver
DIST_DIR := $(DIST_DIR_TOP)/$(SV_TIMESTAMP)

REAL_DIST_DIR_TOP = package/svonedsolver
REAL_DIST_DIR := $(DIST_DIR_TOP)/$(SV_TIMESTAMP)

ifeq ($(CLUSTER), x64_cygwin)
  DIST_DIR_TOP := /cygdrive/q/package/svonedsolver
  DIST_DIR := /cygdrive/q/package/svonedsolver/$(SV_TIMESTAMP)
  DIST_DIR_WIN := "Q:/package/svonedsolver/$(SV_TIMESTAMP)"
else
  DIST_DIR_TOP := $(REAL_DIST_DIR_TOP)
  DIST_DIR := $(REAL_DIST_DIR) 
endif

ifneq ($(CLUSTER),x64_linux)
  all: use_subst windows-distribution msi msi-bundle-svonedsolver
else
  all: linux-distribution taritup taritup-svonedsolver
endif

use_subst:
ifeq ($(CLUSTER),x64_cygwin)
	-subst Q: /D
	-subst Q: $(shell cygpath -w `pwd`)
endif

veryclean: use_subst msiclean-svonedsolver msi-bundle-svonedsolver-clean
	-rm -f tar_files/*.tar
	-rm -f tar_files/*.tar.gz
	-rm -Rf $(DIST_DIR_TOP)
	-rm -Rf package
	-rm -Rf tmp/*

#
#  Create msi installer for simvascular svonedsolver
#

msiclean-svonedsolver:
	-rm -f simvascular-svonedsolver.msi
	-rm -f tmp/simvascular-svonedsolver.wixobj
	-rm -f tmp/simvascular-svonedsolver.wxs
	-rm -f simvascular-svonedsolver.wixpdb
	-rm -f tmp/uuidgen.exe*
	-rm -f tmp/uuidgen.obj*

msi-svonedsolver:	msiclean-svonedsolver
ifeq ($(CXX_COMPILER_VERSION),mingw-gcc)
	x86_64-w64-mingw32-gcc -o tmp/uuidgen.exe windows_msi_helpers/uuidgen.c -lole32
else
	CL /Fetmp/uuidgen.exe /Fotmp/uuidgen.obj windows_msi_helpers/uuidgen.c ole32.lib
endif
	-echo "Create WXS file..."
	-tclsh windows_msi_helpers/makewxs-svonedsolver.tcl $(SV_TIMESTAMP) $(SV_TIMESTAMP_PERIOD) $(REAL_DIST_DIR_TOP) 
	-echo "Run candle..."
	-$(WIXUI_PATH)/candle.exe -ext WixBalExtension -ext WixUIExtension -ext WixUtilExtension -out tmp/simvascular-svonedsolver.wixobj tmp/simvascular-svonedsolver.wxs
	-echo "Run light..."
	-$(WIXUI_PATH)/light.exe -ext WixBalExtension -ext WixUIExtension -ext WixUtilExtension -out simvascular-svonedsolver.msi tmp/simvascular-svonedsolver.wixobj

msi-bundle-svonedsolver-clean:
	-rm -f tmp/bundle.wixobj
	-rm -f simvascular-svonedsolver-bundle.wixpdb
	-rm -f simvascular-svonedsolver-bundle.exe
	-rm -f tmp/wix-simvascular-svonedsolver-bundle.wxs

msi-bundle-svonedsolver: msi-bundle-svonedsolver-clean mksedscript
	-echo "Create WXS file..."
	-sed -f tmp/sedscript windows_installer_packages/wix-simvascular-svonedsolver-bundle.wxs > tmp/wix-simvascular-svonedsolver-bundle.wxs
	-echo "Run candle..."
	-$(WIXUI_PATH)/candle.exe -ext WixBalExtension -ext WixUIExtension -ext WixUtilExtension -out tmp/bundle.wixobj tmp/wix-simvascular-svonedsolver-bundle.wxs
	-echo "Run light..."
	-$(WIXUI_PATH)/light.exe tmp/bundle.wixobj -ext WixBalExtension -ext WixUIExtension -ext WixUtilExtension -out simvascular-svonedsolver-bundle.exe

msi: msiclean-svonedsolver msi-svonedsolver

mksedscript:
	-mkdir -p tmp
	-echo "puts \"s/REPLACE_TIMESTAMP/[file tail $(DIST_DIR)]/g\"" | tclsh | cat > tmp/sedscript
	-echo "puts \"s/REPLACE_SV_TIMESTAMP/[file tail $(DIST_DIR)]/g\"" | tclsh | cat >> tmp/sedscript
	-echo "puts \"s/REPLACE_SV_OS/$(SV_OS)/g\"" | tclsh | cat >> tmp/sedscript
	-echo "puts \"s/REPLACE_SV_VERSION/svonedsolver/g\"" | tclsh | cat >> tmp/sedscript
	-echo "puts \"s/REPLACE_SV_PLATFORM/$(SV_PLATFORM)/g\"" | tclsh | cat >> tmp/sedscript
	-echo "puts \"s/REPLACE_SV_POSTFIX/$(SV_POSTFIX)/g\"" | tclsh | cat >> tmp/sedscript
	-echo "puts \"s/REPLACE_SV_MAJOR_VER_NO/$(SV_MAJOR_VER_NO)/g\"" | tclsh | cat >> tmp/sedscript
	-echo "puts \"s/REPLACE_SV_FULL_VER_NO/$(SV_FULL_VER_NO)/g\"" | tclsh | cat >> tmp/sedscript
	-echo "puts \"s/REPLACE_SV_CLUSTER/$(CLUSTER)/g\"" | tclsh | cat >> tmp/sedscript
	-echo "puts \"s/REPLACE_TCL_VERSION_NO/$(notdir $(TCL_LIBRARY))/g\"" | tclsh | cat >> tmp/sedscript
	-echo "puts \"s/REPLACE_TK_VERSION_NO/$(notdir $(TK_LIBRARY))/g\"" | tclsh | cat >> tmp/sedscript
	-echo "puts \"s/REPLACE_TCLSH/$(notdir $(TCLSH))/g\"" | tclsh | cat >> tmp/sedscript

common-pre-distribution: tmp/sedscript
	-date > tmp/release-date
	-rm -f tar_files/*.tar
	-rm -f tar_files/*.tar.gz
	-rm -Rf $(DIST_DIR_TOP)
	-mkdir -p $(DIST_DIR)
	-date > $(DIST_DIR)/release-date
	-chmod -R a+rwx $(DIST_DIR_TOP)/../../package
	-echo "puts \"[file tail $(DIST_DIR)]\"" | tclsh > tmp/timestamp
	-cp -f License.rtf $(DIST_DIR)

copy-sv-executables:
	-cp -f $(TOP)/svOneDSolver.exe $(DIST_DIR)/svOneDSolver.exe

windows-distribution: mksedscript common-pre-distribution copy-sv-executables 

taritup-svonedsolver:
	-mkdir tar_files
	-rm -Rf package
	-mkdir -p package
	-date > $(DIST_DIR)/release-date
	-chmod -R a+rwx package
	-echo "puts \"[file tail $(DIST_DIR)]\"" | tclsh > tmp/timestamp
	-cp -f tmp/timestamp $(DIST_DIR)/
	-tar --preserve-permissions --transform 's,$(DIST_DIR)/,,g'  -cvf tar_files/simvascular-svonedsolver-$(SV_OS)-$(SV_PLATFORM).`awk '{print $$2"-"$$3"-"$$6}' tmp/release-date`.tar package
	-gzip tar_files/simvascular-svonedsolver-$(SV_OS)-$(SV_PLATFORM).`awk '{print $$2"-"$$3"-"$$6}' tmp/release-date`.tar

taritup:
	-mkdir tar_files
	-tar --preserve-permissions -cvf tar_files/simvascular-$(SV_OS)-$(SV_PLATFORM).`awk '{print $$2"-"$$3"-"$$6}' tmp/release-date`.tar package
	-gzip tar_files/simvascular-$(SV_OS)-$(SV_PLATFORM).`awk '{print $$2"-"$$3"-"$$6}' tmp/release-date`.tar
