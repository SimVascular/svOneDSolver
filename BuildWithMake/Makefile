
UNAME := $(shell uname)
$(info UNAME="$(UNAME)") 

ifeq ($(UNAME),Darwin)
  include Makefile.darwin
endif

ifeq ($(UNAME),CYGWIN_NT-10.0)
  include Makefile.cygwin
endif

src=../Code/Source/*.cxx
obj = $(src:.cpp=.o)

exe = svOneDSolver

$(exe): $(obj)
	$(CXX) -o $(exe) $(obj) $(libs) 

.cpp.o:
	$(CXX) $(copt) $(flags) -c $<

