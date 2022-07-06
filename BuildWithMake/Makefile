
UNAME := $(shell uname)
$(info UNAME="$(UNAME)") 

ifeq ($(UNAME),Darwin)
  include Makefile.darwin
else
  include Makefile.Linux
endif

src=../Code/Source/*.cxx
obj = $(src:.cpp=.o)

exe = svOneDSolver

$(exe): $(obj)
	$(CC) -o $(exe) $(obj) $(libs) 

.cpp.o:
	$(CC) $(copt) $(flags) -c $<

