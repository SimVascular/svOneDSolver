
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

# Set option to use Skyline solver
flags=/DUSE_SKYLINE

exe = svOneDSolver

$(exe): $(obj)
	$(CXX) $(OPT_FLAGS) $(flags) -o $(exe) $(obj) $(libs) 

clean:
	rm *.obj *.exe

