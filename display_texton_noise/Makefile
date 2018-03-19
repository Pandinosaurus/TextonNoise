####### Variables
CPP      = g++
BINDIR   = bin
OBJDIR   = obj
SRCDIR   = src
LIBDIR   = lib
CXXFLAGS = -Wall -O3 -fopenmp -fno-strict-aliasing
INCPATH  = -Ilib/AntTweakBar/include -Iusr/X11R6/include
LDFLAGS  = -lglut -lGL -lGLEW -lX11 -lm -lpng -lGLU
LIBS     = lib/AntTweakBar/lib/libAntTweakBar.a

####### Files
# source files without extension:
TARGET_FILES = $(notdir $(basename $(shell find $(SRCDIR)/ -type f -exec grep -l "int main" {} \;)))
SRC_FILES    = $(patsubst $(SRCDIR)/%.cpp,%,$(shell find $(SRCDIR)/ -name '*.cpp' -type f))
SRC_FILES   += $(patsubst $(SRCDIR)/%.c,%,$(shell find $(SRCDIR)/ -name '*.c' -type f))
OBJ_FILES    = $(addprefix $(OBJDIR)/,$(addsuffix .o, $(filter-out $(TARGET_FILES), $(SRC_FILES))))

# name of the application:
TARGET       = $(addprefix $(BINDIR)/,$(TARGET_FILES))

####### Build rules
all: $(TARGET)

$(BINDIR)/%: $(OBJDIR)/%.o $(OBJ_FILES) $(LIBS)
	@echo "===== Link $@ ====="
	@mkdir -p $(@D)
	$(CXX) -o $@ $^ $(LDFLAGS) $(CXXFLAGS)
	@cp shader/* $(BINDIR)

$(OBJDIR)/%.o: src/%.cpp
	@echo "===== Compile $< ====="
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCPATH) -c $< -o $@

$(OBJDIR)/%.o: src/%.c
	@echo "===== Compile $< ====="
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCPATH) -c $< -o $@

clean:
	@echo "===== Clean $< ====="
	@$(MAKE) clean -C lib/AntTweakBar/src
	@rm -rf $(BINDIR) $(OBJDIR) 

####### Specific rules
lib/AntTweakBar/lib/libAntTweakBar.a:
	@$(MAKE) -C lib/AntTweakBar/src
