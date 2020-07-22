######################################################################################################
#
# Makefile for compiling/mexing all C stuff (code & wrappers) - strictly gcc!
#
# You need the Matlab standalone 'mex' on your (system) path. You must specify MEXEXT (see below).
# Best run from Matlab command prompt (in C subdirectory) as:
#
# >> [status, result] = system(['make MEXEXT=' mexext ' ' <target>],'-echo')
#
#####################################################################################################

# C wrappers

CWRAP =           \
	genvar_mex.c

# The output (mex files). MEXEXT is the platform mex file extension, as returned by 'mexext' command.

BIN = $(patsubst %.c,../mex/%.$(MEXEXT),$(CWRAP))

# Standalone 'mex' must be on your (system) path

MEX = mex

# The -largeArrayDims is probably advisable and -std=c99 too as some C source might need C99
# functionality.libut handles the (undocumented) Matlab mex call 'utIsInterruptPending()' which
# allows trapping of Ctrl-C in mex files. Might not be present on some installations, apparently.
# The libm link is probably rendundant (mex probably does it anyway) - but that's what the 'man'
# for log2() says.

# MFLAGS = -O -largeArrayDims CFLAGS='$$CFLAGS -Wall -Werror -std=c99 -O3' -lut -lm
MFLAGS = -O -largeArrayDims CFLAGS='$$CFLAGS -Wall -Werror -O3'

.PHONY: all clean diag

all: $(BIN)

clean:
	rm -f $(BIN)

diag:
	@echo "Nothing to see here, please move along"
	@echo $(BIN)

$(BIN): ../mex/%.$(MEXEXT): %.c
	$(MEX) $(MFLAGS) $< -o $@
