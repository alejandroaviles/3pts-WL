# ----- Makefile FOR 3pts CODE -----
#

OPTIMIZE =
OPT1    =

opt_home := -std=c99 -Wno-missing-braces -Wno-missing-field-initializers -I/usr/local/include -L/usr/local/lib -lgsl -lfftw3 -lgslcblas -lm -g -O3 -std=gnu99 -ffast-math -funroll-loops

OPTIONS =  $(OPTIMIZE) \
	   $(OPT1)


EXEC = a

all: $(EXEC)

# C compiler:
CC = gcc
# With icc (Intel):
#CC = icc

# Default CFLAGS:
#~ INCLUDES    =  -I ./ -I ./libs/
CFLAGS = -g -O3 $(opt_home)
#
# With OpenMP:
#CFLAGS = -g -O3 -fopenmp $(OPTIONS)
# For icc use instead:
#CFLAGS = -g -O3 -qopenmp $(OPTIONS)


#
# Nothing to do below
#
H_PATH = libs

S_PATH = src


OBJS	= $(S_PATH)/main.o $(S_PATH)/functions.o \
    $(S_PATH)/tests.o $(S_PATH)/procedures.o $(S_PATH)/background.o \
    $(S_PATH)/libs.o $(S_PATH)/zetam.o $(S_PATH)/twobessel.o $(S_PATH)/utils.o
 

INCL	= $(S_PATH)/global.h $(S_PATH)/functions.h \
    $(S_PATH)/procedures.h $(S_PATH)/background.h  \
    $(S_PATH)/twobessel.h $(S_PATH)/utils.h $(S_PATH)/zetam.h

#~ $(EXEC): $(OBJS) 
#~ 	($(CC) $(OBJS) $(LIBS) $(CFLAGS) -o $@ -lm; cp $(EXEC) ../)
$(EXEC): $(OBJS) 
	($(CC) $(OBJS) $(LIBS) $(CFLAGS) -o $@ -lm)

$(OBJS): $(INCL)

clean:
	(rm -f $(OBJS) $(EXEC); rm -fR mgpt.dSYM; rm $(EXEC))
#~ clean:
#~ 	(rm -f $(OBJS) $(EXEC); rm -fR mgpt.dSYM; rm ../$(EXEC))


.PHONY : all clean check install uninstall
