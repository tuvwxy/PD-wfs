all: pc linux mac
NAME=wfs~
PDVER=0.42-5
# ----------------------- NT -----------------------
# For MinGW Compiler
pc: pd_nt
pd_nt: $(NAME).dll

.SUFFIXES: .dll

PDNT = C:/pd
# PDNTCFLAGS = -DPD -DNT -W3 -WX -Wall -Werror -Wno-unused -mms-bitfields \
#     -Wno-parentheses -Wno-switch -O6 -funroll-loops \
#    -fomit-frame-pointer -pedantic
PDNTCFLAGS = -ansi -std=c99 -pedantic -Wall

PDNTINCLUDE = -I. -I$(PDNT)/src
PDNTLDFLAGS = --export_dynamic -shared
PDNTLIB = $(PDNT)/bin/pd.dll

.c.dll:
	gcc $(PDNTCFLAGS) $(PDNTINCLUDE) -g -c $*.c
	gcc $(PDNTLDFLAGS) -o bin/$*.dll $*.o $(PDNTLIB)
	rm -r $*.o

# ----------------------- IRIX 5.x -----------------------

pd_irix5: $(NAME).pd_irix5

.SUFFIXES: .pd_irix5

SGICFLAGS5 = -o32 -DPD -DUNIX -DIRIX -O2

SGIINCLUDE =  -I../../src

.c.pd_irix5:
	$(CC) $(SGICFLAGS5) $(SGIINCLUDE) -o $*.o -c $*.c
	ld -elf -shared -rdata_shared -o $*.pd_irix5 $*.o
	rm $*.o

# ----------------------- IRIX 6.x -----------------------

pd_irix6: $(NAME).pd_irix6

.SUFFIXES: .pd_irix6

SGICFLAGS6 = -n32 -DPD -DUNIX -DIRIX -DN32 -woff 1080,1064,1185 \
	-OPT:roundoff=3 -OPT:IEEE_arithmetic=3 -OPT:cray_ivdep=true \
	-Ofast=ip32

.c.pd_irix6:
	$(CC) $(SGICFLAGS6) $(SGIINCLUDE) -o $*.o -c $*.c
	ld -n32 -IPA -shared -rdata_shared -o $*.pd_irix6 $*.o
	rm $*.o

# ----------------------- LINUX i386 -----------------------
linux: pd_linux
pd_linux: $(NAME).pd_linux

.SUFFIXES: .pd_linux

LINUXCFLAGS = -DPD -O2 -funroll-loops -fomit-frame-pointer -fPIC \
    -Wall -W -Wshadow -Wstrict-prototypes \
    -Wno-unused -Wno-parentheses -Wno-switch $(CFLAGS)

LINUXINCLUDE =  -I.

.c.pd_linux:
	$(CC) $(LINUXCFLAGS) $(LINUXINCLUDE) -o $*.o -c $*.c
	$(CC) -export_dynamic -shared -o bin/$*.pd_linux $*.o -lc -lm
	strip --strip-unneeded bin/$*.pd_linux
	rm -f $*.o

# ----------------------- Mac OSX -----------------------
mac: pd_fat
pd_ppc: $(NAME).d_ppc
pd_fat: $(NAME).d_fat

.SUFFIXES: .d_ppc .d_fat

DARWINCFLAGS = -DPD -O2 -Wall -W -Wshadow -Wstrict-prototypes \
    -Wno-unused -Wno-parentheses -Wno-switch
DARWINPD = /Applications/Pd-$(PDVER).app/Contents/Resources
DARWININCLUDE = -I. -I../../src -I$(DARWINPD)/src


.c.d_ppc:
	$(CC) $(DARWINCFLAGS) $(DARWININCLUDE) -o $*.o -c $*.c
	$(CC) -bundle -undefined suppress -flat_namespace -o bin/$(NAME).pd_darwin $*.o 
	rm -f $*.o
	cp bin/$(NAME).pd_darwin $(DARWINPD)/extra/$(NAME)

.c.d_fat:
	$(CC) -arch i386 -arch ppc $(DARWINCFLAGS) $(DARWININCLUDE) -o $*.o -c $*.c
	$(CC) -arch i386 -arch ppc -bundle -undefined suppress -flat_namespace -o bin/$*.d_fat $*.o 
	rm -f $*.o
	cp bin/$(NAME).d_fat $(DARWINPD)/extra/$(NAME)/

# ----------------------------------------------------------

clean:
	rm -f *.o *.pd_* so_location
