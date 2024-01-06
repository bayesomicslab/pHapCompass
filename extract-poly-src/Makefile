
# part of HAPCUT2 MAKEFILE
default: all

CC=gcc -g -O3 -Wall -D_GNU_SOURCE
CFLAGS=-c -Wall

# DIRECTORIES
B=build
H=hairs-src
X=hapcut2-src
HTSLIB=sub/htslib
SAMTOOLS=sub/samtools
T=test
# below is the path to CUnit directory, change if need be
CUNIT=/usr/include/CUnit

all: $(B)/extractHAIRS 

# if samtools makefile not present, then submodules have not yet been downloaded (init & updated)
# first check if git present, else print error message
# credit to lhunath for the one-line check for git below http://stackoverflow.com/questions/592620/check-if-a-program-exists-from-a-bash-script
$(SAMTOOLS)/Makefile:
	@hash git 2>/dev/null || { echo >&2 "Git not installed (required for auto-download of required submodules). Install git and retry, or manually download samtools 1.2 and htslib 1.2.1 and place the unzipped folders in \"submodules\" directory with names matching their Makefile variables (default \"samtools\" and \"htslib\")"; exit 1; }
	git submodule init
	git submodule update

# just call the samtools makefile target because it downloads both htslib and samtools submodules
$(HTSLIB)/Makefile: $(SAMTOOLS)/Makefile

$(SAMTOOLS)/libbam.a: $(SAMTOOLS)/Makefile
	echo "Building Samtools bam library..."
	make -C $(SAMTOOLS) lib

$(HTSLIB)/libhts.a: $(HTSLIB)/Makefile
	echo "Building HTSlib libraries..."
	make -C $(HTSLIB) lib-static

# BUILD HAIRS

#temporarily removed -O2 flag after -I$(HTSLIB)
$(B)/extractHAIRS: $(B)/bamread.o $(B)/hashtable.o $(B)/readvariant.o $(B)/readfasta.o $(B)/hapfragments.o $(H)/extracthairs.c $(SAMTOOLS)/libbam.a $(HTSLIB)/libhts.a $(H)/parsebamread.c $(H)/realignbamread.c $(H)/nw.c $(H)/realign_pairHMM.c $(H)/estimate_hmm_params.c | $(B)
	$(CC) -I$(SAMTOOLS) -I$(HTSLIB) -g $(B)/bamread.o $(B)/hapfragments.o $(B)/hashtable.o $(B)/readfasta.o $(B)/readvariant.o -o $(B)/extractHAIRS $(H)/extracthairs.c  -L$(SAMTOOLS) -L$(HTSLIB) -pthread -lhts -lbam -lm -lz
#temporarily removed -O2 flag after -I$(HTSLIB)

$(B)/hapfragments.o: $(H)/hapfragments.c $(H)/hapfragments.h $(H)/readvariant.h | $(B)
	$(CC) -c $(H)/hapfragments.c -o $(B)/hapfragments.o

$(B)/readvariant.o: $(H)/readvariant.c $(H)/readvariant.h $(H)/hashtable.h $(H)/hashtable.c | $(B)
	$(CC) -c -I$(HTSLIB) $(H)/readvariant.c -o $(B)/readvariant.o

$(B)/bamread.o: $(H)/bamread.h $(H)/bamread.c $(H)/readfasta.h $(H)/readfasta.c $(SAMTOOLS)/libbam.a $(HTSLIB)/libhts.a | $(B)
	$(CC) -I$(SAMTOOLS) -I$(HTSLIB) -c $(H)/bamread.c -o $(B)/bamread.o

$(B)/hashtable.o: $(H)/hashtable.h $(H)/hashtable.c | $(B)
	$(CC) -c $(H)/hashtable.c -o $(B)/hashtable.o

$(B)/readfasta.o: $(H)/readfasta.c $(H)/readfasta.h | $(B)
	$(CC) -c $(H)/readfasta.c -o $(B)/readfasta.o



# create build directory
$(B):
	mkdir -p $(B)

