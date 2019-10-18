BOOST = /home/yutongq/SQUID/Tools/boost_1_64_0/
BAMTOOLS = /home/yutongq/SQUID/Tools/bamtools/usr/local/
# BAMTOOLS = /home/congm1/ocean/oceancong02/Software/bamtools-2.5.1/bin/
GLPK = /home/yutongq/SQUID/Tools/glpk-4.62/bin

CC = gcc
CXX = g++
INCLUDES = -g -I $(BAMTOOLS)/include/bamtools -I $(GLPK)/include -I $(BOOST)
CXXFLAGS = -std=c++11 $(INCLUDES)
LDADD = $(BAMTOOLS)/lib/libbamtools.a $(GLPK)/lib/libglpk.a
LDLIBS = -lz -lm
# RPATH =  $(BAMTOOLS)/lib/:$(GLPK)/lib/

SRCS = src/main.cpp src/ReadRec.cpp src/SegmentGraph.cpp src/WriteIO.cpp src/Config.cpp
TESTSRCs = src/ReadRec.cpp src/SegmentGraph.cpp src/WriteIO.cpp src/Config.cpp src/testall.cpp

all: bin/diploidsquid

bin/test: $(subst .cpp,.o,$(TESTSRCs))
	$(CXX) -o $@ $^ $(LDADD) $(LDLIBS)

bin/diploidsquid: $(subst .cpp,.o,$(SRCS))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDADD) $(LDLIBS)

clean:
	rm -f bin/squid src/*.o test
