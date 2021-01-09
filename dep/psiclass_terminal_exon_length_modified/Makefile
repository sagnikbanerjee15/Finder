CXX = g++
CXXFLAGS= -Wall -O3 #-g #-std=c++11 #-Wall #-g
#CXXFLAGS= -Wall -g #-std=c++11 #-Wall #-g
LINKPATH= -I./samtools-0.1.19 -L./samtools-0.1.19
LINKFLAGS = -lbam -lz -lm -lpthread 
DEBUG=
OBJECTS = stats.o subexon-graph.o 

all: subexon-info combine-subexons classes vote-transcripts junc grader trust-splice add-genename addXS

subexon-info: subexon-info.o $(OBJECTS)
	if [ ! -f ./samtools-0.1.19/libbam.a ] ; \
	        then \
		                cd samtools-0.1.19 ; make ;\
	fi ;
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) subexon-info.o $(LINKFLAGS)

combine-subexons: combine-subexons.o $(OBJECTS)
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) combine-subexons.o $(LINKFLAGS)

classes: classes.o constraints.o transcript-decider.o $(OBJECTS)
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) constraints.o transcript-decider.o classes.o $(LINKFLAGS)

trust-splice: trust-splice.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) trust-splice.o $(LINKFLAGS)

vote-transcripts: vote-transcripts.o 
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) vote-transcripts.o $(LINKFLAGS)

junc: junc.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) junc.o $(LINKFLAGS)

grader: grader.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) grader.o $(LINKFLAGS)

addXS: addXS.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) addXS.o $(LINKFLAGS)

add-genename: add-genename.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) add-genename.o $(LINKFLAGS)

subexon-info.o: SubexonInfo.cpp alignments.hpp blocks.hpp support.hpp defs.h stats.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
combine-subexons.o: CombineSubexons.cpp alignments.hpp blocks.hpp support.hpp defs.h stats.hpp SubexonGraph.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
stats.o: stats.cpp stats.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
subexon-graph.o: SubexonGraph.cpp SubexonGraph.hpp 
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
constraints.o: Constraints.cpp Constraints.hpp SubexonGraph.hpp alignments.hpp BitTable.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
transcript-decider.o: TranscriptDecider.cpp TranscriptDecider.hpp Constraints.hpp BitTable.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
classes.o: classes.cpp SubexonGraph.hpp SubexonCorrelation.hpp BitTable.hpp Constraints.hpp alignments.hpp TranscriptDecider.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
trust-splice.o: GetTrustedSplice.cpp alignments.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
vote-transcripts.o: Vote.cpp TranscriptDecider.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
junc.o: FindJunction.cpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
grader.o: grader.cpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
addXS.o: AddXS.cpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
add-genename.o: AddGeneName.cpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)

clean:
	rm -f *.o *.gch subexon-info combine-subexons trust-splice vote-transcripts junc grader add-genename addXS
