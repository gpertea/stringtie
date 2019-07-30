/* SuperRead pipeline
 * Copyright (C) 2012  Genome group at University of Maryland.
 * 
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


// For usage, see --help
// #define DEBUG 1
// #define DEBUG120626
#define NEW_STUFF // Put in to get node-to-node connections
// #define KILLED111115
// #define KILL120102
// #define DEBUG120911
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <charb.hpp>
#include <unistd.h>
#include <signal.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <set>
#include <map>
#include <utility>
#include <vector>
#include <list>
#include <stack>
#include <iterator>

#include <jellyfish/err.hpp>
#include <misc.hpp>
#include <heap.hpp>
#include <exp_buffer.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jflib/multiplexed_parser.hpp>
#include <jflib/multiplexed_io.hpp>
#include <src2/joinKUnitigs_v3_cmdline.hpp>
#include <rb_tree.hpp>

namespace err = jellyfish::err;
using jellyfish::thread_exec;

#define FRONT_END 1
#define BACK_END 2

// The command line arguments
cmdline_parse args;

// Global / Constant / parameter. Merge into switches parsing?
static const int maxDiffInsertSizesForPrinting      = 10; // Const / Parameter
static const int maxTotAllowableMissingOnEnds       = 2; // Const / Parameter
static const int default_max_offset_considered_same = 5; // Changed from 5
// static const int max_offset_to_test                 = 10000;


struct overlapDataStruct
{
     int unitig1;
     int unitig2;
     int ahg;
     int bhg;
     char ori;
};

struct unitigLocStruct
{
     int unitig2;
     int frontEdgeOffset;
     char ori;
     bool operator>(const unitigLocStruct &rhs) const {
	  if(frontEdgeOffset > rhs.frontEdgeOffset)
	       return true;
	  if(frontEdgeOffset < rhs.frontEdgeOffset)
	       return false;
	  if (unitig2 > rhs.unitig2)
	       return true;
	  if (unitig2 < rhs.unitig2)
	       return false;
	  if (ori > rhs.ori)
	       return true;
	  if (ori < rhs.ori)
	       return false;
	  return false;
     }
     bool operator<(const unitigLocStruct &rhs) const {
	  if ((frontEdgeOffset == rhs.frontEdgeOffset) &&
	      (unitig2 == rhs.unitig2) &&
	      (ori == rhs.ori))
	       return false;
	  return !operator>(rhs);
     }
};

// In the following, lastOverlappingOffset gives us the offset for the node
//  where we have the largest overlap into the node
// It is set artificially to the length of the first unitig to start
struct abbrevUnitigLocStruct
{
     int            frontEdgeOffset;
     mutable unsigned short pathNum; 
     char           ori;

     bool operator<(const abbrevUnitigLocStruct& rhs) const {
	  if(ori != rhs.ori) return ori < rhs.ori;
	  return frontEdgeOffset < rhs.frontEdgeOffset;
     }
     bool operator==(const abbrevUnitigLocStruct& rhs) const {
	  return ori == rhs.ori && frontEdgeOffset == rhs.frontEdgeOffset;
     }
};

struct kuniToReadMatchStruct
{
     int matchRgnBegin;
     int matchRgnEnd;
     int ahg;
     int bhg;
     int kUnitigMatchBegin;
     int kUnitigMatchEnd;
     int orientedReadMatchBegin;
     int orientedReadMatchEnd;
     int readMatchBegin;
     int readMatchEnd;
     int kUnitigLength;
     int readLength;
     int kUnitigNumber;
     char ori;
};
struct unitigConnectionsForPathStruct
{
     int unitig1;
     int unitig2;
     int frontEdgeOffset1;
     int frontEdgeOffset2;
     char ori1;
     char ori2;
};

struct augmentedUnitigPathPrintStruct
{
     int unitig1;
     int frontEdgeOffset;
     int numOverlapsIn;
     int numOverlapsOut;
     int beginOffset;
     int endOffset;
     char ori;
};

struct kUniBeginOffsetStruct
{
     unsigned int kUnitig;
     int begin;
     int end;
};

// Global variables. Read only after initialization in main
long int                 *startOverlapByUnitig;
int                      *unitigLengths; // RO
struct overlapDataStruct *overlapData; // RO
long int                 *startOverlapIndexByUnitig2, *unitig2OverlapIndex; // RO
std::map<unsigned         int, struct kUniBeginOffsetStruct>    origUnitigToNewUnitigPlacement; // RO

double mean[256][256], stdev[256][256]; // RO. Clarify why 256x256. Convert to ExpBuffer

typedef heap<unitigLocStruct>::min min_heap;
typedef heap<unitigLocStruct>::max max_heap;
typedef std::map<unitigLocStruct, int>::iterator unitigLocMap_iterator;
typedef std::set<std::pair<int, int> >::iterator edge_iterator;
struct nodePair {
     int node1;
     int node2;
};
// RB tree data stuff
typedef std::set<abbrevUnitigLocStruct> unitig_ori_offsets;
typedef std::map<int, unitig_ori_offsets> unitig_to_ori_offsets;

/** Parse the overlap information. Each line contain the information
    about 1 read and mated reads (consecutive even/odd read number) are
    kept together.
*/
class overlap_parser : public multiplexed_parser<charb> {
     std::ifstream input_;

public:
     overlap_parser(const char* input_file, int nb_threads) :
	  multiplexed_parser<charb>(nb_threads, 100),
	  input_(input_file)
	  { start_parsing(); }

     virtual void parser_loop() {
	  charb extra;
	  bool  has_extra = false;

	  while(input_) {
	       elt e(elt_init());
	       size_type& i = e->nb_filled;
	       i = 0;
	       if(has_extra) {
		    e->elements[0].swap(extra);
		    ++i;
	       }
	       for( ; i < group_size(); ++i) {
		    if(!getline(input_, e->elements[i]))
			 break;
	       }
      
	       has_extra = false;
	       if(i > 0) {
		    long read_number = atol(e->elements[i-1] + 2);
		    if(read_number % 2 == 0 && i == group_size()) {
			 // Keep it for next iteration
			 extra.swap(e->elements[i-1]);
			 --i;
			 has_extra = true;
		    }
	       }
	  }
     }
};

// Main worker class
class KUnitigsJoinerThread {
     ExpBuffer<kuniToReadMatchStruct>                 evenReadMatchStructs;
     ExpBuffer<kuniToReadMatchStruct>                 oddReadMatchStructs;
     ExpandingBuffer<unsigned char>                   matchStructIsUsed;
     ExpandingBuffer<unitigConnectionsForPathStruct>  unitigConnectionsForPathData;
     ExpandingBuffer<augmentedUnitigPathPrintStruct>  augmentedUnitigPathPrintData;
     int                                              numUnitigPathPrintRecsOnPath;
     ExpandingBuffer<int>                             fwdStartIndices;
     ExpandingBuffer<int>                             revStartIndices;
     ExpandingBuffer<int>                             fwdNumIndices;
     ExpandingBuffer<int>                             revNumIndices;
     ExpandingBuffer<int>                             newNodeNumsFromOld;
     // unitigLocStruct                                 *unitigLocData1;
     // unitigLocStruct                                 *unitigLocData2;
     int                                              mateUnitig1;
     int                                              mateUnitig2;
     unsigned char                                    mateUnitig1ori;
     unsigned char                                    mateUnitig2ori;
     int                                              beginUnitig;
     int                                              endUnitig;
     unsigned char                                    beginUnitigOri;
     unsigned char                                    endUnitigOri;
     max_heap                                         backward_path_unitigs;
     std::set<unitigLocStruct>                        startingNodes;
     std::set<unitigLocStruct>                        endingNodes;
     int                                              startingNodeNumber;
     std::vector<unitigLocStruct>                     nodeArray;
     std::map<unitigLocStruct, int>                   nodeToIndexMap;
     std::set<std::pair<int, int>>                    edgeList;
     std::set<std::pair<int, int>>                    sortedEdgeList;
     std::vector<int>                                 unitigNodeNumbersForPath;
     std::vector<nodePair>                            fwdConnections;
     std::vector<nodePair>                            revConnections;
     int                                              curPathNum;
     int                                              treeSize;
     charb                                            outputString;
#ifdef KILL120103
     charb                                            stderrOutputString;
#endif
     char                                             rdPrefix[3];
     char                                             rdPrefixHold[3];
     long long                                        readNum;
     long long                                        readNumHold;
     int                                              approxNumPaths;
     double                                           insertLengthMeanBetweenKUnisForInsertGlobal;
     double                                           insertLengthStdevGlobal;
     // The following keeps track of the distance the 2 read mates are from the
     // ends of the k-unitigs at the end
     int                                              lengthAdjustment1;
     int                                              lengthAdjustment2;
     charb                                            superReadName;
     int                                              splitJoinWindowMin;
     int                                              splitJoinWindowMax;
     int                                              numUnitigConnectionsForPathData;
     int                                              numPairsInOneUnitig;
     int                                              numSimplyJoinable;
     int                                              numJoinableAfterRead1Analysis;
     int                                              numJoinableAfterBothReadAnalysis;
     int                                              numJoinableUnresolvedAtEnd;
     int                                              numUnjoinableMissingSequence;
     int                                              numUnjoinableOverMaxNodes;
     int                                              joinCode;
     static const char                               *joinCodeNames[6];
  //     int                                             *treeReinitList;
     int                                              maxPathNumUsed;
     bool                                             tooManyPossibleInsertLengths;

public:
     struct unitigPathPrintStruct
     {
	  int unitig1;
	  int frontEdgeOffset;
	  mutable int numOverlapsIn;
	  mutable int numOverlapsOut;
	  char ori;
  
	  bool operator<(const unitigPathPrintStruct& rhs) const {
	       if(frontEdgeOffset != rhs.frontEdgeOffset) return frontEdgeOffset < rhs.frontEdgeOffset;
	       if(unitig1 != rhs.unitig1) return unitig1 < rhs.unitig1;
	       return ori < rhs.ori;
	  }
     };
     typedef std::set<unitigPathPrintStruct> unitig_print_path;

private:
     unitig_to_ori_offsets                                   treeArr;
     unitig_print_path                                       treeArr2;

     friend struct unitigPathPrintStruct;

public:
     KUnitigsJoinerThread() :
	  mateUnitig1ori('F'), mateUnitig2ori('R'),
	  numPairsInOneUnitig(0), numSimplyJoinable(0), numJoinableAfterRead1Analysis(0), numJoinableAfterBothReadAnalysis(0),
	  numJoinableUnresolvedAtEnd(0), numUnjoinableMissingSequence(0), numUnjoinableOverMaxNodes(0)
	  {
	       rdPrefix[2] = rdPrefixHold[2] = '\0';
	  }

     int process_input_file(overlap_parser& ovp_parser, 
			    jflib::omstream& m_out, int index);
  
private:
     void generateSuperReadPlacementLinesForJoinedMates ();
     void updateMatchRecords(int readNum, ExpBuffer<char*>& flds);
     int joinKUnitigsFromMates (int insertLengthMean, int insertLengthStdev);
     void printIfGood (struct abbrevUnitigLocStruct *ptr);
     template<typename T>
     void printPathNode (const T& ptr);
     template<typename T>
     void completePathPrint (const T& ptr);
     int setSuperReadNameFromAugmentedPath (void);
     int getSuperReadLength(void);
     template<typename T>
     void funcToGetTreeSize (const T& ptr); // Adds 1 to treeSize each time
     void findSingleReadSuperReads(char *readName, jflib::omstream& m_out);
     void getSuperReadsForInsert (jflib::omstream& m_out);
     int processKUnitigVsReadMatches (overlap_parser::stream& ovp_stream, jflib::omstream& m_out);
};
const char* KUnitigsJoinerThread::joinCodeNames[6] = { "NJ", "SU", "J", "A", "MS", "TMN" };


class KUnitigsJoiner : public thread_exec {
     overlap_parser       ovp_parser;
     std::ofstream        out;
     jflib::o_multiplexer multiplexer;

public:
     KUnitigsJoiner(const char* input_file, const char* output_file, int nb_threads) :
	  ovp_parser(input_file, nb_threads), out(output_file), 
	  multiplexer(&out, 3 * nb_threads, 4096)
	  { 
	       if(!out.good())
                   throw std::runtime_error(err::msg() << "Failed to open '" << output_file << "'" << err::no);
	  }

     virtual void start(int thid) {
	  KUnitigsJoinerThread joiner;
	  jflib::omstream      m_out(multiplexer);
	  int ret = joiner.process_input_file(ovp_parser, m_out, thid);
#ifndef DEBUG120626
	  std::cerr << "joiner thread " << thid << " returned " << ret << std::endl;
#else
	  ret = ret;
#endif
     }
};

bool firstNodeSort (struct nodePair val1, struct nodePair val2);
bool secondNodeSort (struct nodePair val1, struct nodePair val2);

bool unitigLocStructCompare (struct unitigLocStruct ptr1,
			     struct unitigLocStruct ptr2);
bool unitigLocStructCompareReversed (struct unitigLocStruct ptr1,
				     struct unitigLocStruct ptr2);
void loadKUnitigTranslationTable(const char *fileName);
FILE *Fopen (const char *fn, const char *mode);
FILE *Popen (const char *fn, const char *mode);
// The following returns the overlap length if it is greater than the
// existing largest overlap on the end. It returns -1 if not.
int getOvlLenFromOvlIndicesPlus (long int maxOvlIndex, long int j, int maxOvlLen, int whichEnd);
long int findOtherOverlapIndex (long int ovlIndex1);

int getInt (const char *fname);


// TODO: merge and template the following two functions
unitig_ori_offsets::iterator find_within(unitig_ori_offsets& tree,
					 abbrevUnitigLocStruct x, int delta) {
     x.frontEdgeOffset -= delta;
     auto res = tree.lower_bound(x);
     if(res == tree.end())
	  return res;
     if(res->frontEdgeOffset > x.frontEdgeOffset + 2 * delta)
	  return tree.end();
     return res;
}

KUnitigsJoinerThread::unitig_print_path::iterator find_within(KUnitigsJoinerThread::unitig_print_path& tree,
							      KUnitigsJoinerThread::unitigPathPrintStruct x, 
							      int delta) {
     x.frontEdgeOffset -= delta;
     auto res = tree.lower_bound(x);
     while (1) {
	  if(res == tree.end())
	       return res;
	  if(res->frontEdgeOffset > x.frontEdgeOffset + 2 * delta)
	       return tree.end();
	  if ((res->unitig1 != x.unitig1) || (res->ori != x.ori)) {
	       ++res;
	       continue; }
	  break;
     }
     return res;
}

#ifndef mallocOrDie
#define mallocOrDie(name, num, type) fprintf (stderr, "Allocating %lu bytes for %s.\n", (unsigned long) ((num) * sizeof ( type )), #name); \
     name = (type *) calloc (num, sizeof ( type ));			\
     if (name == NULL) { fprintf (stderr, "Couldn't allocate space for '%s'\nBye!\n", #name ); exit (-1); }
#endif

// #define DEBUG
int KUnitigsJoinerThread::process_input_file(overlap_parser& ovp_parser,
                                             jflib::omstream& m_out, int index)
{
     overlap_parser::stream ovp_stream(ovp_parser);
     int ret = processKUnitigVsReadMatches(ovp_stream, m_out);
#ifndef DEBUG120626
     fprintf (stderr, 
	      "Num pairs with both reads in same unitig: %d\n"
	      "Num pairs uniquely joinable: %d\n"
	      "Num pairs after disambiguation to beginning of insert: %d\n"
	      "Num pairs after disambiguation to end of insert: %d\n"
	      "Num still joinable but not uniquely joinable: %d\n"
	      "Num pairs unjoinable due to missing sequence: %d\n"
	      "Num pairs unjoinable because there are too many nodes: %d\n",
	      numPairsInOneUnitig, numSimplyJoinable, numJoinableAfterRead1Analysis, numJoinableAfterBothReadAnalysis, numJoinableUnresolvedAtEnd, numUnjoinableMissingSequence, numUnjoinableOverMaxNodes);
     fflush (stderr);
#endif
     return ret;
}

int main (int argc, char **argv)
{
     args.parse(argc, argv);   // Parse arguments
     FILE *infile;
     charb line(2000);
     int unitig1, unitig2;
     long int overlapCount = 0;
     int unitigNum, numUnitigs;
     static const int firstUnitigNum = 0; // Is it really needed?
     int numFlds;
     ExpBuffer<char*> flds;

     // Read in library information
     infile = Fopen (args.mean_and_stdev_by_prefix_file_arg, "r");
     while (fgets (line, 2000, infile)) {
	  getFldsFromLine(line, flds);
	  mean[(int)flds[0][0]][(int)flds[0][1]] = atof (flds[1]);
	  stdev[(int)flds[0][0]][(int)flds[0][1]] = atof (flds[2]);
     }
     fclose (infile);

     // Get the number of unitigs
     numUnitigs = getInt (args.num_kunitigs_file_arg) + 1;
     mallocOrDie (startOverlapByUnitig, numUnitigs + 1 + firstUnitigNum,long int);
     mallocOrDie (startOverlapIndexByUnitig2, numUnitigs + 1 + firstUnitigNum,long int);

     // Getting the unitig lengths
     mallocOrDie (unitigLengths, numUnitigs + 1 + firstUnitigNum, int);
     // Here we read in the unitig lengths, allowing for either type of length
     // format
     infile = Fopen (args.unitig_lengths_file_arg, "r");
     if (! fgets (line, 2000, infile))
       err::die(err::msg() << "File '" << args.unitig_lengths_file_arg << "' is of length 0 or can't be read");
     numFlds = getFldsFromLine (line, flds);
     rewind (infile);
     if (numFlds == 1) {
	  for (int i = 0, *iptr = unitigLengths + firstUnitigNum; i < numUnitigs; i++, iptr++) {
	       if(fscanf (infile, "%d\n", iptr) != 1)
                   err::die(err::msg() << "Failed to parse " << i << "th integer in: " << (unitigLengths + firstUnitigNum));
	  }
     }
     else {
	  while (fgets (line, 2000, infile)) {
	       getFldsFromLine (line, flds);
	       unitigLengths[atoi(flds[0])] = atoi (flds[1]);
	  }
     }
     fclose (infile);
     // End section getting unitig lengths

     if (args.kunitigs_translation_file_given)
	  loadKUnitigTranslationTable(args.kunitigs_translation_file_arg);

     // Map the overlaps file
     int fd = open(args.overlaps_file_arg, O_RDONLY);
     if(fd == -1) {
	  perror("open failed");
	  exit(1);
     }
     struct stat stat_buf;
     if(fstat(fd, &stat_buf) == -1) {
	  perror("stat failed");
	  exit(1);
     }
     off_t fileSize = stat_buf.st_size;
     if (fileSize == 0) // So that the 'mmap' doesn't fail
	  fileSize = 1;
     overlapData = (struct overlapDataStruct *)mmap(0, fileSize, PROT_READ, MAP_SHARED, fd, 0);
     if(overlapData == MAP_FAILED) {
	  perror("mmap failed");
	  exit(1);
     }
     close(fd);

     // Force it in memory by touching every page
     char *end = (char*)overlapData + stat_buf.st_size;
     char whatever = 0;
     for(char *ptr = (char*)overlapData; ptr < end; ptr += getpagesize())
	  whatever ^= *ptr;

     overlapCount=(long int)((double)stat_buf.st_size/(double)sizeof(struct overlapDataStruct)+.01);
   
     // Loading the overlap data (startOverlapByUnitig, overlapData)
     for(long int j=0;j<overlapCount;j++)
     {
	  int unitig1=overlapData[j].unitig1;
	  int unitig2=overlapData[j].unitig2;
#if 0
	  if (unitig1 > unitig2) 
	       continue;
	  else if ((unitig1 == unitig2) && (overlapData[j].ahg < 0))
	       continue;
#endif
	  if (unitig1 >= unitig2)
	       continue;
	  startOverlapByUnitig[unitig1]++;
	  startOverlapByUnitig[unitig2]++;
     }
     mallocOrDie (unitig2OverlapIndex, overlapCount,long int);
     for (unitigNum = 1; unitigNum < numUnitigs + 1 + firstUnitigNum; unitigNum++)
	  startOverlapByUnitig[unitigNum] += startOverlapByUnitig[unitigNum - 1];
     for (unitigNum = 0; unitigNum < numUnitigs + 1 + firstUnitigNum; unitigNum++)
	  startOverlapIndexByUnitig2[unitigNum] = startOverlapByUnitig[unitigNum];

// Unitig in the overlaps file
     for(long int j=0;j<overlapCount;j++)
     {
	  unitig1=overlapData[j].unitig1;
	  unitig2=overlapData[j].unitig2;

	  startOverlapByUnitig[unitig1]--;
	  startOverlapIndexByUnitig2[unitig2]--;
	  unitig2OverlapIndex[startOverlapIndexByUnitig2[unitig2]] = j;
     }

     free(startOverlapIndexByUnitig2);
     // Done loading the overlap data.
     // overlapData for indices in [startOverlapByUnitig[uni],startOverlapByUnitig[uni+1])
     // has the overlaps for unitig 'uni'

     KUnitigsJoiner joiners(args.input_file_arg, args.output_arg, args.threads_arg);
     joiners.exec_join(args.threads_arg);
     
     munmap(overlapData, fileSize);
     
     return (0);
}

void loadKUnitigTranslationTable(const char *fileName)
{
     FILE *infile = Fopen (fileName, "r");
     struct kUniBeginOffsetStruct kubos;
     charb line;
     ExpBuffer<char*> flds;
     unsigned long origUnitigNumber;

     while (fgets (line, 2000, infile)) {
	  getFldsFromLine(line, flds);
	  origUnitigNumber = atoll(flds[0]);
	  kubos.kUnitig = atoll(flds[1]);
	  kubos.begin = atoi(flds[2]);
	  kubos.end = atoi(flds[3]);
	  unitigLengths[origUnitigNumber] = -(abs(kubos.end - kubos.begin)+1);
	  origUnitigToNewUnitigPlacement[origUnitigNumber] = kubos;
     }
     fclose (infile);

     return;
}

// This loads a line of matches of reads to k-unitigs (output of findMatchesBetweenKUnitigsAndReads)
void KUnitigsJoinerThread::updateMatchRecords(int readNum, ExpBuffer<char*>& flds)
{
     int readLength = atoi(flds[1]);
     ExpBuffer<struct kuniToReadMatchStruct> *structs;
     if (readNum % 2 == 0) 
	  structs = &evenReadMatchStructs;
     else
	  structs = &oddReadMatchStructs;

     for(size_t i = 2; i < flds.size(); i += 3) {
	  structs->push_back(kuniToReadMatchStruct());
	  struct kuniToReadMatchStruct &kUTRMS = structs->back();
	  kUTRMS.ori           = *(flds[i+2]);
	  kUTRMS.ahg           = atoi (flds[i+1]);
	  kUTRMS.readLength    = readLength;
	  kUTRMS.kUnitigNumber = atoi (flds[i]);
	  kUTRMS.kUnitigLength = abs (unitigLengths[kUTRMS.kUnitigNumber]);
	  kUTRMS.bhg = kUTRMS.ahg + kUTRMS.readLength - kUTRMS.kUnitigLength;
	  if (kUTRMS.ahg <= 0)
	       kUTRMS.kUnitigMatchBegin = 0;
	  else
	       kUTRMS.kUnitigMatchBegin = kUTRMS.ahg;
	  if (kUTRMS.bhg > 0)
	       kUTRMS.kUnitigMatchEnd = kUTRMS.kUnitigLength;
	  else
	       kUTRMS.kUnitigMatchEnd = kUTRMS.kUnitigLength + kUTRMS.bhg;
	  // orientedReadMatchBegin is 0-based for now
	  kUTRMS.orientedReadMatchBegin = kUTRMS.kUnitigMatchBegin - kUTRMS.ahg;
	  kUTRMS.orientedReadMatchEnd = kUTRMS.kUnitigMatchEnd - kUTRMS.ahg;
	  if (kUTRMS.ori == 'F') {
	       kUTRMS.readMatchBegin = kUTRMS.orientedReadMatchBegin;
	       kUTRMS.readMatchEnd = kUTRMS.orientedReadMatchEnd; }
	  else {
	       kUTRMS.orientedReadMatchBegin = kUTRMS.readLength - kUTRMS.orientedReadMatchBegin;
	       kUTRMS.orientedReadMatchEnd = kUTRMS.readLength - kUTRMS.orientedReadMatchEnd;
	       kUTRMS.readMatchBegin = kUTRMS.orientedReadMatchEnd;
	       kUTRMS.readMatchEnd = kUTRMS.orientedReadMatchBegin; }
#ifdef DEBUG120626
	  fprintf (stderr, "readNum = %d, readMatchBegin = %d, readMatchEnd = %d, orientedReadMatchBegin = %d, orientedReadMatchEnd = %d, ahg = %d, bhg = %d, kUnitigNumber = %d\n", (int) readNum, kUTRMS.readMatchBegin, kUTRMS.readMatchEnd, kUTRMS.orientedReadMatchBegin, kUTRMS.orientedReadMatchEnd, kUTRMS.ahg, kUTRMS.bhg, kUTRMS.kUnitigNumber);
#endif
//	  printf ("Loop before k-unitig translations, i = %d, flds.size() = %d\n", (int)i, (int)flds.size());
	  if (! args.kunitigs_translation_file_given)
	       continue;
	  if (unitigLengths[kUTRMS.kUnitigNumber] > 0)
	       continue;
	  auto map_it = origUnitigToNewUnitigPlacement.find(kUTRMS.kUnitigNumber);
	  // The following should not be necessary, but just in case...
	  if (map_it == origUnitigToNewUnitigPlacement.end()) {
	       fprintf (stderr, "Error at 100, readNum = %d, kUnitig = %d\n", readNum, kUTRMS.kUnitigNumber);
	       return; }

	  int amtKUniBeforeOldKUni, amtKUniAfterOldKUni;
	  int newUnitigLength = unitigLengths[map_it->second.kUnitig];
	  if (map_it->second.begin < map_it->second.end) { // The k-uni lies in the new one in F dir.
	       amtKUniBeforeOldKUni = map_it->second.begin-1;
	       amtKUniAfterOldKUni = newUnitigLength - map_it->second.end;
	       kUTRMS.ahg += amtKUniBeforeOldKUni;
	       kUTRMS.bhg -= amtKUniAfterOldKUni;
	  }
	  else {
	       amtKUniBeforeOldKUni = map_it->second.end-1;
	       amtKUniAfterOldKUni =  newUnitigLength - map_it->second.begin;
	       int ahgHold = kUTRMS.ahg;
	       kUTRMS.ahg = amtKUniBeforeOldKUni - kUTRMS.bhg;
	       kUTRMS.bhg = - amtKUniAfterOldKUni - ahgHold;
	       // Change the orientation of the read in the k-unitig
	       kUTRMS.ori = (kUTRMS.ori == 'F') ? 'R' : 'F';
	  }
	  kUTRMS.kUnitigMatchBegin = (kUTRMS.ahg < 0) ? 0 : kUTRMS.ahg;
	  int tempValueAdjustment = (kUTRMS.bhg < 0) ? kUTRMS.bhg : 0;
	  kUTRMS.kUnitigMatchEnd = newUnitigLength + tempValueAdjustment;
	  int amtOfReadBeforeKUni = (kUTRMS.ahg > 0) ? 0 : -kUTRMS.ahg;
	  int amtOfReadAfterKUni  = (kUTRMS.bhg > 0) ? kUTRMS.bhg : 0;
	  if (kUTRMS.ori == 'F') {
	       kUTRMS.readMatchBegin = amtOfReadBeforeKUni;
	       kUTRMS.readMatchEnd = kUTRMS.readLength - amtOfReadAfterKUni;
	       kUTRMS.orientedReadMatchBegin = kUTRMS.readMatchBegin;
	       kUTRMS.orientedReadMatchEnd   = kUTRMS.readMatchEnd; }
	  else {
	       kUTRMS.readMatchBegin = amtOfReadAfterKUni;
	       kUTRMS.readMatchEnd = kUTRMS.readLength - amtOfReadBeforeKUni;
	       kUTRMS.orientedReadMatchBegin = kUTRMS.readMatchEnd;
	       kUTRMS.orientedReadMatchEnd   = kUTRMS.readMatchBegin; }
	  kUTRMS.kUnitigNumber = map_it->second.kUnitig;
	  kUTRMS.kUnitigLength = unitigLengths[kUTRMS.kUnitigNumber];
//	  printf ("At end of loop, i = %d, flds.size() = %d\n", (int)i, (int)flds.size());
	  //     fprintf (stderr, "After change: readMatchBegin = %d, readMatchEnd = %d, orientedReadMatchBegin = %d, orientedReadMatchEnd = %d, ahg = %d, bhg = %d, kUnitigNumber = %d\n", kUTRMS.readMatchBegin, kUTRMS.readMatchEnd, kUTRMS.orientedReadMatchBegin, kUTRMS.orientedReadMatchEnd, kUTRMS.ahg, kUTRMS.bhg, kUTRMS.kUnitigNumber);
     }
}

int KUnitigsJoinerThread::processKUnitigVsReadMatches (overlap_parser::stream& ovp_stream,
                                                       jflib::omstream& m_out)
{
     ExpBuffer<char*> flds;

     rdPrefixHold[0] = rdPrefixHold[1] = 0;
     evenReadMatchStructs.clear();
     oddReadMatchStructs.clear();
     
     // For each overlap line
     for( ; ovp_stream; ++ovp_stream) {
	  getFldsFromLine(*ovp_stream, flds);
	  // Parse read info
	  rdPrefix[0] = flds[0][0];
	  rdPrefix[1] = flds[0][1];
	  readNum = atoll (flds[0] + 2);

	  // Compute SuperRead if lonely read or got mate pair
	  if ((strcmp (rdPrefix, rdPrefixHold) != 0) ||
	      (readNum != readNumHold+1) ||
	      (readNum % 2 == 0)) {
	       if(rdPrefixHold[0] != '\0') {
		    approxNumPaths = 0;
		    joinCode = 0;
		    getSuperReadsForInsert(m_out);
		    m_out << jflib::endr;
		    // Set up and load the new data
		    evenReadMatchStructs.clear();
		    oddReadMatchStructs.clear();
	       }
	       strcpy (rdPrefixHold, rdPrefix);
	  }
	  updateMatchRecords(readNum, flds);
	  readNumHold = readNum;
     }
     if(!evenReadMatchStructs.empty() || !oddReadMatchStructs.empty()) {
	  joinCode = 0;
	  getSuperReadsForInsert(m_out); }
     return (0);
}
     
// returns 1 if successful, 0 if too many nodes (so failure)
// The only routine that sets treeArr
int KUnitigsJoinerThread::joinKUnitigsFromMates (int insertLengthMean, int insertLengthStdev)
{
     
     int lastOffsetToTest = 6000, lastOffsetToTestIfNotMate2, maxOffsetToAllow;
     long int j;
     struct unitigLocStruct unitigLocVal;
     struct abbrevUnitigLocStruct abbrevUnitigLocVal;
     size_t maxNodes;
     int unitig1, unitig2;
     char ori; 
     int offset;
     int overlapLength;
     int ahg, bhg;
     int forcedStop;
     min_heap forward_path_unitigs;

     lastOffsetToTest = insertLengthMean+args.num_stdevs_allowed_arg*insertLengthStdev;
//     if (lastOffsetToTest > max_offset_to_test)
//	  lastOffsetToTest = max_offset_to_test;
     // The following assumes that all the overlaps are of length
     // minOverlapLength
     lastOffsetToTestIfNotMate2 = lastOffsetToTest - (unitigLengths[mateUnitig2]-args.min_overlap_length_arg);
     // Adjust overlaps for mateUnitig1 if ori not 'F' and for
     //    mateUnitig2 if ori not 'R' to what they would be if they had
     //    the desired orientation. Note that mateUnitig1 is only a
     //    source and mateUnitig2 is only a sink
     unitigLocVal.unitig2 = mateUnitig1;
     unitigLocVal.frontEdgeOffset = unitigLengths[mateUnitig1];
     unitigLocVal.ori = mateUnitig1ori;
     abbrevUnitigLocVal.frontEdgeOffset = unitigLocVal.frontEdgeOffset;
     abbrevUnitigLocVal.ori = unitigLocVal.ori;
     abbrevUnitigLocVal.pathNum = 0;

     treeArr.clear();
     treeArr[mateUnitig1].insert(abbrevUnitigLocVal);
     assert(treeArr.find(mateUnitig1) != treeArr.end());
     unitig2 = mateUnitig1; // Initialized to make the compiler happy
#if 0
     fprintf (stderr, "Inserting at 1 in the RB tree at %d: fEO = %d, pN = %u ori = %c\n", mateUnitig1, abbrevUnitigLocVal.frontEdgeOffset, abbrevUnitigLocVal.pathNum, abbrevUnitigLocVal.ori);
#endif     
//     puts ("At 101\n"); fflush (stdout);
     forward_path_unitigs.clear();
     int numUnitigsPushed = 0;
     forward_path_unitigs.push(unitigLocVal);
     ++numUnitigsPushed;
     
//     puts ("At 102\n"); fflush (stdout);
     startingNodes.clear();
     startingNodes.insert(unitigLocVal);
     nodeArray.clear();
     nodeToIndexMap.clear();
//     puts ("At 103\n"); fflush (stdout);
//     fwdEdgeList.clear();
//     revEdgeList.clear();
     nodeArray.push_back(unitigLocVal);
     nodeToIndexMap.insert (std::pair<unitigLocStruct, int> (unitigLocVal, nodeArray.size()-1) );
     maxNodes = 1;
     forcedStop = 0;
//     puts ("At 104\n"); fflush (stdout);
     while (!forward_path_unitigs.empty())
     {
	  if (forward_path_unitigs.size() > maxNodes)
	       maxNodes = forward_path_unitigs.size();
#if DEBUG
	  fprintf (stderr, "Starting new offset: "); fflush (stderr);
	  for(min_heap::iterator it = forward_path_unitigs.begin(); it != forward_path_unitigs.end(); ++it)
	       fprintf(stderr, "%d ", it->frontEdgeOffset); fflush(stderr);
	  fprintf (stderr, "\n"); fflush (stderr);
#endif
          unitigLocVal = forward_path_unitigs.pop();
	  unitig1 = unitigLocVal.unitig2;
#if DEBUG
	  fprintf (stderr, "unitig1 = %d, ", unitig1); fflush (stdout);
#endif
	  ori = unitigLocVal.ori;
	  offset = unitigLocVal.frontEdgeOffset;
	  abbrevUnitigLocVal.frontEdgeOffset = offset;
#if DEBUG
	  fprintf (stderr, "offset = %d, ", offset); fflush (stdout);
#endif
	  abbrevUnitigLocVal.ori = ori;
          
	  for (j = startOverlapByUnitig[unitig1];
	       j < startOverlapByUnitig[unitig1 + 1]; j++)
	  {
#if DEBUG
	       fprintf (stderr, "Starting an overlap:");
#endif
	       unitig2 = overlapData[j].unitig2;
	       if (overlapData[j].ahg >= 0)
		    overlapLength = unitigLengths[unitig1] - overlapData[j].ahg;
	       else
		    overlapLength = unitigLengths[unitig2] + overlapData[j].ahg;
	       if (overlapLength > unitigLengths[unitig1])
		    overlapLength = unitigLengths[unitig1];
	       if (overlapLength > unitigLengths[unitig2])
		    overlapLength = unitigLengths[unitig2];
#if DEBUG
	       fprintf (stderr, "; ovl len = %d", overlapLength);
	       fprintf (stderr, "; unitig2 = %d, ori = %c, ahg = %d, bhg = %d\n", unitig2,
			overlapData[j].ori, overlapData[j].ahg,
			overlapData[j].bhg);
#endif
	       if (ori == 'F')
	       {
		    if (overlapData[j].bhg <= 0) 
			 continue;
		    bhg = overlapData[j].bhg;
		    if (overlapData[j].ori == 'N')
			 abbrevUnitigLocVal.ori = 'F';
		    else
			 abbrevUnitigLocVal.ori = 'R';
		    abbrevUnitigLocVal.frontEdgeOffset = offset + bhg;
	       }
	       else
	       {
		    if (overlapData[j].ahg >= 0)
			 continue;
		    ahg = overlapData[j].ahg;
		    if (overlapData[j].ori == 'N')
			 abbrevUnitigLocVal.ori = 'R';
		    else
			 abbrevUnitigLocVal.ori = 'F';
		    abbrevUnitigLocVal.frontEdgeOffset = offset - ahg;
	       }
	       // Skip if the offset is too large
#if DEBUG
	       fprintf (stderr, "frontEdgeOffset = %d, lastOffsetToTest = %d\n", abbrevUnitigLocVal.frontEdgeOffset, lastOffsetToTest);
#endif
	       if (unitig2 == mateUnitig2)
		    maxOffsetToAllow = lastOffsetToTest;
	       else
		    maxOffsetToAllow = lastOffsetToTestIfNotMate2;
	       if (abbrevUnitigLocVal.frontEdgeOffset > maxOffsetToAllow)
		    continue;
#if DEBUG
	       fprintf (stderr, "cur front = %d\n", abbrevUnitigLocVal.frontEdgeOffset);
#endif
	       // Skip if abbrevUnitigLocVal al unitig y seen for unitig2
               auto unitig2_tree = treeArr.find(unitig2);
               if(unitig2_tree != treeArr.end()) {
		    auto element = find_within(unitig2_tree->second, abbrevUnitigLocVal,
					       default_max_offset_considered_same);
		    if(element != unitig2_tree->second.end())
			 continue;
               }
	       
#if DEBUG
	       fprintf (stderr, "Adding to the tree, unitig2 = %d\n", unitig2);
#endif
	       // Insert this value in the priority queue
	       unitigLocVal.unitig2 = unitig2;
	       unitigLocVal.frontEdgeOffset = abbrevUnitigLocVal.frontEdgeOffset;
	       unitigLocVal.ori = abbrevUnitigLocVal.ori;
               forward_path_unitigs.push(unitigLocVal);
	       ++numUnitigsPushed;
	       abbrevUnitigLocVal.pathNum = 0;
	       // Add offset to list for tree
               treeArr[unitig2].insert(abbrevUnitigLocVal);
	       if ((unitig2 == mateUnitig2) && (unitigLocVal.ori == mateUnitig2ori) && (args.join_aggressive_arg > 0))
		    break;
#if DEBUG
	       if (unitig2 == mateUnitig2)
		    fprintf (stderr, "Adding distance %d for the (rev oriented) mate, unitig %d\n",
			     abbrevUnitigLocVal.frontEdgeOffset, mateUnitig2);
#endif		
	       //   Make sure the root of the tree is updated (if needed)
	  }			// End of going through overlaps for unitig
//	  if (maxNodes > MAX_NODES_ALLOWED)
//	  if (forward_path_unitigs.size() > MAX_NODES_ALLOWED) {
	  if (numUnitigsPushed > args.max_nodes_allowed_arg) {
	       ++numUnjoinableOverMaxNodes;
	       joinCode = 5;
//          if (treeArr[unitig2].size() > MAX_NODES_ALLOWED) {
	       forcedStop = 1;
	       break; }
     }			// Ends !forward_path_unitigs.empty() line
     forward_path_unitigs.clear();

//     if (maxNodes > MAX_NODES_ALLOWED)
     if (forcedStop)
	  return (0);
     else
	  return (1);
}

// The following returns the overlap length if it is greater than the
// existing largest overlap on the end. It returns -1 if not.
int getOvlLenFromOvlIndicesPlus (long int maxOvlIndex,long int j, int maxOvlLen, int whichEnd)
{
     int thisUnitig, otherUnitig, overlapLength;
     otherUnitig = overlapData[maxOvlIndex].unitig2;
     thisUnitig = overlapData[j].unitig2;
     if (whichEnd == FRONT_END)
	  overlapLength = unitigLengths[thisUnitig] - overlapData[j].ahg;
     else
	  overlapLength = unitigLengths[thisUnitig] + overlapData[j].bhg;
     if (overlapLength < maxOvlLen) return (-1);
     if (overlapLength > maxOvlLen) return (overlapLength);
     if (unitigLengths[otherUnitig] > unitigLengths[thisUnitig]) return (-1);
     if (unitigLengths[otherUnitig] < unitigLengths[thisUnitig]) return (overlapLength);
     if (otherUnitig < thisUnitig) return (-1);
     if (otherUnitig > thisUnitig) return (overlapLength);
     if (overlapData[maxOvlIndex].ori == 'N') return (-1);
     return (overlapLength);
}

long int findOtherOverlapIndex (long int ovlIndex1)
{
     int unitig1, unitig2;
     long int itemp;
     unitig1 = overlapData[ovlIndex1].unitig1;
     unitig2 = overlapData[ovlIndex1].unitig2;
     for (itemp=startOverlapByUnitig[unitig2]; itemp<startOverlapByUnitig[unitig2+1]; itemp++) {
	  if (overlapData[itemp].unitig2 != unitig1) continue;
	  if (overlapData[itemp].ori != overlapData[ovlIndex1].ori) continue;
	  if (overlapData[itemp].ori == 'N') {
	       if (overlapData[itemp].ahg + overlapData[ovlIndex1].ahg != 0)
		    continue;
	  }
	  else {
	       if (overlapData[itemp].ahg != overlapData[ovlIndex1].bhg)
		    continue;
	  }
	  // If we get here they agree
	  break;
     }
     return (itemp);
}

void KUnitigsJoinerThread::printIfGood (struct abbrevUnitigLocStruct *ptr)
{
     int val;
     if (ptr->ori == 'R') {
	  val = ptr->frontEdgeOffset;
	  if (mateUnitig1ori != 'F') val -= unitigLengths[mateUnitig1];
	  if (mateUnitig2ori != 'R') val -= unitigLengths[mateUnitig2];
#ifndef NO_OUTPUT
	  fprintf (stderr, "%d\n", val);
#endif
     }
}

template<typename T>
void KUnitigsJoinerThread::printPathNode (const T& ptr) // take a ptr/iterator to a unitigPathPrintStruct
{
     int beginOffset, endOffset;
     if (ptr->ori == 'F') {
	  endOffset = ptr->frontEdgeOffset;
	  beginOffset = endOffset - unitigLengths[ptr->unitig1];
     }
     else {
	  beginOffset = ptr->frontEdgeOffset;
	  endOffset = beginOffset - unitigLengths[ptr->unitig1];
     }
//     fprintf (stderr, "In printPathNode:\n");
//     fprintf (stderr, "uni = %d, offset = %d, ori = %c, beginOffset = %d, endOffset = %d, numOvlsIn = %d, numOvlsOut = %d\n", ptr->unitig1, ptr->frontEdgeOffset, ptr->ori, beginOffset, endOffset, ptr->numOverlapsIn, ptr->numOverlapsOut);
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].unitig1 = ptr->unitig1;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].frontEdgeOffset = ptr->frontEdgeOffset;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].ori = ptr->ori;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].beginOffset = beginOffset;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].endOffset = endOffset;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsIn = ptr->numOverlapsIn;
     augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsOut = ptr->numOverlapsOut;
     ++numUnitigPathPrintRecsOnPath;
     if (ptr->numOverlapsIn == 0)
	  approxNumPaths += ptr->numOverlapsOut;
     else if (ptr->numOverlapsOut > 0)
	  approxNumPaths += (ptr->numOverlapsOut - 1);
}

template<typename T> // A pointer/iterator to a struct abbrevUnitigLocStruct
void KUnitigsJoinerThread::completePathPrint (const T& ptr)
{
     struct abbrevUnitigLocStruct abbrevUnitigLocVal;
     struct unitigLocStruct unitigLocVal;
     struct unitigPathPrintStruct unitigPathPrintVal;
     struct unitigConnectionsForPathStruct unitigConnectionsForPathRec;
     char ori;
#ifdef NEW_STUFF
     char tempOri1, tempOri2;
#endif
     int isSpecialCase, finalOffset, minConnectingOffset;
     long int i, index;
     int unitig1, unitig2, offset, overlapLength;
     int minConnectingUnitig=0;
     //     long int minConnectingOverlapIndex; // Unused!
     char minConnectingOri=' ';
#ifdef KILLED111115
     double numStdevsFromMean;
#endif
     // In the following we assume we move from left to right when moving from
     // beginUnitig to endUnitig.
     ++curPathNum;
     if (curPathNum > maxPathNumUsed)
	  maxPathNumUsed = curPathNum;
#if 0
     fprintf (stderr, "curPathNum = %d, nodePathNum = %d; ", curPathNum, ptr->pathNum);
     fprintf (stderr, "frontEdgeOffset = %d, ori = %c\n", ptr->frontEdgeOffset, ptr->ori);
#endif
     ori = ptr->ori;
     if (ori != endUnitigOri) {
	  tooManyPossibleInsertLengths = false;
	  return; }
     isSpecialCase = 1;
     finalOffset = ptr->frontEdgeOffset;
     minConnectingOffset = finalOffset + 1000000;
#ifdef KILLED111115
     numStdevsFromMean = (finalOffset - insertLengthMeanBetweenKUnisForInsertGlobal)/insertLengthStdevGlobal;
     fprintf (stderr, "%d %f\n", finalOffset, numStdevsFromMean);
#endif
     for (i=startOverlapByUnitig[endUnitig]; i<startOverlapByUnitig[endUnitig+1]; i++) {
	  index = unitig2OverlapIndex[i];
	  unitig1 = overlapData[index].unitig1;
	  if (overlapData[index].ahg >= 0)
	       overlapLength = unitigLengths[unitig1] - overlapData[index].ahg;
	  else
	       overlapLength = unitigLengths[endUnitig] + overlapData[index].ahg;
	  if (overlapLength > unitigLengths[unitig1])
	       overlapLength = unitigLengths[unitig1];
	  if (overlapLength > unitigLengths[endUnitig])
	       overlapLength = unitigLengths[endUnitig];
	  if (overlapData[index].ori == 'N')
	       abbrevUnitigLocVal.ori = ori;
	  else {
	       if (ori == 'F') abbrevUnitigLocVal.ori = 'R';
	       else abbrevUnitigLocVal.ori = 'F';
	  }
	  if (abbrevUnitigLocVal.ori == 'F')
	       abbrevUnitigLocVal.frontEdgeOffset = finalOffset - overlapData[index].bhg;
	  else
	       abbrevUnitigLocVal.frontEdgeOffset = finalOffset + overlapData[index].ahg;
          auto unitig_tree = treeArr.find(unitig1);
          if(unitig_tree == treeArr.end())
	       continue;
          auto element = find_within(unitig_tree->second, abbrevUnitigLocVal, 
                                     default_max_offset_considered_same);
          if(element == unitig_tree->second.end())
	       continue;
#if 0
	  fprintf (stderr, "frontEdgeOffset = %d\n", element->frontEdgeOffset);
#endif
	  if (element->frontEdgeOffset < finalOffset) {
	       isSpecialCase = 0;
	       break;
	  }
	  else if (element->frontEdgeOffset < minConnectingOffset) {
	       minConnectingOffset = element->frontEdgeOffset;
	       minConnectingUnitig = unitig1;
	       minConnectingOri = abbrevUnitigLocVal.ori;
               //	       minConnectingOverlapIndex = index;
	       if (args.join_aggressive_arg > 0) {
		    isSpecialCase = 0;
		    break; }
	  }
     }
     if (isSpecialCase) {
	  fprintf (stdout, "We shouldn't get here\n");
	  unitigLocVal.unitig2 = minConnectingUnitig;
	  unitigLocVal.frontEdgeOffset = minConnectingOffset;
	  unitigLocVal.ori = minConnectingOri;
	  // ..and for the path
	  unitigPathPrintVal.unitig1 = endUnitig;
	  unitigPathPrintVal.numOverlapsIn = 1;
	  unitigPathPrintVal.numOverlapsOut = 0;
	  unitigPathPrintVal.ori = 'R'; // Forced; may be adjusted later
	  unitigPathPrintVal.frontEdgeOffset = finalOffset;
          
          treeArr2.insert(unitigPathPrintVal);
	  unitigPathPrintVal.unitig1 = minConnectingUnitig;
	  unitigPathPrintVal.frontEdgeOffset = minConnectingOffset;	
	  unitigPathPrintVal.numOverlapsIn = 0;
	  unitigPathPrintVal.numOverlapsOut = 1;
	  unitigPathPrintVal.ori = minConnectingOri;
          treeArr2.insert(unitigPathPrintVal);
     }
     else {
	  unitigLocVal.unitig2 = endUnitig;
	  unitigLocVal.frontEdgeOffset = finalOffset;
	  unitigLocVal.ori = endUnitigOri;
	  unitigPathPrintVal.unitig1 = endUnitig;
	  unitigPathPrintVal.frontEdgeOffset = finalOffset;
	  unitigPathPrintVal.numOverlapsIn = 0;
	  unitigPathPrintVal.numOverlapsOut = 0;
	  unitigPathPrintVal.ori = endUnitigOri;
#if 0
	  fprintf (stderr, "Inserting unitig1 = %d, offset = %d, ori = %c\n", unitigPathPrintVal.unitig1, unitigPathPrintVal.frontEdgeOffset, unitigPathPrintVal.ori);
#endif
          treeArr2.insert(unitigPathPrintVal);
     }
#if 0
     fprintf (stderr, "isSpecialCase = %d, unitigLocVal = %d, %d, %c\n", isSpecialCase, endUnitig, finalOffset, unitigLocVal.ori);
#endif
     backward_path_unitigs.clear();
     backward_path_unitigs.push(unitigLocVal);
     while (!backward_path_unitigs.empty()) {
          unitigLocVal = backward_path_unitigs.pop();
	  unitig2 = unitigLocVal.unitig2;
	  offset = unitigLocVal.frontEdgeOffset;
	  ori = unitigLocVal.ori;
	  unitigPathPrintVal.unitig1 = unitig2;
	  unitigPathPrintVal.frontEdgeOffset = offset;
	  unitigPathPrintVal.ori = ori;
          auto front_unitig = find_within(treeArr2, unitigPathPrintVal,
                                          default_max_offset_considered_same);
#if 0
	  fprintf (stderr, "unitig2 = %d, offset = %d, ori = %c, front_unitigs's unitig = %d\n", unitig2, offset, ori, front_unitig->unitig1);
#endif
	  for (i=startOverlapByUnitig[unitig2]; i<startOverlapByUnitig[unitig2+1]; i++) {
	       index = unitig2OverlapIndex[i];
	       unitig1 = overlapData[index].unitig1;
#if 0
	       fprintf (stderr, "unitig1 = %d\n", unitig1);
#endif
	       if (overlapData[index].ahg >= 0)
		    overlapLength = unitigLengths[unitig1] - overlapData[index].ahg;
	       else 
		    overlapLength = unitigLengths[unitig2] + overlapData[index].ahg;
	       if (overlapLength > unitigLengths[unitig1])
		    overlapLength = unitigLengths[unitig1];
	       if (overlapLength > unitigLengths[unitig2])
		    overlapLength = unitigLengths[unitig2];
	       if (overlapData[index].ori == 'N')
		    abbrevUnitigLocVal.ori = ori;
	       else {
		    if (ori == 'F')
			 abbrevUnitigLocVal.ori = 'R';
		    else
			 abbrevUnitigLocVal.ori = 'F';
	       }
	       if (abbrevUnitigLocVal.ori == 'F')
		    abbrevUnitigLocVal.frontEdgeOffset = offset - overlapData[index].bhg;
	       else
		    abbrevUnitigLocVal.frontEdgeOffset = offset + overlapData[index].ahg;
               auto unitig1_tree = treeArr.find(unitig1);
               if(unitig1_tree == treeArr.end())
		    continue;
               auto element = find_within(unitig1_tree->second, abbrevUnitigLocVal, 
                                          default_max_offset_considered_same);
               if(element == unitig1_tree->second.end())
		    continue;
	       if (element->frontEdgeOffset >= offset) continue;
	       // It hasn't been seen in the retrace, so put on the queue
	       if (element->pathNum < curPathNum) {
		    element->pathNum = curPathNum;
		    unitigLocVal.unitig2 = unitig1;
		    unitigLocVal.frontEdgeOffset = element->frontEdgeOffset;
		    unitigLocVal.ori = element->ori;
#if DEBUG
		    fprintf (stderr, "Adding node: unitig2 = %d, offset = %d, ori = %c, curPathNum = %d\n", unitigLocVal.unitig2, unitigLocVal.frontEdgeOffset, unitigLocVal.ori, (int) curPathNum);
#endif
                    backward_path_unitigs.push(unitigLocVal);
		    unitigPathPrintVal.unitig1 = unitig1;
		    unitigPathPrintVal.frontEdgeOffset = unitigLocVal.frontEdgeOffset;
		    unitigPathPrintVal.ori = unitigLocVal.ori;
		    unitigPathPrintVal.numOverlapsIn = 0;
		    unitigPathPrintVal.numOverlapsOut = 0;
                    treeArr2.insert(unitigPathPrintVal);
	       }
	       unitigPathPrintVal.unitig1 = unitig1;
	       unitigPathPrintVal.frontEdgeOffset = element->frontEdgeOffset;
	       unitigPathPrintVal.ori = element->ori;
               auto rear_unitig = find_within(treeArr2, unitigPathPrintVal, 
                                              default_max_offset_considered_same);
	       int frontEdgeOffset1 = rear_unitig->frontEdgeOffset;
	       int frontEdgeOffset2 = front_unitig->frontEdgeOffset;
//	       fprintf (stderr, "At 31411: front_unitig = %d, rear_unitig = %d, frontEdgeOffset2 = %d, frontEdgeOffset1 = %d, firstUnitigLen = %d, minOvlLen = %d\n", (int) front_unitig->unitig1, (int) rear_unitig->unitig1, (int) frontEdgeOffset2, (int) frontEdgeOffset1, (int) unitigLengths[front_unitig->unitig1], (int) args.min_overlap_length_arg);
//	       fprintf (stderr, "frontEdgeOffset1 = %d, frontEdgeOffset2 = %d\n", frontEdgeOffset1, frontEdgeOffset2); fflush (stderr);
	       int diff = (frontEdgeOffset2 - frontEdgeOffset1) - (unitigLengths[front_unitig->unitig1] - args.min_overlap_length_arg);
	       if (diff < 0)
		    diff = - diff;
//	       if (frontEdgeOffset2 - frontEdgeOffset1 != unitigLengths[front_unitig->unitig1] - args.min_overlap_length_arg)
	       if (diff > default_max_offset_considered_same)
		    continue;
//	       fprintf (stderr, "At 31415: front_unitig = %d, rear_unitig = %d\n", (int) front_unitig->unitig1, (int) rear_unitig->unitig1);
	       ++(rear_unitig->numOverlapsOut);
	       if (rear_unitig->numOverlapsOut > 1) {
		    if (rear_unitig->frontEdgeOffset < splitJoinWindowMin)
			 splitJoinWindowMin = rear_unitig->frontEdgeOffset; }
	       ++(front_unitig->numOverlapsIn);
	       if (front_unitig->numOverlapsIn > 1) {
		    if (front_unitig->frontEdgeOffset > splitJoinWindowMax)
			 splitJoinWindowMax = front_unitig->frontEdgeOffset; }
#ifdef NEW_STUFF
	       tempOri1 = front_unitig->ori;
	       if (front_unitig->unitig1 == beginUnitig) tempOri1 = beginUnitigOri;
	       if (front_unitig->unitig1 == endUnitig) tempOri1 = endUnitigOri;
	       tempOri2 = rear_unitig->ori;
	       if (rear_unitig->unitig1 == beginUnitig) tempOri2 = beginUnitigOri;
	       if (rear_unitig->unitig1 == endUnitig) tempOri2 = endUnitigOri;
	       unitigConnectionsForPathRec.unitig1 = rear_unitig->unitig1;
	       unitigConnectionsForPathRec.unitig2 = front_unitig->unitig1;
	       unitigConnectionsForPathRec.frontEdgeOffset1 = rear_unitig->frontEdgeOffset;
	       unitigConnectionsForPathRec.frontEdgeOffset2 = front_unitig->frontEdgeOffset;
	       unitigConnectionsForPathRec.ori1 = tempOri2;
	       unitigConnectionsForPathRec.ori2 = tempOri1;
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].unitig1 = unitigConnectionsForPathRec.unitig1;
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].unitig2 = unitigConnectionsForPathRec.unitig2;
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].frontEdgeOffset1 = unitigConnectionsForPathRec.frontEdgeOffset1;
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].frontEdgeOffset2 = unitigConnectionsForPathRec.frontEdgeOffset2;
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].ori1 = unitigConnectionsForPathRec.ori1;
	       unitigConnectionsForPathData[numUnitigConnectionsForPathData].ori2 = unitigConnectionsForPathRec.ori2;
	       ++numUnitigConnectionsForPathData;
#endif
	       if(args.join_aggressive_arg > 0)
		    break;
	       
	  }
     }
// #ifdef KILLED111115
#ifdef KILL120102
     fprintf (stderr, "Entering loop in completePathPrint\n");
#endif
     for (i=0; i<(int) unitigConnectionsForPathData.size(); i++) {
	  unitigLocStruct uLS;
//	  int index1, index2;
	  int firstIndex, secondIndex;
	  unitigLocMap_iterator it;
	  std::list<int>::iterator l_it;
          //	  struct nodePair nodeValue; // Unused!
//	  fprintf (stderr, "(uni1,front1,ori1) = (%d,%d,%c); (uni2,front2,ori2) = (%d,%d,%c)\n", unitigConnectionsForPathData[i].unitig1, unitigConnectionsForPathData[i].frontEdgeOffset1, unitigConnectionsForPathData[i].ori1, unitigConnectionsForPathData[i].unitig2, unitigConnectionsForPathData[i].frontEdgeOffset2, unitigConnectionsForPathData[i].ori2);
	  uLS.unitig2 = unitigConnectionsForPathData[i].unitig1;
	  uLS.frontEdgeOffset = unitigConnectionsForPathData[i].frontEdgeOffset1;
	  uLS.ori = unitigConnectionsForPathData[i].ori1;
	  it = nodeToIndexMap.find (uLS);
	  if (it == nodeToIndexMap.end()) {
	       nodeArray.push_back(uLS);
	       firstIndex = nodeArray.size()-1;
	       nodeToIndexMap.insert(std::pair<unitigLocStruct, int> (uLS, firstIndex ) ); }
	  else {
	       firstIndex = it->second; }
	  uLS.unitig2 = unitigConnectionsForPathData[i].unitig2;
	  uLS.frontEdgeOffset = unitigConnectionsForPathData[i].frontEdgeOffset2;
	  uLS.ori = unitigConnectionsForPathData[i].ori2;
	  it = nodeToIndexMap.find (uLS);
	  if (it == nodeToIndexMap.end()) {
	       nodeArray.push_back(uLS);
	       secondIndex = nodeArray.size()-1;
	       nodeToIndexMap.insert(std::pair<unitigLocStruct, int> (uLS, secondIndex ) ); }
	  else {
	       secondIndex = it->second; }
	  // nodeValue.node1 = firstIndex;
	  // nodeValue.node2 = secondIndex;
	  edgeList.insert(std::pair<int, int> (firstIndex, secondIndex));
	  
#ifdef KILL120102	  
	  fprintf (stderr, "%s%lld Node (%d, %d, %c) -> (%d, %d, %c)\n", rdPrefixHold, readNumHold, unitigConnectionsForPathData[i].unitig1, unitigConnectionsForPathData[i].frontEdgeOffset1, unitigConnectionsForPathData[i].ori1, unitigConnectionsForPathData[i].unitig2, unitigConnectionsForPathData[i].frontEdgeOffset2, unitigConnectionsForPathData[i].ori2);
#endif
     }
// #endif

     numUnitigPathPrintRecsOnPath = 0;
     if (treeSize <= maxDiffInsertSizesForPrinting)
	  for(auto it = treeArr2.begin(); it != treeArr2.end(); ++it)
	       printPathNode(it);
     else
	  tooManyPossibleInsertLengths = true;
#if DEBUG
     fprintf (stderr, "treeSize = %d, maxDiffInsertSizesForPrinting = %d\n", (int) treeSize, (int) maxDiffInsertSizesForPrinting);
#endif
#ifdef KILLED111115
     for (i=0; i<numUnitigPathPrintRecsOnPath; i++)
	  fprintf (stderr, "uni = %d, offset = %d, ori = %c, beginOffset = %d, endOffset = %d, numOvlsIn = %d, numOvlsOut = %d\n", augmentedUnitigPathPrintData[i].unitig1, augmentedUnitigPathPrintData[i].frontEdgeOffset, augmentedUnitigPathPrintData[i].ori, augmentedUnitigPathPrintData[i].beginOffset, augmentedUnitigPathPrintData[i].endOffset, augmentedUnitigPathPrintData[i].numOverlapsIn, augmentedUnitigPathPrintData[i].numOverlapsOut);
#endif
     if (approxNumPaths == 1) {
	  joinCode = 2;
	  generateSuperReadPlacementLinesForJoinedMates(); }
#if 0
//     printf ("final offset = %d, arraySize = %d\n", finalOffset, dataArr2.arraySize);
     fprintf (stderr, "final offset = %d\n", finalOffset);
#endif
     // TODO: delete
     // treeArr2[0].root = TREE_NIL;
     // dataArr2.arraySize = 0;
     treeArr2.clear();
}

void KUnitigsJoinerThread::generateSuperReadPlacementLinesForJoinedMates ()
{
     int isReversed, superReadLength;
     // the following uses augmentedUnitigPathPrintData
     superReadLength = getSuperReadLength ();
     isReversed = setSuperReadNameFromAugmentedPath ();
     sprintf (outputString,"%s%lld %s ", rdPrefixHold, readNumHold-1, (char*)superReadName);
     if (! isReversed)
	  sprintf_append (outputString,"%d F ", lengthAdjustment1);
     else
	  sprintf_append (outputString,"%d R ", superReadLength - lengthAdjustment1);
     sprintf_append (outputString, "%s\n", joinCodeNames[joinCode]);
     sprintf_append (outputString,"%s%lld %s ", rdPrefixHold, readNumHold, (char*)superReadName);
     if (! isReversed)
	  sprintf_append (outputString,"%d R ", superReadLength - lengthAdjustment2);
     else
	  sprintf_append (outputString,"%d F ", lengthAdjustment2);
     sprintf_append (outputString, "%s\n", joinCodeNames[joinCode]);
}

int KUnitigsJoinerThread::setSuperReadNameFromAugmentedPath (void)
{
     int isReversed=0, i;
     for (i=0; i<numUnitigPathPrintRecsOnPath/2; i++) {
	  if (augmentedUnitigPathPrintData[i].unitig1 != augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath-i-1].unitig1) {
	       if (augmentedUnitigPathPrintData[i].unitig1 < augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath-i-1].unitig1)
		    isReversed = 0;
	       else
		    isReversed = 1;
	       break;
	  }
	  if (augmentedUnitigPathPrintData[i].ori == augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath-i-1].ori) {
	       if (augmentedUnitigPathPrintData[i].ori == 'F')
		    isReversed = 0;
	       else
		    isReversed = 1;
	       break;
	  }
     }
     superReadName.clear();
     if (isReversed == 0) {
	  sprintf_append (superReadName, "%d%c", augmentedUnitigPathPrintData[0].unitig1, augmentedUnitigPathPrintData[0].ori);
	  for (i=1; i<numUnitigPathPrintRecsOnPath; i++) {
	       sprintf_append (superReadName, "_%d%c", augmentedUnitigPathPrintData[i].unitig1, augmentedUnitigPathPrintData[i].ori);
	  }
     }
     else {
	  sprintf_append (superReadName, "%d%c", augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath-1].unitig1, (augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath-1].ori == 'F') ? 'R' : 'F');
	  for (i=numUnitigPathPrintRecsOnPath-2; i>=0; i--) {
	       sprintf_append (superReadName, "_%d%c", augmentedUnitigPathPrintData[i].unitig1, (augmentedUnitigPathPrintData[i].ori == 'F') ? 'R' : 'F');
	  }
     }
     return (isReversed);
}

int KUnitigsJoinerThread::getSuperReadLength(void)
{
     int totLen, i;
     totLen = unitigLengths[augmentedUnitigPathPrintData[0].unitig1];
     for (i=1; i<numUnitigPathPrintRecsOnPath; i++)
	  totLen += (unitigLengths[augmentedUnitigPathPrintData[i].unitig1] - args.min_overlap_length_arg);
     
     return (totLen);
}

template<typename T> // A ptr/iterator to a abbrevUnitigLocStruct
void KUnitigsJoinerThread::funcToGetTreeSize (const T& ptr)
{
     struct unitigLocStruct localUnitigLoc;
     if (ptr->ori == endUnitigOri) {
	  ++treeSize;
	  localUnitigLoc.unitig2 = endUnitig;
	  localUnitigLoc.frontEdgeOffset = ptr->frontEdgeOffset;
	  localUnitigLoc.ori = endUnitigOri;
	  endingNodes.insert (localUnitigLoc); 
	  nodeArray.push_back(localUnitigLoc);
	  nodeToIndexMap.insert (std::pair<unitigLocStruct, int> (localUnitigLoc, nodeArray.size()-1) ); }
}

bool unitigLocStructCompare (struct unitigLocStruct uLS1,
                             struct unitigLocStruct uLS2)
{
     if (uLS1.frontEdgeOffset != uLS2.frontEdgeOffset)
	  return (uLS1.frontEdgeOffset < uLS2.frontEdgeOffset);
     if (uLS1.unitig2 != uLS2.unitig2)
	  return (uLS1.unitig2 < uLS2.unitig2);
     return (uLS1.ori < uLS2.ori);
}

bool unitigLocStructCompareReversed (struct unitigLocStruct uLS1,
                                     struct unitigLocStruct uLS2)
{
     if (uLS1.frontEdgeOffset < uLS2.frontEdgeOffset)
	  return (-1);
     if (uLS1.frontEdgeOffset > uLS2.frontEdgeOffset)
	  return (1);
     if (uLS1.unitig2 < uLS2.unitig2)
	  return (-1);
     if (uLS1.unitig2 > uLS2.unitig2)
	  return (1);
     if (uLS1.ori < uLS2.ori)
	  return (-1);
     if (uLS1.ori > uLS2.ori)
	  return (1);
     return (0);
}

FILE *Fopen (const char *fn, const char *mode)
{
     FILE *result;
     result = fopen (fn, mode);
     if (result == NULL)
     {
	  fprintf (stderr, "Couldn't open file '%s' for ", fn);
	  switch (mode[0])
	  {
	  case 'r':
	       fprintf (stderr, "reading");
	       break;
	  case 'w':
	       fprintf (stderr, "writing");
	       break;
	  case 'a':
	       fprintf (stderr, "appending");
	       break;
	  default:
	       fprintf (stderr, "unknown operation code '%c'", mode[0]);
	       break;
	  }
	  fprintf (stderr, ". Bye!\n");
	  exit (-1);
     }

     return (result);
}

FILE *Popen (const char *fn, const char *mode)
{
     FILE *result;
     result = popen (fn, mode);
     if (result == NULL)
     {
	  fprintf (stderr, "Couldn't open file '%s' for ", fn);
	  switch (mode[0])
	  {
	  case 'r':
	       fprintf (stderr, "reading");
	       break;
	  case 'w':
	       fprintf (stderr, "writing");
	       break;
	  case 'a':
	       fprintf (stderr, "appending");
	       break;
	  default:
	       fprintf (stderr, "unknown operation code '%c'", mode[0]);
	       break;
	  }
	  fprintf (stderr, ". Bye!\n");
	  exit (-1);
     }

     return (result);
}

int getInt (const char *fname)
{
     FILE *infile;
     int tval;

     infile = Fopen (fname, "r");
     if (! fscanf (infile, "%d\n", &tval)) {
	  fprintf (stderr, "Couldn't read file %s. Bye!\n", fname);
	  exit (1);
     }
     fclose (infile);
     return (tval);
}

void KUnitigsJoinerThread::findSingleReadSuperReads(char *readName, jflib::omstream& m_out)
{
     long long tempInt;
     char *cptr=readName+2;
     int countOfMatchingKUnitigs, offsetOfReadInSuperRead;
     int minReadOffset, maxReadOffset, minReadOffsetSeen, maxReadOffsetSeen;
     int i, j, recNumToUse=0;
     int isReversed=0;
     struct kuniToReadMatchStruct *kUTRMSptr;

     tempInt = atoll(cptr);
#ifdef KILLED111115
     printf ("findSingleReadSuperReads\n");
#endif
     if (tempInt % 2 == 0) {
          countOfMatchingKUnitigs = evenReadMatchStructs.size();
	  kUTRMSptr = &(evenReadMatchStructs[0]);
     }
     else {
          countOfMatchingKUnitigs = oddReadMatchStructs.size();
	  kUTRMSptr = &(oddReadMatchStructs[0]);
     }
     
//     printf ("countOfMatchingKUnitigs = %d\n", countOfMatchingKUnitigs);
     i = 0;
     minReadOffsetSeen = kUTRMSptr[i].readMatchBegin;
     maxReadOffsetSeen = kUTRMSptr[i].readMatchEnd;
     matchStructIsUsed[i] = 1;
     for (i=1; i<countOfMatchingKUnitigs; i++) {
	  matchStructIsUsed[i] = 0;
	  if (kUTRMSptr[i].readMatchEnd <= maxReadOffsetSeen)
	       continue;
	  if (kUTRMSptr[i].readMatchBegin < maxReadOffsetSeen-args.min_overlap_length_arg)
	       continue; // Otherwise the k-unitigs overlap too much
	  if (kUTRMSptr[i].readMatchBegin > maxReadOffsetSeen)
	       return; // Part of the middle of the read is uncovered by k-unis
	  maxReadOffsetSeen = kUTRMSptr[i].readMatchEnd;
	  matchStructIsUsed[i] = 1;
     }
     if (minReadOffsetSeen + (kUTRMSptr[0].readLength - maxReadOffsetSeen) > maxTotAllowableMissingOnEnds)
	  return;
     
     i=-1; j=countOfMatchingKUnitigs;
     isReversed = 0;
     while (1) {
	  ++i; --j;
	  while (!matchStructIsUsed[i])
	       ++i;
	  while (!matchStructIsUsed[j])
	       --j;
	  if (kUTRMSptr[i].kUnitigNumber != kUTRMSptr[j].kUnitigNumber) {
	       if (kUTRMSptr[i].kUnitigNumber < kUTRMSptr[j].kUnitigNumber)
		    isReversed = 0;
	       else
		    isReversed = 1;
	       break;
	  }
	  if (kUTRMSptr[i].ori == kUTRMSptr[j].ori) {
	       if (kUTRMSptr[i].ori == 'F')
		    isReversed = 0;
	       else
		    isReversed = 1;
	       break;
	  }
	  if (j<=i)
	       break;
     }
     superReadName.clear();
     if (isReversed == 0) {
	  for (i=0; 1; i++)
	       if (matchStructIsUsed[i])
		    break;
	  recNumToUse = i;
	  sprintf_append (superReadName, "%d%c", kUTRMSptr[i].kUnitigNumber, kUTRMSptr[i].ori);
	  maxReadOffset = kUTRMSptr[i].readMatchEnd;
	  for (++i; i<countOfMatchingKUnitigs; i++) {
	       if (! matchStructIsUsed[i])
		    continue;
	       // The next is the overlap amount between k-unitigs, which we now require to be minOverlapLength
	       if (maxReadOffset-kUTRMSptr[i].readMatchBegin != args.min_overlap_length_arg)
		    return;
	       sprintf_append(superReadName, "_%d%c", kUTRMSptr[i].kUnitigNumber, kUTRMSptr[i].ori);
	       maxReadOffset = kUTRMSptr[i].readMatchEnd;
	  }
	  // Must do the output here
	  if (kUTRMSptr[recNumToUse].ori == 'F')
	       offsetOfReadInSuperRead = kUTRMSptr[recNumToUse].kUnitigMatchBegin - kUTRMSptr[recNumToUse].readMatchBegin;
	  else
	       offsetOfReadInSuperRead = unitigLengths[kUTRMSptr[recNumToUse].kUnitigNumber] - kUTRMSptr[recNumToUse].kUnitigMatchEnd - kUTRMSptr[recNumToUse].readMatchBegin;
	  m_out << readName << " " << superReadName << " " << offsetOfReadInSuperRead << " F " << joinCodeNames[joinCode] << "\n";
	  //	  fprintf (outputFile, "%s %s %d %c\n", readName, (char*)superReadName, offsetOfReadInSuperRead, 'F');	  
     }
     else { // The k-unitigs are reversed from those reported
	  for (i=countOfMatchingKUnitigs-1; 1; i--)
	       if (matchStructIsUsed[i])
		    break;
	  recNumToUse = i;
	  sprintf_append(superReadName, "%d%c", kUTRMSptr[i].kUnitigNumber, (kUTRMSptr[i].ori == 'F') ? 'R' : 'F');
	  minReadOffset = kUTRMSptr[i].readMatchBegin;
	  for (--i; i>=0; i--) {
	       if (! matchStructIsUsed[i])
		    continue;
	       // The next is the overlap amount between k-unitigs, which we now require to be minOverlapLength
	       if (kUTRMSptr[i].readMatchEnd-minReadOffset != args.min_overlap_length_arg)
		    return;
	       sprintf_append(superReadName, "_%d%c", kUTRMSptr[i].kUnitigNumber, (kUTRMSptr[i].ori == 'F') ? 'R' : 'F');
	       minReadOffset = kUTRMSptr[i].readMatchBegin;
	  }
	  // Must do the output here
	  if (kUTRMSptr[recNumToUse].ori == 'F')
	       offsetOfReadInSuperRead = (unitigLengths[kUTRMSptr[recNumToUse].kUnitigNumber]-kUTRMSptr[recNumToUse].kUnitigMatchBegin) + kUTRMSptr[recNumToUse].readMatchBegin;
	  else
	       offsetOfReadInSuperRead = kUTRMSptr[recNumToUse].kUnitigMatchEnd + kUTRMSptr[recNumToUse].readMatchBegin;
	  // The k-unitigs are reversed from those reported
	  m_out << readName << " " << superReadName << " " << offsetOfReadInSuperRead << " R " << joinCodeNames[joinCode] << "\n";
	  //	  fprintf (outputFile, "%s %s %d %c\n", readName, (char*)superReadName, offsetOfReadInSuperRead, 'R');
     }
//     printf ("At 50\n");
}

void KUnitigsJoinerThread::getSuperReadsForInsert (jflib::omstream& m_out)
{
     charb readNameSpace;
     int insertLengthMean;
     int successCode;
     struct abbrevUnitigLocStruct abbULS1, abbULS1Hold;
     struct unitigLocStruct tempULS;
     int numPossibleLengths=0;
     int startValue;
     int pathNum=0;
     std::set<unitigLocStruct>::iterator it1;
     unitigLocMap_iterator it2;
     ExpandingBuffer<int> pathNumArray;
     std::stack<int> nodeIntArray;
     long int localNodeNumber = 0, overlapMatchIndexHold = 0;
     int localUnitigNumber = 0, localNodeNumberHold = 0;
     int lastGoodNodeNumber;
     int localFrontEdgeOffset = 0, localSuperReadLength = 0;
     int doMinimalWorkHere;
     int distFromEndOfSuperRead = 0;
     bool last_element_is_nil = false;
     bool isJoinable = false, wasJoined = false, wasDeclaredUnresolvedJoinable = false;
     unitig_to_ori_offsets::iterator end_tree;
//     int minPathNumInForwardDirection = 0;
     
//     fprintf (stderr, "Entering getSuperReadsForInsert\n");
     
     // Make sure it is initialized
     abbULS1.frontEdgeOffset = 0;
     abbULS1.ori = 'F';
     maxPathNumUsed = 0;
     // Make the compiler happy
     abbULS1Hold.frontEdgeOffset = 0;
     abbULS1Hold.ori = 'F';
     
     
     // Output the stuff for the old pair
#ifdef KILL120103
     stderrOutputString.clear();
#endif
#ifdef KILLED111115
     fprintf (stderr, "%s%lld %ld %ld\n", rdPrefixHold, readNumHold, evenReadMatchStructs.size(), oddReadMatchStructs.size());
#endif
     tooManyPossibleInsertLengths = false;
     sprintf (readNameSpace, "%s%lld", rdPrefixHold, readNumHold);
//     printf ("numEvenReadMatchStructs = %d, numOddReadMatchStructs = %d\n", (int) evenReadMatchStructs.size(), (int) oddReadMatchStructs.size());
     if (evenReadMatchStructs.empty() || oddReadMatchStructs.empty()) {
	  findSingleReadSuperReads(readNameSpace, m_out);
	  return; }
     if (args.max_nodes_allowed_arg == 0) {
	  sprintf (readNameSpace, "%s%lld", rdPrefixHold, readNumHold-1);
	  findSingleReadSuperReads (readNameSpace, m_out);
	  sprintf (readNameSpace, "%s%lld", rdPrefixHold, readNumHold);
	  findSingleReadSuperReads (readNameSpace, m_out);
	  return; }
//     puts ("Got to 1\n"); fflush (stdout);
     // If we get here both the even read and the odd read have
     // matches to k-unitigs
     // The next takes care of the case where both the source and
     // destination k-unitig are the same (We don't join in this case)
     mateUnitig1 = evenReadMatchStructs[0].kUnitigNumber;
     mateUnitig2 = oddReadMatchStructs[0].kUnitigNumber;
     if (mateUnitig1 == mateUnitig2) {
	  joinCode = 1;
	  sprintf (readNameSpace, "%s%lld", rdPrefixHold, readNumHold-1);
//	  puts ("Entering findSingleReadSuperReads from 1\n"); fflush (stdout);
	  findSingleReadSuperReads (readNameSpace, m_out);
	  sprintf (readNameSpace, "%s%lld", rdPrefixHold, readNumHold);
//	  puts ("Entering findSingleReadSuperReads from 2\n"); fflush (stdout);
	  findSingleReadSuperReads (readNameSpace, m_out);
	  ++numPairsInOneUnitig;
	  return;
     }
//     puts ("Got to 2\n"); fflush (stdout);
     // evenReadMatchStructs keeps track of the k-unitig matches with the even-numbered read of a mate pair
     // oddReadMatchStructs keeps track of the k-unitig matches with the odd-numbered read of the pair
     // We allow up to maxTotAllowableMissingOnEnds bases (total) to be uncovered at the ends of the insert
     if ((evenReadMatchStructs[0].readMatchBegin + oddReadMatchStructs[0].readMatchBegin <= maxTotAllowableMissingOnEnds)) {
	  mateUnitig1ori = evenReadMatchStructs[0].ori;
	  // Since we assume that the second read of the mate pair is reversed, we reverse the orientation
	  // of the k-unitig matching the beginning of the second read when we determine the "goal" k-unitig
	  // to match
	  if (oddReadMatchStructs[0].ori == 'F')
	       mateUnitig2ori = 'R';
	  else
	       mateUnitig2ori = 'F';
	  // We have given distances between the ends of the reads; since we are connecting k-unitigs and are
	  // using distances between k-unitigs, we must adjust the allowed distances between the ending k-unitigs
	  // to account for the amount they extend beyond the ends of the reads
	  if (mateUnitig1ori == 'F')
	       lengthAdjustment1 = evenReadMatchStructs[0].ahg;
	  else
	       lengthAdjustment1 = - evenReadMatchStructs[0].bhg;
	  if (mateUnitig2ori == 'R')
	       lengthAdjustment2 = oddReadMatchStructs[0].ahg;
	  else
	       lengthAdjustment2 = - oddReadMatchStructs[0].bhg;
	  insertLengthMean = mean[(int)rdPrefixHold[0]][(int)rdPrefixHold[1]] + (lengthAdjustment1 + lengthAdjustment2);
	  
#ifdef KILLED111115
	  fprintf (stderr, "joinKUnitigsFromMates for pair %s%lld %s%lld using mean %d\n", rdPrefixHold, readNumHold-1, rdPrefixHold, readNumHold, insertLengthMean);
#endif
	  // The following to check what we're doing
	  insertLengthMeanBetweenKUnisForInsertGlobal = insertLengthMean;
	  insertLengthStdevGlobal = stdev[(int)rdPrefixHold[0]][(int)rdPrefixHold[1]];
//	  puts ("Entering joinKUnitigsFromMates from 3\n"); fflush (stdout);
	  // joinKUnitigsFromMates starts from the first k-unitig and finds all k-unitigs (with offsets) that can obtained
	  // by using overlaps between k-unitigs from a k-unitig and offset that have already been added to the list
	  // and only considering overlaps which continue away from the starting k-unitig.
	  // For each k-unitig it creates a tree of all the offsets where it has been observed (usually 1).
	  // It returns 0 if the mate k-unitig is never encountered with the appropriate offset and within
	  // the appropriate distance.
	  successCode = joinKUnitigsFromMates (insertLengthMean, stdev[(int)rdPrefixHold[0]][(int)rdPrefixHold[1]]);
#if DEBUG
	  fprintf (stderr, "Leaving joinKUnitigsFromMates at A with successCode = %d\n", successCode);
#endif
//	  puts ("Leaving joinKUnitigsFromMates\n"); fflush (stdout);
	  if (! successCode)
//	       goto afterSuperRead;
	       goto outputTheReadsIndividually;
	  // Now doing the back trace
//	  puts ("Got to 5\n"); fflush (stdout);
	  // We now start from the target k-unitig and work backwards, only using k-unitig
	  // overlaps that join to nodes we saw on the outward path
	  curPathNum = 0;
	  approxNumPaths = 0;
	  beginUnitig = mateUnitig1; beginUnitigOri = mateUnitig1ori;
	  endUnitig = mateUnitig2; endUnitigOri = mateUnitig2ori;
	  // If the target k-unitig was not found then don't try to join
	  // (treeArr is the array of distances at which occur the k-unitig in question, in this case mateUnitig2)
          if(treeArr.find(mateUnitig2) == treeArr.end())
	       goto afterSuperRead;
	  treeSize = 0;
	  edgeList.clear();
	  endingNodes.clear();
	  fwdConnections.clear();
	  revConnections.clear();
          end_tree = treeArr.find(endUnitig);
          for(auto it = end_tree->second.begin(); it != end_tree->second.end(); ++it)
	       funcToGetTreeSize(it);
#ifdef KILLED111115
	  fprintf (stderr, "treeSize = %d\n", treeSize);
#ifndef NO_OUTPUT
	  fprintf (stderr, "endUnitig = %d\n", endUnitig); // This prints
#endif
#endif
	  // splitJoinWindownMin and Max are used to disambiguate when there are multiple paths
	  // They are the offsets for the minimum offset for a split and the maximal offset of
	  // a join when there are multiple paths. To be able to disambiguate multiple paths,
	  // we must be able to select some k-unitigs over others in the region between the
	  // splitJoinWindowMin and Max. Outside of this region we are dealing with a unique path.
	  splitJoinWindowMin = INT_MAX;
	  splitJoinWindowMax = INT_MIN;
	  // Don't bother if treeSize > 1; we must disambiguate
	  if (treeSize > 1) { // treeSize indicates the number of possible offsets for the target unitig
	       splitJoinWindowMin = INT_MIN;
	       splitJoinWindowMax = INT_MAX; }
	  
	  unitigConnectionsForPathData.clear();
	  numUnitigConnectionsForPathData = 0;
#if DEBUG
	  fprintf (stderr, "Got to end of simply joinable section\n");
#endif
	  // This does the back trace from the target k-unitig
	  for(auto it = end_tree->second.begin(); it != end_tree->second.end(); ++it) {
	       completePathPrint(it);
	       if (args.join_aggressive_arg > 0)
		    if (approxNumPaths > 0)
			 break;
	  }
//	  minPathNumInForwardDirection = curPathNum;
#if 0
	  fprintf (stderr, "approxNumPaths = %d\n", (int) approxNumPaths);
#endif
	  // If the number of paths was 0, the mates were unjoinable (possibly due to too many nodes)
	  // If the number of paths was 1, the mates were uniquely joinable
	  // If the number of paths was > 1, the mates were joinable but not uniquely so, although the
	  // actual correct number of paths should be similar to but not necessarily equal to the 
	  // number of paths reported
	  if (approxNumPaths >= 1)
	       isJoinable = true;
	  if (approxNumPaths == 1) {
	       wasJoined = true;
	       ++numSimplyJoinable;
	       joinCode = 2; }
	  if (approxNumPaths <= 1)
	       goto afterSuperRead;
	  // approxNumPaths > 1, joinable but not uniquely so
#ifdef DEBUG
	  fprintf (stderr, "Trying to advance the unitigs in the first read to join them uniquely\n");
#endif
#if 1
	  // nodeArray has unitigNumber and frontEdgeOffset (the latter is binned)
	  // This sorts them in frontEdgeOffset order.
	  sort (nodeArray.begin(), nodeArray.end(), unitigLocStructCompare);
	  // Create the reverse map
	  for (unsigned int i=0; i<nodeArray.size(); i++) {
	       unitigLocMap_iterator it;
	       it = nodeToIndexMap.find (nodeArray[i]);
	       int j = it->second;
	       it->second = i;
	       newNodeNumsFromOld[j] = i; }
	  // replace the node names in edgeList with the new (reduced) node names from newNodeNumsFromOld
	  // with this we are working only with the nodes that appear on the reverse path instead of all k-unitigs
	  // we have found in our searches. this helps us determine the min and max for the window as described above.
	  // this will also allow us to output a restricted set of paths rather than one or nothing
	  // which we may be able to use in the future, as well as allowing us to much more effectively debug the program.
	  sortedEdgeList.clear();
	  for (edge_iterator it3141=edgeList.begin(); it3141 != edgeList.end(); it3141++) {
	       sortedEdgeList.insert(std::pair<int, int> (newNodeNumsFromOld[it3141->first], newNodeNumsFromOld[it3141->second]));
	  }
	  edgeList = sortedEdgeList;
#ifdef KILL120102
	  fprintf (stderr, "Sorted node list:\n");
	  for (unsigned int i=0; i<nodeArray.size(); i++)
	       fprintf (stderr, "i = %u, uni = %d, frontEdgeOffset = %d, ori = %c\n", i, (int) nodeArray[i].unitig2, (int) nodeArray[i].frontEdgeOffset, nodeArray[i].ori);
#endif
#ifdef KILL120103
	  sprintf (stderrOutputString, "%s%lld\n", rdPrefixHold, readNumHold);
	  for (edge_iterator it3141=edgeList.begin(); it3141 != edgeList.end(); it3141++) {
	       sprintf_append(stderrOutputString, "Node %d %d %c -> %d %d %c\n", (int) nodeArray[it3141->first].unitig2, (int) nodeArray[it3141->first].frontEdgeOffset, (char) nodeArray[it3141->first].ori, (int) nodeArray[it3141->second].unitig2, (int) nodeArray[it3141->second].frontEdgeOffset, (char) nodeArray[it3141->second].ori);
	  }
#endif
#endif
	  // creating the joins between this restricted node set
	  struct nodePair tNodePair;
	  for (edge_iterator it3141=edgeList.begin(); it3141 != edgeList.end(); it3141++) {
	       tNodePair.node1 = it3141->first;
	       tNodePair.node2 = it3141->second;
	       fwdConnections.push_back(tNodePair);
	       revConnections.push_back(tNodePair);
	  }
	  sort (fwdConnections.begin(), fwdConnections.end(), firstNodeSort);
	  sort (revConnections.begin(), revConnections.end(), secondNodeSort);
	  fwdStartIndices.clear();
	  revStartIndices.clear();
	  fwdNumIndices.clear();
	  revNumIndices.clear();
	  pathNumArray.clear();
	  for (unsigned int j=0; j<nodeArray.size(); j++) {
	       fwdNumIndices[j] = revNumIndices[j] = 0;
	       fwdStartIndices[j] = revStartIndices[j] = 0; }
	  ++pathNum;
	  // path num array keeps track of the maximal-numbered path on which a node is to be found
	  // (a node is a (k-unitig, offset) pair
#if DEBUG
	  fprintf (stderr, "Reporting pathNumArray before overwriting all with %d:\n", (int) pathNum);
	  for (unsigned int j=0; j<nodeArray.size(); j++)
	       fprintf (stderr, "pathNumArray[%d] = %d\n", (int) j, (int) pathNumArray[j]);
#endif
	  // Re-setting pathNumArray since we will (hopefully) be killing off some possible paths
	  // Loading the reduced set of joins, forward and reverse
	  for (unsigned int j=0; j<nodeArray.size(); j++)
	       pathNumArray[j] = pathNum;
	  for (unsigned int j=0; j<fwdConnections.size(); j++) {
	       fwdStartIndices[fwdConnections[j].node1] = j;
	       ++fwdNumIndices[fwdConnections[j].node1];
	  }
#ifdef KILL120102
	  fprintf (stderr, "At 3 revStartIndices[1] = %d\n", (int) revStartIndices[1]);
#endif
	  for (unsigned int j=0; j<revConnections.size(); j++) {
	       revStartIndices[revConnections[j].node2] = j;
	       ++revNumIndices[revConnections[j].node2]; }
#ifdef KILL120102
	  for (unsigned int j=0; j<nodeArray.size(); j++)
	       fprintf (stderr, "At 4 revStartIndices[1] = %d\n", (int) revStartIndices[1]);
#endif
	  for (unsigned int j=0; j<fwdConnections.size(); j++) {
	       if (fwdNumIndices[j] > 0)
		    fwdStartIndices[j] -= (fwdNumIndices[j]-1);
	       else
		    fwdStartIndices[j] = 0; }
#ifdef KILL120102
	  fprintf (stderr, "At 5 revStartIndices[1] = %d\n", (int) revStartIndices[1]);
#endif
	  for (unsigned int j=0; j<revConnections.size(); j++) {
	       if (revNumIndices[j] > 0)
		    revStartIndices[j] -= (revNumIndices[j]-1);
	       else
		    revStartIndices[j] = 0; }

#ifdef KILL120102
	  fprintf (stderr, "revStartIndices[1] = %d\n", (int) revStartIndices[1]);
#endif
//     mustSplit1:
//	  fprintf (stderr, "minPathNumInForwardDirection = %d\n", (int) minPathNumInForwardDirection);
	  // if doMinimalWorkHere is set, it means that we were not able to eliminate paths by advancing
	  // along the first read
	  doMinimalWorkHere = 0;
	  // The first read has only one match to a k-unitig, no advancement possible
	  if (evenReadMatchStructs.size() == 1) {
	       doMinimalWorkHere = 1; }
	  if (evenReadMatchStructs[0].ori == 'F')
	       startValue = evenReadMatchStructs[0].ahg;
	  else
	       startValue = - evenReadMatchStructs[0].bhg;
	  startValue += evenReadMatchStructs[0].readLength;
#ifdef KILL120102
	  fprintf (stderr, "Starting nodes:\n");
	  for (it1=startingNodes.begin(); it1!= startingNodes.end(); it1++)
	       fprintf (stderr, "%d %d %c\n", it1->unitig2, it1->frontEdgeOffset, it1->ori);
	  fprintf (stderr, "Ending nodes:\n");
	  for (it1=endingNodes.begin(); it1!= endingNodes.end(); it1++)
	       fprintf (stderr, "%d %d %c\n", it1->unitig2, it1->frontEdgeOffset, it1->ori);
#endif
	  // Starting from the last k-unitig matching to the first read and working backwards,
	  // we find the k-unitig with the greatest offset into the first read which appears
	  // on a path and whose offset is between splitJoinWindowMin and Max, and break out of
	  // the loop at that k-unitig
	  for (long int i=evenReadMatchStructs.size()-1; i>=0; i--) {
	       abbULS1.ori = evenReadMatchStructs[i].ori;
	       // Find the starting offset for the k-unitig
	       if (evenReadMatchStructs[i].ori == 'F')
		    abbULS1.frontEdgeOffset = startValue - evenReadMatchStructs[i].bhg;
	       else
		    abbULS1.frontEdgeOffset = startValue + evenReadMatchStructs[i].ahg;
	       // Did we find this k-unitig on the forward path?
               auto match_tree = treeArr.find(evenReadMatchStructs[i].kUnitigNumber);
	       // If not, don't continue with this k-unitig
               if(match_tree == treeArr.end())
		    continue;
	       // If we get here we have a unitig on the path
               auto element = find_within(match_tree->second, abbULS1,
                                          default_max_offset_considered_same);
               pathNum = element->pathNum;
#ifdef KILL120102
	       fprintf (stderr, "%s %d %d %d SUCCESS %d %d %d\n", (char *) readNameSpace, treeSize, (int) i, abbULS1.frontEdgeOffset, splitJoinWindowMin, splitJoinWindowMax, pathNum);
#endif
	       if (pathNum == 0) // If this (k-unitig, offset) pair were not found on the reverse traversal, skip
		    continue;
	       if (abbULS1.frontEdgeOffset >= splitJoinWindowMax) // If the offset of the k-unitig puts it beyond the repeat region, skip
		    continue;
	       // Checking for useless result
	       // If the offset puts it before the repeat region, since all subsequent k-unitig matches with the first read will have
	       // offsets less than this one, we indicate that we can get no information by tring to split using the k-unitigs matching
	       // the first read.
	       if (abbULS1.frontEdgeOffset <= splitJoinWindowMin) { 
		    doMinimalWorkHere = 1;
		    break; }
	       // If we get here we've passed all the tests and we can use it
	       localUnitigNumber = evenReadMatchStructs[i].kUnitigNumber;
	       overlapMatchIndexHold = i;
	       break;
	  }
	  
	  // Nothing is gained by moving the first k-unitig forward, so we just the same first
	  // k-unitig as before
	  if (doMinimalWorkHere) {
	       localUnitigNumber = evenReadMatchStructs[0].kUnitigNumber;
	       abbULS1.frontEdgeOffset = unitigLengths[localUnitigNumber];
	       abbULS1.ori = evenReadMatchStructs[0].ori; }
	  // tempULS is the k-unitig from which we will now be trying to join
	  tempULS.unitig2 = localUnitigNumber;
	  tempULS.frontEdgeOffset = abbULS1.frontEdgeOffset;
	  tempULS.ori = abbULS1.ori;
#ifdef KILL120102
	  fprintf (stderr, "unitig2 = %d, frontEdgeOffset = %d, ori = %c\n", tempULS.unitig2, tempULS.frontEdgeOffset, tempULS.ori);
#endif
	  localNodeNumber = localNodeNumberHold = nodeToIndexMap[tempULS]; // local node number for the (new) starting k-unitig
#ifdef KILL120102
	  fprintf (stderr, "localNodeNumber = %d\n", (int) localNodeNumber);
#endif
	  // Change to a pathNum larger than what we have used before for this pair so we can keep track of both
	  // (a) if the node was on a path before, and (b) if the node is on the new (restricted) path
	  ++pathNum;
#ifdef DEBUG120626
	  fprintf (stderr, "Changing path numbers for nodes at 1\n");
	  fprintf (stderr, "localNodeNumber = %d, pathNum = %d\n", (int) localNodeNumber, (int) pathNum);
#endif
	  pathNumArray[localNodeNumber] = pathNum; // put the local node number as the starting k-unitig and mark all the k-unitigs on the
	  // path from here to the ending k-unitig as having the new pathNum; in this way we know which nodes are achievable from the new starting node
	  nodeIntArray.push(localNodeNumber);
	  while (! nodeIntArray.empty()) {
	       int localLoopNodeNumber = nodeIntArray.top();
	       nodeIntArray.pop();
	       // for each forward overlap with the current node, add the new node to the possible nodes, since if we could join this node to
	       // target node before we certainly still can
	       for (int j=fwdStartIndices[localLoopNodeNumber]; j<fwdStartIndices[localLoopNodeNumber]+fwdNumIndices[localLoopNodeNumber]; j++) {
		    localNodeNumber = fwdConnections[j].node2;
		    if (pathNumArray[localNodeNumber] < pathNum) {
			 nodeIntArray.push (localNodeNumber);
#ifdef DEBUG120626
			 fprintf (stderr, "localNodeNumber = %d, pathNum = %d, old value = %d (in loop 1)\n", (int) localNodeNumber, (int) pathNum, (int) pathNumArray[localNodeNumber]);
#endif
			 pathNumArray[localNodeNumber] = pathNum; }
	       }
	  }
	  if (doMinimalWorkHere) { // Nothing more to do for this section; none of the changes has disambiguated anything, so we go on to trying to split using the second read's k-unitigs
	       for (it1=startingNodes.begin(); it1!= startingNodes.end(); it1++) { // Must set start node
		    tempULS.unitig2 = it1->unitig2;
		    tempULS.frontEdgeOffset = it1->frontEdgeOffset;
		    tempULS.ori = it1->ori;
		    startingNodeNumber = nodeToIndexMap[tempULS]; }
	       goto mustSplit2;
	  }
#ifdef DEBUG120626
	  fprintf (stderr, "overlapMatchIndexHold = %d\n", (int) overlapMatchIndexHold);
#endif

	  // Now collecting the prefix to the path generated immediately above
	  // For the overlaps involving k-unitigs before the one with which we started
	  for (int i=overlapMatchIndexHold-1; i>=0; i--) {
	       unitigLocMap_iterator it;
	       int isGood = 0;
	       tempULS.ori = evenReadMatchStructs[i].ori;
	       // startValue = offset of k-unitig to which we had advanced in the first read
	       if (evenReadMatchStructs[i].ori == 'F')
		    tempULS.frontEdgeOffset = startValue - evenReadMatchStructs[i].bhg;
	       else
		    tempULS.frontEdgeOffset = startValue + evenReadMatchStructs[i].ahg;
	       tempULS.unitig2 = evenReadMatchStructs[i].kUnitigNumber;
	       it = nodeToIndexMap.find (tempULS);
	       if (it == nodeToIndexMap.end()) // If the k-unitig not found at the offset
		    break;
	       int nodeNum = it->second;
	       if (pathNumArray[nodeNum] == 0)
		    break;
	       isGood = 0;
	       for (int j=revStartIndices[localNodeNumberHold]; j<revStartIndices[localNodeNumberHold]+revNumIndices[localNodeNumberHold]; j++) {
		    localNodeNumber = revConnections[j].node1;
		    if (localNodeNumber == nodeNum) {
#ifdef DEBUG120626
			 fprintf (stderr, "nodeNum = %d, pathNum = %d, old value = %d (in loop 2)\n", (int) nodeNum, (int) pathNum, (int) pathNumArray[nodeNum]);
#endif
			 pathNumArray[nodeNum] = pathNum;
			 localNodeNumberHold = i;
			 isGood = 1;
			 break; }
	       }
	       if (! isGood)
		    break;
	       localNodeNumberHold = nodeNum;
	  }
	       
	  nodeIntArray.push(localNodeNumberHold);
	  while (! nodeIntArray.empty()) {
	       int localLoopNodeNumber = nodeIntArray.top();
	       nodeIntArray.pop();
	       for (int j=revStartIndices[localLoopNodeNumber]; j<revStartIndices[localLoopNodeNumber]+revNumIndices[localLoopNodeNumber]; j++) {
		    localNodeNumber = revConnections[j].node1;
		    if (pathNumArray[localNodeNumber] < pathNum) {
			 nodeIntArray.push (localNodeNumber);
#ifdef DEBUG120626
			 fprintf (stderr, "localNodeNumber = %d, pathNum = %d, old value = %d (in loop 3)\n", (int) localNodeNumber, (int) pathNum, (int) pathNumArray[localNodeNumber]);
#endif
			 pathNumArray[localNodeNumber] = pathNum; }
	       }
	  }
	  lastGoodNodeNumber = -1;
	  for (it1=startingNodes.begin(); it1!= startingNodes.end(); it1++) {
	       tempULS.unitig2 = it1->unitig2;
	       tempULS.frontEdgeOffset = it1->frontEdgeOffset;
	       tempULS.ori = it1->ori;
	       startingNodeNumber = nodeToIndexMap[tempULS]; }
	  if (nodeIntArray.size() > 0)
	       fprintf (stderr, "ERROR in nodeIntArray: size should be 0\n");
	  unitigNodeNumbersForPath.clear();
	  nodeIntArray.push (startingNodeNumber);
	  while (! nodeIntArray.empty()) {
	       int localLoopNodeNumber = nodeIntArray.top();
	       unitigNodeNumbersForPath.push_back(localLoopNodeNumber);
	       nodeIntArray.pop();
	       for (int j=fwdStartIndices[localLoopNodeNumber]; j<fwdStartIndices[localLoopNodeNumber]+fwdNumIndices[localLoopNodeNumber]; j++) {
                    localNodeNumber = fwdConnections[j].node2;
                    if (pathNumArray[localNodeNumber] == pathNum) {
                         nodeIntArray.push (localNodeNumber);
#ifdef DEBUG120626
			 fprintf (stderr, "localNodeNumber = %d , pathNum = %d, old value = %d (in loop 4)\n", (int) localNodeNumber, (int) pathNum, (int) pathNumArray[localNodeNumber]);
#endif
                         pathNumArray[localNodeNumber] = pathNum; }
               }
	       if (nodeIntArray.size() > 1) {
		    lastGoodNodeNumber = localLoopNodeNumber;
		    while (! nodeIntArray.empty())
			 nodeIntArray.pop();
		    break; }
          }
	  if (lastGoodNodeNumber < 0) {
	       for (numUnitigPathPrintRecsOnPath=0; numUnitigPathPrintRecsOnPath<(int) unitigNodeNumbersForPath.size(); numUnitigPathPrintRecsOnPath++) {
		    augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].unitig1 = nodeArray[unitigNodeNumbersForPath[numUnitigPathPrintRecsOnPath]].unitig2;
		    augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].frontEdgeOffset = nodeArray[unitigNodeNumbersForPath[numUnitigPathPrintRecsOnPath]].frontEdgeOffset;
		    if (numUnitigPathPrintRecsOnPath > 0)
			 augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsIn = 1;
		    else
			 augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsIn = 0;
		    if (numUnitigPathPrintRecsOnPath<(int)unitigNodeNumbersForPath.size() - 1)
			 augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsOut = 1;
		    else
			 augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsOut = 0;
		    augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].ori = nodeArray[unitigNodeNumbersForPath[numUnitigPathPrintRecsOnPath]].ori;
	       }
	       joinCode = 2;
	       generateSuperReadPlacementLinesForJoinedMates();
	       approxNumPaths = 1;
	       ++numJoinableAfterRead1Analysis;
	       wasJoined = true;
#ifdef KILL120102
	       fprintf (stderr, "We have successfully joined the mates after advancing the k-unitig in the first read!\n");
	       for (unsigned int i=0; i<unitigNodeNumbersForPath.size(); i++) {
		    fprintf (stderr, "%d %c ", nodeArray[unitigNodeNumbersForPath[i]].unitig2, nodeArray[unitigNodeNumbersForPath[i]].ori);
		    fprintf (stderr, "\n"); }
#endif
	       goto afterSuperRead; }
#ifdef KILL120102
	  fprintf (stderr, "The node before the split is uni = %d offset = %d ori = %c\n", nodeArray[lastGoodNodeNumber].unitig2, nodeArray[lastGoodNodeNumber].frontEdgeOffset, nodeArray[lastGoodNodeNumber].ori);
#endif	  
	  
     mustSplit2:
#ifdef DEBUG
	  fprintf (stderr, "Trying to advance the unitigs in the second read to join them uniquely\n");
#endif

	  if (oddReadMatchStructs.size() == 1)
	       goto afterSuperRead;
	  // startValue is the number of bases in the super-read beyond
	  // the last base in the read
          if (oddReadMatchStructs[0].ori == 'F')
               startValue = oddReadMatchStructs[0].ahg;
          else
               startValue = - oddReadMatchStructs[0].bhg;

	  // Here we are measuring the distance of the frontEdgeOffset
	  // from the end of the super-read (using a positive distance)
          last_element_is_nil = false;
          for (int i=oddReadMatchStructs.size()-1; i>=0; i--) {
	       int distFromEndOfSuperRead = 0;
	       if (oddReadMatchStructs[i].ori == 'F')
		    abbULS1.ori = 'R';
	       else
		    abbULS1.ori = 'F';
	       if (oddReadMatchStructs[i].ori == 'F')
		    distFromEndOfSuperRead = startValue - oddReadMatchStructs[i].ahg;
	       else
		    distFromEndOfSuperRead = startValue + oddReadMatchStructs[i].bhg;
	       // We now have to modify the offset so as to be in relation to the
	       // left end of the super-read. We must consider all possible
	       // distances to do this.
	       numPossibleLengths = 0;
	       for (it1=endingNodes.begin(); it1!=endingNodes.end(); it1++) {
		    abbULS1.frontEdgeOffset = it1->frontEdgeOffset - distFromEndOfSuperRead;
                    auto match_tree = treeArr.find(oddReadMatchStructs[i].kUnitigNumber);
                    if(match_tree == treeArr.end()) {
			 last_element_is_nil = true;
			 continue;
                    }
                    auto element = find_within(match_tree->second, abbULS1,
                                               default_max_offset_considered_same);
                    if(element == match_tree->second.end()) {
			 last_element_is_nil = true;
			 continue;
                    }
                    last_element_is_nil = false;
		    tempULS.unitig2 = oddReadMatchStructs[i].kUnitigNumber;
		    tempULS.frontEdgeOffset = abbULS1.frontEdgeOffset;
		    tempULS.ori = abbULS1.ori;
#if DEBUG
		    fprintf (stderr, "At line 2090: unitig = %d, frontEdgeOffset = %d, ori = %c\n", (int) tempULS.unitig2, (int) tempULS.frontEdgeOffset, (char) tempULS.ori);
#endif
		    it2 = nodeToIndexMap.find (tempULS);
		    if (it2 == nodeToIndexMap.end())
			 continue;
#ifdef DEBUG120911
		    fprintf (stderr, "it2->second = %d, pathNumArray[%d] = %d, pathNum = %d\n", (int) it2->second,  (int) it2->second, (int) pathNumArray[(int) it2->second], (int) pathNum);
#endif
		    if (pathNumArray[it2->second] != pathNum)
			 continue;
#ifdef DEBUG120911
		    fprintf (stderr, "At line 2127\n");
#endif
		    localSuperReadLength = it1->frontEdgeOffset;
		    abbULS1Hold = abbULS1;
		    ++numPossibleLengths;
	       }
	       if (numPossibleLengths != 1) {
#ifdef KILL120102
		    if (numPossibleLengths == 0)
			 fprintf (stderr, "No good read2 unitig match %s %d %d %d FAIL %d %d\n", (char *) readNameSpace, treeSize, i, abbULS1.frontEdgeOffset, splitJoinWindowMin, splitJoinWindowMax);
		    else
			 fprintf (stderr, "Too many good read2 unitig matches. Last is %s %d %d %d FAIL %d %d\n", (char *) readNameSpace, treeSize, i, abbULS1.frontEdgeOffset, splitJoinWindowMin, splitJoinWindowMax);
#endif
		    continue;
	       }
               // If we get here we have a unitig on the path
#ifdef KILL120102
	       fprintf (stderr, "At 2145 %s %d %d %d SUCCESS %d %d %d\n", (char *) readNameSpace, treeSize, i, abbULS1Hold.frontEdgeOffset, splitJoinWindowMin, splitJoinWindowMax, pathNum);
#endif
               if (abbULS1Hold.frontEdgeOffset <= splitJoinWindowMin)
		    continue;
#ifdef DEBUG120911
	       fprintf (stderr, "At 2149\n");
#endif
               // Checking for useless result
#ifdef DEBUG120911
	       fprintf (stderr, "frontEdgeOffset = %d, splitJoinWindowMax = %d\n", (int) abbULS1Hold.frontEdgeOffset, (int) splitJoinWindowMax);
#endif
               if (abbULS1Hold.frontEdgeOffset >= splitJoinWindowMax)
                    goto afterSuperRead;
               // If we get here we've passed all the tests and we can use it
#ifdef DEBUG120911
	       fprintf (stderr, "At 2157\n");
#endif
               localUnitigNumber = oddReadMatchStructs[i].kUnitigNumber;
	       localFrontEdgeOffset = abbULS1Hold.frontEdgeOffset;
               overlapMatchIndexHold = i;
	       break;
          }
          if(last_element_is_nil)
	       goto afterSuperRead;
	       
          tempULS.unitig2 = localUnitigNumber;
          tempULS.ori = abbULS1Hold.ori;
          tempULS.frontEdgeOffset = localFrontEdgeOffset;
#ifdef KILL120102
          fprintf (stderr, "unitig2 = %d, frontEdgeOffset = %d, ori = %c\n", tempULS.unitig2, tempULS.frontEdgeOffset, tempULS.ori);
#endif
          localNodeNumber = localNodeNumberHold = nodeToIndexMap[tempULS];
#ifdef KILL120102
	  fprintf (stderr, "localNodeNumber = %d\n", (int) localNodeNumber);
#endif
          ++pathNum;
#ifdef DEBUG120626
          fprintf (stderr, "Changing path numbers for nodes at 5\n");
          fprintf (stderr, "localNodeNumber = %d, pathNum = %d\n", (int) localNodeNumber, (int) pathNum);
#endif
          pathNumArray[localNodeNumber] = pathNum;
          nodeIntArray.push(localNodeNumber);
          while (! nodeIntArray.empty()) {
               int localLoopNodeNumber = nodeIntArray.top();
               nodeIntArray.pop();
               for (int j=revStartIndices[localLoopNodeNumber]; j<revStartIndices[localLoopNodeNumber]+revNumIndices[localLoopNodeNumber]; j++) {
                    localNodeNumber = revConnections[j].node1;
                    if (pathNumArray[localNodeNumber] == pathNum-1) {
                         nodeIntArray.push (localNodeNumber);
#ifdef DEBUG120626
                         fprintf (stderr, "localNodeNumber = %d, pathNum = %d, old value = %d (in loop 5)\n", (int) localNodeNumber, (int) pathNum, (int) pathNumArray[localNodeNumber]);
#endif
                         pathNumArray[localNodeNumber] = pathNum; }
               }
          }
#ifdef DEBUG120626
          fprintf (stderr, "overlapMatchIndexHold = %d\n", (int) overlapMatchIndexHold);
#endif

          for (int i=overlapMatchIndexHold-1; i>=0; i--) {
               unitigLocMap_iterator it;
               int isGood = 0;
#ifdef DEBUG120626
	       fprintf (stderr, "i = %d\n", (int) i);
#endif
	       if (oddReadMatchStructs[i].ori == 'F')
		    tempULS.ori = 'R';
	       else
		    tempULS.ori = 'F';
               if (oddReadMatchStructs[i].ori == 'F')
                    distFromEndOfSuperRead = startValue - oddReadMatchStructs[i].ahg;
               else
		    distFromEndOfSuperRead = startValue + oddReadMatchStructs[i].bhg;
	       tempULS.frontEdgeOffset = localSuperReadLength - distFromEndOfSuperRead;
#ifdef DEBUG120626
	       fprintf (stderr, "localSuperReadLength = %d, kUnitig = %d, distFromEndOfSuperRead = %d, startValue = %d, ahg = %d, bhg = %d, ori = %c\n", (int) localSuperReadLength, (int) oddReadMatchStructs[i].kUnitigNumber, (int) distFromEndOfSuperRead, (int) startValue, (int) oddReadMatchStructs[i].ahg, (int) oddReadMatchStructs[i].bhg, (char) oddReadMatchStructs[i].ori);
#endif
               tempULS.unitig2 = oddReadMatchStructs[i].kUnitigNumber;
               it = nodeToIndexMap.find (tempULS);
#ifdef DEBUG120626
	       fprintf (stderr, "At 10101\n");
	       fprintf (stderr, "tempULS (uni, offset, ori) = (%d, %d, %c)\n", (int) tempULS.unitig2, (int) tempULS.frontEdgeOffset, (char) tempULS.ori);
#endif
               if (it == nodeToIndexMap.end())
                    break;
#ifdef DEBUG120626
	       fprintf (stderr, "At 10102\n");
#endif
               int nodeNum = it->second;
	       if (pathNumArray[nodeNum] < pathNum-1)
                    break;
               isGood = 0;
#ifdef DEBUG120626
	       fprintf (stderr, "localNodeNumberHold = %d, fwdStartIndices[%d] = %d, fwdStartIndices[%d] = %d\n", (int) localNodeNumberHold, localNodeNumberHold, fwdStartIndices[localNodeNumberHold],  localNodeNumberHold+1, fwdStartIndices[localNodeNumberHold+1]);
#endif
               for (int j=fwdStartIndices[localNodeNumberHold]; j<fwdStartIndices[localNodeNumberHold]+fwdNumIndices[localNodeNumberHold]; j++) {
                    localNodeNumber = fwdConnections[j].node2;
                    if (localNodeNumber == nodeNum) {
#ifdef DEBUG120626
                         fprintf (stderr, "nodeNum = %d, pathNum = %d, old value = %d (in loop 6)\n", (int) nodeNum, (int) pathNum, (int) pathNumArray[nodeNum]);
#endif
                         pathNumArray[nodeNum] = pathNum;
                         localNodeNumberHold = i;
                         isGood = 1;
                         break; }
               }
               if (! isGood)
                    break;
               localNodeNumberHold = nodeNum;
          }

          nodeIntArray.push(localNodeNumberHold);

          while (! nodeIntArray.empty()) {
               int localLoopNodeNumber = nodeIntArray.top();
               nodeIntArray.pop();
               for (int j=fwdStartIndices[localLoopNodeNumber]; j<fwdStartIndices[localLoopNodeNumber]+fwdNumIndices[localLoopNodeNumber]; j++) {
                    localNodeNumber = fwdConnections[j].node2;
                    if (pathNumArray[localNodeNumber] == pathNum-1) {
                         nodeIntArray.push (localNodeNumber);
#ifdef DEBUG120626
                         fprintf (stderr, "localNodeNumber = %d, pathNum = %d, old value = %d (in loop 7)\n", (int) localNodeNumber, (int) pathNum, (int) pathNumArray[localNodeNumber]);
#endif
                         pathNumArray[localNodeNumber] = pathNum; }
               }
          }
               
	  // We now analyze the path we have to see if it's unique
          lastGoodNodeNumber = -1;
          if (nodeIntArray.size() > 0)
               fprintf (stderr, "ERROR in nodeIntArray: size should be 0\n");
          unitigNodeNumbersForPath.clear();
          nodeIntArray.push (startingNodeNumber);
          while (! nodeIntArray.empty()) {
               int localLoopNodeNumber = nodeIntArray.top();
#if DEBUG
	       fprintf (stderr, "Starting loop for node %d, unitig = %d\n", localLoopNodeNumber, (int) nodeArray[localLoopNodeNumber].unitig2);
#endif
               unitigNodeNumbersForPath.push_back(localLoopNodeNumber);
               nodeIntArray.pop();
               for (int j=fwdStartIndices[localLoopNodeNumber]; j<fwdStartIndices[localLoopNodeNumber]+fwdNumIndices[localLoopNodeNumber]; j++) {
                    localNodeNumber = fwdConnections[j].node2;
                    if (pathNumArray[localNodeNumber] == pathNum) {
                         nodeIntArray.push (localNodeNumber);
#ifdef DEBUG120626
                         fprintf (stderr, "localNodeNumber = %d , pathNum = %d, old value = %d (in loop 8)\n", (int) localNodeNumber, (int) pathNum, (int) pathNumArray[localNodeNumber]);
#endif
#if DEBUG
			 fprintf (stderr, "Pushing node %d (unitig %d) in the loop\n", (int) localNodeNumber, nodeArray[localNodeNumber].unitig2);
#endif
                         pathNumArray[localNodeNumber] = pathNum; }
               }
#if DEBUG
	       fprintf (stderr, "unitig = %d, numPaths = %d\n", nodeArray[unitigNodeNumbersForPath[numUnitigPathPrintRecsOnPath]].unitig2, (int) nodeIntArray.size());
#endif
               if (nodeIntArray.size() > 1) {
		    lastGoodNodeNumber = localLoopNodeNumber;
		    if (args.join_aggressive_arg) {
			 ++numSimplyJoinable;
			 joinCode = 2;
			 wasDeclaredUnresolvedJoinable = false; }
		    else {
			 ++numJoinableUnresolvedAtEnd;
			 joinCode = 3;
			 wasDeclaredUnresolvedJoinable = true;
		    }
                    while (! nodeIntArray.empty())
                         nodeIntArray.pop();
                    break; }
          }
          if (lastGoodNodeNumber < 0) {
               for (numUnitigPathPrintRecsOnPath=0; numUnitigPathPrintRecsOnPath<(int) unitigNodeNumbersForPath.size(); numUnitigPathPrintRecsOnPath++) {
                    augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].unitig1 = nodeArray[unitigNodeNumbersForPath[numUnitigPathPrintRecsOnPath]].unitig2;
                    augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].frontEdgeOffset = nodeArray[unitigNodeNumbersForPath[numUnitigPathPrintRecsOnPath]].frontEdgeOffset;
                    if (numUnitigPathPrintRecsOnPath > 0)
                         augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsIn = 1;
                    else
                         augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsIn = 0;
                    if (numUnitigPathPrintRecsOnPath<(int)unitigNodeNumbersForPath.size() - 1)
                         augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsOut = 1;
                    else
                         augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].numOverlapsOut = 0;
                    augmentedUnitigPathPrintData[numUnitigPathPrintRecsOnPath].ori = nodeArray[unitigNodeNumbersForPath[numUnitigPathPrintRecsOnPath]].ori;
               }
	       joinCode = 2;
               generateSuperReadPlacementLinesForJoinedMates();
               approxNumPaths = 1;
	       ++numJoinableAfterBothReadAnalysis;
	       wasJoined = true;
#ifdef KILL120102
	       fprintf (stderr, "We have successfully joined the mates!\n");
	       for (unsigned int i=0; i<unitigNodeNumbersForPath.size(); i++) {
		    fprintf (stderr, "%d %c ", nodeArray[unitigNodeNumbersForPath[i]].unitig2, nodeArray[unitigNodeNumbersForPath[i]].ori);
                    fprintf (stderr, "\n"); }
#endif
               goto afterSuperRead; }
#ifdef KILL120102
          fprintf (stderr, "The node before the split is uni = %d offset = %d ori = %c\n", nodeArray[lastGoodNodeNumber].unitig2, nodeArray[lastGoodNodeNumber].frontEdgeOffset, nodeArray[lastGoodNodeNumber].ori);
#endif
          
     afterSuperRead:
	  
#ifdef KILLED111115
	  fprintf (stderr, "Approx num paths returned = %d\n", approxNumPaths);
#endif
	  
	  // Doing the output (if possible)
	  
	  if (approxNumPaths == 1) {
	       m_out << outputString;
	       if (isJoinable && (! wasJoined) && (! wasDeclaredUnresolvedJoinable)) {
		    ++numJoinableUnresolvedAtEnd;
		    joinCode = 3; }
	       else if (successCode == 0) {
		    ++numUnjoinableMissingSequence;
		    joinCode = 4; }
	       //	       fputs (outputString, outputFile);
	       return; }
#ifdef KILL120103
	  if (stderrOutputString[0] != 0)
	       fputs (stderrOutputString, stderr);
#endif
	  if (approxNumPaths < 1) {
#if DEBUG
	       if (tooManyPossibleInsertLengths)
		    fprintf (stderr, "tooManyPossibleInsertLengths = true\n");
#endif
	       if (tooManyPossibleInsertLengths) {
		    if (args.join_aggressive_arg) {
			 ++numSimplyJoinable;
			 joinCode = 2;
			 wasDeclaredUnresolvedJoinable = false; }
		    else {
			 ++numJoinableUnresolvedAtEnd;
			 joinCode = 3;
			 wasDeclaredUnresolvedJoinable = true; }
	       }
	       else {
		    ++numUnjoinableMissingSequence; 
		    joinCode = 4; }
	       goto outputTheReadsIndividually; }
	  // Now we move up the unitig on read1 and try again
	  
     outputTheReadsIndividually:
//	  fprintf (stderr, "readNumHold = %d\n", (int) readNumHold);
	  sprintf (readNameSpace, "%s%lld", rdPrefixHold, readNumHold-1);
	  findSingleReadSuperReads(readNameSpace, m_out);
	  sprintf (readNameSpace, "%s%lld", rdPrefixHold, readNumHold);
	  findSingleReadSuperReads(readNameSpace, m_out);
	  if (isJoinable && (! wasJoined) && (! wasDeclaredUnresolvedJoinable)) {
	       ++numJoinableUnresolvedAtEnd;
	       joinCode = 3;
	  }

     }

    return;
}

bool firstNodeSort (struct nodePair val1, struct nodePair val2)
{
     return (val1.node1 < val2.node1);
}
 
bool secondNodeSort (struct nodePair val1, struct nodePair val2)
{
     return (val1.node2 < val2.node2);
}

