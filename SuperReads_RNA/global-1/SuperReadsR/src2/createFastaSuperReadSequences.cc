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


/*
I)   Allocate space for the k-unitigs
        A) Use the length of guillaumeKUnitigsAtLeast32bases_all.fasta
II)  Read in the file with the k-unitigs
        A) guillaumeKUnitigsAtLeast32bases_all.fasta
III) Allocate space for the pointers to the k-unitigs
        A) Find the largest number assigned to a k-unitig
              1) Use the files guillaumeKUnitigsAtLeast32bases_*.fa
IV)  Process the file (in memory) with the k-unitigs
        A) Set the pointers to the k-unitigs
        B) Eliminate spaces from the k-unitigs
V)   For each super-read
        A) Create the header
              1) Use the line with the order of the k-unitigs, their oris, and
                   their overlaps
                    a) No space between the k-unitig number and its orientation
                    b) Underscores between all other fields
        B) Generate the sequence
              1) Create the reverse-complement of a k-unitig if necessary
                    a) Must generate space if necessary ahead of time
              2) Can use strcat and strcpy
              3) Keep track of the output length
        C) Output the fastq format for the read
        D) Output quals 'a' for the length of the read
VI)  Make sure you allow for an appropriate set of params
        A) Working directory (where to find k-unitig sequence files)
        B) File listing the super-reads (with path)
              1) allData/superReadGroups.onePerLine.withReadInfoIncluded.txt
        C) -seqdiffmax # allows one to specify the number of differences one
	     allows between the sequences of overlapping k-unitigs in a
	     super-read; otherwise it fails (default is 0).
	D) -nosequence just outputs the names of the passing super-reads
	     instead of the name and sequence in fasta format
	E) -error-filename filename sends the error output to the specified
	     file
	F) -min-ovl-len # allows one to specify the minimum overlap amount
	     between k-unitigs in the file
        G) -minreadsinsuperread # allows one to specify the minimum number
	     of inserts to allow the super-read to exist
	H) -kunitigsfile filename specifies the input k-unitig filename instead
	     of using the default:
	     $workingDir/guillaumeKUnitigsAtLeast32bases_all.fasta
	I) -maxunitignumberfile filename specifies the file which contains
	     the largest k-unitig number. The default is
	     maxKUnitigNumber.txt within the working directory
	J) -good-sequence-output-file filename sends the output super-reads
	     to 'filename' instead of to stdout
	H) -super-read-name-and-lengths-file filename creates the file
	     'filename' which has the super-read name and length for
	     each good super-read.

*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <cctype>
#include <iostream>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>
#include <string>
#include <charb.hpp>
#include <exp_vector.hpp>
#include <misc.hpp>
using namespace std;

#define NUM_KUNITIGS_FILENAME "maxKUnitigNumber.txt"
#define KUNITIG_FILE_COMPLETE "guillaumeKUnitigsAtLeast32bases_all.fasta"
#define DEFAULT_SUPER_READ_LIST_FILE "superReadGroups.onePerLine.withReadInfoIncluded.txt"
#define AFTER_NEWLINE 1
#define AFTER_HEADER_NEWLINE 2
#define IN_SEQUENCE 3
#define MAX_READ_LEN 1

char *kUnitigSpace, **kUnitigSeq;
int *kUnitigLengths;
charb line(1000000);
charb reverseComplementSpace(1000000), outputSeqSpace(3000000);
ExpBuffer<char*> flds;
bool isJumpLibrary;


void reverseTheString (char *str);
void generateReverseComplement (char *seq, int seqLen);
FILE *Fopen (const char *fn, const char *mode);

#define processTheChar if (! isspace (kUnitigSpace[i64])) { \
		           if (i64 != j64) \
			       kUnitigSpace[j64] = kUnitigSpace[i64]; \
		               ++j64; } \
	               state = IN_SEQUENCE; \
	               continue;

#define mallocOrDie(name, num, type) name = (type *) calloc (num, sizeof ( type )); \
if (name == NULL) { fprintf (stderr, "Couldn't allocate space for '%s'\nBye!\n", #name ); exit (-1); }

int main (int argc, char **argv)
{
     char *workingDir;
     charb superReadListFile(512), fname(512), kUnitigFilename(512), maxUnitigsFilename(512);
     charb goodSequenceOutputFilename(512), superReadNameAndLengthsFilename(512);
     struct stat statbuf;
     uint64_t kUnitigSeqFileSize, fsize;
     uint64_t i64, j64=0;
     charb goodFilename(500);
     FILE *infile, *goodFile = stdout; // goodFile initted makes compiler happy
     int lastKUnitigNumber, kUnitigNumber=0, kUnitigNumberHold, i, argNum;
     int state;
     char *cptr, *cptr2;
     char* superReadName;
     char ori, oriHold;
     int overlap = 30;
     int outputSeqLen = 0;
     int seqDiffMax;
     int fail;
     charb errorMessage(2000), errorMessageLine(2000);
     int numReads=0;
     char pluralStr[2];
     int noSequence = 0;
     int minReadsInSuperRead = 2;
     int minOvlLen = 30;
     int minReadLength = 64, maxReadLength = 2047;
     bool renameSuperReads = false;
     int maxSuperReadLength = 0;
     exp_vector<int> maxSuperReadNumbersUsed;
     std::map<uint64_t, string> oldSuperReadName;
     std::map<uint64_t, string> superReadSequences;

     // (VI) above
     workingDir = (char *) ".";
     sprintf (superReadListFile, "%s/%s", workingDir, DEFAULT_SUPER_READ_LIST_FILE);
     seqDiffMax = 0;
     goodFilename[0] = 0;
     argNum = 0;
     strcpy (kUnitigFilename, "");
     strcpy (maxUnitigsFilename, "");
     strcpy (goodSequenceOutputFilename, "");
     strcpy (superReadNameAndLengthsFilename, "");

     isJumpLibrary = false;
     for (i=1; i<argc; i++) {
	  if (strcmp (argv[i], "-nosequence") == 0) {
	       noSequence = 1;
	       continue; }
	  if (strcmp (argv[i], "-seqdiffmax") == 0) {
	       ++i;
	       seqDiffMax = atoi (argv[i]);
	       continue; }
	  if (strcmp (argv[i], "-good-sr-filename") == 0) {
	       ++i;
	       strcpy (goodFilename, argv[i]);
	       continue; }
	  if (strcmp (argv[i], "-min-ovl-len") == 0) {
	       ++i;
	       overlap = minOvlLen = atoi (argv[i]);
	       continue; }
	  if (strcmp (argv[i], "-minreadsinsuperread") == 0) {
	       ++i;
	       minReadsInSuperRead = atoi (argv[i]);
	       continue; }
	  if (strcmp (argv[i], "-kunitigsfile") == 0) {
	       ++i;
	       strcpy (kUnitigFilename, argv[i]);
	       continue; }
	  if (strcmp (argv[i], "-maxunitignumberfile") == 0) {
	       ++i;
	       strcpy (maxUnitigsFilename, argv[i]);
	       continue; }
	  if (strcmp (argv[i], "-good-sequence-output-file") == 0) {
	       ++i;
	       strcpy (goodSequenceOutputFilename, argv[i]);
	       continue; }
	  if (strcmp (argv[i], "-super-read-name-and-lengths-file") == 0) {
	       ++i;
	       strcpy (superReadNameAndLengthsFilename, argv[i]);
	       continue; }
	  if (strcmp (argv[i], "-min-read-length") == 0) {
	       ++i;
	       minReadLength = atoi(argv[i]);
	       continue; }
	  if (strcmp (argv[i], "-max-read-length") == 0) {
	       ++i;
	       maxReadLength = atoi(argv[i]);
	       continue; }
	  if (strcmp (argv[i], "-jump-library") == 0) {
	       isJumpLibrary = true;
	       continue; }
	  if (strcmp (argv[i], "-rename-super-reads") == 0) {
	       renameSuperReads = true;
	       continue; }
	  if (argNum == 0)
	       workingDir = argv[i];
	  else if (argNum == 1)
	       strcpy (superReadListFile, argv[i]);
	  else {
	       fprintf (stderr, "The program %s is called incorrectly. Too many args. The extra arg is '%s'. Bye!\n", argv[0], argv[i]);
	       exit (1); }
	  ++argNum; }
     
     if (strlen (goodSequenceOutputFilename) != 0)
	  if (noSequence != 0) {
	       fprintf (stderr, "-good-sequence-output-file arg used, as well as -nosequence flag.\n-nosequence flag is ignored.\n");
	       noSequence = 0; }
     if (strlen (goodFilename) == 0)
	  sprintf (goodFilename, "%s/superReadNames.txt", workingDir);

     if (kUnitigFilename[0] == 0)
	  sprintf (kUnitigFilename, "%s/%s", workingDir, KUNITIG_FILE_COMPLETE);
     stat (kUnitigFilename, &statbuf);
     fsize = statbuf.st_size;
     kUnitigSeqFileSize = fsize;
     // (I) above
     mallocOrDie (kUnitigSpace, fsize, char);
     // (II) above
     infile = Fopen (kUnitigFilename, "r");
     size_t bytes_read = fread (kUnitigSpace, 1, kUnitigSeqFileSize, infile);
     if(bytes_read != kUnitigSeqFileSize) {
       fprintf(stderr, "Failed to read the entire file '%s'. Bye!\n", (char*)fname);
       exit(2);
     }
     fclose (infile);

     // (III) above
     // Find out the last kUnitig number
     if (maxUnitigsFilename[0] == 0)
	  sprintf (maxUnitigsFilename, "%s/%s", workingDir, NUM_KUNITIGS_FILENAME);
     infile = Fopen (maxUnitigsFilename, "r");
     int fields_read = fscanf (infile, "%d\n", &lastKUnitigNumber);
     if(fields_read != 1) {
       fprintf(stderr, "Failed to read one int from '%s'. Bye!\n", (char*)maxUnitigsFilename);
       exit(2);
     }
     fclose (infile);
     mallocOrDie (kUnitigSeq, lastKUnitigNumber+1, char *);
     mallocOrDie (kUnitigLengths, lastKUnitigNumber+1, int);

     // (IV) above
     state = AFTER_NEWLINE;
     
     for (i64=0; i64<kUnitigSeqFileSize; i64++) {
	  if (kUnitigSpace[i64] == '\n') {
	       kUnitigSpace[j64] = 0;
	       state = AFTER_NEWLINE;
	       continue; }
	  else if (state == IN_SEQUENCE) {
	       processTheChar; }
	  else if (state == AFTER_NEWLINE) {
	       if (kUnitigSpace[i64] == '>') {
		    ++i64;
		    kUnitigNumber = atoi (kUnitigSpace + i64);
		    while (kUnitigSpace[i64] != '\n')
			 ++i64;
		    state = AFTER_HEADER_NEWLINE;
		    continue; }
	       else { // It's after a regular newline
		    processTheChar; } }
	  else { // state == AFTER_HEADER_NEWLINE
	       j64 = i64;
	       kUnitigSeq[kUnitigNumber] = kUnitigSpace + i64;
	       processTheChar; }
     }

     for (i=0; i<=lastKUnitigNumber; i++)
	  if (kUnitigSeq[i] != NULL)
	       kUnitigLengths[i] = strlen (kUnitigSeq[i]);

#if 0
     for (i=0; i<=lastKUnitigNumber; i++)
	  if (kUnitigSeq[i] != NULL)
	       printf ("i = %d ; len = %d ; seq = %s\n", i, kUnitigLengths[i], kUnitigSeq[i]);
#endif

     FILE *outfile;
     if (strlen (goodSequenceOutputFilename) != 0)
	  outfile = Fopen (goodSequenceOutputFilename, "w");
     else
	  outfile = stdout;
     // (V) above
     strcpy (fname, superReadListFile);
     infile = Fopen (fname, "r");
     if (! isJumpLibrary)
	  goodFile = Fopen (goodFilename, "w");
     std::vector<std::string> emptyStringVector;
     char *outputSeqPtr;
     while (fgets (line, MAX_READ_LEN, infile)) {
	  int readStart=0, readOri='F';
	  char *readName = (char *)0;
	  getFldsFromLine (line, flds);
	  superReadName = flds[1];
	  cptr = superReadName;
	  if (isJumpLibrary) {
	       readStart = atoi(flds[2]);
	       if (readStart < 0)
		    readStart = 0;
	       readOri = *(flds[3]);
	       readName = flds[0]; }
	  else {
	       numReads = atoi (flds[0]);
	       if (numReads < minReadsInSuperRead) {
		    if (numReads == 1)
			 strcpy (pluralStr, "");
		    else
			 strcpy (pluralStr, "s");
		    sprintf (errorMessageLine, "%s has only %d insert%s, which is less than %d. Skipping.\n", superReadName, numReads, pluralStr, minReadsInSuperRead);
		    fputs (errorMessageLine, stderr);
		    continue;
	       }
	  }
	  kUnitigNumber = atoi (cptr);
	  while (isdigit (*cptr)) ++cptr;
	  ori = *cptr;
	  kUnitigNumberHold = kUnitigNumber;
	  oriHold = ori;
	  ++cptr;
	  if (ori == 'F')
	       strcpy (outputSeqSpace, kUnitigSeq[kUnitigNumber]);
	  else {
	       generateReverseComplement (kUnitigSeq[kUnitigNumber], kUnitigLengths[kUnitigNumber]);
	       strcpy (outputSeqSpace, reverseComplementSpace); }
	  fail = 0;
	  errorMessage[0] = 0;
	  while (1) {
	       char *cptr3, *cptr4;
	       int numdiffs;
	       if (*cptr != '_')
		    break;
	       ++cptr;
	       kUnitigNumber = atoi (cptr);
	       while (isdigit (*cptr)) ++cptr;
	       ori = *cptr;
	       ++cptr;
	       if (ori == 'F') {
		    cptr3 = kUnitigSeq[kUnitigNumber];
		    cptr2 = kUnitigSeq[kUnitigNumber] + overlap; }
	       else {
		    generateReverseComplement (kUnitigSeq[kUnitigNumber], kUnitigLengths[kUnitigNumber]);
		    cptr3 = reverseComplementSpace;
		    cptr2 = reverseComplementSpace + overlap; }
	       cptr4 = outputSeqSpace + (strlen(outputSeqSpace)-overlap);
	       numdiffs = 0;
	       for (i=0; i<overlap; i++)
		    if (cptr3[i] != cptr4[i])
			 ++numdiffs;
	       if (numdiffs > seqDiffMax) {
		    fail = 1;
		    if (numdiffs == 1)
			 strcpy (pluralStr, "");
		    else 
			 strcpy (pluralStr, "s");
		    sprintf (errorMessageLine, "The %d-base overlap between %d%c and %d%c has %d difference%s.\n", overlap, kUnitigNumberHold, oriHold, kUnitigNumber, ori, numdiffs, pluralStr);
                    strcat (errorMessage, errorMessageLine); 
		    }
	       strcat (outputSeqSpace, cptr2);
	       kUnitigNumberHold = kUnitigNumber;
	       oriHold = ori; }
	  outputSeqLen = strlen (outputSeqSpace);
#if 0
	  printf ("superReadName = %s ; seq = %s\n", superReadName, outputSeqSpace);
#endif
	  if (fail) {
	       if (numReads == 1)
		    strcpy (pluralStr, "");
	       else
		    strcpy (pluralStr, "s");
	       if (isJumpLibrary)
		    sprintf (errorMessageLine, "super-read %s for read %s fails\n", (char *)superReadName, readName);
	       else
		    sprintf (errorMessageLine, "%s (with %d read%s) fails\n", (char *)superReadName, numReads, pluralStr);
	       fputs (errorMessageLine, stderr);
	       fputs (errorMessage, stderr); 
	       continue; }
	  if (isJumpLibrary) {
	       int amtOfSequenceToOutput;
	       if (readOri == 'F') {
		    amtOfSequenceToOutput = outputSeqLen - readStart;
		    if (amtOfSequenceToOutput < minReadLength) {
			 readStart = outputSeqLen - minReadLength;
			 if (readStart < 0)
			      readStart = 0;
			 amtOfSequenceToOutput = outputSeqLen - readStart; }
	       }
	       else {
		    if (readStart < minReadLength)
			 readStart = minReadLength;
		    if (readStart > outputSeqLen)
			 readStart = outputSeqLen;
		    amtOfSequenceToOutput = readStart;
	       }
	       if (amtOfSequenceToOutput < minReadLength) {
//		    sprintf (errorMessageLine, "super-read for read %s is too short\n", readName);
//		    fputs (errorMessageLine, stderr);
		    continue; }
	       if (amtOfSequenceToOutput > maxReadLength)
		    amtOfSequenceToOutput = maxReadLength;
	       if (readOri == 'F') {
		    outputSeqSpace[readStart+amtOfSequenceToOutput] = 0;
		    outputSeqPtr = outputSeqSpace+readStart;
	       }
	       else {
		    outputSeqSpace[readStart] = 0;
		    outputSeqPtr = outputSeqSpace + (readStart - amtOfSequenceToOutput);
		    reverseTheString(outputSeqPtr); }
	       fputc ('>', outfile);
	       fputs (readName, outfile); fputc('\n', outfile);
	       fputs (outputSeqPtr, outfile); fputc ('\n', outfile);
	       continue;
	  }
			 
	  // If we get here it's a regular super-reads run,
	  // not a jumping library
	  uint64_t tempVar;
	  if ((strlen(superReadNameAndLengthsFilename) != 0) ||
	      renameSuperReads) {
	       std::string superReadNameString = std::string (superReadName);
	       int superReadLen = strlen (outputSeqSpace);
	       if (maxSuperReadLength < superReadLen) {
		    maxSuperReadLength = superReadLen;
		    maxSuperReadNumbersUsed[superReadLen] = 0; }
	       tempVar = 1000000000;
	       tempVar *= superReadLen;
	       tempVar += maxSuperReadNumbersUsed[superReadLen];
	       oldSuperReadName[tempVar] = superReadNameString;
	       ++maxSuperReadNumbersUsed[superReadLen]; }

	  if (! renameSuperReads) {
	       if (! noSequence)
		    fputc ('>', outfile);
	       fputs (superReadName, outfile); fputc ('\n', outfile);
	       fputs (superReadName,goodFile); fputc ('\n',goodFile);
	       if (! noSequence) {
		    fputs (outputSeqSpace, outfile); fputc ('\n', outfile); }
	  }
	  else {
	       superReadSequences[tempVar] = std::string (outputSeqSpace);
	  }
     }

     fclose (infile);

     if (renameSuperReads && (! isJumpLibrary)) {
	  int fauxSuperReadNumber = 0;
	  for (int i2=maxSuperReadLength; i2>=0; i2--) {
	       for (int j=0; j<maxSuperReadNumbersUsed[i2]; j++) {
		    uint64_t tempSuperReadName = 1000000000;
		    tempSuperReadName *= i2;
		    tempSuperReadName += j;
		    if (! noSequence)
			 fputc('>', outfile);
		    fprintf (outfile, "%d\n", fauxSuperReadNumber);
		    fputs (oldSuperReadName[tempSuperReadName].c_str(), goodFile);
		    fputc ('\n', goodFile);
		    if (! noSequence) {
			 fputs (superReadSequences[tempSuperReadName].c_str(), outfile);
			 fputc ('\n', outfile); }
		    ++fauxSuperReadNumber;
	       }
	  }
     }

     if (! isJumpLibrary)
	  fclose (goodFile);
     if (strlen(goodSequenceOutputFilename) != 0)
	  fclose (outfile);

     if (strlen(superReadNameAndLengthsFilename) != 0) {
	  outfile = Fopen (superReadNameAndLengthsFilename, "w");

	  for (int i2=maxSuperReadLength; i2>=0; i2--) {
	       for (int j=0; j<maxSuperReadNumbersUsed[i2]; j++) {
		    uint64_t tempSuperReadName = 1000000000;
		    tempSuperReadName *= i2;
		    tempSuperReadName += j;
		    fputs (oldSuperReadName[tempSuperReadName].c_str(), outfile);
		    fputc (' ', outfile);
		    fprintf (outfile, "%d\n", i2);
	       }
	  }
	  fclose (outfile);
     }
     
     return (0);
}

void reverseTheString (char *str)
{
     int len = strlen(str);
//     printf ("Begin reverseTheString, string is %s\n", str);
     for (int i=0, j=len-1; i<=j; i++, j--) {
	  char tmp = str[i];
	  switch (str[j]) {
	  case 'a': str[i] = 't'; break;
	  case 'c': str[i] = 'g'; break;
	  case 'g': str[i] = 'c'; break;
	  case 't': str[i] = 'a'; break;
	  case 'A': str[i] = 'T'; break;
	  case 'C': str[i] = 'G'; break;
	  case 'G': str[i] = 'C'; break;
	  case 'T': str[i] = 'A'; break;
	  default: str[i] = str[j]; break;
	  }
	  if (i == j)
	       break;
	  switch (tmp) {
	  case 'a': str[j] = 't'; break;
	  case 'c': str[j] = 'g'; break;
	  case 'g': str[j] = 'c'; break;
	  case 't': str[j] = 'a'; break;
	  case 'A': str[j] = 'T'; break;
	  case 'C': str[j] = 'G'; break;
	  case 'G': str[j] = 'C'; break;
	  case 'T': str[j] = 'A'; break;
	  default: str[j] = tmp; break;
	  }
     }
//     printf ("After reverseTheString, string is %s\n", str);
}

void generateReverseComplement (char *seq, int seqLen)
{
     int j, k;
     for (j=0, k=seqLen-1; k>=0; ++j, --k) {
	  switch (seq[j]) {
	  case 'a': reverseComplementSpace[k] = 't'; break;
	  case 'A': reverseComplementSpace[k] = 'T'; break;
	  case 'c': reverseComplementSpace[k] = 'g'; break;
	  case 'C': reverseComplementSpace[k] = 'G'; break;
	  case 'g': reverseComplementSpace[k] = 'c'; break;
	  case 'G': reverseComplementSpace[k] = 'C'; break;
	  case 't': reverseComplementSpace[k] = 'a'; break;
	  case 'T': reverseComplementSpace[k] = 'A'; break;
	  default: reverseComplementSpace[k] = seq[j]; break; }
     }
     reverseComplementSpace[seqLen] = 0;
}

FILE *Fopen (const char *fn, const char *mode)
{
     FILE *result;
     result = fopen (fn, mode);
     if (result == NULL)
     {
          fprintf (stderr, "Couldn't open file '%s' for ", fn);
          switch (mode[0]) {
          case 'r': fprintf (stderr, "reading"); break;
          case 'w': fprintf (stderr, "writing"); break;
          case 'a': fprintf (stderr, "appending"); break;
          default: fprintf (stderr, "unknown operation code '%c'", mode[0]);
               break;
          }
          fprintf (stderr, ". Bye!\n");
          exit (-1);
     }

     return (result);
}

