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
  The program outputs the coords and overlaps files for a given set of
  k-unitigs that overlap by exactly (k-1) bases. We work, by default,
  with k-unitigs generated using a k-mer size of 31, but if another
  value was used, use the flag
  -kmervalue kmerLen
  to specify the k-mer length used when generating the k-unitigs.

  The first non-flag arg is the prefix used for the k-unitigs file. (But see
  NOTE1 below.)
  The first flag assumes that the files are named *_#.fa, where
  * is the prefix specified and the #s start from 0 and continue until
  the last one.
  NOTE1: You may also use this arg to specify a filename. The filename must
  have the k-unitigs in numeric order.

  The second non-flag arg is the prefix used for the output files.
  The program will generate the files
  prefix.coords and prefix.overlaps.

  So the final syntax is
  createKUnitigMaxOverlaps [flags] inputPrefix outputPrefix OR
  createKUnitigMaxOverlaps [flags] inputFile   outputPrefix
  where the possible flags are
  -h: help and exit
  -kmervalue kmerLen
  -largest-kunitig-number largestKUnitigNumber

  NOTE2: This program assumes the k-unitig numbers are in increasing order
  in each file (unless the -largest-kunitig-number flag is used).
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <cassert>
#include <string>
#include <vector>

#include <jellyfish/err.hpp>
#include <charb.hpp>
#include <jellyfish/mer_dna.hpp>

#define KMER_LENGTH 31
#define EST_OVLS_PER_KUNITIG 5

namespace err = jellyfish::err;

using jellyfish::mer_dna;
class endKUnitigKmerStruct {
public:
     mer_dna kMerValue;  // (kmerLen as an arg?)
     int kUnitigNumber;
     unsigned char kUnitigEnd; // 0 or 1
     unsigned char ori; // 0 or 1
     endKUnitigKmerStruct (unsigned int k) : kMerValue(k) { }
} **ptrsToEndKUnitigKmerStructs;

struct overlapDataStruct
{
     int kUni1;
     int kUni2; 
     int ahg;
     int bhg;
     char netOri;
};

ExpandingBuffer< struct overlapDataStruct > overlapData(1000);

charb line(2000);
char **kUnitigSequences;
unsigned char **kUnitigSequenceCounts;
unsigned char *endIsDone;
int *kUnitigLengths;
uint64_t largestKUnitigNumber;
char *inputPrefix, *outputPrefix;
uint64_t *startOverlapByUnitig;
uint64_t *overlapDataToSave;
bool createCoordsFile;

void reportKUnitigEndMatches (void);
int kmerStructCompare (const endKUnitigKmerStruct **ptr1, const endKUnitigKmerStruct **ptr2);
void loadKUnitigEndingKMerValues (void);
int getLargestKUnitigNumber (char *prefix, int numInputFiles);
int getNumInputFiles (char *prefix);
void loadKUnitigSequences (char *prefix, int numInputFiles);
FILE *Fopen (const char *fn, const char *mode);
void giveUsageAndExit (void);
void processArgs (int argc, char **argv);

#define mallocOrDie(name, num, type) name = (type *) calloc (num, sizeof ( type )); \
if (name == NULL) { fprintf (stderr, "Couldn't allocate space for '%s'\nBye!\n", #name ); exit (-1); }

int main (int argc, char *argv[])
{
     int numInputFiles;
     uint64_t llval;

     processArgs (argc, argv);
     numInputFiles = getNumInputFiles (inputPrefix);
     if (numInputFiles >= 1)
	  printf ("num input files = %d\n", numInputFiles);
     else
	  printf ("input file = %s\n", inputPrefix);

     if (largestKUnitigNumber == 0)
	  largestKUnitigNumber = getLargestKUnitigNumber (inputPrefix, numInputFiles);
     printf ("largestKUnitigNumber = %llu\n", (long long unsigned int) largestKUnitigNumber);

     mallocOrDie (kUnitigSequences, largestKUnitigNumber+1, char *);
     mallocOrDie (kUnitigLengths, largestKUnitigNumber+1, int);
     mallocOrDie (startOverlapByUnitig, largestKUnitigNumber+2, uint64_t);
     loadKUnitigSequences (inputPrefix, numInputFiles);

     llval = largestKUnitigNumber+1; llval *= 4;
     mallocOrDie (ptrsToEndKUnitigKmerStructs, llval, endKUnitigKmerStruct *);
     loadKUnitigEndingKMerValues ();
     qsort (ptrsToEndKUnitigKmerStructs, 4*(largestKUnitigNumber+1), sizeof (endKUnitigKmerStruct *), (int (*)(const void*, const void*)) kmerStructCompare);

     mallocOrDie (endIsDone, largestKUnitigNumber+1, unsigned char);

     reportKUnitigEndMatches ();

     return (0);
}

void reportKUnitigEndMatches (void)
{
     uint64_t beginIndex, endIndex;
     endKUnitigKmerStruct *ptr1, *ptr2;
     uint64_t i, j;
     int totKUniStartSep;
     int kUni1, kUni2, ahg, bhg, begin1, end1, begin2, end2;
     int isGoodOverlap, skipThis;
     uint64_t numOvlsOutput=0;
     char netOri;
     char filename[500];
     FILE *coordsFile=NULL, *overlapsFile, *selfOverlapsFile;

     if (createCoordsFile) {
	  sprintf (filename, "%s.coords", outputPrefix);
	  coordsFile = Fopen (filename, "w"); }
     sprintf (filename, "%s.overlaps", outputPrefix);
     overlapsFile = Fopen (filename, "wb");
     sprintf (filename, "%s.selfOverlaps.txt", outputPrefix);
     selfOverlapsFile = Fopen (filename, "w");

     beginIndex=0;
     while (beginIndex<4*(largestKUnitigNumber+1)) {
	  ptr1 = ptrsToEndKUnitigKmerStructs[beginIndex];
	  for (endIndex=beginIndex+1; (endIndex<4*(largestKUnitigNumber+1)) && (ptrsToEndKUnitigKmerStructs[endIndex]->kMerValue == ptr1->kMerValue); endIndex++) {
	  }
	  skipThis = 0;
	  if (ptrsToEndKUnitigKmerStructs[beginIndex]->kUnitigEnd == 0) {
	       if (endIsDone[ptrsToEndKUnitigKmerStructs[beginIndex]->kUnitigNumber] >= 2)
		    skipThis = 1;
	  }
	  else {
	       if (endIsDone[ptrsToEndKUnitigKmerStructs[beginIndex]->kUnitigNumber] % 2 == 1)
		    skipThis = 1;
	  }
	  if (skipThis)
	       goto endOfLoop;
	  for (i=beginIndex; i<endIndex; i++) {
	       ptr1 = ptrsToEndKUnitigKmerStructs[i];
	       kUni1 = ptr1->kUnitigNumber;
	       if (ptr1->kUnitigEnd == 0)
		    endIsDone[kUni1] += 2;
	       else
		    endIsDone[kUni1] += 1;
	       if (kUnitigLengths[kUni1] == 0)
		    continue;
	       
	       for (j=beginIndex; j<endIndex; j++) {
		    if (i == j)
			 continue;
		    ptr2 = ptrsToEndKUnitigKmerStructs[j];
		    kUni2 = ptr2->kUnitigNumber;
		    if (kUnitigLengths[kUni2] == 0)
			 continue;
		    if (ptr1->ori == ptr2->ori)
			 netOri = 'N';
		    else
			 netOri = 'I';
		    if ((ptr1->kUnitigEnd + ptr2->kUnitigEnd + ptr1->ori + ptr2->ori) % 2 == 0)
			 isGoodOverlap = 0;
		    else
			 isGoodOverlap = 1;
		    if (isGoodOverlap) {
			 if (ptr1->kUnitigEnd == 0) {
			      ahg = (mer_dna::k()) - kUnitigLengths[kUni2];
			      bhg = (mer_dna::k()) - kUnitigLengths[kUni1]; }
			 else {
			      ahg = kUnitigLengths[kUni1] - (mer_dna::k());
			      bhg = kUnitigLengths[kUni2] - (mer_dna::k()); }
		    }
		    else {
			 if (ptr1->kUnitigEnd == 0) {
			      ahg = 0;
			      bhg = kUnitigLengths[kUni2] - kUnitigLengths[kUni1]; }
			 else {
			      ahg = kUnitigLengths[kUni1] - kUnitigLengths[kUni2];
			      bhg = 0; }
		    }
		    if (ptr1->kUnitigEnd == 0) {
			 begin1 = 1;
			 end1 = mer_dna::k(); }
		    else {
			 end1 = kUnitigLengths[kUni1];
			 begin1 = end1 - (mer_dna::k()) + 1; }
		    if (netOri == 'N') {
			 begin2 = begin1 - ahg;
			 end2 = end1 - ahg; }
		    else {
			 totKUniStartSep = kUnitigLengths[kUni1] + bhg;
			 begin2 = totKUniStartSep - (begin1-1);
			 end2 = (totKUniStartSep - end1) + 1;
		    }
		    if (createCoordsFile)
			 fprintf (coordsFile, "%d %d %d %d 100.00 %d %d %d %d\n", begin1, end1, begin2, end2, kUnitigLengths[kUni1], kUnitigLengths[kUni2], kUni1, kUni2);
		    if (isGoodOverlap) {
			 if(kUni1 != kUni2){   //temporary dirty fix by Aleksey
			      overlapData[numOvlsOutput].kUni1 = kUni1;
			      overlapData[numOvlsOutput].kUni2 = kUni2;
			      overlapData[numOvlsOutput].ahg = ahg;
			      overlapData[numOvlsOutput].bhg = bhg;
			      overlapData[numOvlsOutput].netOri = netOri;
			      ++numOvlsOutput; }
			 else {
			      fprintf (selfOverlapsFile, "%d %d %c %d %d 0.0 0.0\n", kUni1, kUni2, netOri, ahg, bhg);
//			      fprintf (overlapsFile, "%d %d %c %d %d 0.0 0.0\n", kUni1, kUni2, netOri, ahg, bhg);
			 }
		    }
	       }
	  }
     endOfLoop:
	  beginIndex = endIndex;
     }
     free(kUnitigLengths);
     free(kUnitigSequences);
     free(ptrsToEndKUnitigKmerStructs);
     //mallocOrDie (overlapDataToSave, numOvlsOutput, struct overlapDataStruct);
     mallocOrDie (overlapDataToSave, numOvlsOutput, uint64_t);	
     for (uint64_t j=0; j<numOvlsOutput; j++)
     {
	  overlapDataToSave[j]=0;
	  int unitig1 = overlapData[j].kUni1, unitig2 = overlapData[j].kUni2;
	  if (unitig1 >= unitig2)
	       continue;
	  startOverlapByUnitig[unitig1]++;
	  startOverlapByUnitig[unitig2]++;
     }
     for (uint64_t unitigNum = 1; unitigNum < largestKUnitigNumber + 2; unitigNum++)
	  startOverlapByUnitig[unitigNum] += startOverlapByUnitig[unitigNum - 1];
     
     for (uint64_t j=0; j<numOvlsOutput; j++)
	  overlapDataToSave[--startOverlapByUnitig[overlapData[j].kUni1]]=j;

     size_t written = 0;
     while(written < numOvlsOutput) {
       size_t res;
       res = fwrite (overlapData+overlapDataToSave[written], sizeof (struct overlapDataStruct),
                            1, overlapsFile);
       if(res == 0)
         throw std::runtime_error(err::msg() << "Failed to write overlaps to file: " << err::no);
       written +=res;
     }
     
     if (createCoordsFile)
	  fclose (coordsFile);
     fclose (overlapsFile);
     fclose (selfOverlapsFile);
}
	  
int kmerStructCompare (const endKUnitigKmerStruct **ptr1, const endKUnitigKmerStruct **ptr2)
{
     if ((*ptr1)->kMerValue < (*ptr2)->kMerValue) return (-1);
     if ((*ptr2)->kMerValue < (*ptr1)->kMerValue) return (1);
     if ((*ptr1)->kUnitigNumber < (*ptr2)->kUnitigNumber) return (-1);
     if ((*ptr1)->kUnitigNumber > (*ptr2)->kUnitigNumber) return (1);
     if ((*ptr1)->kUnitigEnd < (*ptr2)->kUnitigEnd) return (-1);
     if ((*ptr1)->kUnitigEnd > (*ptr2)->kUnitigEnd) return (1);
     if ((*ptr1)->ori < (*ptr2)->ori) return (-1);
     if ((*ptr1)->ori > (*ptr2)->ori) return (1);
     return (0);
}

void loadKUnitigEndingKMerValues (void)
{
     uint64_t kUnitigNumber, index;
     mer_dna maxVal;
     maxVal.polyA();
     maxVal.reverse_complement();
     
     for (kUnitigNumber=0; kUnitigNumber<=largestKUnitigNumber; kUnitigNumber++) {
	  index = 4 * kUnitigNumber;
	  if (kUnitigLengths[kUnitigNumber] == 0) {
	       ptrsToEndKUnitigKmerStructs[index] = new endKUnitigKmerStruct(mer_dna::k());
	       ptrsToEndKUnitigKmerStructs[index+1] = new endKUnitigKmerStruct(mer_dna::k());
	       ptrsToEndKUnitigKmerStructs[index+2] = new endKUnitigKmerStruct(mer_dna::k());
	       ptrsToEndKUnitigKmerStructs[index+3] = new endKUnitigKmerStruct(mer_dna::k());
	       ptrsToEndKUnitigKmerStructs[index]->kMerValue = ptrsToEndKUnitigKmerStructs[index+1]->kMerValue = ptrsToEndKUnitigKmerStructs[index+2]->kMerValue = ptrsToEndKUnitigKmerStructs[index+3]->kMerValue = maxVal;
	       ptrsToEndKUnitigKmerStructs[index]->kUnitigNumber = ptrsToEndKUnitigKmerStructs[index+1]->kUnitigNumber = ptrsToEndKUnitigKmerStructs[index+2]->kUnitigNumber = ptrsToEndKUnitigKmerStructs[index+3]->kUnitigNumber = kUnitigNumber;
	       ptrsToEndKUnitigKmerStructs[index]->kUnitigEnd = ptrsToEndKUnitigKmerStructs[index+1]->kUnitigEnd = 0;
	       ptrsToEndKUnitigKmerStructs[index+2]->kUnitigEnd = ptrsToEndKUnitigKmerStructs[index+3]->kUnitigEnd = 1;
	       ptrsToEndKUnitigKmerStructs[index]->ori = ptrsToEndKUnitigKmerStructs[index+2]->ori = 0;
	       ptrsToEndKUnitigKmerStructs[index+1]->ori = ptrsToEndKUnitigKmerStructs[index+3]->ori = 1;
	       continue;
	  }
	  // k-mer at beginning of k-unitig, forward ori
	  std::string kUnitigSequence = std::string (kUnitigSequences[kUnitigNumber]);
	  ptrsToEndKUnitigKmerStructs[index] = new endKUnitigKmerStruct(mer_dna::k());
	  ptrsToEndKUnitigKmerStructs[index]->kMerValue = kUnitigSequence.substr(0, mer_dna::k());
	  ptrsToEndKUnitigKmerStructs[index]->kUnitigNumber = kUnitigNumber;
	  ptrsToEndKUnitigKmerStructs[index]->kUnitigEnd = 0;
	  ptrsToEndKUnitigKmerStructs[index]->ori = 0;
	  ++index;
	  // k-mer at beginning of k-unitig, reverse ori
	  ptrsToEndKUnitigKmerStructs[index] = new endKUnitigKmerStruct(mer_dna::k());
	  ptrsToEndKUnitigKmerStructs[index]->kMerValue = (ptrsToEndKUnitigKmerStructs[index-1]->kMerValue).get_reverse_complement();
	  ptrsToEndKUnitigKmerStructs[index]->kUnitigNumber = kUnitigNumber;
	  ptrsToEndKUnitigKmerStructs[index]->kUnitigEnd = 0;
	  ptrsToEndKUnitigKmerStructs[index]->ori = 1;
	  ++index;
	  // k-mer at end of k-unitig, forward ori
	  ptrsToEndKUnitigKmerStructs[index] = new endKUnitigKmerStruct(mer_dna::k());
	  ptrsToEndKUnitigKmerStructs[index]->kMerValue = kUnitigSequence.substr(kUnitigLengths[kUnitigNumber]-mer_dna::k(), mer_dna::k());
	  ptrsToEndKUnitigKmerStructs[index]->kUnitigNumber = kUnitigNumber;
	  ptrsToEndKUnitigKmerStructs[index]->kUnitigEnd = 1;
	  ptrsToEndKUnitigKmerStructs[index]->ori = 0;
	  ++index;
	  // k-mer at end of k-unitig, reverse ori
	  ptrsToEndKUnitigKmerStructs[index] = new endKUnitigKmerStruct(mer_dna::k());
	  ptrsToEndKUnitigKmerStructs[index]->kMerValue = (ptrsToEndKUnitigKmerStructs[index-1]->kMerValue).get_reverse_complement();
	  ptrsToEndKUnitigKmerStructs[index]->kUnitigNumber = kUnitigNumber;
	  ptrsToEndKUnitigKmerStructs[index]->kUnitigEnd = 1;
	  ptrsToEndKUnitigKmerStructs[index]->ori = 1;
     }
}	       

int getLargestKUnitigNumber (char *prefix, int numInputFiles)
{
     struct stat statbuf;
     FILE *infile;
     int i;
     char fname[500];
     off_t fileOffset;
     unsigned int maxKUnitig = 0;
     unsigned int currentKUnitig;
     int numInputFilesTemp;

     numInputFilesTemp = (numInputFiles == 0) ? 1 : numInputFiles;

     for (i=0; i<numInputFilesTemp; i++) {
	  if (numInputFiles > 0)
	       sprintf (fname, "%s_%d.fa", prefix, i);
	  else
	       strcpy (fname, prefix);
	  stat (fname, &statbuf);
	  fileOffset = statbuf.st_size;
	  if (fileOffset >= 1000000)
	       fileOffset -= 1000000;
	  else
	       fileOffset = 0;
	  infile = Fopen (fname, "r");
	  fseek (infile, fileOffset, SEEK_SET);
	  if(!fgets (line, infile))
            err::die(err::msg() << "Error reading file '" << fname << "'" << err::no);
	  while (fgets (line, 1000000, infile)) {
	       if (line[0] != '>')
		    continue;
	       currentKUnitig = atoi (line+1);
	       if (currentKUnitig > maxKUnitig)
		    maxKUnitig = currentKUnitig;
	  }
	  fclose (infile);
     }
     return (maxKUnitig);
}

int getNumInputFiles (char *prefix)
{
     struct stat statbuf;
     int i;
     char fname[500];
     for (i=0; 1; i++) {
	  sprintf (fname, "%s_%d.fa", prefix, i);
	  if (stat (fname, &statbuf) == 0)
	       continue;
	  return (i);
     }
     return (0);
}

void loadKUnitigSequences (char *prefix, int numInputFiles)
{
     FILE *infile;
     int i, kUnitigNumber, length;
     char fname[500];
     int numInputFilesTemp;
     charb header;

     numInputFilesTemp = (numInputFiles == 0) ? 1 : numInputFiles;

     for (i=0; i<numInputFilesTemp; i++) {
	  if (numInputFiles > 0)
	       sprintf (fname, "%s_%d.fa", prefix, i);
	  else
	       strcpy (fname, prefix);
	  infile = Fopen (fname, "r");
	  if(!infile)
	    err::die(err::msg() << "Can't open file '" << fname << "'" << err::no);
	  
	  int next_char = fgetc(infile);
	  if(next_char != '>')
	    err::die(err::msg() << "Badly formatted fasta file '" << fname << "'. Missing header");

	  while(fgets (header, infile)) {
	    kUnitigNumber = atoi (header);
	    next_char = fgetc(infile);
	    line.clear();
	    while(next_char != EOF && next_char != '>') {
	      ungetc(next_char, infile);
	      fgets_append (line, infile);
	      fflush(stdout);
	      line.chomp();
	      next_char = fgetc(infile);
	    }
	    length = line.len();
	    mallocOrDie (kUnitigSequences[kUnitigNumber], length+1, char);
	    kUnitigLengths[kUnitigNumber] = length;
	    strcpy (kUnitigSequences[kUnitigNumber], line);
	  }
	  fclose (infile);
     }
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

void giveUsageAndExit (void)
{
     fprintf (stderr, 
              "This program outputs the overlaps file as well as (optionally) the coords file for a given set of\n" 
              " k-unitigs generated by k-mers of length K that overlap by exactly\n"
              " (K-1) bases. We work, by default, with k-unitigs generated using a\n"
              " k-mer size of 31, but if another k-mer size was used, use the flag\n"
              "-kmervalue kMerSize\n"
              " to specify the k-mer size used when generating the k-unitigs.\n"
              "\n"
              "The first non-flag arg is the prefix used for the k-unitigs files.\n"
              " It assumes that the files are named *_#.fa, where\n"
              " * is the prefix specified and the #s start from 0 and continue until\n"
              " the last one. This arg may also be used to specify the complete\n"
              " filename. Note that all input files are assumed to have k-unitig\n"
              " numbers in ascending order.\n"
              "\n"
              "The second non-flag arg is the prefix used for the output files.\n"
              " The program will generate the files\n"
              " prefix.coords   and   prefix.overlaps.\n"
              "\n"
              " So the final syntax is\n"
              "\n"
              "createKUnitigMaxOverlaps [flags] inputPrefix outputPrefix\n"
              " where the possible flags are\n"
              "   -h: help and exit\n"
              "   -kmervalue kMerSize\n"
              "   -create-coords-file to output the coords file as well as the overlaps file\n"
              "   -largest-kunitig-number largestKUnitigNumber (in this case the\n"
              "       k-unitigs don't have to be in numeric order in the files.)\n");
     exit (0);
}

void processArgs (int argc, char **argv)
{
     int numArgsSeen, i;
     inputPrefix = outputPrefix = NULL;
     mer_dna::k(KMER_LENGTH-1);
     numArgsSeen = 0;
     largestKUnitigNumber = 0;
     createCoordsFile = false;
     for (i=1; i<argc; i++) {
	  if (strcmp (argv[i], "-h") == 0)
	       giveUsageAndExit();
	  if (strcmp (argv[i], "-kmervalue") == 0) {
	       ++i;
	       mer_dna::k(atoi(argv[i]) - 1);
	       continue; }
	  if (strcmp (argv[i], "-largest-kunitig-number") == 0) {
	       ++i;
	       largestKUnitigNumber = atoi (argv[i]);
	       continue; }
	  if (strcmp (argv[i], "-create-coords-file") == 0) {
	       createCoordsFile = true;
	       continue; }
	  if (argv[i][0] == '-') {
	       fprintf (stderr, "\nUnrecognized flag: %s\n\n", argv[i]);
	       giveUsageAndExit(); }
	  if (numArgsSeen == 0)
	       inputPrefix = argv[i];
	  else if (numArgsSeen == 1)
	       outputPrefix = argv[i];
	  else {
	       fprintf (stderr, "\nToo many args supplied. 2 are expected.\n\n");
	       giveUsageAndExit(); }
	  ++numArgsSeen;
     }
     if (numArgsSeen < 2) {
	  fprintf (stderr, "\nToo few args supplied. 2 are expected.\n\n");
	  giveUsageAndExit();
     }
}

