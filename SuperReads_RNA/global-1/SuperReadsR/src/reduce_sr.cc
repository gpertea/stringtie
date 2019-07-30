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


#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<algorithm>
#include<charb.hpp>
#include<src/reduce_sr_cmdline.hpp>
// #define DEBUG 0

typedef ExpandingBuffer<int> int_buf;
typedef ExpandingBuffer<int_buf, remaper<int_buf> > twoD_int_buf;
typedef ExpandingBuffer<charb, remaper<charb> > twoD_charb_buf;


FILE *Fopen (const char *fn, const char *mode);

static int int_compare(const void *p1, const void *p2){
     if(*(const int *) p1 < *(const int *)p2)
	  return(-1);
     else if(*(const int *) p1 > *(const int *)p2)
	  return(1);
     else
	  return(0);
}


void reverse_sr(const char * original, char * reversed){
     //will write the reversed sr into *reversed. Assume reversed is at
     //least as large as original.
     charb original_local=original;
     char * token, *savespace;
     int offset=strlen(original);
  
     reversed[offset]='\0';

     token=strtok_r(original_local,"_",&savespace);
  
     if(token[strlen(token)-1]=='F')
	  token[strlen(token)-1]='R';
     else
	  token[strlen(token)-1]='F';

     offset-=strlen(token);
     for(size_t i=0;i<strlen(token);i++)
	  reversed[offset+i]=token[i];

     while(1){
	  token=strtok_r(NULL,"_",&savespace);
	  if(token==NULL)
	       break;
	  if(token[strlen(token)-1]=='F')
	       token[strlen(token)-1]='R';
	  else if(token[strlen(token)-1]=='R')
	       token[strlen(token)-1]='F';
	  reversed[--offset]='_';
	  offset-=strlen(token);
	  for(size_t i=0;i<strlen(token);i++)
	       reversed[offset+i]=token[i];
     }
}

// cptr is the candidate string ending at end.
int compute_lengthAdjustment(char ori, char* cptr, char* end, int_buf& kUnitigLengths, int overlapLength) {
  int   res                = 0;

  while (cptr < end) {
    int unitigNumber  = atoi (cptr);
    res              += (kUnitigLengths[unitigNumber] - overlapLength);
    while (cptr < end && isdigit(*cptr)) ++cptr;
    while (cptr < end && ! isdigit(*cptr)) ++cptr;
    }
    if(ori == 'R')
      res += overlapLength;

  return res;
}

int main(int argc,char *argv[]){
     cmdline_parse args(argc, argv);
     int i,j,k,l=0,lastKUnitigIndex,irreducibleSuperReadIndex=0;
     time_t time_start=time(NULL);
     int_buf candidates(100),kUnitigsInSuperRead(1000),kUnitigLengths(args.largestkunitig_arg + 1);
     twoD_charb_buf irreducibleSuperReadNames(10000000);
     charb superReadName_reverse(1000),superReadName(1000),superReadName_save(1000),line(1000000);
     char *token, *saveptr;
     char *tptr;
     char relativeOri='F';
     twoD_int_buf superReadIndicesForKUnitig(args.largestkunitig_arg);
     int overlapLength = args.kmerlen_arg - 1;
     int lengthAdjustment = 0;
#ifndef DEBUG
     FILE* output = Fopen(args.output_arg, "w");
#endif
     FILE *infile = Fopen (args.kunitigLengthsFile_arg, "r");
     FILE *sr_sizes = Fopen (args.SuperReads_sizes_arg, "r");


     while (fgets (line, infile)) {
	  int uniNum, uniLen;
	  sscanf (line, "%d %d", &uniNum, &uniLen);
	  kUnitigLengths[uniNum] = uniLen; }
     fclose (infile);

     while(fgets(line,sr_sizes)){
	  l++;	
	  if(l%500000==0){
	       fprintf(stderr,"Processed %d super reads, irreducible %d, processing %d super reads per second\n",l,irreducibleSuperReadIndex,(int)floor(500000/difftime(time(NULL),time_start)));
	       time_start=time(NULL);
#if DEBUG
		fprintf(stderr,"USAGE:\n");
                int mem=0;
		for (size_t k=0;k<irreducibleSuperReadNames.size();k++)
			mem+=irreducibleSuperReadNames[k].size();
		fprintf(stderr,"irreducibleSuperReadNames takes %d bytes\n",mem);
		mem=0;
                for (size_t k=0;k<args.largestkunitig_arg;k++)
                        mem+=superReadIndicesForKUnitig[k].size()*4;
                fprintf(stderr,"superReadIndicesForKUnitig takes %d bytes\n",mem);

                fprintf(stderr,"CAPACITY:\n");
                mem=0;
                for (size_t k=0;k<irreducibleSuperReadNames.capacity();k++)
                        mem+=irreducibleSuperReadNames[k].capacity();
                fprintf(stderr,"irreducibleSuperReadNames takes %d bytes\n",mem);
                mem=0;
                for (size_t k=0;k<args.largestkunitig_arg;k++)
                        mem+=superReadIndicesForKUnitig[k].capacity()*4;
                fprintf(stderr,"superReadIndicesForKUnitig takes %d bytes\n",mem);
#endif
	  }
	  //first we parse the super read line, space separated, to get the name
	  token=strtok_r(line," ",&saveptr);
	  if(!token) {
	       fprintf(stderr, "Invalid empty line number %d, expected a SuperRead name\n", l);
	       continue;
	  }
#if DEBUG
	  printf("Found super read %s\n",token);
#endif
	  //we need to save the super read name
	  superReadName=token;
	  superReadName_save=token;
	  //here we need to examine first and last k-unitig in the super read name
	  j=0;
	  token=strtok_r(superReadName,"FR_",&saveptr);
	  for(i=0; token != NULL;i++){
	       kUnitigsInSuperRead[j++]=atoi(token);
	       token=strtok_r(NULL,"FR_",&saveptr); 
	  }
	  lastKUnitigIndex=j-1;
#if DEBUG
	  printf("first k_unitig: %d last k_unitig %d, super read %s candidates for first %d candidates for last %d\n",kUnitigsInSuperRead[0],kUnitigsInSuperRead[lastKUnitigIndex],(char*)superReadName_save,(int)superReadIndicesForKUnitig[kUnitigsInSuperRead[0]].size(),(int)superReadIndicesForKUnitig[kUnitigsInSuperRead[lastKUnitigIndex]].size());
#endif



	  //now we try to reduce
	  k = 0;
	  const int_buf& first_srs = superReadIndicesForKUnitig[kUnitigsInSuperRead[0]];
	  const int_buf& last_srs = superReadIndicesForKUnitig[kUnitigsInSuperRead[lastKUnitigIndex]];
	  if(first_srs.size() > 0 && last_srs.size() > 0) {
	       int       max_first_index = 0;
	       const int max_k           = args.maximum_search_depth_arg;
	       const int max_2k          = max_k*2;
      
	       max_first_index = std::max(first_srs[0], last_srs[0]);
      
	       for(i=0; i< (int)first_srs.size() && k < max_k; ++i)
		    if(first_srs[i] >= max_first_index)
			 candidates[k++] = first_srs[i];

	       for(i=0; i< (int)last_srs.size() && k < max_2k; ++i)
		    if(last_srs[i] >= max_first_index)
			 candidates[k++] = last_srs[i];
      
#if DEBUG
	       printf("Found %d candidates\n",k);
#endif

#if DEBUG
	       for(i = 0; i < k; ++i)
		    printf("Before qsort:Candidate %d, super read %s\n",candidates[i], (char*)irreducibleSuperReadNames[candidates[i]]);
#endif

	       qsort(candidates,k,sizeof(int),int_compare);

	       superReadName_reverse=superReadName_save;

	       reverse_sr(superReadName_save,superReadName_reverse);
	       //now we go through the sorted candidates and figure out which one matches
#if DEBUG
	       for(i = 0; i < k; i++)
		    printf("Candidate %d, super read %s\n",candidates[i], (char*)irreducibleSuperReadNames[candidates[i]]);
#endif

	       //now we go through the sorted candidates and figure out which one matches, we only look at the candidate if it is encountered twice
	       int last_candidate  = -1;
	       int candidate_count = 1;
	       for(i = 0; i < k; last_candidate = candidates[i], ++i) {
                 if(candidates[i] != last_candidate || candidate_count != 1) {
                   candidate_count=1;
                   continue;
                 }      
#if DEBUG
                 printf("Checking candidate %s\n",(char*)irreducibleSuperReadNames[candidates[i]]);
#endif
                 tptr = strstr(irreducibleSuperReadNames[candidates[i]],superReadName_save);
                 if(tptr!=NULL) {
                   relativeOri = 'F';
                   lengthAdjustment = compute_lengthAdjustment(relativeOri, irreducibleSuperReadNames[candidates[i]], tptr,
                                                               kUnitigLengths, overlapLength);
                   break;
                 }
                 tptr = strstr(irreducibleSuperReadNames[candidates[i]],superReadName_reverse);
                 if(tptr!=NULL) {
                   relativeOri = 'R';
                   lengthAdjustment = compute_lengthAdjustment(relativeOri, irreducibleSuperReadNames[candidates[i]],
                                                               tptr + strlen(superReadName_reverse),
                                                               kUnitigLengths, overlapLength);
                   break;
                 }
                 candidate_count++;
	       }

#if DEBUG
	       if(i<k) {
		    printf("Reduced %s to %s\n",(char*)superReadName_save,(char*)irreducibleSuperReadNames[candidates[i]]); 	
		    continue;
	       }
#else
	       if(i < k) {
                 fprintf(output, "%s %s %c %d\n",
                         (char*)superReadName_save,(char*)irreducibleSuperReadNames[candidates[i]], relativeOri, lengthAdjustment);
                  
		    continue;
	       }
#endif
	  }
#if DEBUG
	  printf("Irreducible %s, index %d\n",(char*)superReadName_save,irreducibleSuperReadIndex);
#endif

	  //if we got here, then the super read is irreducible :(
	  //here is what we do with an irreducible super read
	  irreducibleSuperReadNames[irreducibleSuperReadIndex]=superReadName_save;
	  for(i = 0; i <= lastKUnitigIndex; ++i)
            superReadIndicesForKUnitig[kUnitigsInSuperRead[i]].push_back(irreducibleSuperReadIndex);
	  irreducibleSuperReadIndex++;
     }
     fclose(sr_sizes);
     return(0);
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

