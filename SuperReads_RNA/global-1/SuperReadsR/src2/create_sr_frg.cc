// Replaces create_sr_frg.pl
#include <cstdio>
#include <cstring>
#include <vector>
#include <string>
#include <charb.hpp>
#include <misc.hpp>
void reverse_complement (charb &inSeq, charb &outSeq);

int main (int argc, char **argv)
{
     charb readName(200);
     int shootingIndex = 0;
     charb line(1000), seq(1000), revCompSeq(1000);
     std::vector<char *> flds;
    
     while (fgets (line, 1000, stdin)) {
//	  printf ("line = %s\n", (char *)line);
	  if (line[0] == '>') {
//	       printf ("readNameLen = %d\n", readName.len());
	       if (readName.len() > 0) {
		    int len = seq.len();
#if 0
		    reverse_complement (seq, revCompSeq); 
		    int oriCode = 0;
		    if (strcmp (revCompSeq, seq) < 0) {
			 seq = revCompSeq;
			 oriCode = 1; }
#endif
		    if (len < 2048)
			 printf ("%s\n%s\n", (char *) readName, (char *)seq);
		    else {
			 int k = 0;
			 int offset = 1550;
			 char bufferHold[2048];
			 while (1) {
			      bool atEndOfString = false;
			      printf ("%s.%d\n", (char *) readName, k * offset);
			      if (len-k*offset >= 2047) {
				   strncpy (bufferHold, ((char *)seq)+(k*offset), 2047);
				   bufferHold[2047] = 0; }
			      else {
				   strcpy (bufferHold, ((char *)seq) + (k*offset));
				   atEndOfString = true; }
			      ++k;
			      printf ("%s\n", bufferHold);
			      if (atEndOfString)
				   break; }
		    }
	       }
	       getFldsFromLine (line, flds);
	       if (strlen (flds[0]) > 100) {
		    fprintf (stderr, "SR%d %s\n", shootingIndex, flds[0]+1);
		    sprintf (readName, "SR%d:super-read", shootingIndex);
		    ++shootingIndex; }
	       else
		    sprintf (readName, "%s:super-read", flds[0]);
	       seq.clear(); }
	  else {
	       line.chomp();
	       strcat (seq, line);
	  }
     }
     if (readName.len() > 0) {
	  int len = seq.len();
#if 0
	  reverse_complement (seq, revCompSeq); 
	  int oriCode = 0;
	  if (strcmp (revCompSeq, seq) < 0) {
	       seq = revCompSeq;
	       oriCode = 1; }
#endif
	  if (len < 2048)
	       printf ("%s\n%s\n", (char *)readName, (char *)seq);
	  else {
	       int k = 0;
	       int offset = 1550;
	       char bufferHold[2048];
	       while (1) {
		    bool atEndOfString = false;
		    printf ("%s.%d\n", (char *)readName, k * offset);
		    if (len-k*offset >= 2047) {
			 strncpy (bufferHold, ((char *)seq)+(k*offset), 2047);
			 bufferHold[2047] = 0; }
		    else {
			 strcpy (bufferHold, ((char *)seq) + (k*offset));
			 atEndOfString = true; }
		    ++k;
		    printf ("%s\n", bufferHold);
		    if (atEndOfString)
			 break; }
	  }
     }
     
     return (0);
}

void reverse_complement (charb &inSeq, charb &outSeq)
{
     outSeq = inSeq;
     for (int i=inSeq.len()-1, j=0; i>=0; --i, ++j) {
	  switch (inSeq[i]) {
	  case 'a': outSeq[j] = 't'; break;
	  case 'A': outSeq[j] = 'T'; break;
	  case 'c': outSeq[j] = 'g'; break;
	  case 'C': outSeq[j] = 'G'; break;
	  case 'g': outSeq[j] = 'c'; break;
	  case 'G': outSeq[j] = 'C'; break;
	  case 't': outSeq[j] = 'a'; break;
	  case 'T': outSeq[j] = 'A'; break;
	  default:
	       outSeq[j] = inSeq[i]; break; }
     }
}

	       
