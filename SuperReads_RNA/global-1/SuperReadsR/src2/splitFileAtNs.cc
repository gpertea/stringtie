/* Give this program a scaffold sequence file and it will generate a contig sequence
   file, where the contigs are regions of the scaffold without (>= numNsToBreakScaffold)
   N's together. It also modifies the scaffold name. The auxiliary files generated are:
   scaffNameTranslations.txt: The correspondence between input and output scaffold names
   genome.posmap.ctgscf: A genome.posmap.ctgscf file relating contigs and scaffs
   genome.asm: A faux but syntactically correct asm file giving seps and stdevs for gaps
  
   The arguments are as follows:
   (1) Input filename
   (2) (Optional) The minimum number of N's that force a break in contigs (see DEFAULT_NS_TO_BREAK_SCAFFOLD)
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <vector>
#include <misc.hpp>
#include <charb.hpp>

void outputAsmFile (FILE *outfile, char *scaff, char *scaffHold, char *contigName, char *contigNameHold, int gapSize);

const double STDEV_DECIMAL = 0.2;
const int MIN_STD_SIZE = 100;
const int DEFAULT_NS_TO_BREAK_SCAFFOLD = 2;

int main (int argc, char **argv)
{
     char *fn = argv[1];
     int numNsToBreakScaffold;
     if (argc <= 2)
	  numNsToBreakScaffold = DEFAULT_NS_TO_BREAK_SCAFFOLD;
     else
	  numNsToBreakScaffold = atoi (argv[2]);
     if (numNsToBreakScaffold <= 0)
	  numNsToBreakScaffold = 1; // We may want to set this to DEFAULT_NS_TO_BREAK_SCAFFOLD
     int contigNum = 0;
     FILE *infile = fopen (fn, "r");
     FILE *outfile = fopen ("genome.posmap.ctgscf", "w");
     FILE *outfile2 = fopen ("genome.asm", "w");
     FILE *scaffNameTranslationFile = fopen ("scaffNameTranslations.txt", "w");
     int scaffoldNumber = 0;
     bool isFirstLine = true;
     charb line(100);
     std::vector<char *> flds;
     int numBadInARow = -1;
     charb badSequence(50);
     charb scaff(30), scaffHold, scaffName;
     charb contigName(30), contigNameHold(30);
     int curScaffOffset = 0;
     int beginContigOffsetInScaff = 0, endContigOffsetInScaff = 0, endContigOffsetInScaffHold = 0;
     bool contigIsStarted = true;
     strcpy (scaffHold, ""); strcpy (contigNameHold, "");
     while (fgets (line, 100, infile)) {
	  if (line[0] == '>') {
	       if (! isFirstLine) {
		    fputc ('\n', stdout);
		    ++endContigOffsetInScaff;
		    fprintf (outfile, "%s %s %d %d f\n", (char *)(contigName+3), (char *)(scaff+3), beginContigOffsetInScaff, endContigOffsetInScaff);
		    // Output .asm file
		    outputAsmFile (outfile2, (char *)scaff, (char *)scaffHold, (char *)contigName, (char *)contigNameHold, beginContigOffsetInScaff - endContigOffsetInScaffHold);
	       }
	       else
		    isFirstLine = false;
	       // set to create .asm file
	       strcpy (contigNameHold, contigName);
	       endContigOffsetInScaffHold = endContigOffsetInScaff;
	       strcpy (scaffHold, scaff);

	       contigIsStarted = false;
	       getFldsFromLine (line, flds);
	       char *cptr = flds[0]+1;
	       if (*cptr)
		    strcpy (scaff, cptr);
	       else
		    strcpy (scaff, "scaff");
	       // Output scaff and new scaff name
	       fprintf (scaffNameTranslationFile, "%s ", (char *) scaff);
	       sprintf (scaff, "jcf7190%09d", scaffoldNumber);
	       fprintf (scaffNameTranslationFile, "%s\n", (char *) scaff);
	       ++scaffoldNumber;
	       sprintf (contigName, "ctg7180%09d", contigNum);
	       printf (">%s\n", (char *)contigName);
	       ++contigNum;
	       numBadInARow = -1; // Negative until the first good base
	       curScaffOffset = -1;
	       continue; }
	  for (char *cptr = line; *cptr; ++cptr) {
	       if (! isspace (*cptr))
		    ++curScaffOffset;
	       switch (*cptr) {
	       case 'A': case 'a':
	       case 'C': case 'c':
	       case 'G': case 'g':
	       case 'T': case 't':
		    if (! contigIsStarted) {
			 beginContigOffsetInScaff = curScaffOffset;
			 contigIsStarted = true; }
		    if (numBadInARow > 0) {
			 if (numBadInARow >= numNsToBreakScaffold) {
			      fputc ('\n', stdout);
			      ++endContigOffsetInScaff;
			      fprintf (outfile, "%s %s %d %d f\n", (char *)(contigName+3), (char *)(scaff+3), beginContigOffsetInScaff, endContigOffsetInScaff);
			      // Output .asm file
			      outputAsmFile (outfile2, (char *)scaff, (char *)scaffHold, (char *)contigName, (char *)contigNameHold, beginContigOffsetInScaff - endContigOffsetInScaffHold);
			      // .. and re-set
			      strcpy (contigNameHold, contigName);
			      endContigOffsetInScaffHold = endContigOffsetInScaff;
			      strcpy (scaffHold, scaff);
			      
			      sprintf (contigName, "ctg7180%09d", contigNum);
			      printf (">%s\n", (char *)contigName);
			      beginContigOffsetInScaff = curScaffOffset;
			      ++contigNum; }
			 else 
			      fputs ((char *)badSequence, stdout);
			 badSequence.clear(); }
		    endContigOffsetInScaff = curScaffOffset;
		    numBadInARow = 0;
		    fputc (*cptr, stdout);
		    break;
	       default:
		    if (isspace (*cptr))
			 break;
		    if (numBadInARow < 0)
			 break;
		    ++numBadInARow;
		    if (numBadInARow < numNsToBreakScaffold)
			 strcat (badSequence, "N"); // Replacing any form of bad sequence with N's
		    break;
	       }
	  } // Matches for (char *cptr = line...
     } // Matches fgets line

     if (! isFirstLine) {
	  fputc ('\n', stdout);
	  ++endContigOffsetInScaff;
	  fprintf (outfile, "%s %s %d %d f\n", (char *)(contigName+3), (char *)(scaff+3), beginContigOffsetInScaff, endContigOffsetInScaff);
	  // Output .asm file
	  outputAsmFile (outfile2, (char *)scaff, (char *)scaffHold, (char *)contigName, (char *)contigNameHold, beginContigOffsetInScaff - endContigOffsetInScaffHold);
	  // .. and re-set
	  fprintf (outfile2, "}\n");
     }

     fclose (infile);
     fclose (outfile);
     fclose (outfile2);
     fclose (scaffNameTranslationFile);

     return (0);
}

void outputAsmFile (FILE *outfile, char *scaff, char *scaffHold, char *contigName, char *contigNameHold, int gapSize)
{
     if (strlen (scaffHold) == 0)
	  fprintf (outfile, "{SCF\nacc:(%s,0)\n", scaff+3);
     else {
	  if (strcmp (scaff, scaffHold) == 0) {
	       int std = STDEV_DECIMAL * gapSize + .50000001;
	       if (std < MIN_STD_SIZE)
		    std = MIN_STD_SIZE;
	       fprintf (outfile, "{CTP\nct1:%s\nct2:%s\nmea:%d\nstd:%d\nori:N\n}\n", contigNameHold+3, contigName+3, gapSize, std); }
	  else
	       fprintf (outfile, "}\n{SCF\nacc:(%s,0)\n", (char *)(scaff+3));
     }
}
