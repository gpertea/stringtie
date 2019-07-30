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


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <set>
#include <charb.hpp>
#include <misc.hpp>
using namespace std;

charb line(1000000), line2;
vector<char *> flds;
set<string> readnames;

FILE *Fopen (const char *fn, const char *mode);

int main(int argc, char **argv)
{
     char *file1=NULL, *file2=NULL;
     int fieldNum=0, fieldNum0=0;
     int argNum = 0;
     for (int i=1; i<argc; i++) {
	  if (strcmp(argv[i], "--fld-num") == 0) {
	       ++i;
	       fieldNum0 = atoi(argv[i]);
	       continue; }
	  if (argNum == 0)
	       file1 = argv[i];
	  else if (argNum == 1)
	       file2 = argv[i];
	  else if (argNum == 2)
	       fieldNum = atoi (argv[i]);
	  else {
	       fprintf (stderr, "There are too many arguments to '%s'. Bye!\n", argv[0]);
	       exit (1); }
	  ++argNum; }
     if (argNum < 3) {
	  fprintf (stderr, "There are only %d args to '%s', 3 must be provided. Bye!\n", argNum, argv[0]);
	  exit (1); }
     FILE *infile = Fopen (file1, "r");
     while (fgets (line, 1000000, infile)) {
	  getFldsFromLine (line, flds);
	  readnames.insert (string(flds[fieldNum0])); }
     fclose (infile);

     infile = Fopen (file2, "r");
     int isOn = 0;
     while (fgets (line, 1000000, infile)) {
	  if (line[0] == '>') {
	       strcpy (line2, line); // getFldsFromLine will overwrite line
	       getFldsFromLine (line2, flds);
	       // Special for the first field of the line
	       ++flds[0];
	       if (readnames.find(string(flds[fieldNum])) == readnames.end())
		    isOn = 1;
	       else
		    isOn = 0; }
	  if (isOn)
	       fputs (line, stdout);
     }
     fclose (infile);

     return (0);
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

