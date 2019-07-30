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


#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include <charb.hpp>

int main(int argc, char *argv[])
{
  ssize_t bytes_read,i,j;
  charb my_string(100000),rev_sequence(100000);
  
  while(1)
    {
      if(fgets (my_string, 1, stdin)==NULL)
	break;
      else
	{
	  if(my_string[0]=='>')
	    {
	      printf("%s",(char*)my_string);
	    }
	  else
	    {
             bytes_read=strlen(my_string);
	      j=0;
	      for(i=bytes_read-2;i>=0;i--)
		{
		  switch (my_string[i])
		    {
		    case 'A': rev_sequence[j]='T';break;
		    case 'C': rev_sequence[j]='G';break;
		    case 'G': rev_sequence[j]='C';break;
		    case 'T': rev_sequence[j]='A';break;
		    case 'a': rev_sequence[j]='t';break;
		    case 'c': rev_sequence[j]='g';break;
		    case 'g': rev_sequence[j]='c';break;
		    case 't': rev_sequence[j]='a';break;
		    default:  rev_sequence[j]='N';
		    }
		  j++;
		}
              rev_sequence[j]='\0';
	      printf("%s\n",(char*)rev_sequence);
	    }
	}
    }
return(0);
}

