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
#include <charb.hpp>
int main(int argc, char *argv[])
{
  ssize_t i,sequence_len;
  charb my_string(100),sequence(100);
  ssize_t *index_array;
  ssize_t max_val=0,max_ind=0;
  double num_cg=0;
  char last_letter;
  index_array = (ssize_t*)malloc(sizeof(ssize_t)*1000000);
  
  if(argc==1)
    {
      printf("Need one agrument -- minimum peak height to trim!\n");
      return(1);
    }

  while(1)
    {
      if(!fgets(my_string,stdin)){ 
	break;
	}
      else
	{
	  if(my_string[0]=='>')
	    {
	      max_val=0;
	      max_ind=0;
	      num_cg=0;
              if(!fgets(sequence, stdin))
                {
                return(1);
                }
              sequence_len=strlen(sequence);
//              printf("sequence_len=%d\n",sequence_len);
	      last_letter=sequence[sequence_len-2];
//	      printf("last_letter=%c\n",last_letter);
   	      for(i=sequence_len-3;i>=0;i--)
		    {
		    if(sequence[i] != last_letter)
			{
			num_cg--;
			last_letter=sequence[i];
			}
		    else
			{
			num_cg++;
			}
			index_array[sequence_len-3-i]=num_cg;
		    }
	      for(i=0;i<sequence_len-2;i++)
		{
		if(index_array[i]>max_val)
			{
			max_val=index_array[i];
			max_ind=i+1;
			}
		}
//                printf("max_val=%d max_ind=%d\n",max_val,max_ind);
		if(max_val<(ssize_t)atoi(argv[1]))
			max_ind=0;
		sequence[sequence_len-max_ind-1]='\0';
		printf("%s%s\n",(char*)my_string,(char *)sequence);
	}
	}
}
return(0);
}
