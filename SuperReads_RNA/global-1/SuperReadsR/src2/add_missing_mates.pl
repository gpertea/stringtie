#!/usr/bin/env perl
# SuperRead pipeline
# Copyright (C) 2012  Genome group at University of Maryland.
# 
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#
#This program adds missing mates to the read fasta files "on the fly"

use strict;
use warnings;

my $readnumberHold = -1;
my ($editlineHold, $prefixHold, $sequenceHold);
while(my $line=<STDIN>){
    chomp($line);
    my ($first, $editline) = split(' ', $line, 2);
    my $prefix             = substr($first, 1, 2);
    my $readnumber         = int(substr($first, 3));

    if(($readnumber % 2) == 0){#if the read is even we simply remember it
	if($readnumberHold!=-1){
          print(">", $prefixHold, $readnumberHold, " ", $editlineHold, "\n", $sequenceHold,
                "\n>", $prefixHold, $readnumberHold+1, "\nN\n");
	}
	$prefixHold     = $prefix;
	$readnumberHold = $readnumber;
        $editlineHold   = $editline;
        $sequenceHold   = <STDIN>;
        chomp($sequenceHold);
    }
    elsif($readnumberHold==-1){#the previous even read is missing
      print(">", $prefix, $readnumber-1,"\nN\n", $line, "\n", scalar(<STDIN>));
    }
    elsif($readnumber-1!=$readnumberHold){#previous mate is missing odd and current is missing even
      print(">", $prefixHold, $readnumberHold, " ", $editlineHold, "\n", $sequenceHold, 
            "\n>", $prefixHold, $readnumberHold+1, "\nN\n>", $prefix, $readnumber-1, "\nN\n",
            $line, "\n", scalar(<STDIN>));
      $readnumberHold = -1;
    }
    else { # if($readnumber-1==$readnumberHold)
      print(">", $prefixHold, $readnumberHold, " ", $editlineHold, "\n", $sequenceHold, "\n", $line, "\n", scalar(<STDIN>));
      $readnumberHold = -1;
    }
}

if($readnumberHold != -1) {
  print(">", $prefixHold, $readnumberHold, " ", $editlineHold, "\n", $sequenceHold, 
        "\n>", $prefixHold, $readnumberHold+1, "\nN\n");
}

