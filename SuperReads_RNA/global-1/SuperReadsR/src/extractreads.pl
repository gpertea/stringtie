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

$infile1 = $ARGV[0]; #file with read names
$infile2 = $ARGV[1]; #fasta file
$fieldnum = $ARGV[2]-1; #field number

open (FILE1, $infile1);
open (FILE2, $infile2);
my %readnames;
while ($line = <FILE1>){
    chomp($line);
    $readnames{$line}=1;
  }
close(FILE1);

my $sequence="";
my $readname="";
open (OUTFILE, ">$outfile");
while ($line = <FILE2>){
    if ($line =~ /^>/){
	if(defined $readnames{$readname}){
	print ">$readname\n$sequence\n";
	}
	chomp($line);
        @f=split(/\s+/,substr($line,1));
        $readname=$f[$fieldnum];
        $sequence="";
    }else{
     chomp($line);
     $sequence.=$line;
    }
}
        if(defined $readnames{$readname})
        {
        print ">$readname\n$sequence\n";
        }

close (FILE2);

