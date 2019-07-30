#!/usr/bin/env perl
# Pass the name of a fasta file as an arg and this gives the number of
# bases in each read of the fasta file
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

if ($#ARGV != 0) {
    open (FILE, $0);
    while ($line = <FILE>) {
	last unless ($line =~ /^\#/);
	print $line;
    }
    close (FILE);
    exit;
}
$file = $ARGV[0];
$cmd = "zcat -f $file |";
open (FILE, $cmd);
$isFirstRead = 1;
while ($line = <FILE>) {
    if ($line =~ /^>/) {
	if (! $isFirstRead) { print "$readLen\n"; }
	$readLen = 0;
	$isFirstRead = 0;
    }
    else {
	$len = length ($line)-1;
	$readLen += $len;
    }
}
if (! $isFirstRead) { print "$readLen\n"; }
