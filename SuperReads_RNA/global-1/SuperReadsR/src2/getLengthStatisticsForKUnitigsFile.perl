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

# Cat the k-unitig fasta file through
$workingDir = ".";
$outputPrefix = "";
for ($i=0; $i<=$#ARGV; $i++) {
    $arg = $ARGV[$i];
    if ($arg eq "-output-prefix") {
	++$i;
	$outputPrefix = $ARGV[$i];
	if ($outputPrefix !~ /\.$/) { $outputPrefix .= "."; }
	next; }
    $workingDir = $arg;
}

$numKUnitigsFile = $outputPrefix . "numKUnitigs.txt";
$maxKUnitigNumberFile = $outputPrefix . "maxKUnitigNumber.txt";
$totBasesInKUnitigsFile = $outputPrefix . "totBasesInKUnitigs.txt";
$isFirstRead = 1;
while ($line = <STDIN>) {
    if ($line =~ /^>/) {
	if (! $isFirstRead) { $kUnitigLengths[$kUnitig] = $kUnitigLength; }
	$kUnitigLength = 0;
	$isFirstRead = 0;
	($kUnitig) = ($line =~ /^.(\S+)\s/);
    }
    else {
	$len = length ($line)-1;
	$kUnitigLength += $len;
    }
}
if (! $isFirstRead) { $kUnitigLengths[$kUnitig] = $kUnitigLength; }

for ($i=0; $i<=$#kUnitigLengths; $i++) {
    $length = $kUnitigLengths[$i];
    $totBasesInKUnitigs += $length;
    if (! $length) {
	$length = 0; }
    else {
	++$numKUnitigs; }
    print "$i $length\n";
}
open (OUTFILE, ">", "$workingDir/$numKUnitigsFile");
print OUTFILE "$numKUnitigs\n";
close (OUTFILE);
open (OUTFILE, ">", "$workingDir/$maxKUnitigNumberFile");
$arraySizeForKUnitigData = $#kUnitigLengths+1;
print OUTFILE "$arraySizeForKUnitigData\n";
close (OUTFILE);
open (OUTFILE, ">", "$workingDir/$totBasesInKUnitigsFile");
print OUTFILE "$totBasesInKUnitigs\n";
close (OUTFILE);
