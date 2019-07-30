#!/usr/bin/env perl
# Take the input from the first arg
# Takes an optional second arg: The max number of times to output a k-mer
#    (default is 2)
$maxTimesToRepeat = 2;
$inputFile = $ARGV[0];
if ($#ARGV >= 1) {
    $maxTimesToRepeat = $ARGV[1]; }

print ">0\n";
$cmd = "jellyfish dump -c $inputFile |";
open (FILE, $cmd);
while ($line = <FILE>) {
    chomp ($line);
    ($kmer, $count) = split (" ", $line);
    if ($count> $maxTimesToRepeat) {
	$countOut = $maxTimesToRepeat; }
    else {
	$countOut = $count; }
    for ($i=0; $i<$countOut; $i++) {
	print $kmer,"N"; }
}
close (FILE);
print "\n";
