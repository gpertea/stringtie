#! /usr/bin/env perl
# Pass 2 args: the name of the files containing
# (1) The original super-read names, and
# (2) The name of the input reduce file (e.g. reduce.tmp)
$superReadCount = 0;
open (FILE, $ARGV[0]);
while ($line = <FILE>) {
    chomp ($line);
    $newSuperReadName{$line} = $superReadCount;
    ++$superReadCount; }
close (FILE);

open (FILE, $ARGV[1]);
while ($line = <FILE>) {
    chomp ($line);
    @flds = split (" ", $line);
    $flds[0] = $newSuperReadName{$flds[0]};
    $flds[1] = $newSuperReadName{$flds[1]};
    print "@flds\n";
}
close (FILE);

