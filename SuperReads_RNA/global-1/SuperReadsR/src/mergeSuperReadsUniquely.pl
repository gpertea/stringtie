#!/usr/bin/env perl

$workingDir=$ARGV[0];

open(FILE,"$workingDir/reduce.tmp");
while($line=<FILE>){
	chomp($line);
	@f=split(/\s+/,$line);
	$reduce{$f[0]}=1;
	}
open(FILE,"$workingDir/superReadNames.txt");
open(OUTFILE,">$workingDir/superReadNames.reduced.txt");
while($line=<FILE>){
	chomp($line);
	print OUTFILE $line,"\n" if(!defined($reduce{$line}));
}

open(FILE,"$workingDir/superReadNames.reduced.txt");
while($line=<FILE>){
	chomp($line);
	if($line =~ /_/){
		@f=split(/_/,$line);
		for($i=1;$i<$#f;$i++){
			$mid_sr{substr($f[$i],0,length($f[$i])-1)}=1;
			}
	}else{
		$mid_sr{substr($line,0,length($line)-1)}=1;
	}
}
close(FILE);

open(FILE,"$workingDir/superReadNames.reduced.txt");
while($line=<FILE>){
	chomp($line);
	@f=split(/_/,$line);
	if(not(defined($mid_sr{substr($f[0],0,length($f[0])-1)}))){
		$merge_k_u{substr($f[0],0,length($f[0])-1)}++;
	}
	if(not(defined($mid_sr{substr($f[-1],0,length($f[-1])-1)}))){
		$merge_k_u{substr($f[-1],0,length($f[-1])-1)}++;
	}
}

%osr=();
open(FILE,"$workingDir/superReadNames.reduced.txt");
while($line=<FILE>){
	chomp($line);
	@f=split(/_/,$line); 
	($eku1,$ori1,$eku2,$ori2)=($line=~/(^\d+)(F|R).*_(\d+)(F|R)$/);
	if($merge_k_u{$eku1}==2){
		if($ori1 eq "F"){
                $$osr{$eku1}[1]=$line;
		}else{
		$$osr{$eku1}[0]=reverse_sr($line);
		}
	}
        if($merge_k_u{$eku2}==2){
                if($ori2 eq "F"){
		$$osr{$eku2}[0]=$line;
                }else{
		$$osr{$eku2}[1]=reverse_sr($line);
                }
        }
} 

%usedKUnitigs=();
#now we do the actual merging.  pick up k-unitig k, and merge the two reads AB , then examine the outer k-unitigs of the two reads and put these on the list, etc.  The super read A specifiessss orientation on the whole thing.  We first build a space-separated list of the super reads for each super merge, then we assign all super reads to the super merge and then we do the actual super merge
#

foreach $k(keys %merge_k_u){
next if($merge_k_u{$k}!=2);
next if(not(defined($$osr{$k}[0]) || not(defined($$osr{$k}[1]))));
$valid_merges{$k}=1;
}

foreach $k(keys %valid_merges){
next if($usedKUnitigs{$k});
#merge loop
$merged_sr="$$osr{$k}[0] $$osr{$k}[1]";
$usedKUnitigs{$k}=1;
#now we examine the end k-unitigs of the $merged_sr and see if we can extend it
$cont=1;
while($cont){
$cont=0;
($eku1,$ori1,$eku2,$ori2)=($merged_sr=~/(^\d+)(F|R).*_(\d+)(F|R)$/);
@f=split(" ",$merged_sr);
if(defined($valid_merges{$eku1}) && !$usedKUnitigs{$eku1}){
#means we have another merge on eku1
$usedKUnitigs{$eku1}=1;
if($ori1 eq "F"){
$merged_sr="$$osr{$eku1}[0] $merged_sr";
}else{
$merged_sr=reverse_sr($$osr{$eku1}[1])." $merged_sr";
}
$cont++;
}
if(defined($valid_merges{$eku2}) && !$usedKUnitigs{$eku2}){
#means we have another merge on eku2
$usedKUnitigs{$eku2}=1;
if($ori2 eq "F"){
$merged_sr="$merged_sr $$osr{$eku2}[1]";
}else{
$merged_sr="$merged_sr ".reverse_sr($$osr{$eku2}[0]);
}
$cont++;
}
}
#as we exit the loop, we should have our preliminary merged sr on $merged_sr, and all super reads in merged_sr pointing to it
#we mangled the super read orientations.  The rule of output -- first   k-unitig is larger than the last
#here we do output
#print "DEBUG final $merged_sr\n";
#now we have our merged_sr, let's record which reads reduced to it and their orientation
@f=split(" ",$merged_sr);
$merged_sr_final=$f[0];
for($i=1;$i<=$#f;$i++){
($ttt)=($f[$i]=~/^\d+.(_.*)/);
$merged_sr_final.=$ttt;
}

if($eku1>$eku2){
$merged_sr_final=reverse_sr($merged_sr_final);
$reversed=1;
}
push(@newmerged,$merged_sr_final);
#print "DEBUG final $merged_sr_final\n";
foreach $sr(@f){
($eku1,$ori1,$eku2,$ori2)=($sr=~/(^\d+)(F|R).*_(\d+)(F|R)$/);
if($eku1<$eku2){
#then this is the original
if($reversed){
$merge{$sr}="$merged_sr_final 9999 R";
print "$sr $merged_sr_final 9999 R\n";
}else{
$merge{$sr}="$merged_sr_final 0 F";
print "$sr $merged_sr_final 0 F\n";
}
}else{
if($reversed){
$merge{reverse_sr($sr)}="$merged_sr_final 0 F";
print reverse_sr($sr)," $merged_sr_final 0 F\n";
}else{
$merge{reverse_sr($sr)}="$merged_sr_final 9999 R";
print reverse_sr($sr)," $merged_sr_final 9999 R\n";
}
}
}
}
#now 

open(OUTFILE,">$workingDir/superReadNames.newmerged.txt");
foreach $k(@newmerged){
	print OUTFILE $k,"\n";
	}
close(OUTFILE);

open(FILE,"$workingDir/reduce.tmp");
open(OUTFILE,">$workingDir/reduce.merged.tmp");
while($line=<FILE>){
	chomp($line);
	@f=split(/\s/,$line);
	if(defined($merge{$f[1]})){
	print OUTFILE "$f[0] $merge{$f[1]}\n";
	}else{
	print OUTFILE "$line\n";
	}
}
open(FILE,"$workingDir/superReadNames.reduced.txt");
while($line=<FILE>){
        chomp($line);
        if(defined($merge{$line})){
        	print OUTFILE "$line $merge{$line}\n";
	}
}
close(OUTFILE);


sub reverse_sr
{
    my $old_sr=$_[0];
    my $new_sr;
    my @f=split(/_/,$old_sr);
    $new_sr=$f[-1];
    $new_sr=~tr/FR/RF/;
    for(my $i=$#f-1;$i>=0;$i--){
	$f[$i]=~tr/FR/RF/;
	$new_sr.="_".$f[$i];
    }
    return($new_sr);
}

