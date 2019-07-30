package MasurcaSuperReads;

use strict;
use warnings;
use MasurcaConf qw(fail);
use MasurcaCommon qw($rerun_pe $rerun_sj $list_frg_files);


sub filter_jump {
  #globals that can be changed here: $rerun_pe, $rerun_sj
  #const globals:  $$config{NUM_THREADS}, $$config{EXTEND_JUMP_READS} LIMIT_JUMP_COVERAGE 
  my ($out, %config) = @_;

  #creating super reads for filtering
  if($rerun_pe==1 || $rerun_sj==1 || not(-e "work2")){
    print $out "rm -rf work2\n";
    $rerun_sj=1;
  }
  print $out "createSuperReadsForDirectory.perl  -maxnodes 2000 -minreadsinsuperread 1 -l \$KMER_J -join-aggressive 1 -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.sj.txt -kunitigsfile guillaumeKUnitigsAtLeast32bases_all.jump.fasta -t $config{NUM_THREADS} -mikedebug work2 sj.cor.fa 1> super2.err 2>&1\n";

  #check if the super reads pipeline finished successfully
  print $out "if [[ ! -e work2/superReads.success ]];then\n";
  print $out "fail Super reads failed, check super2.err and files in ./work2/\n";
  print $out "fi\n";

  #now, using read positions in super reads, we find out which mates got joined -- these are the ones that do not have the biotin in the middle, call them chimeric
  if(not(-e "chimeric_sj.txt")||$rerun_pe==1||$rerun_sj==1){
    print $out "filter_alt.pl outtie < work2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt >  chimeric_sj.txt \n";
    $rerun_sj=1;
  }

  #we also do initial redundancy filtering here, based on positions of reads in suoer reads
  if(not(-e "redundant_sj.txt")||$rerun_pe==1||$rerun_sj==1){
    print $out "filter_redundancy.pl 2 < work2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt > redundant_sj.txt\n";
    $rerun_sj=1;
  } 
  print $out "echo 'Chimeric/Redundant jump reads:';wc -l  chimeric_sj.txt redundant_sj.txt;\n";

  #remove all chimeric and all redundant reads from sj.cor.fa
  if(not(-e "sj.cor.clean.rev.fa")||$rerun_pe==1||$rerun_sj==1){
    print $out "extractreads.pl <(cat chimeric_sj.txt redundant_sj.txt | perl -e '{
		while(\$line=<STDIN>){
		chomp(\$line);
		\$h{\$line}=1
		}
		open(FILE,\$ARGV[0]);
		while(\$line=<FILE>){
		chomp(\$line);
		print \$line,\"\\n\" if(not(defined(\$h{\$line})));
		}
		}' <(awk '{
			prefix=substr(\$1,1,2); 
			readnumber=int(substr(\$1,3));  
			if(readnumber\%2==0){
				last_readnumber=readnumber; 
				last_prefix=prefix;
			}else{
				if(last_readnumber==readnumber-1 && last_prefix==prefix){
					print prefix\"\"last_readnumber\"\\n\"prefix\"\"readnumber;
				}
			}
			}' work2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt)) sj.cor.fa 1 | putReadsIntoGroupsBasedOnSuperReads --super-read-sequence-file work2/superReadSequences.fasta --read-placements-file work2/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt > sj.cor.clean.fa\n";

    #here we perform another round of filtering bad mates
    print $out "findReversePointingJumpingReads_bigGenomes.perl --jellyfish-hash-size \$JF_SIZE --kmer-step-size 5 --reads-file sj.cor.clean.fa --reads-for-kunitigs-file pe.cor.fa --reads-for-kunitigs-file sj.cor.fa --dir-to-change-to  work2.1 --dir-for-kunitigs work2.1 --min-kmer-len 41 --max-kmer-len 81 -t $config{NUM_THREADS} --maxnodes 1000 --kmer-step-size 10 1>findReversePointingJumpingReads.err 2>&1 \n";
        print $out "extractreads_not.pl work2.1/readsToExclude.txt sj.cor.clean.fa 1 > sj.cor.clean2.fa\n";
    print $out "echo Found extra chimeric mates: \n";
    print $out "wc -l work2.1/readsToExclude.txt\n";
    print $out "rm -f sj.cor.clean.rev.fa\n";
    foreach my $v (@{$config{JUMP_INFO}}) {
      my @f        = @$v;
      my $if_innie = "";
      $if_innie    = " | reverse_complement " if($f[1]>0);
      print $out "grep --text -A 1 '^>$f[0]' sj.cor.clean2.fa | grep --text -v '^\\-\\-' $if_innie >> sj.cor.clean.rev.fa\n";
    }
    $rerun_sj=1;
  }

  #here we extend the jumping library reads if they are too short
  if(not($config{SOAP_ASSEMBLY})){
    if(not(-e "sj.cor.ext.fa")||$rerun_pe==1||$rerun_sj==1){      
      $rerun_sj=1;
      if($config{EXTEND_JUMP_READS} == 1){
        print $out "createSuperReadsForDirectory.perl -jumplibraryreads -minreadsinsuperread 1 -l \$KMER_J -mean-and-stdev-by-prefix-file meanAndStdevByPrefix.sj.txt -kunitigsfile guillaumeKUnitigsAtLeast32bases_all.jump.fasta -t $config{NUM_THREADS} -mikedebug work3 sj.cor.clean.rev.fa 1> super2.err 2>&1\n";

        #check if the super reads pipeline finished successfully
        print $out "if [[ ! -e work3/superReads.success ]];then\n";
        print $out "fail Super reads failed, check super2.err and files in ./work2/\n";
        print $out "fi\n";

        print $out "ln -sf work3/superReadSequences.jumpLibrary.fasta sj.cor.ext.fa\n";
      }else{
	print $out "ln -sf sj.cor.clean.rev.fa sj.cor.ext.fa\n";
      }
    }

    #here we create the frg files for CA from the jump libraries: each jump library will contribute one jump frg file and one additional frg file of linking information from "chimers"
    print $out "log 'Creating FRG files'\n";
    print $out "rm -rf compute_jump_coverage.txt\n";
    
    foreach my $v (@{$config{JUMP_INFO}}) {
      my @f = @$v;
      print $out "echo -n \"",abs($f[1])," \" >> compute_jump_coverage.txt\n";
      print $out "grep --text -A 1 '^>$f[0]' sj.cor.ext.fa | grep --text -v '^\\-\\-' > $f[0].tmp\n";
      print $out "error_corrected2frg $f[0] ",abs($f[1])," $f[2] 2000000000 $f[0].tmp | grep --text '^{LKG' |wc -l >> compute_jump_coverage.txt\n";
    }
    print $out "JUMP_BASES_COVERED=`awk 'BEGIN{b=0}{b+=\$1*\$2;}END{print b}' compute_jump_coverage.txt`\n";
    print $out "save JUMP_BASES_COVERED\n";

    #here we reduce jump library coverage: we know the genome size (from k-unitigs) and JUMP_BASES_COVERED contains total jump library coverage :)
    foreach my $v (@{$config{JUMP_INFO}}) {
      my @f = @$v;
      $list_frg_files .= "$f[0].cor.clean.frg ";
      print $out "grep --text -A 1 '^>$f[0]' sj.cor.ext.fa | grep --text -v '^\\-\\-' | sample_mate_pairs.pl $config{LIMIT_JUMP_COVERAGE} `perl -e 'print int('\$JUMP_BASES_COVERED'/'\$ESTIMATED_GENOME_SIZE'/1.6)'` 1 > $f[0].tmp\n";
      print $out "error_corrected2frg $f[0] ",abs($f[1])," $f[2] 2000000000 $f[0].tmp > $f[0].cor.clean.frg\n";
      print $out "rm -f $f[0].tmp\n";
    }
  }
}

sub create_pe_linking_mates {
  my ($out, %config) = @_;
  if(not(-e "pe.linking.fa")||$rerun_pe==1){ 
  print $out "ufasta extract -f <( awk 'BEGIN{last_readnumber=-1;last_super_read=\"\"}{readnumber=int(substr(\$1,3));if(readnumber%2>0){readnumber--}super_read=\$2;if(readnumber==last_readnumber){if(super_read!=last_super_read){print read;print \$1;}}else{read=\$1;last_super_read=\$2}last_readnumber=readnumber}' work1/readPlacementsInSuperReads.final.read.superRead.offset.ori.txt )  pe.cor.fa > pe.linking.fa.tmp && mv pe.linking.fa.tmp pe.linking.fa\n";
  }

  print $out "NUM_LINKING_MATES=`wc -l pe.linking.fa | perl -ane '{print int(\$F[0]/2)}'`\n";
  if(@{$config{JUMP_INFO}}){
    print $out "MAX_LINKING_MATES=`perl -e '{\$g=int('\$ESTIMATED_GENOME_SIZE'/25);\$g=100000000 if(\$g>100000000);print \$g}'`\n";
  } else {
    print $out "MAX_LINKING_MATES=`perl -e '{\$g=int('\$ESTIMATED_GENOME_SIZE'/2);\$g=100000000 if(\$g>100000000);print \$g}'`\n";
  }

  foreach my $v (@{$config{PE_INFO}}){
    my @f = @$v;
    $list_frg_files .= "$f[0].linking.frg ";
    if(not(-e "$f[0].linking.frg")||$rerun_pe==1){
      print $out <<"EOS";      
grep --text -A 1 '^>$f[0]' pe.linking.fa | grep --text -v '^\\-\\-' | sample_mate_pairs.pl \$MAX_LINKING_MATES \$NUM_LINKING_MATES 1 > $f[0].tmp && error_corrected2frg $f[0] $f[1] $f[2] 2000000000 $f[0].tmp > $f[0].linking.frg.tmp && rm $f[0].tmp && mv $f[0].linking.frg.tmp $f[0].linking.frg
EOS
    }
  }
  print $out "echo \"Using linking mates\"\n";
}


1;
