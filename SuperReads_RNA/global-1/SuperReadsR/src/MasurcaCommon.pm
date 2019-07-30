package MasurcaCommon;

use strict;
use warnings;

use MasurcaConf qw(fail);

our (@ISA, @EXPORT, @EXPORT_OK);

BEGIN {
  require Exporter;
  @ISA = qw(Exporter);
  @EXPORT_OK = qw($rerun_pe $rerun_sj $list_frg_files);
}

our ($rerun_pe, $rerun_sj) = (0, 0);

sub rename_reads {
  my ($out, $type, @read_info) = @_;
  my @res;
  my $mean_file = "meanAndStdevByPrefix.$type.txt";
  my $do_wait = 0;

  print $out "log 'Processing $type library reads'\n";
  print $out "rm -rf $mean_file\n";
  foreach my $v (@read_info){
    my ($name, $mean, $stdev, $fr, $rr) = @$v;
    my $renamed = "$name.renamed.fastq";
    ($mean, $stdev) = (500, 100) if($type eq "sj");
    push(@res, $renamed);
    print $out "echo '$name $mean $stdev' >> $mean_file\n";
    next if(-e $renamed);
    $do_wait = 1;
    $rerun_pe = 1;
    $rerun_sj = 1;
    print($out "run_bg rename_filter_fastq '$name' <(exec expand_fastq '$fr' | awk '{if(length(\$0>200)) print substr(\$0,1,200); else print \$0;}') ",
           ($fr eq $rr) ? "''" : "<(exec expand_fastq '$rr' | awk '{if(length(\$0>200)) print substr(\$0,1,200); else print \$0;}' )",
           " > '$renamed'\n");
  }
  return @res;
}

sub estimate_optimal_kmer{
  my ($out, $filelist, $name) = @_;
  my $max_kmer=127;
  
  if($name eq "KMER_J"){
    print $out "$name=31\n";
  } else {
    print $out "$name=`for f in $filelist;do head -n 80000 \$f |tail -n 40000;done | perl -e 'while(\$line=<STDIN>){\$line=<STDIN>;chomp(\$line);push(\@lines,\$line);\$line=<STDIN>;\$line=<STDIN>}\$min_len=100000;\$base_count=0;foreach \$l(\@lines){\$base_count+=length(\$l);push(\@lengths,length(\$l));\@f=split(\"\",\$l);foreach \$base(\@f){if(uc(\$base) eq \"G\" || uc(\$base) eq \"C\"){\$gc_count++}}} \@lengths =sort {\$b <=> \$a} \@lengths; \$min_len=\$lengths[int(\$\#lengths*.75)];  \$gc_ratio=\$gc_count/\$base_count;\$kmer=0;if(\$gc_ratio<0.5){\$kmer=int(\$min_len*.7);}elsif(\$gc_ratio>=0.5 && \$gc_ratio<0.6){\$kmer=int(\$min_len*.5);}else{\$kmer=int(\$min_len*.33);} \$kmer++ if(\$kmer\%2==0); \$kmer=31 if(\$kmer<31); \$kmer=$max_kmer if(\$kmer>$max_kmer); print \$kmer'`\n";
  }
}

sub get_MIN_Q_CHAR{
  my ($out, $filelist) = @_;
  my @fns=split(/\s+/, $filelist);
  print $out <<"EOS";
MIN_Q_CHAR=`head -n 50000 $fns[0] | awk 'BEGIN{flag=0}{if(\$0 ~ /^\\+/){flag=1}else if(flag==1){print \$0;flag=0}}'  | perl -ne 'BEGIN{\$q0_char=\"\@\";}{chomp;\@f=split \"\";foreach \$v(\@f){if(ord(\$v)<ord(\$q0_char)){\$q0_char=\$v;}}}END{\$ans=ord(\$q0_char);if(\$ans<64){print \"33\\n\"}else{print \"64\\n\"}}'`
save MIN_Q_CHAR
echo MIN_Q_CHAR: \$MIN_Q_CHAR
EOS
}

sub count_kmers_for_quorum {
  my ($out, $filelist, $NUM_THREADS, $thresh) = @_;
  print $out <<"EOS";
quorum_create_database -t $NUM_THREADS -s \$JF_SIZE -b 7 -m 24 -q \$((MIN_Q_CHAR + $thresh)) -o quorum_mer_db.jf.tmp $filelist && mv quorum_mer_db.jf.tmp quorum_mer_db.jf
if [ $? != 0 ]; then
  fail Increase JF_SIZE in config file, the recommendation is to set this to genome_size*coverage/2
fi
EOS
}
