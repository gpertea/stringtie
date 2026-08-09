#!/usr/bin/env perl
use strict;
use warnings;

sub usage {
    die "Usage: compare_merge_outputs.pl guide.gtf old.gtf new.gtf input1.gtf [input2.gtf ...]\n";
}

sub transcript_id {
    my ($attrs) = @_;
    return $1 if $attrs =~ /(?:^|;\s*)transcript_id\s+"([^"]+)"/;
    return '';
}

sub parse_gtf {
    my ($path) = @_;
    open(my $fh, '<', $path) or die "Error opening $path: $!\n";
    my @records;
    my $current;
    while (my $line = <$fh>) {
        next if $line =~ /^\s*#/;
        chomp $line;
        my @f = split(/\t/, $line, 9);
        next unless @f == 9;
        my $tid = transcript_id($f[8]);
        next if $tid eq '';
        if ($f[2] eq 'transcript' || $f[2] eq 'mRNA') {
            $current = {
                ref => $f[0], start => 0 + $f[3], end => 0 + $f[4],
                strand => $f[6], tid => $tid, exons => []
            };
            push(@records, $current);
        }
        elsif ($f[2] eq 'exon') {
            next unless defined($current) && $current->{tid} eq $tid;
            push(@{$current->{exons}}, [0 + $f[3], 0 + $f[4]]);
        }
    }
    close($fh);

    ## derive exact-structure and intron-chain keys after all exons are loaded
    for my $r (@records) {
        if (!@{$r->{exons}}) {
            push(@{$r->{exons}}, [$r->{start}, $r->{end}]);
        }
        my $coords = join(',', map { $_->[0] . '-' . $_->[1] } @{$r->{exons}});
        $r->{key} = join('|', $r->{ref}, $r->{strand}, $coords);
        my @introns;
        for (my $i = 1; $i < @{$r->{exons}}; $i++) {
            push(@introns, $r->{exons}->[$i - 1]->[1] . '-' . $r->{exons}->[$i]->[0]);
        }
        $r->{chain} = join('|', $r->{ref}, $r->{strand}, @introns);
    }
    return \@records;
}

sub make_set {
    my ($records, $field) = @_;
    my %set;
    $set{$_->{$field}} = 1 for @$records;
    return \%set;
}

sub count_intersection {
    my ($a, $b) = @_;
    my $count = 0;
    $count++ for grep { exists($b->{$_}) } keys(%$a);
    return $count;
}

sub summarize {
    my ($records, $guide_by_id, $source_structures) = @_;
    my (%ids, %structures, %structure_records, %met_chains, %guide_records);
    my ($sets, $mets, $source_records) = (0, 0, 0);
    for my $r (@$records) {
        $ids{$r->{tid}}++;
        $structures{$r->{key}}++;
        push(@{$structure_records{$r->{key}}}, $r);
        if (@{$r->{exons}} == 1) {
            $sets++;
        }
        else {
            $mets++;
            $met_chains{$r->{chain}} = 1;
        }
        push(@{$guide_records{$r->{tid}}}, $r) if exists($guide_by_id->{$r->{tid}});
        $source_records++ if exists($source_structures->{$r->{key}});
    }

    my ($duplicate_ids, $extra_ids) = (0, 0);
    for my $count (values(%ids)) {
        if ($count > 1) {
            $duplicate_ids++;
            $extra_ids += $count - 1;
        }
    }
    my ($duplicate_structures, $extra_structures, $guide_duplicate_groups,
        $generated_duplicate_groups, $generated_duplicate_records) = (0, 0, 0, 0, 0);
    for my $key (keys(%structures)) {
        my $count = $structures{$key};
        if ($count > 1) {
            $duplicate_structures++;
            $extra_structures += $count - 1;
            my $generated = grep { !exists($guide_by_id->{$_->{tid}}) } @{$structure_records{$key}};
            if ($generated) {
                $generated_duplicate_groups++;
                $generated_duplicate_records += $generated;
            }
            else {
                $guide_duplicate_groups++;
            }
        }
    }

    my ($guide_exact, $guide_changed, $guide_missing, $guide_outer, $guide_intron) = (0, 0, 0, 0, 0);
    for my $tid (keys(%$guide_by_id)) {
        my $guide = $guide_by_id->{$tid};
        my $outs = $guide_records{$tid};
        if (!$outs) {
            $guide_missing++;
            next;
        }
        if (grep { $_->{key} eq $guide->{key} } @$outs) {
            $guide_exact++;
            next;
        }
        $guide_changed++;
        my $out = $outs->[0];
        $guide_outer++ if $out->{start} != $guide->{start} || $out->{end} != $guide->{end};
        $guide_intron++ if $out->{chain} ne $guide->{chain};
    }

    my $source_unique = 0;
    $source_unique++ for grep { exists($source_structures->{$_}) } keys(%structures);
    return {
        transcript_records => scalar(@$records),
        unique_transcript_ids => scalar(keys(%ids)),
        duplicate_id_values => $duplicate_ids,
        duplicate_id_extra_records => $extra_ids,
        set_records => $sets,
        met_records => $mets,
        unique_structures => scalar(keys(%structures)),
        duplicate_structure_groups => $duplicate_structures,
        duplicate_structure_extra_records => $extra_structures,
        guide_only_duplicate_structure_groups => $guide_duplicate_groups,
        duplicate_groups_with_generated_records => $generated_duplicate_groups,
        generated_records_in_duplicate_groups => $generated_duplicate_records,
        unique_met_intron_chains => scalar(keys(%met_chains)),
        guide_id_records => scalar(map { @$_ } values(%guide_records)),
        generated_id_records => scalar(@$records) - scalar(map { @$_ } values(%guide_records)),
        guide_ids_exact => $guide_exact,
        guide_ids_changed => $guide_changed,
        guide_ids_missing => $guide_missing,
        guide_ids_outer_changed => $guide_outer,
        guide_ids_intron_changed => $guide_intron,
        source_structure_records => $source_records,
        unique_source_structures => $source_unique,
        synthesized_structure_records => scalar(@$records) - $source_records,
        unique_synthesized_structures => scalar(keys(%structures)) - $source_unique,
        structure_set => \%structures,
    };
}

@ARGV >= 4 or usage();
my ($guide_path, $old_path, $new_path, @input_paths) = @ARGV;
my $guides = parse_gtf($guide_path);
my $old = parse_gtf($old_path);
my $new = parse_gtf($new_path);

my %guide_by_id;
$guide_by_id{$_->{tid}} = $_ for @$guides;
my %source_structures = %{make_set($guides, 'key')};
for my $path (@input_paths) {
    my $records = parse_gtf($path);
    $source_structures{$_->{key}} = 1 for @$records;
}

my $old_summary = summarize($old, \%guide_by_id, \%source_structures);
my $new_summary = summarize($new, \%guide_by_id, \%source_structures);
my @metrics = qw(
    transcript_records unique_transcript_ids duplicate_id_values duplicate_id_extra_records
    set_records met_records unique_structures duplicate_structure_groups
    duplicate_structure_extra_records guide_only_duplicate_structure_groups
    duplicate_groups_with_generated_records generated_records_in_duplicate_groups
    unique_met_intron_chains guide_id_records generated_id_records guide_ids_exact
    guide_ids_changed guide_ids_missing guide_ids_outer_changed guide_ids_intron_changed
    source_structure_records unique_source_structures synthesized_structure_records
    unique_synthesized_structures
);

print "metric\tv3.0.3\tnew_default\tdelta\tpercent_change\n";
for my $metric (@metrics) {
    my $old_value = $old_summary->{$metric};
    my $new_value = $new_summary->{$metric};
    my $delta = $new_value - $old_value;
    my $percent = $old_value ? sprintf('%.2f', 100 * $delta / $old_value) : 'NA';
    print join("\t", $metric, $old_value, $new_value, $delta, $percent), "\n";
}

my $old_set = $old_summary->{structure_set};
my $new_set = $new_summary->{structure_set};
my $common = count_intersection($old_set, $new_set);
my $old_only = scalar(keys(%$old_set)) - $common;
my $new_only = scalar(keys(%$new_set)) - $common;
my $union = $common + $old_only + $new_only;
print "comparison\tvalue\n";
print "source_unique_structures\t", scalar(keys(%source_structures)), "\n";
print "common_unique_structures\t$common\n";
print "v3.0.3_only_unique_structures\t$old_only\n";
print "new_only_unique_structures\t$new_only\n";
print "structure_jaccard_percent\t", sprintf('%.2f', $union ? 100 * $common / $union : 0), "\n";
