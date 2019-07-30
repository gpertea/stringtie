#!/usr/bin/env perl

# lightweight no-dependency version of pkg-config. This will work on any machine
# with Perl installed.

# Copyright (C) 2012 M. Nunberg.
# You may use and distribute this software under the same terms and conditions
# as Perl itself.

package PkgConfig::Vars;
# this is a namespace for .pc files to hold their variables without
# relying on lexical scope.

package PkgConfig::UDefs;
# This namespace provides user-defined variables which are to override any
# declarations within the .pc file itself.

package PkgConfig;

#First two digits are Perl version, second two are pkg-config version
our $VERSION = '0.05020';

require 5.005;

use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Class::Struct; #in core since 5.004
our $UseDebugging;

use Data::Dumper;

################################################################################
### Check for Log::Fu                                                        ###
################################################################################
BEGIN {
    my $ret = eval {
        require Log::Fu;
        Log::Fu->import({level => "warn"}); 1;
    };
    
    if(!$ret) {
        my $log_base = sub {
            my (@args) = @_;
            print STDERR "[DEBUG] ", join(' ', @args);
            print STDERR "\n";
        };
        *log_debug = *log_debugf = sub { return unless $UseDebugging; goto &$log_base };
        *log_err = *log_errf = *log_warn = *log_warnf = *log_info = *log_infof =
            $log_base;
        
    }
}

our $VarClassSerial = 0;

################################################################################
### Sane Defaults                                                            ###
################################################################################
our @DEFAULT_SEARCH_PATH = qw(
    /usr/lib/pkgconfig /usr/share/pkgconfig
    /usr/local/lib/pkgconfig /usr/local/share/pkgconfig

);

my @ENV_SEARCH_PATH = split(/:/, $ENV{PKG_CONFIG_PATH} || "");

push @DEFAULT_SEARCH_PATH, @ENV_SEARCH_PATH;

our @DEFAULT_EXCLUDE_CFLAGS = qw(-I/usr/include -I/usr/local/include);
# don't include default link/search paths!
our @DEFAULT_EXCLUDE_LFLAGS = qw(
    -L/usr/lib -L/lib -L/lib64 -L/lib32
    -L/usr/lib32 -L/usr/lib64
    -L/usr/local/lib
    
    -R/lib -R/usr/lib -R/usr/lib64 -R/lib32 -R/lib64 -R/usr/local/lib
);


my $LD_OUTPUT_RE = qr/
    SEARCH_DIR\("
    ([^"]+)
    "\)
/x;

sub GuessPaths {
    my $pkg = shift;
    local $ENV{LD_LIBRARY_PATH} = "";
    local $ENV{C_INCLUDE_PATH} = "";
    local $ENV{LD_RUN_PATH} = "";
    
    my $ld = $ENV{LD} || 'ld';
    my $ld_output = qx(ld -verbose);
    my @defl_search_dirs = ($ld_output =~ m/$LD_OUTPUT_RE/g);
    
    @DEFAULT_EXCLUDE_LFLAGS = ();
    foreach my $path (@defl_search_dirs) {
        push @DEFAULT_EXCLUDE_LFLAGS, (map { "$_".$path }
            (qw(-R -L -rpath= -rpath-link= -rpath -rpath-link))); 
    }
    log_debug("Determined exclude LDFLAGS", @DEFAULT_EXCLUDE_LFLAGS);
    
    #now get the include paths:
    my @cpp_output = qx(cpp --verbose 2>&1 < /dev/null);
    @cpp_output = map  { chomp $_; $_ } @cpp_output;
    #log_info(join("!", @cpp_output));
    while (my $cpp_line = shift @cpp_output) {
        chomp($cpp_line);
        if($cpp_line =~ /\s*#include\s*<.+search starts here/) {
            last;
        }
    }
    #log_info(@cpp_output);
    my @include_paths;
    while (my $path = shift @cpp_output) {
        if($path =~ /\s*End of search list/) {
            last;
        }
        push @include_paths, $path;
    }
    @DEFAULT_EXCLUDE_CFLAGS = map { "-I$_" } @include_paths;
    log_debug("Determine exclude CFLAGS", @DEFAULT_EXCLUDE_CFLAGS);
}


################################################################################
### Define our fields                                                        ###
################################################################################
struct(
    __PACKAGE__,
    [
     # .pc search paths, defaults to PKG_CONFIG_PATH in environment
     'search_path' => '@',

     # whether to also spit out static dependencies
     'static' => '$',
     
     # whether we replace references to -L and friends with -Wl,-rpath, etc.
     'rpath' => '$',
     
     # build rpath-search,
     
     # no recursion. set if we just want a version, or to see if the
     # package exists.
     'no_recurse' => '$',

     #list of cflags and ldflags to exclude
     'exclude_ldflags' => '@',
     'exclude_cflags' => '@',

     # what level of recursion we're at
     'recursion' => '$',

     # hash of libraries, keyed by recursion levels. Lower recursion numbers
     # will be listed first
     'libs_deplist' => '*%',

     # cummulative cflags and ldflags
     'ldflags'   => '*@',
     'cflags'    => '*@',

     # whether we print the c/ldflags
     'print_cflags' => '$',
     'print_ldflags' => '$',

     # information about our top-level package
     'pkg'  => '$',
     'pkg_exists' => '$',
     'pkg_version' => '$',
     'pkg_url', => '$',
     'pkg_description' => '$',
     'errmsg'   => '$',
     
     # classes used for storing persistent data
     'varclass' => '$',
     'udefclass' => '$',
     'filevars' => '*%',
     'uservars' => '*%',
     
     # options for printing variables
     'print_variables' => '$',
     'print_values' => '$',
     'defined_variables' => '*%',
    ]
);

################################################################################
################################################################################
### Variable Storage                                                         ###
################################################################################
################################################################################

sub _get_pc_varname {
    my ($self,$vname_base) = @_;
    return $self->varclass . "::" . $vname_base;
}

sub _get_pc_udefname {
    my ($self,$vname_base) = @_;
    return $self->udefclass . "::" . $vname_base;
}

sub _pc_var {
    my ($self,$vname) = @_;
    $vname =~ s,\.,DOT,g;
    no strict 'refs';
    $vname = $self->_get_pc_varname($vname);
    my $glob = *{$vname};
    return unless $glob;
    
    return $$glob;
}

sub assign_var {
    my ($self,$field,$value) = @_;
    no strict 'refs';
    
    # if the user has provided a definition, use that.
    if(exists ${$self->udefclass."::"}{$field}) {
        log_debug("Prefix already defined by user");
        return;
    }
    my $evalstr = sprintf('$%s = %s',
                    $self->_get_pc_varname($field), $value);
    
    log_debug("EVAL", $evalstr);
    eval $evalstr;
    if($@) {
        log_err($@);
    }
}

sub prepare_vars {
    my $self = shift;
    my $varclass = $self->varclass;
    no strict 'refs';
    
    %{$varclass . "::"} = ();
    
    while (my ($name,$glob) = each %{$self->udefclass."::"}) {
        my $ref = *$glob{SCALAR};
        next unless defined $ref;
        ${"$varclass\::$name"} = $$ref;
    }
}

################################################################################
################################################################################
### Initializer                                                              ###
################################################################################
################################################################################
sub find {
    my ($cls,$library,%options) = @_;
    my @uspecs = (
        ['search_path', \@DEFAULT_SEARCH_PATH],
        ['exclude_ldflags', \@DEFAULT_EXCLUDE_LFLAGS],
        ['exclude_cflags', \@DEFAULT_EXCLUDE_CFLAGS]
    );
    
    foreach (@uspecs) {
        my ($basekey,$default) = @$_;
        my $list = [ @{$options{$basekey} ||= [] } ];
        if($options{$basekey . "_override"}) {
            @$list = @{ delete $options{$basekey."_override"} };
        } else {
            push @$list, @$default;
        }

        $options{$basekey} = $list;
        #print "$basekey: " . Dumper($list);
    }
    
    $VarClassSerial++;
    $options{varclass} = sprintf("PkgConfig::Vars::SERIAL_%d", $VarClassSerial);
    $options{udefclass} = sprintf("PkgConfig::UDefs::SERIAL_%d", $VarClassSerial);
    
    
    my $udefs = delete $options{VARS} || {};
    
    while (my ($k,$v) = each %$udefs) {
        no strict 'refs';
        my $vname = join('::', $options{udefclass}, $k);
        ${$vname} = $v;
    }
    
    my $o = $cls->new(%options);
    
    my @libraries;
    if(ref $library eq 'ARRAY') {
        @libraries = @$library;
    } else {
        @libraries = ($library);
    }
    
    foreach my $lib (@libraries) {
        $o->recursion(0);
        $o->find_pcfile($lib);
    }
    
    return $o;
}

################################################################################
################################################################################
### Modify our flags stack                                                   ###
################################################################################
################################################################################
sub append_ldflags {
    my ($self,@flags) = @_;
    my @ld_flags = _split_flags(@flags);
    
    foreach my $ldflag (@ld_flags) {
        next unless $ldflag =~ /^-Wl/;

        my (@wlflags) = split(/,/, $ldflag);
        shift @wlflags; #first is -Wl,
        filter_omit(\@wlflags, $self->exclude_ldflags);
        
        if(!@wlflags) {
            $ldflag = "";
            next;
        }
        
        $ldflag = join(",", '-Wl', @wlflags);
    }
    
    @ld_flags = grep $_, @ld_flags;
    return unless @ld_flags;
    
    push @{($self->libs_deplist->{$self->recursion} ||=[])},
        @ld_flags;
}

# notify us about extra compiler flags
sub append_cflags {
    my ($self,@flags) = @_;
    push @{$self->cflags}, _split_flags(@flags);
}


################################################################################
################################################################################
### All sorts of parsing is here                                             ###
################################################################################
################################################################################
sub get_requires {
    my ($self,$requires) = @_;
    return () unless $requires;

    my @reqlist = split(/[\s,]+/, $requires);
    my @ret;
    while (defined (my $req = shift @reqlist) ) {
        my $reqlet = [ $req ];
        push @ret, $reqlet;
        last unless @reqlist;
        #check if we need some version scanning:

        my $cmp_op;
        my $want;

        GT_PARSE_REQ:
        {
            #all in one word:
            ($cmp_op) = ($req =~ /([<>=]+)/);
            if($cmp_op) {
                if($req =~ /[<>=]+$/) {
                    log_debug("comparison operator spaced ($cmp_op)");
                    ($want) = ($req =~ /([^<>=]+$)/);
                    $want ||= shift @reqlist;
                } else {
                    $want = shift @reqlist;
                }
                push @$reqlet, ($cmp_op, $want);
            } elsif ($reqlist[0] =~ /[<>=]+/) {
                $req = shift @reqlist;
                goto GT_PARSE_REQ;
            }
        }
    }
    #log_debug(@ret);
    return @ret;
}


sub parse_line {
    my ($self,$line,$evals) = @_;
    no strict 'vars';
    
    $line =~ s/#[^#]+$//g; # strip comments
    return unless $line;
    
    my ($tok) = ($line =~ /([=:])/);
    
    my ($field,$value) = split(/[=:]/, $line, 2);
    return unless defined $value;
    
    if($tok eq '=') {
        $self->defined_variables->{$field} = $value;
    }
    
    #strip trailing/leading whitespace:
    $field =~ s/(^\s+)|(\s+)$//msg;
    
    #remove trailing/leading whitespace from value
    $value =~ s/(^\s+)|(\s+$)//msg;

    log_debugf("Field %s, Value %s", $field, $value);
    
    #perl variables can't have '.' in them:
    $field =~ s/\./DOT/g;
    
    
    $field = lc($field);
    
    #remove quoutes from field names
    $field =~ s/['"]//g;
    

    # pkg-config escapes a '$' with a '$$'. This won't go in perl:
    $value =~ s/[^\\]\$\$/\\\$/g;
    $value =~ s/([@%&])/\$1/g;
    
    
    # append our pseudo-package for persistence.
    my $varclass = $self->varclass;
    $value =~ s/(\$\{[^}]+\})/lc($1)/ge;
    
    $value =~ s/\$\{/\$\{$varclass\::/g;
    
    #quoute the value string, unless quouted already
    $value = "\"$value\"" unless $value =~ /^["']/;
    
    #get existent variables from our hash:
    
    
    $value =~ s/'/"/g; #allow for interpolation
    
    $self->assign_var($field, $value);
    
}

sub parse_pcfile {
    my ($self,$pcfile,$wantversion) = @_;
    #log_warn("Requesting $pcfile");
    open my $fh, "<", $pcfile or die "$pcfile: $!";
    
    $self->prepare_vars();
    
    my @lines = (<$fh>);
    close($fh);
    
    my $text = join("", @lines);
    $text =~ s,\\[\r\n],,g;
    @lines = split(/[\r\n]/, $text);
    
    my @eval_strings;
    
    #Fold lines:
    
    foreach my $line (@lines) {
        $self->parse_line($line, \@eval_strings);
    }
    
    #now that we have eval strings, evaluate them all within the same
    #lexical scope:
    

    $self->append_cflags(  $self->_pc_var('cflags') );
    $self->append_ldflags( $self->_pc_var('libs') );
    if($self->static) {
        $self->append_ldflags( $self->_pc_var('libs.private') );
    }

    my @deps;
    my @deps_dynamic = $self->get_requires( $self->_pc_var('requires'));
    my @deps_static = $self->get_requires( $self->_pc_var('requires.private') );
    @deps = @deps_dynamic;


    if($self->static) {
        push @deps, @deps_static;
    }

    if($self->recursion == 1 && (!$self->pkg_exists())) {
        $self->pkg_version( $self->_pc_var('version') );
        $self->pkg_url( $self->_pc_var('url') );
        $self->pkg_description( $self->_pc_var('description') );
        $self->pkg_exists(1);        
    }
    
    unless ($self->no_recurse) {
        foreach (@deps) {
            my ($dep,$cmp_op,$version) = @$_;
            $self->find_pcfile($dep);
        }
    }
}


################################################################################
################################################################################
### Locate and process a .pc file                                            ###
################################################################################
################################################################################
sub find_pcfile {
    my ($self,$libname,$version) = @_;

    $self->recursion($self->recursion + 1);

    my $pcfile = "$libname.pc";
    my $found = 0;
    my @found_paths = (grep {
        -e File::Spec->catfile($_, $pcfile)
        } @{$self->search_path});

    if(!@found_paths) {
        my @search_paths = @{$self->search_path};
        $self->errmsg(
            join("\n",
                 "Can't find $pcfile in any of @search_paths",
                 "use the PKG_CONFIG_PATH environment variable, or",
                 "specify extra search paths via 'search_paths'",
                 ""
                )
        ) unless $self->errmsg();
        return;
    }

    $pcfile = File::Spec->catfile($found_paths[0], $pcfile);

    $self->parse_pcfile($pcfile);

    $self->recursion($self->recursion - 1);
}

################################################################################
################################################################################
### Public Getters                                                           ###
################################################################################
################################################################################

sub get_cflags {
    my $self = shift;
    my @cflags = @{$self->cflags};

    filter_omit(\@cflags, $self->exclude_cflags);
    filter_dups(\@cflags);
    return @cflags;
}

sub get_ldflags {
    my $self = shift;
    my @ordered_libs;
    my @lib_levels = sort keys %{$self->libs_deplist};
    my @ret;

    @ordered_libs = @{$self->libs_deplist}{@lib_levels};
    foreach my $liblist (@ordered_libs) {
        my $lcopy = [ @$liblist ];
        filter_dups($lcopy);
        filter_omit($lcopy, $self->exclude_ldflags);
        push @ret, @$lcopy;
    }

    @ret = reverse @ret;
    filter_dups(\@ret);
    @ret = reverse(@ret);
    return @ret;
}



################################################################################
################################################################################
### Utility functions                                                        ###
################################################################################
################################################################################

#split a list of tokens by spaces
sub _split_flags {
    my @flags = @_;
    if(!@flags) {
        return @flags;
    }
    if(@flags == 1) {
        my $str = shift @flags;
        return () if !$str;
        @flags = split(/\s+/, $str);
    }
    @flags = grep $_, @flags;
    return @flags;
}



sub filter_dups {
    my $array = shift;
    my @ret;
    my %seen_hash;
    #@$array = reverse @$array;
    foreach my $elem (@$array) {
        if(exists $seen_hash{$elem}) {
            next;
        }
        $seen_hash{$elem} = 1;
        push @ret, $elem;
    }
    #print Dumper(\%seen_hash);
    @$array = @ret;
}

sub filter_omit {
    my ($have,$exclude) = @_;
    my @ret;
    #print Dumper($have);
    foreach my $elem (@$have) {
        #log_warn("Checking '$elem'");
        if(grep { $_ eq $elem } @$exclude) {
            #log_warn("Found illegal flag '$elem'");
            next;
        }
        push @ret, $elem;
    }
    @$have = @ret;
}

sub version_2_array {
    my $string = shift;
    my @chunks = split(/\./, $string);
    my @ret;
    my $chunk;
    while( ($chunk = pop @chunks)
        && $chunk =~ /^\d+$/) {
        push @ret, $chunk;
    }
    return @ret;
}


sub version_check {
    my ($want,$have) = @_;
    my @a_want = version_2_array($want);
    my @a_have = version_2_array($have);

    my $max_elem = scalar @a_want > scalar @a_have
        ? scalar @a_have
        : scalar @a_want;

    for(my $i = 0; $i < $max_elem; $i++) {
        if($a_want[$i] > $a_have[$i]) {
            return 0;
        }
    }
    return 1;
}


if(caller) {
    return 1;
}

################################################################################
################################################################################
################################################################################
################################################################################
### Script-Only stuff                                                        ###
################################################################################
################################################################################
################################################################################
################################################################################
package PkgConfig::Script;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my $quiet_errors = 1;
my @ARGV_PRESERVE = @ARGV;

my @POD_USAGE_SECTIONS = (
    "NAME",
    'DESCRIPTION/SCRIPT OPTIONS/USAGE',
    "DESCRIPTION/SCRIPT OPTIONS/ARGUMENTS|ENVIRONMENT",
    "AUTHOR & COPYRIGHT"
);

my @POD_USAGE_OPTIONS = (
    -verbose => 99,
    -sections => \@POD_USAGE_SECTIONS
);

GetOptions(
    'libs' => \my $PrintLibs,
    'static' => \my $UseStatic,
    'cflags' => \my $PrintCflags,
    'exists' => \my $PrintExists,
    
    'silence-errors' => \my $SilenceErrors,
    'print-errors' => \my $PrintErrors,
    
    'define-variable=s', => \my %UserVariables,
    
    'print-variables' => \my $PrintVariables,
    'print-values'  => \my $PrintValues,
    'variable=s',   => \my %OutputVariableValue,
    
    'modversion'    => \my $PrintVersion,
    'version',      => \my $PrintAPIversion,
    'real-version' => \my $PrintRealVersion,
    'atleast-pkgconfig-version=s' => \my $AtLeastVersion,
    
    'debug'         => \my $Debug,
    'with-path=s',    => \my @ExtraPaths,
    'env-only',     => \my $EnvOnly,
    'guess-paths',  => \my $GuessPaths,
    
    'h|help|?'      => \my $WantHelp
) or pod2usage(@POD_USAGE_OPTIONS);


if($WantHelp) {
    pod2usage(@POD_USAGE_OPTIONS, -exitval => 0);
}

if($Debug) {
    eval {
    Log::Fu::set_log_level('PkgConfig', 'DEBUG');
    };
    $PkgConfig::UseDebugging = 1;
}

if($GuessPaths) {
    PkgConfig->GuessPaths();
}

if($PrintAPIversion) {
    print "0.20\n";
    exit(0);
}

if($PrintRealVersion) {

    printf STDOUT ("pkg-config.pl - cruftless pkg-config\n" .
            "Version: %s\n", $PkgConfig::VERSION);
    exit(0);
}

if(defined($AtLeastVersion)) {
  my ($major, $minor, $bug) = split(/\./, $AtLeastVersion);
  exit(0) if !defined($major) || $major eq "";
  exit(1) if $major > 0;
  exit(0) if !defined($minor) || $minor eq "";
  exit(1) if $major == 0 && $minor > 20;
  exit(1) if defined($bug) && $bug ne "" && $minor == 20;
  exit(0);
}

my @FINDLIBS = @ARGV or die "Must specify at least one library";

if($PrintErrors) {
    $quiet_errors = 0;
}
if($SilenceErrors) {
    $quiet_errors = 1;
}

my $WantFlags = ($PrintCflags || $PrintLibs || $PrintVersion);

if($WantFlags) {
    $quiet_errors = 0 unless $SilenceErrors;
}

my %pc_options;
if($PrintExists) {
    $pc_options{no_recurse} = 1;
}


$pc_options{static} = $UseStatic;
$pc_options{search_path} = \@ExtraPaths;

if($EnvOnly) {
    delete $pc_options{search_path};
    $pc_options{search_path_override} = [ @ExtraPaths, @ENV_SEARCH_PATH];
}

$pc_options{print_variables} = $PrintVariables;
$pc_options{print_values} = $PrintValues;
$pc_options{VARS} = \%UserVariables;


my $o = PkgConfig->find(\@FINDLIBS, %pc_options);

if($o->errmsg) {
    print STDERR $o->errmsg unless $quiet_errors;
    exit(1);
}

if($o->print_variables) {
    while (my ($k,$v) = each %{$o->defined_variables}) {
        print $k;
        if($o->print_values) {
            print "=$v";
        } else {
            print "\n";
        }
    }
}

if(!$WantFlags) {
    exit(0);
}

if($PrintVersion) {
    print $o->pkg_version . "\n";
    exit(0);
}

if($PrintCflags) {
    print join(" ", $o->get_cflags) . " ";
}

if($PrintLibs) {
    print join(" ", $o->get_ldflags) . " ";
}

print "\n";
exit(0);

__END__

=head1 NAME

PkgConfig - Pure-Perl Core-Only replacement for C<pkg-config>

=head1 NOTE

The script is not actually installed yet (i haven't settled on a good name), but
will decide based on input in a future version.

Additionally, some dependencies are superficially included for debugging, and
will be sanitized in a future 'release/stable' version.

=head1 SYNOPSIS

=head2 As a replacement for C<pkg-config>

    $ pkg-config.pl --libs --cflags --static gio-2.0

    #outputs (lines artifically broken up for readability):
    # -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include
    # -pthread -lgio-2.0 -lz -lresolv -lgobject-2.0
    # -lgmodule-2.0 -ldl -lgthread-2.0 -pthread -lrt -lglib-2.0


Compare to:
    $ pkg-config --libs --cflags --static gio-2.0

    #outputs ( "" ):
    # -pthread -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include
    # -pthread -lgio-2.0 -lz -lresolv -lgobject-2.0 -lgmodule-2.0
    # -ldl -lgthread-2.0 -lrt -lglib-2.0


=head2 From another Perl module

    use PkgConfig;

    my $o = PkgConfig->find('gio');
    if($o->errmsg) {
        #handle error
    } else {
        my @cflags = $o->get_cflags;
        my @ldflags = $o->get_ldflags;
    }

=head1 DESCRIPTION

C<PkgConfig> provides a pure-perl, core-only replacement for the C<pkg-config>
utility.

This is not a description of the uses of C<pkg-config> but rather a description
of the differences between the C version and the Perl one.

While C<pkg-config> is a compiled binary linked with glib, the pure-perl version
has no such requirement, and will run wherever Perl ( >= 5.04 ) does.

The main supported options are the common C<--libs>, C<--cflags>,
C<--static>, C<--exists> and C<--modversion>.

=head2 SCRIPT OPTIONS

=head3 USAGE

    <packagename1 pkgname2..> [ --options ]

=head3 ARGUMENTS

By default, a library name must be supplied unless one of L<--version>,
or L<--real-version> is specified.

The output should normally be suitable for passing to your favorite compiler.

=head4 --libs

(Also) print linker flags. Dependencies are traverse in order. Top-level dependencies
will appear earlier in the command line than bottom-level dependencies.

=head4 --cflags

(Also) print compiler and C preprocessor flags.

=head4 --static

Use extra dependencies and libraries if linking against a static version of the
requested library

=head4 --exists

Return success (0) if the package exists in the search path.

=head4 --with-path=PATH

Prepend C<PATH> to the list of search paths containing C<.pc> files.

This option can be specified multiple times with different paths, and they will
all be added.

=head4 --env-only

Using this option, B<only> paths specified in C<PKG_CONFIG_PATH> are recognized
and any hard-coded defaults are ignored.

=head4 --guess-paths

Invoke C<gcc> and C<ld> to determine default linker and include paths. Default
paths will be excluded from explicit -L and -I flags.

=head4 --define-variable=VARIABLE=VALUE

Define a variable, overriding any such variable definition in the .pc file, and
allowing your value to interpolate with subsequent uses.

=head4 --print-variables

Print all defined variables found in the .pc files.



=head4 --version

The target version of C<pkg-config> emulated by this script

=head4 --real-version

The actual version of this script

=head4 --debug

Print debugging information

=head4 --silence-errors

Turn off errors. This is the default for non-libs/cflag/modversion
arguments

=head4 --print-errors

This makes all errors noisy and takes precedence over
C<--silence-errors>



=head3 ENVIRONMENT

the C<PKG_CONFIG_PATH> variable is honored and used as a colon-delimited list
of directories with contain C<.pc> files.

=head2 MODULE OPTIONS

=head4 I<< PkgConfig->find >>

    my $result = PkgConfig->find($libary, %options);

Find a library and return a result object.
C<$library> can be either a single name of a library, or a reference to an
array of library names

The options are in the form of hash keys and values, and the following are
recognized:

=over

=item C<search_path>

=item C<search_path_override>

Prepend search paths in addition to the paths specified in C<$ENV{PKG_CONFIG_PATH}>
The value is an array reference.

the C<_override> variant ignores defaults (like c<PKG_CONFIG_PATH).

=item C<exclude_cflags>

=item C<exclude_ldflags>

=item C<exclude_cflags_override>

=item C<exclude_ldflags_override>


Some C<.pc> files specify default compiler and linker search paths, e.g.
C<-I/usr/include -L/usr/lib>. Specifying them on the command line can be
problematic as it drastically changes the search order.

The above options will either append or replace the options which are excluded
and filtered.

The default excluded linker and compiler options can be obtained via
C<@PkgConfig::DEFAULT_EXCLUDE_LFLAGS> and C<@PkgConfig::DEFAULT_EXCLUDE_CFLAGS>,
respectively.

=item C<static>

Also specify static libraries.

=item C<no_recurse>

Do not recurse dependencies. This is useful for just doing version checks.

=item C<VARS>

Define a hashref of variables to override any variable definitions within
the .pc files. This is equivalent to the C<--define-variable> command-line
option.

=back

A C<PkgConfig> object is returned and may be queried about the results:

=head4 I<< $o->errmsg >>

An error message, if any. This is a string and indicates an error.

=head4 I<< $o->pkg_exists >>

Boolean value, true if the package exists.

=head4 I<< $o->pkg_version >>

The version of the package

=head4 I<< $o->get_cflags >>

=head4 I<< $o->get_ldflags >>

Returns a list of compiler and linker flags, respectively.

=head4 I<< PkgConfig->Guess >>

This is a class method, and will replace the hard-coded default linker and include
paths with those discovered by invoking L<ld(1)> and L<cpp(1)>.

Currently this only works with GCC-supplied C<ld> and GNU C<ld>.

=head2 BUGS

The order of the flags is not exactly matching to that of C<pkg-config>. From my
own observation, it seems this module does a better job, but I might be wrong.

Version checking is not yet implemented.

Unlike C<pkg-config>, the scripts C<--exists> function will return nonzero if
a package B<or> any of its dependencies are missing. This differs from the
behavior of C<pkg-config> which will just check for the definition of the
package itself (without dependencies).

=head1 SEE ALSO

L<ExtUtils::PkgConfig>, a wrapper around the C<pkg-config> binary

L<pkg-config|http://www.freedesktop.org/wiki/Software/pkg-config>

=head1 AUTHOR & COPYRIGHT

Copyright (C) 2012 M. Nunberg

You may use and distribute this software under the same terms and conditions as
Perl itself.
