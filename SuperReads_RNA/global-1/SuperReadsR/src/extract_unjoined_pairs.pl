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

my $flag=0;

while($line=<STDIN>)
{
chomp($line);
@F=split(/\s+/,$line);
if($flag==0)
{
$sr_pair=$F[1];
if(int(substr($F[0],2))%2==0)
{
$insert=$F[0];
}
else
{
$insert=substr($F[0],0,2).int(substr($F[0],2))-1;
}
$flag=1;
}
elsif($flag==1)
{
$srh{"$sr_pair $F[1]"}.="$insert ";
$srh{"$F[1] $sr_pair"}.="$insert ";
push(@pairs,"$sr_pair $F[1]");
$flag=0;
}
}

foreach $v (@pairs)
{
@f=split(/ /,$v);
if(defined($srh{"$f[0] $f[1]"}) && defined($srh{"$f[1] $f[0]"}))
{
print "$v : ",$srh{"$f[0] $f[1]"},"\n";
}
}


