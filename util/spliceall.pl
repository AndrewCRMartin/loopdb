#!/usr/bin/perl 
use strict;

my $bindir     = $ENV{'HOME'} . "/bin";

my $splicepdb    = "$bindir/splicepdb";
my $renumabloop  = "$bindir/renumabloop";
my $pdbprep      = '/data/pdb/pdb';
my $pdbext       = '.ent';

my $framework = shift(@ARGV);

while(<>)
{
    chomp;
    my @fields = split;
    my $pdb = $fields[0];
    my $start = $fields[1];
    my $stop = $fields[2];
    my $deviation = $fields[14];
    my $pdbfile = $pdbprep . $pdb . $pdbext;
    my $splicedfile = "${framework}_${pdb}";

    if( -e $pdbfile)
    {
        `$splicepdb -l $start $stop $pdbfile H95 H102 $framework $splicedfile`;
        `$renumabloop $splicedfile $splicedfile.num`;

        print STDERR "$splicedfile.num\n";

        last;  ################## BREAK OUT #####################
    }
    else
    {
        print STDERR "PDB file does not exist (skipped): $pdbfile\n";
    }
}
