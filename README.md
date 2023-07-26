loopdb
======

loopdb is a set of programs for creating and screening a database of
loops extracted from all files in the Protein Databank. This is done
by having a set of distances between loop takeoff residues (the three
residues either side of a loop) and finding all fragments in the PDB
that satisfy these distance constraints.

By default, the database building creates a database of loops that
have takeoff positions suitable for replacing CDR-H3 in an
antibody. This can be overridden by providing a distance table such as
the example `distanceTable.txt`

The loop scanning program then takes a required length for the loop
and a PDB file onto which you wish to fit a loop. It scans the
database to find loops of the correct length and sorts them by the
RMSD over the takeoff residues.

INSTALLATION
------------

### Build the programs

To update the distances information, simply remove the file src/distances.h
If this file is missing, the PDB files for antibodies will be downloaded and
the mean and standard deviation distances calculated and written to 
src/distances.h

    cd src
    make
    make install

The executables will be installed in the bin directory

### Update the loop database

Assuming you have the PDB installed in /data/pdb, do

    ./bin/buildloopdb /data/pdb >data/loops.db

### Search the database

To search the database for a loop of a given length that fits onto a PDB
file of interest:

    ./bin/scanloopdb -l looplen data/loops.db file.pdb > file.hits

(where `looplen` is the length of the loop of interest).

DOCUMENTATION
-------------

### buildloopdb.c


Builds the loop database `./buildloopdb -h` for help. By default this
will build a loopm database for antibody CDR-H3. This can be changed
by generating a different version of distances.h (by modifying
finddist.c to specify a different set of residues). Alternatively you
can provide a file to the program via the `-t` flag that will override
the default distances.

Note that the resulting loop database contains the following fields:

1.  The PDB code
2.  The first residue (including the 3 residue overlap with framework)
3.  The last residue (including the 3 residue overlap with framework)
4.  The length of the loop (excluding the overlap with framework)
5.  n0-c0
6.  n0-c1
7.  n0-c2
8.  n1-c0 
9.  n1-c1 
10. n1-c2
11. n2-c0
12. n2-c1
13. n2-c2

### scanloopdb.c

Scans a PDB file against the loop database to find potential loops
`./scanloopdb -h` for help. Again the default is to find loops that
fit onto the takeoff residues for antibody CDR-H3, but this can be
overridden using the `-r` flag.

The following programs do not need to be run by the end user
------------------------------------------------------------

### makedistances.pl

Script to build distances.h from data in abdb directory. This is all
run and data downloaded to run it by doing `(cd src; rm distances.h;
make distances.h)`

### finddist.c

Calculates the distance matrix for an individual antibody file. Used
by makedistances.pl. Not designed to be called by the end user

