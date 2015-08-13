# Unpack the Chothia numbered antibody files into directory and build
# the distances.h include file
make distances

# Build the programs
make

# Build a loop database
./buildloopdb /data/pdb >loops.db

# Scan a structure against the loop database
./scanloopdb -l looplen loops.db file.pdb > file.hits


# NOTES
# -----
# makedistances.pl - Script to build distances.h from data in abdb
#                    directory. This is all run and data downloaded to
#                    run it by doing 'make distances'
# finddist.c       - Calculates the distance matrix for an individual
#                    antibody file. Used by makedistances.pl. Not
#                     designed to be called by the end user
# buildloopdb.c    - Builds the loop database
#                    './buildloopdb -h' for help
# scanloopdb.c     - Scans a PDB file against the loop database to
#                    find potential loops
#                    './scanloopdb -h' for help




