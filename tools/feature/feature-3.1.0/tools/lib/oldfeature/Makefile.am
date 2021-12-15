## Process this file with automake to produce Makefile.in

# the script files will get installed in a directory like /usr/local/feature
instdir = ${prefix}/feature

# these are the executables
inst_SCRIPTS = train_model.py scan_model.py

# these are other files that need to be installed
inst_DATA = callable.py featureexts.py fetchpdb.py \
    feature.py fetchdssp.py fetchproteinfile.py  \
    calciums.site feature.ini params.ini README COPYING

# this ensure the files get into the tarball
EXTRA_DIST = ${inst_SCRIPTS} ${inst_DATA}
