
import sys,os
from pymol import cmd
import gzip
import pdbutils

def gzload(id):
    # if filename doesn't exist, try it as a PDBID
    if not os.path.exists(id):
        filename = pdbutils.pdbfilename(id)
        label = id
    else:
        filename = id
        label = os.path.basename(filename)

    # Abort if file doesn't exist
    if not filename or not os.path.exists(filename):
        print "Error: id=%(id)s Can not find %(filename)s" % vars()
        return

    # read the data
    fh = gzip.open(filename)
    data = fh.read()
    fh.close()

    # load it into pymol
    cmd.read_pdbstr(data,label)

if __name__ == 'pymol' or __name__ == 'feature.gui.pymol.gzload':
    cmd.extend('gzload',gzload)
else:
    print __name__
