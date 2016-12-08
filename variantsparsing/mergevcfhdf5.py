#!/usr/bin/env python

import h5py
import os
import scipy as sp
import pdb
import utilities.hdf5 as hdf5
import fnmatch
import sys

if __name__ == "__main__":
    tempvarfiles = sys.argv[1]
    fn_output = sys.argv[2]

    files = []
    for rts,dirs,fs in os.walk(tempvarfiles):
        fs = fnmatch.filter(fs, '*.hdf5')
        files.extend([os.path.join(rts, f) for f in fs])


    OUT   = h5py.File(fn_output, 'w')

    for i,f in enumerate(files):

        print "Processed %i/%i files" % (i, len(files))

        print f
        IN = h5py.File(os.path.join(tempvarfiles, f) , 'r')

        

        if not 'gt' in IN.keys(): ### empty vcf and so hdf5 file
            continue
    
    
        if not 'gtid' in OUT.keys(): ### this is the first one
            OUT.create_dataset(name = 'gtid', data = IN['gtid'][:], chunks= True, compression = 'gzip')
            OUT.create_dataset(name = 'allele_alt', data = IN['allele_alt'][:], chunks = True, compression = 'gzip', maxshape= (None,))
            OUT.create_dataset(name = 'allele_ref', data = IN['allele_ref'][:], chunks =True, compression= 'gzip', maxshape = (None,))
            OUT.create_dataset(name = 'pos', data = IN['pos'][:], chunks = True, compression = 'gzip', maxshape = (None, 2))
            OUT.create_dataset(name = 'gt', data = IN['gt'][:], compression = 'gzip', maxshape = (None, IN['gt'].shape[1]), chunks = (5000, IN['gt'].shape[1]))
        else:
            hdf5.appendToHDF5(OUT, IN['allele_alt'][:], 'allele_alt')
            hdf5.appendToHDF5(OUT, IN['allele_ref'][:], 'allele_ref')
            hdf5.appendToHDF5(OUT, IN['pos'][:], 'pos')
            hdf5.appendToHDF5(OUT, IN['gt'][:], 'gt')
        IN.close()
    OUT.close()


