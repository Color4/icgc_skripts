#!/usr/bin/env python2.7


import sys
import scipy as sp
import os
import pdb
import gzip
import h5py
import libs.gatk as gatk
if __name__ == "__main__":
    fn = sys.argv[1]
    path = sys.argv[2]
    fn_out = fn.split('/')[-1] + '.hdf5'
    fn_out = os.path.join(path, fn_out)
    OUT    = h5py.File(fn_out, 'w')
    gatk.readVCF_pyvcf(fn, f_transgt = OUT)
    OUT.close()
    sys.exit(0)
