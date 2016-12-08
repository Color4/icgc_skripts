#!/usr/bin/env python
import scipy as sp
import h5py
import pdb
import utilities.preprocess as prep
import oop_libs.genotype as geno
import pdb
import sys
import progressbar
offset = 10000

if __name__ == "__main__":
    fn_gt = sys.argv[1]
    fn_kinship = sys.argv[2]
    SNP = h5py.File(fn_gt, 'r')
    Genotype = geno.genotype()

    print "CREATING GENOTYPE OBJECT"
    bar = progressbar.ProgressBar(max_value = SNP['gt'].shape[0] / offset)
    for i in xrange(0, SNP['gt'].shape[0], offset):
        bar.update(i / offset)
        gt = SNP['gt'][i:min(i+offset, SNP['gt'].shape[0]),:]
        pos = SNP['pos'][i:min(i+offset, SNP['gt'].shape[0]),:]
        gt = gt[::20,:]
        pos = pos[::20,:]
        Genotype.addgt(gt, pos, maf = 0.1, msf = 0.05)

    Genotype.getibd()
    OUT = h5py.File(fn_kinship, 'w')
    OUT.create_dataset(name = 'K', data = Genotype.K)
    OUT.close()
        
