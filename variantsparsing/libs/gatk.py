import scipy as sp
import numpy as np
import scipy.stats as spst
import gzip as gz
import pdb
import os
import fnmatch
import warnings
import time
import sys
import utilities.hdf5 as hdf5
import time
gt_dict = {'0/0':0, '1/1':2, '0/1':1, '1/0':1, './.':sp.nan}
gt_dict_indel = {'0/0':1, '1/1':1, '0/1':1, '1/0':1, './.':sp.nan}




def readVCF_pyvcf(filename, delimiter = '\t', missing = './.', quality = 20, useInDel = False, f_transgt=None):
    import vcf
    ### chrTidx converter

    ### check if we got an empty file
    if os.stat(filename)[6] == 0:
        return dict()


    vcf_reader  = vcf.Reader(filename= filename)
    vcf_samples = vcf_reader.samples
    f_transgt.create_dataset(name  = 'gtid', data = sp.array(vcf_samples))
    first       = True

    for cnter,record in enumerate(vcf_reader):#.fetch()):#chrm,int(start),int(start)+int(offset))):

        if record.INFO['AF'][0]  < 0.001: ### set allele frequency filter to at least 0.1 %
            continue
        if len(record.FILTER) > 0:
            continue

        if len(record.ALT) > 1: ### only bi-allelic
            continue
        if record.QUAL < quality: ## require quality of at least 100 (really weak threshold considering)
            continue



        ### get position
        if record.CHROM.isdigit():
            pos = sp.array([int(record.CHROM),record.POS])
        else:
            pos = sp.array([chrTidx(record.CHROM),record.POS])

        ### get call
        if record.is_snp:
            gt        = sp.array([gt_dict[x['GT']] for x in record.samples], dtype = 'float')
        elif record.is_indel:
            gt        = sp.array([gt_dict_indel[x['GT']] for x in record.samples], dtype = 'float')
            
        else:
            continue
        if sp.unique(gt[~sp.isnan(gt)]).shape[0] <= 1:
            continue

        ra = record.REF
        aa = record.ALT
        aa = ",".join(sp.array(aa).astype('string'))
        maf = record.INFO['AF'][0]
        if first:
            f_transgt.create_dataset(name = 'pos', data = pos[sp.newaxis,:], compression = 'gzip', chunks = True, maxshape = (None, 2))
            f_transgt.create_dataset(name  = 'gt', data = gt[sp.newaxis, :], compression = 'gzip', chunks = True, maxshape = (None, gt.shape[0]))
            f_transgt.create_dataset(name = 'allele_ref', data = [ra], chunks = True,compression = 'gzip', maxshape = (None,))
            f_transgt.create_dataset(name = 'allele_alt', data = [aa], chunks = True, compression = 'gzip', maxshape = (None,))
            f_transgt.create_dataset(name = 'maf', data = [maf], chunks = True, compression = 'gzip', maxshape = (None,))
            
            first = False
        else:
            hdf5.appendToHDF5(f_transgt, pos[sp.newaxis,:], 'pos')
            hdf5.appendToHDF5(f_transgt, gt[sp.newaxis,:], 'gt')
            hdf5.appendToHDF5(f_transgt, sp.array([ra]), 'allele_ref')
            hdf5.appendToHDF5(f_transgt, sp.array([aa]), 'allele_alt')
            hdf5.appendToHDF5(f_transgt, sp.array([maf]), 'maf')
