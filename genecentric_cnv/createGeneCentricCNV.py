import scipy as sp

import pdb

import sys
import os

import fnmatch

from intervaltree import Interval, IntervalTree

def readcnv(fn, intervalltree):
    cnvdata = sp.loadtxt(fn, delimiter = '\t', dtype = 'string')
    cnvdata = cnvdata[1:,:]
    dict_results = {}
    for i,c in enumerate(cnvdata):
        if int(c[-1]) <= 1: ###
            continue
        if int(c[3]) == 2: ### normal CNV
            continue
        gene_region =  intervalltree[int(c[0])][ int(c[1]): int(c[2])]
        if len(gene_region) == 0:
            continue
        cnvgenes = sp.array([x.data[1].split(' ')[-1] for x in gene_region])
        cnvgenes = sp.array([x.strip('\'').strip('\"') for x in cnvgenes])
        for g in cnvgenes:
            if g in dict_results:
                dict_results[g].append(int(c[3]))
            else:
                dict_results[g] = [int(c[3])]

    gid = []
    cnv = []

    for g in dict_results.keys():
        gid.append(g)
        cnv.append(sp.mean(dict_results[g]))
    return gid, cnv
#            if r.data[1] in dict_results:
                


if __name__ == "__main__":
    ### gencode file in gtf format
    fn_anno = sys.argv[1]
    ### aliquot list in tsv (claudia format)
    fn_aliquotlist =sys.argv[2] 
    ### directory of cnv calls (one aliquot per file)
    dir_cnv = sys.argv[3] 
    ### out file
    fn_out = sys.argv[4] 


    ### read annotation
    anno_data  = sp.loadtxt(fn_anno, delimiter = '\t', dtype = 'string', usecols=[0,2,3,4,8])
    # filter to gene annotation only
    anno_data = anno_data[anno_data[:,1] == 'gene',:]
    # remove X/Y/MT (since we do not care in QTL)
    iAutosome = sp.array([x[0].isdigit() for x in anno_data])
    anno_data = anno_data[iAutosome,:]
    
    ## create dictionary of interval trees
    intervalltree={}
    for i in xrange(1,23):
        anno_data_chr = anno_data[anno_data[:,0] == str(i),:]
        ttree = IntervalTree()
        for gn in anno_data_chr:
            ttree[int(gn[2]):int(gn[3])] = [gn[0], gn[4].split(';')[0]]
        intervalltree[i] = ttree



    ### load aliquot data
    aq_data = sp.loadtxt(fn_aliquotlist, delimiter = '\t', dtype = 'string', usecols = [2])
    aq_data = aq_data[1:]
    

    ### get aliquot file list
    cnv_files = sp.array(os.listdir(dir_cnv))

    allgids = sp.array([x.split(';')[0].split(' ')[1].strip('\"') for x in anno_data[:,4]])
    allgids = sp.sort(allgids)
    cnvdata = sp.zeros(( aq_data.shape[0], anno_data.shape[0] ))
    cnvdata[:] = sp.nan
    for i,aq in enumerate(aq_data): ### read each aliquot
        print aq
        ix_f = sp.array([x.startswith(aq) for x in cnv_files])
        if not ix_f.any():
            continue

        gid, cnv = readcnv(os.path.join(dir_cnv, cnv_files[ix_f][0]), intervalltree)
        gid = sp.array(gid)
        cnv = sp.array(cnv)
    
        sidx = sp.argsort(gid)
        gid  = gid[sidx]
        cnv = cnv[sidx]
        
        midx = sp.in1d( allgids, gid)
        cnvdata[i,midx] = cnv 

    ### OUTPUT statement
    OUT = h5py.File(fn_out, 'w')
    OUT.create_dataset(name = 'cnv', data = cnvdata.T)
    OUT.create_dataset(name = 'gene_name', data = allgids)
    OUT.create_dataset(name = 'wgs_aliquot_id', data = aq_data)
    OUT.close()

