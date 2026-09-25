# cython: language_level=3
"""Optional screening without Python alignment objects (pysam build required)."""
from libc.stdint cimport uint8_t
from libc.string cimport memcmp
from cpython.exc cimport PyErr_CheckSignals
from pysam.libcalignmentfile cimport AlignmentFile
from pysam.libchtslib cimport bam1_t, bam_get_seq, bam_init1, bam_destroy1, sam_read1, sam_write1

cdef extern from *:
    const char *TELOMERECAT_PYSAM_VERSION

import pysam
if pysam.__version__ != (<bytes>TELOMERECAT_PYSAM_VERSION).decode('ascii'):
    raise ImportError('The optional HTS screening extension must be rebuilt for this pysam version')

cdef int matches(bam1_t *read) except -1:
    cdef int length=read.core.l_qseq
    cdef int i, last
    cdef uint8_t *seq
    if length==0:
        raise TypeError("argument of type 'NoneType' is not iterable")
    if length<12:
        return 0
    seq=bam_get_seq(read)
    last=(length-12)//2
    for i in range(last+1):
        if seq[i]==0x88 and memcmp(seq+i,b'\x88\x14\x44\x88\x14\x44',6)==0:
            return 1
        if seq[i]==0x22 and memcmp(seq+i,b'\x22\x28\x11\x22\x28\x11',6)==0:
            return 1
        if 2*i+13<=length:
            if seq[i]&15==8 and seq[i+6]>>4==4 and memcmp(seq+i+1,b'\x81\x44\x48\x81\x44',5)==0:
                return 1
            if seq[i]&15==2 and seq[i+6]>>4==1 and memcmp(seq+i+1,b'\x22\x81\x12\x22\x81',5)==0:
                return 1
    return 0

def pairs_to_telbam(AlignmentFile src, AlignmentFile dst):
    if src.closed:
        raise ValueError('I/O operation on closed file')
    cdef bam1_t *first=bam_init1()
    cdef bam1_t *second=bam_init1()
    cdef int status
    cdef unsigned int since_signal_check=0
    if first==NULL or second==NULL:
        bam_destroy1(first);bam_destroy1(second)
        raise MemoryError()
    try:
        while True:
            status=sam_read1(src.htsfile,src.header.ptr,first)
            if status==-1:break
            if status < -1:raise OSError('Reading first mate failed')
            status=sam_read1(src.htsfile,src.header.ptr,second)
            if status==-1:raise StopIteration()
            if status < -1:raise OSError('Reading second mate failed')
            if matches(first) or matches(second):
                if dst.closed:
                    raise ValueError('I/O operation on closed file')
                if sam_write1(dst.htsfile,dst.header.ptr,first)<0:raise OSError('Writing first mate failed')
                if sam_write1(dst.htsfile,dst.header.ptr,second)<0:raise OSError('Writing second mate failed')
            since_signal_check+=1
            if since_signal_check==65536:
                PyErr_CheckSignals()
                since_signal_check=0
    finally:
        bam_destroy1(first);bam_destroy1(second)
