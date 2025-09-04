#!/usr/bin/env python
import os
import optparse
from collections import OrderedDict
from pyBioInfo.IO.File import FastqFile
from fbilr.reader import MatrixReader


def load_fastq(path):
    assert path.endswith(".gz")
    with FastqFile(path) as f:
        for read in f:
            yield read
           
            
def load_matrix(path):
    assert path.endswith(".gz")
    for record in MatrixReader.open(path):
        assert len(record.hits) == 2
        yield record
           
            
def write_fastq(fw, read):
    fw.write("@%s\n%s\n+\n%s\n" % (read.name, read.sequence, read.quality))
           
           
def main():
    parser = optparse.OptionParser(usage="%prog [options] input.fastq.gz matrix.tsv config.tsv outdir")
    parser.add_option("-l", "--min-length", dest="min_length", type="int", default=400, metavar="INT", 
                      help="Minimum length. [%default]")
    parser.add_option("-e", "--max-edit-distance", dest="max_edit_distance", type="int", default=5, metavar="INT", 
                      help="Maximum edit distance. [%default]")
    parser.add_option("-t", "--test", dest="test", action="store_true", default=False, 
                      help="For test running. [%default]")
    parser.add_option("-f", "--output-failed", dest="output_failed", action="store_true", default=False, 
                      help="Output failed reads. [%default]")
    options, args = parser.parse_args()

    f_fastq, f_matrix, f_config, outdir = args
    succeed_fastq_dir = os.path.join(outdir, "succeed")
    failed_fastq_dir = os.path.join(outdir, "failed")
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    if not os.path.exists(succeed_fastq_dir):
        os.mkdir(succeed_fastq_dir)
    if not os.path.exists(failed_fastq_dir):
        os.mkdir(failed_fastq_dir)
        
    barcode2cell = dict()
    fastq_handles = dict()
    read_counter = OrderedDict()
    
    with open(f_config) as f:
        for line in f:
            cell, barcode = line.strip("\n").split("\t")
            barcode2cell[barcode] = cell
            fw = open(os.path.join(succeed_fastq_dir, "%s.fastq" % cell), "w+")
            fastq_handles[barcode] = fw
            read_counter[barcode] = 0   
                 
    max_edit_distance = options.max_edit_distance
    min_read_length = options.min_length
    n_total = 0
    n_too_short = 0
    n_unknown = 0
    n_conflict = 0
    n_assigned = 0
    
    trim_barcode = True # trim barcode and outer
    
    output_failed = options.output_failed
    f_too_short = os.path.join(failed_fastq_dir, "too_short.fastq")
    f_unknown = os.path.join(failed_fastq_dir, "unknown.fastq")
    f_conflict = os.path.join(failed_fastq_dir, "conflict.fastq")
    h_too_short = None
    h_unknown = None
    h_conflict = None
    
    debug = options.test # DEBUG
    debug_reads = 100000
            
    if output_failed:
        h_too_short = open(f_too_short, "w+")
        h_unknown = open(f_unknown, "w+")
        h_conflict = open(f_conflict, "w+")
    
    reads = load_fastq(f_fastq)
    records = load_matrix(f_matrix)
    
    for read, record in zip(reads, records):
        if debug and n_total >= debug_reads:
            break
        
        n_total += 1
        
        if record.length < min_read_length:
            n_too_short += 1
            if output_failed:
                write_fastq(h_too_short, read)
            continue
        
        head, tail = record.hits[0], record.hits[1]
        if head.ed <= max_edit_distance and head.location == "H" and head.direction == "F" \
            and tail.ed <= max_edit_distance and tail.location == "T" and tail.direction == "R":
                
            if head.name == tail.name:
                n_assigned += 1
                if trim_barcode:
                    start = head.end
                    end = tail.start
                    read.name = read.name.split()[0]
                    read.sequence = read.sequence[start:end]
                    read.quality = read.quality[start:end]
                write_fastq(fastq_handles[head.name], read)
                read_counter[head.name] += 1
            else:
                n_conflict += 1
                if output_failed:
                    write_fastq(h_conflict, read)
        else:
            n_unknown += 1
            if output_failed:
                write_fastq(h_unknown, read)

    if h_too_short:
        h_too_short.close()
    if h_unknown:
        h_unknown.close()
    if h_conflict:
        h_conflict.close()
    for h in fastq_handles.values():
        h.close()
        
    # Report
    
    with open(os.path.join(outdir, "reads.tsv"), "w+") as fw:
        fw.write("Barcode\tCell\tReads\tRatio\n")
        for barcode, count in read_counter.items():
            cell = barcode2cell[barcode]
            r = count / n_total
            fw.write("%s\t%s\t%d\t%f\n" % (barcode, cell, count, r))
    
    with open(os.path.join(outdir, "stats.tsv"), "w+") as fw:
        r = n_assigned / n_total
        fw.write("File\tTotal\tTooShort\tUnknown\tConflict\tAssigned\tAssignedRatio\n")
        fw.write("%s\t%d\t%d\t%d\t%d\t%d\t%f\n" % (
            os.path.basename(f_fastq), n_total, n_too_short, 
            n_unknown, n_conflict, n_assigned, r))

    
if __name__ == '__main__':
    main()
    