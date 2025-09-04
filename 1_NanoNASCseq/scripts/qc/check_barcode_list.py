#!/usr/bin/env python
import sys
import os
import subprocess
import gzip
import glob
import shutil


def main():
    f_fastq, f_fasta, threads, barcodes, f_out = sys.argv[1:]
    threads = int(threads)
    barcodes = barcodes.split(",")
    
    tmpdir = f_out + ".TMP"
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
        
    f_fastq_downsample = os.path.join(tmpdir, "downsample.fastq")        
    with gzip.open(f_fastq, "rt") as f, open(f_fastq_downsample, "w+") as fw:
        for i, line in enumerate(f):
            if i >= 400000: # 100k reads
                break
            fw.write(line)
                
    f_nanoplexer_output = os.path.join(tmpdir, "nanoplexer_output")
    cmd = "nanoplexer -t %d -b %s -p %s %s" % (threads, f_fasta, f_nanoplexer_output, f_fastq_downsample)
    subprocess.check_call(cmd, shell=True)
    with open(f_out, "w+") as fw:
        fw.write("File\tBarcode\tReads\tStatus\n")
        for path in glob.glob(os.path.join(f_nanoplexer_output, "*.fastq")):
            barcode = os.path.basename(path)[:-6]
            if barcode == "unclassified":
                continue
            n = 0
            with open(path) as f:
                for line in f:
                    n += 1
            reads = int(n / 4)
            if reads < 100:
                if barcode in barcodes:
                    t = "LowReads"
                else:
                    t = "OK"
            elif reads >= 500:
                if barcode in barcodes:
                    t = "OK"
                else:
                    t = "NotInBarcodeList"
            else:
                if barcode in barcodes:
                    t = "Warn-LowReads"
                else:
                    t = "Warn-NotInBarcodeList"
            fw.write("%s\t%s\t%s\t%s\n" % (os.path.basename(f_fastq), barcode, reads, t))
    
    shutil.rmtree(tmpdir)

    
if __name__ == '__main__':
    main()
    