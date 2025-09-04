#!/usr/bin/env runsnakemake
configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
SPECIES = config["SPECIES"]

def get_seqname_pattern(species):
    return config["%s_SEQNAME_PATTERN" % species.upper()]

def get_fasta(species):
    return config["%s_FASTA" % species.upper()]
    
def get_gtf(species):
    return config["%s_GTF" % species.upper()]

def get_transcript_bed(species):
    return config["%s_TRANSCRIPT_BED" % species.upper()]

def get_annotation_csv(species):
    return config["%s_ANNOTATION_CSV" % species.upper()]