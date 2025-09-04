configfile: "config.yaml"

def get_species(sample):
    if "K562" in sample:
        return "Human"
    elif "Mix" in sample:
        return "Mix"
    else:
        return "Mouse"

def get_genome_fasta(sample):
    return config["%s_GENOME_FASTA" % get_species(sample).upper()]

def get_annotation_gtf(sample):
    return config["%s_ANNOTATION_GTF" % get_species(sample).upper()]

def get_star_genome(sample):
    return config["%s_STAR_GENOME" % get_species(sample).upper()]
