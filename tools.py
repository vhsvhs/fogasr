def read_fasta(fpath):
    """Returns a hash: key = taxon name, value = sequence from the fasta file.
        Ignores redundant taxa."""
    fin = open(fpath, "r")    
    taxanames = []
    taxa_seq = {}
    last_taxa = None
    last_seq = ""
    okay = True
    for l in fin.xreadlines():
        l = l.strip()
        if l.__len__() <= 2:
            pass
        elif l.startswith(">"):
            okay = True
            taxaname = re.sub(">", "", l.split()[0] )
            if taxaname in taxa_seq:
                okay = False
            else:
                if last_taxa != None:
                    taxa_seq[last_taxa] = last_seq
                last_taxa = taxaname
                last_seq = ""
                
        elif okay == True:
            last_seq += l
    if last_taxa != None:
        taxa_seq[last_taxa] = last_seq
    
    fin.close()
    
    return taxa_seq