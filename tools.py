import re
import sqlite3 as lite
from fogasrdb_api import *

def read_fastalines(lines):
    """Returns a hash: key = taxon name, value = sequence from the fasta file.
        Ignores redundant taxa."""
    taxanames = []
    taxa_seq = {}
    last_taxa = None
    last_seq = ""
    okay = True
    for l in lines:
        l = l.strip()
        if l.__len__() <= 2:
            pass
        elif l.startswith(">"):
            okay = True
            taxaname = re.sub(">", "", l.split()[0] )
            """Remove underscores"""
            taxaname = re.sub("_", "", taxaname)

            if taxaname in taxa_seq:
                okay = False
            else:
                if last_taxa != None:
                    taxa_seq[last_taxa] = last_seq
                last_taxa = taxaname
                last_seq = ""
                
        elif okay == True:
            """Remove stop codons"""
            l = re.sub("\*", "", l )
            last_seq += l
    if last_taxa != None:
        taxa_seq[last_taxa] = last_seq
    
    return taxa_seq

def read_fasta(fpath):
    """Returns a hash: key = taxon name, value = sequence from the fasta file.
        Ignores redundant taxa."""
    fin = open(fpath, "r")    
    taxanames = []
    taxa_seq = {}
    last_taxa = None
    last_seq = ""
    okay = True
    lines = fin.readlines()
    fin.close()
    return read_fastalines( lines )


def get_taxon_urls():                          
    sock = urllib.urlopen("http://www.broadinstitute.org/regev/orthogroups/sources.html") 
    htmlSource = sock.read()                            
    sock.close()                                        
    
    taxon_aaurl = {}
    taxon_nturl = {}
    
    for l in htmlSource.split("\n"):
        taxon = ""
        aaurl = ""
        nturl = ""
        if l.__contains__("AA</a>"):
            """Taxon Name"""
            x = l.find("<i>")
            if x != -1:
                start = x+3
                for c in l[start:]:
                    if c == "<":
                        break
                    else:
                        taxon += c
                        
            """AA URL"""
            x = l.find("href=aa/")
            if x != -1:
                start = x+5
                for c in l[start:]:
                    if c == ">":
                        break
                    else:
                        aaurl += c
                                    
            """NT URL"""
            x = l.find("href=nt/")
            if x != -1:
                start = x+5
                for c in l[start:]:
                    if c == ">":
                        break
                    else:
                        nturl += c             
            print taxon, aaurl, nturl
            if taxon == "C. lusitaniae":
                taxon_aaurl[taxon] = "aa/Clus.fasta"
                taxon_nturl[taxon] = "nt/Clus.fasta"
            elif taxon == "A. nidulans":
                taxon_aaurl[taxon] = "aa/Anid.fasta"
                taxon_nturl[taxon] = "nt/Anid.fasta"                
            elif taxon == "N. crassa":
                taxon_aaurl[taxon] = "aa/Ncra.fasta"
                taxon_nturl[taxon] = "nt/Ncra.fasta" 
            elif taxon.__len__() > 1 and aaurl.__len__() > 1 and nturl.__len__() > 1:
                taxon_aaurl[taxon] = aaurl
                taxon_nturl[taxon] = nturl
    return (taxon_aaurl, taxon_nturl)
        
def get_fasta_source(url, con):
    sock = urllib.urlopen("http://www.broadinstitute.org/regev/orthogroups/" + url) 
    htmlSource = sock.read()                            
    sock.close() 
    taxa_seq = read_fastalines( htmlSource.split() )
    write_log(con, "I found " + taxa_seq.__len__().__str__() + " sequences in " + url)
    return taxa_seq

def import_sequences(con):
    cur = con.cursor()
    cur.execute("delete from species")
    con.commit()
    cur.execute("delete from seqnames")
    con.commit()
    cur.execute("delete from aaseqs")
    con.commit()
    cur.execute("delete from ntseqs")
    con.commit()
    cur.execute("delete from nt_aa_check")
    con.commit()    
    
    (taxon_aaurl, taxon_nturl) = get_taxon_urls()
    print taxon_aaurl
    print taxon_nturl
    for taxon in taxon_aaurl:
        sql = "select count(*) from species where name='" + taxon + "'"
        cur.execute(sql)
        count = cur.fetchone()[0]
        if count == 0:
            sql = "insert into species (name) values('" + taxon + "')"
            cur.execute(sql)
            con.commit()    
        sql = "select id from species where name='" + taxon + "'"
        cur.execute(sql)
        speciesid = cur.fetchone()[0]
        
        print taxon, taxon_aaurl[taxon]
        
        """aa"""
        taxa_aaseq = get_fasta_source( taxon_aaurl[taxon], con)
        total_count = taxa_aaseq.__len__()
        pcount = 0
        for t in taxa_aaseq:
            pcount += 1
            if pcount%10 == 0:
                sys.stdout.write("\r" + taxon + " aa:    --> %.1f%%" % (100*pcount/float(total_count)) )
                sys.stdout.flush()
            
            sql = "select count(*) from seqnames where speciesid=" + speciesid.__str__() + " and name='" + t + "'"
            cur.execute(sql)
            count = cur.fetchone()[0]
            if count == 0:
                sql = "insert into seqnames (speciesid, name) values(" + speciesid.__str__() + ",'" + t + "')"
                cur.execute(sql)
                con.commit()
            else:
                write_log("Warning: I found a previous copy of an AA sequence named " + t)
            
            sql = "select id from seqnames where speciesid=" + speciesid.__str__() + " and name='" + t + "'"
            cur.execute(sql)
            seqid = cur.fetchone()[0]
            
            sql = "insert into aaseqs (speciesid, sequence, seqid) values("
            sql += speciesid.__str__() + ",'" + taxa_aaseq[t] + "'," + seqid.__str__() + ")"
            cur.execute(sql)
        
        con.commit()
        print ""
        
        """nt"""
        taxa_ntseq = get_fasta_source( taxon_nturl[taxon], con)
        total_count = taxa_ntseq.__len__()
        pcount = 0
        for t in taxa_ntseq:
            this_seq = taxa_ntseq[t]
            
            pcount += 1
            if pcount%10 == 0:
                sys.stdout.write("\r" + taxon + " nt:    --> %.1f%%" % (100*pcount/float(total_count)) )
                sys.stdout.flush()
            
            
            sql = "select count(*) from seqnames where name='" + t + "' and speciesid=" + speciesid.__str__()
            cur.execute(sql)
            count = cur.fetchone()[0]
            if count == 0:            
                """We're expecting to already have seen this sequence name in the AA data"""
                write_log(con, "I can't find the sequence " + t + ".")
                sql = "delete from seqnames where name='" + t + "' and speciesid=" + speciesid.__str__()
                cur.execute(sql)
                con.commit()                
                continue
            
            """Does the aa sequence align with the nt sequence?
                --> we'll import regardless of the answer, but we'll remember which sequences are good."""
            codoncheck = check_nt_vs_aa(con, taxa_aaseq[t], this_seq)
            check = 0
            if codoncheck != None:
                check = 1
            else:
                from Bio.Seq import Seq
                from Bio.Alphabet import generic_dna, generic_protein
                trans = Seq(this_seq, generic_dna).translate()
                print "Error 200:", t
                print "aa:", taxa_aaseq[t]
                print "trans:", trans
                print "nt:", this_seq
                
            sql = "insert into nt_aa_check (seqid, checkval) values(" + seqid.__str__() + "," + check.__str__() + ")"
            cur.execute(sql) 
            
            
            sql = "select id from seqnames where speciesid=" + speciesid.__str__() + " and name='" + t + "'"
            cur.execute(sql)
            seqid = cur.fetchone()[0]
            
            sql = "insert into ntseqs (speciesid, sequence, seqid) values("
            sql += speciesid.__str__() + ",'" + this_seq + "'," + seqid.__str__() + ")"
            cur.execute(sql)            
        con.commit()
        print ""
        
        """Ensure all the aa sequences have a companion nt sequence."""        
        sql = "delete from aaseqs where speciesid=" + speciesid.__str__() + " and seqid not in (select seqid from ntseqs where speciesid=" + speciesid.__str__() + ")"
        cur.execute(sql)
        con.commit()
        
        """Remove names of sequences which lack an AA and/or NT sequence."""
        sql = "delete from seqnames where id not in (select seqid from ntseqs) or id not in (select seqid from aaseqs)"
        cur.execute(sql)
        con.commit()
        
        """Summarize"""
        sql = "select count(*) from seqnames where speciesid=" + speciesid.__str__()
        cur.execute(sql)
        countall = cur.fetchone()[0]        
        
        sql = "select count(*) from seqnames where speciesid=" + speciesid.__str__() + " and id in (select seqid from nt_aa_check where checkval=1)"
        cur.execute(sql)
        countver = cur.fetchone()[0]
        
        write_log(con, "Species " + taxon + " has " + countall.__str__() + " gene sequences, of which " + countver.__str__() + " can be matched between nt and aa.")

def verify_sequences(con):
    cur = con.cursor()
    
    print "\n. Verifying all imported sequences..."
    
    sql = "select id from species"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        speciesid = ii[0]
        speciesname = get_species_name(con, speciesid)

        sql = "select count(*) from seqnames where speciesid=" + speciesid.__str__()
        cur.execute(sql)
        countall = cur.fetchone()[0]        
        
        sql = "select count(*) from seqnames where speciesid=" + speciesid.__str__() + " and id in (select seqid from nt_aa_check where checkval=1)"
        cur.execute(sql)
        countver = cur.fetchone()[0]
        
        write_log(con, "Species " + speciesname + " has " + countall.__str__() + " gene sequences, of which " + countver.__str__() + " can be matched between nt and aa.")


def check_nt_vs_aa(con, aaseq, codonseq):
    """Maps the codon sequence to the aligned (may contain indels) aa seq."""
        
    """Quick sanity check: do we have exactly 3x more nucleotides than amino acids?"""
    aa_no_indels = re.sub("-", "", aaseq)
    nt_no_indels = re.sub("-", "", codonseq)
    
    """Remove stop codon in the nt sequence."""
    stop_codons = ["TAG", "TAA", "TGA", "tag", "taa", "tga"]
    for sc in stop_codons:
        if nt_no_indels.endswith(sc):
            nt_no_indels = nt_no_indels[0:  nt_no_indels.__len__()-3 ]
    
    if float( aa_no_indels.__len__() ) != float(nt_no_indels.__len__())/3.0:
        print aa_no_indels.__len__().__str__(), " versus ", nt_no_indels.__len__().__str__()
        
        return None
        
    
    return True     

def import_orthogroups(con):
    cur = con.cursor()
    sql = "delete from orthogroups"
    cur.execute(sql)
    con.commit()
    sql = "delete from group_seq"
    cur.execute(sql)
    con.commit()
    
    """Download the orthogroups data from the Broad"""
    all_data_url = "http://www.broadinstitute.org/regev/orthogroups/all-output.txt"
    print "\n. Downloading all the data from ", all_data_url
    sock = urllib.urlopen(all_data_url) 
    htmlSource = sock.read()                            
    sock.close()  
    
    lines = htmlSource.split("\n")
    prog_count = 0
    total_count = lines.__len__()
    
    lcounter = 0
    for l in lines:
        if l.__contains__(">"):
            lcounter = 0
        else:
            lcounter += 1
        
        prog_count += 1
        if prog_count%10 == 0:
            sys.stdout.write("\r    --> %.1f%%" % (100*prog_count/float(total_count)) )
            sys.stdout.flush()
        
        
        """Process the lines that contain orthogroup definitions."""
        if lcounter == 1:
            tokens = l.split()
            groupnumber = int( tokens[0] )
            items = []
            items = tokens[3:]
            
            """Ignore singleton genes"""
            if items.__len__() < 2:
                continue
            
            #print "231: items:", items
            
            sql = "insert into orthogroups (name) values('" + groupnumber.__str__() + "')"
            cur.execute(sql)
            con.commit()
            sql = "select id from orthogroups where name='" + groupnumber.__str__() + "'"
            cur.execute(sql)
            groupid = cur.fetchone()[0]
            
            
            for item in items:
                """Each item is formatted as Clus|GeneName, for example"""
                species = item.split("|")[0]

                """Check the first part before the | char."""
                species_frag = species[0].upper() + ". " + species[1:]
                speciesid = get_species_id(con, species_frag)
                if speciesid == None:
                    write_error("I cannot find the species ID for " + species_frag)
                    exit()
                
                """Check the second part after the | char"""                
                seqname = item.split("|")[1]
                seqid = get_seqid(con, seqname, speciesid)                
                if seqid == None:
                    write_log(con, "I cannot find the gene ID for " + seqname + " in species " + species_frag)
                else:
                    sql = "insert into group_seq (groupid, seqid) values(" + groupid.__str__() + "," + seqid.__str__() + ")"
                    cur.execute(sql)
            con.commit()
    
    """How many groups did we find?"""
    sql = "select count(*) from orthogroups"
    cur.execute(sql)
    count = cur.fetchone()[0]
    write_log(con, "I found " + count.__str__() + " non-singleton orthogroups.")

def filter_orthogroups(con):
    cur = con.cursor()
    
    sql = "delete from wgd_groups"
    cur.execute(sql)
    
    groupids = []
    sql = "select id from orthogroups"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        groupids.append( ii[0] )
    
    """Get a list of IDs for species in the cerevisiae clade."""
    
    scer_names = ["S. cerevisiae", "S.paradoxus", "S. mikatae", "S. bayanus", "C. glabrata", "S. castellii"]
    klac_names = ["K. waltii", "K. lactis", "S. kluyveri", "A. gossypii"]
    outgroup_names = ["Y. lipolytica, C. lusitaniae", "D. hansenii", "C. guilliermondii", "C. tropicalis", "C. albicans", "C. parapsilosis", "L. elongisporus"]
    
    sql = "select id from species where name="
    frag = "' or name='".join( scer_names )
    frag = "'" + frag + "'"
    sql += frag
    cur.execute(sql)
    x = cur.fetchall()
    scer_ids = []
    for ii in x:
        scer_ids.append( ii[0] )
    
    sql = "select id from species where name="
    frag = "' or name='".join( klac_names )
    frag = "'" + frag + "'"
    sql += frag
    cur.execute(sql)
    x = cur.fetchall()
    klac_ids = []
    for ii in x:
        klac_ids.append( ii[0] )
    
    sql = "select id from species where name="
    frag = "' or name='".join( outgroup_names )
    frag = "'" + frag + "'"
    sql += frag
    cur.execute(sql)
    x = cur.fetchall()
    outgroup_ids = []
    for ii in x:
        outgroup_ids.append( ii[0] )   
        
    
    for gid in groupids:
        """Get the sequences in this group for which we have both nt and aa versions."""
        sql = "select seqid from group_seq where groupid=" + gid.__str__() + " and seqid in (select seqid from nt_aa_check where checkval>0)"
        cur.execute(sql)
        x = cur.fetchall()
        seqids = [] # seqids is IDs for sequences that can be aligned between nt and aa.
        for ii in x:
            seqids.append( ii[0] )
        
        """Get the species' IDs of these sequences."""
        speciesids = []
        for sid in seqids:
            sql = "select speciesid from seqnames where id=" + sid.__str__()
            cur.execute(sql)
            id = cur.fetchone()[0]
            speciesids.append( id )
        
        """Does this group contain at least entry for the S.cer clade?"""
        foundscer = 0
        for ii in scer_ids:
            if ii in speciesids:
                foundscer += 1
    
        foundklac = 0
        for ii in klac_ids:
            if ii in speciesids:
                foundklac += 1

        foundout = 0
        for ii in outgroup_ids:
            if ii in speciesids:
                foundout += 1
        
        if foundscer>3 and foundklac>3 and foundout>3:
            sql = "insert into wgd_groups (groupid, depth) values(" + gid.__str__() + ", 3)"
            cur.execute(sql)
        elif foundscer>2 and foundklac>2 and foundout>2:
            sql = "insert into wgd_groups (groupid, depth) values(" + gid.__str__() + ", 2)"
            cur.execute(sql)
        elif foundscer>0 and foundklac>0 and foundout>0:
            sql = "insert into wgd_groups (groupid, depth) values(" + gid.__str__() + ", 1)"
            cur.execute(sql)

    con.commit()
    sql = "select count(*) from wgd_groups"
    cur.execute(sql)
    count = cur.fetchone()[0]
    write_log(con, "I found " + count.__str__() + " orthogroups that are useful for studying the whole-genome duplication.")

def setup_all_asr(con):
    cur = con.cursor()
    sql = "select groupid from wgd_groups"
    cur.execute(sql)
    x = cur.fetchall()
    os.system("mkdir data")
    for ii in x:
        print "\n. Building ASR working folder for " + ii[0].__str__()
        setup_asr_analysis(con, ii[0])

def get_gene_family_name(con, orthogroupid):
    cur = con.cursor()
    scer_names = ["S. cerevisiae", "S.paradoxus", "S. mikatae", "S. bayanus", "C. glabrata", "S. castellii"]
    
    """Get a list of species IDs for Scer clade species."""
    sql = "select id from species where name="
    frag = "' or name='".join( scer_names )
    frag = "'" + frag + "'"
    sql += frag
    cur.execute(sql)
    x = cur.fetchall()
    scer_ids = []
    for ii in x:
        scer_ids.append( ii[0] )

    """Get the seqids for this orthogroup"""        
    sql = "select seqid from group_seq where groupid=" + orthogroupid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        for scerid in scer_ids:
            seqid = ii[0]
            sql = "select name from seqnames where id=" + seqid.__str__() + " and speciesid=" + scerid.__str__()
            cur.execute(sql) 
            y = cur.fetchall()
            if y.__len__() > 0:
                return y[0][0]
            
def get_outgroup_string(con, orthogroupid):
    cur = con.cursor()
    names = ["Y. lipolytica, C. lusitaniae", "D. hansenii", "C. guilliermondii", "C. tropicalis", "C. albicans", "C. parapsilosis", "L. elongisporus"]
    return get_group_string(con, orthogroupid, names)

def get_predup_anc_string(con, orthogroupid):
    cur = con.cursor()
    names = ["S. cerevisiae", "S.paradoxus", "S. mikatae", "S. bayanus", "C. glabrata", "S. castellii","K. waltii", "K. lactis", "S. kluyveri", "A. gossypii"]
    return get_group_string(con, orthogroupid, names)

def get_postdup_anc_string(con, orthogroupid):
    cur = con.cursor()
    names = ["S. cerevisiae", "S.paradoxus", "S. mikatae", "S. bayanus", "C. glabrata", "S. castellii"]
    return get_group_string(con, orthogroupid, names)

def get_group_string(con, orthogroupid, names):
    """A helper method."""
    cur = con.cursor()
    """Get a list of species IDs for Scer clade species."""
    sql = "select id from species where name="
    frag = "' or name='".join( names )
    frag = "'" + frag + "'"
    sql += frag
    cur.execute(sql)
    x = cur.fetchall()
    ids = []
    for ii in x:
        ids.append( ii[0] )

    """Get the seqids for this orthogroup"""        
    sql = "select seqid from group_seq where groupid=" + orthogroupid.__str__()
    cur.execute(sql)
    items = []
    x = cur.fetchall()
    for ii in x:
        for id in ids:
            seqid = ii[0]
            sql = "select name from seqnames where id=" + seqid.__str__() + " and speciesid=" + id.__str__()
            cur.execute(sql) 
            y = cur.fetchall()
            if y.__len__() > 0:
                items.append( y[0][0] )
    return "[" + ",".join( items ) + "]"


def setup_asr_analysis(con, orthogroupid):
    """Setup a folder to do the ASR pipeline analysis for one orthogroup."""    
    cur = con.cursor()
    
    """Get a list of sequence IDs that are in this orthogroup, AND which can be validated across their aa and nt sequences."""
    sql = "select seqid from group_seq where groupid=" + orthogroupid.__str__() + " and seqid in (select seqid from nt_aa_check where checkval>0)"
    cur.execute(sql)
    x = cur.fetchall()
    seqids = []
    for ii in x:
        seqids.append( ii[0] )

    os.system("mkdir data/" + orthogroupid.__str__())

    """Write nt FASTA"""
    fout = open("data/" + orthogroupid.__str__() + "/" + orthogroupid.__str__() + ".nt.fasta", "w")    
    for seqid in seqids:
        sql = "select sequence from ntseqs where seqid=" + seqid.__str__()
        cur.execute(sql)
        sequence = cur.fetchone()[0]
        name = get_seqname(con, seqid)
        fout.write(">" + name + "\n")
        fout.write(sequence + "\n")
    fout.close()

    """Write aa FASTA"""
    fout = open("data/" + orthogroupid.__str__() + "/" + orthogroupid.__str__() + ".aa.fasta", "w")    
    for seqid in seqids:
        sql = "select sequence from aaseqs where seqid=" + seqid.__str__()
        cur.execute(sql)
        sequence = cur.fetchone()[0]
        name = get_seqname(con, seqid)
        fout.write(">" + name + "\n")
        fout.write(sequence + "\n")
    fout.close()   
    
    """Write ASR pipeline configuration file."""
    fout = open("data/" + orthogroupid.__str__() + "/" + orthogroupid.__str__() + ".config", "w")    
    fout.write("GENE_ID = " + orthogroupid.__str__() + "\n")
    project_name = get_gene_family_name(con, orthogroupid)
    fout.write("PROJECT_TITLE = " + project_name.__str__() + "\n")
    
    fout.write("SEQUENCES = " + orthogroupid.__str__() + ".aa.fasta\n")
    fout.write("NTFASTA = " + orthogroupid.__str__() + ".nt.fasta\n")
    
    fout.write("RAXML = raxml -T 2\n")
    fout.write("PHYML = phyml\n")
    fout.write("LAZARUS = python /common/lazarus/lazarus.py\n") 
    fout.write("MARKOV_MODEL_FOLDER = /common/lazarus/paml/dat/\n")
    fout.write("ANCCOMP = python ~/Documents/SourceCode/anccomp/compare_ancs.py\n")
    
    fout.write("FASTTREE = FastTree\n")
    
    fout.write("THRESHOLDS_ZORRO = 0.1 0.25 0.5 1.0\n")

    fout.write("USE_MPI = False\n")
    
    # This is how the O.S. should launch shell scripts:
    fout.write("RUN = sh\n")
    
    # In what directories should I expect to find multiple sequence alignments? 
    fout.write("ALIGNMENT_ALGORITHMS = muscle mafft\n")
    
    fout.write("MSAPROBS = msaprobs\n")
    fout.write("MUSCLE = muscle\n")
    fout.write("MAFFT = mafft\n")
    fout.write("ZORRO = zorro\n")
    
    fout.write("MODELS_RAXML = PROTGAMMALG PROTGAMMAWAG PROTCATLG PROTCATWAG\n")
    
    fout.write("SEED_MOTIF_TAXA = " + project_name + "\n")
    
    fout.write("N_BAYES_SAMPLES = 10\n")
    
    ogs = get_outgroup_string(con, orthogroupid)
    fout.write("OUTGROUP = " + ogs + "\n")
    
    fout.write("ANCESTORS = AncPredup AncPostdup\n")
    ancs = get_predup_anc_string(con, orthogroupid)
    fout.write("INGROUP AncPredup " + ancs + "\n")
    ancs = get_postdup_anc_string(con, orthogroupid)
    fout.write("INGROUP AncPostdup " + ancs + "\n")
    
    fout.write("ASRSEED AncPredup " + project_name + "\n")
    fout.write("ASRSEED AncPostdup " + project_name + "\n")
    
    fout.write("COMPARE AncPredup AncPostdup\n")
    
    fout.close()