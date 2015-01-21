import os, re
import sqlite3 as lite
from fogasrdb_api import *
from plotutils import *

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
            if codoncheck:
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
        #con.commit()
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
        return False
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
    #for ii in x:
    #write_log(con, "Restricting the analysis to the first 30 orthogroups only. Checkpoint 486.")    
    for cc in range(0, x.__len__()):
        ii = x[cc] 
        print "\n. Building ASR working folder for " + ii[0].__str__()
        setup_asr_analysis(con, ii[0])
        validate_asr_setup(con, ii[0])

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
    """A helper method for 'get_postdup_anc_string', 'get_predup_anc_string', and 'get_outgroup_string'"""
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

# sql = "select seqid from group_seq where groupid=" + orthogroupid.__str__() + " and seqid in (select seqid from nt_aa_check where checkval>0)"

    """Get the seqids for this orthogroup"""        
    sql = "select seqid from group_seq where groupid=" + orthogroupid.__str__() + " and seqid in (select seqid from nt_aa_check where checkval>0)"
    #sql = "select seqid from group_seq where groupid=" + orthogroupid.__str__()
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

    """Now check that all the amino acid and codon sequences match."""
    ntseqs = {}
    aaseqs = {}
    for seqid in seqids:
        sql = "select sequence from ntseqs where seqid=" + seqid.__str__()
        cur.execute(sql)
        ntsequence = cur.fetchone()[0]
        sql = "select sequence from aaseqs where seqid=" + seqid.__str__()
        cur.execute(sql)
        aasequence = cur.fetchone()[0]
        name = get_seqname(con, seqid)
        ntseqs[name] = ntsequence
        aaseqs[name] = aasequence
        
    bad_taxa = []
    for taxon in aaseqs:
        check = check_nt_vs_aa(con, aaseqs[taxon], ntseqs[taxon])
        if check == False:
            bad_taxa.append( taxon )


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
    fout.write("LAZARUS = python /common/REPOSITORY/lazarus/lazarus.py\n") 
    fout.write("MARKOV_MODEL_FOLDER = /common/REPOSITORY/paml44/dat\n")
    fout.write("ANCCOMP = python /Network/Servers/udp015817uds.ucsf.edu/Users/victor/Documents/EclipseWork/anccomp/compare_ancs.py\n")
    
    fout.write("FASTTREE = FastTree\n")
    
    fout.write("THRESHOLDS_ZORRO = 0.1 0.25 0.5 1.0\n")

    fout.write("USE_MPI = False\n")
    
    # This is how the O.S. should launch shell scripts:
    fout.write("RUN = sh\n")
    
    # In what directories should I expect to find multiple sequence alignments? 
    fout.write("ALIGNMENT_ALGORITHMS = mafft\n")
    
    fout.write("MSAPROBS = msaprobs\n")
    fout.write("MUSCLE = muscle\n")
    fout.write("MAFFT = mafft\n")
    fout.write("ZORRO = zorro\n")
    fout.write("CODEML = codeml\n")
    
    fout.write("MODELS_RAXML = PROTGAMMALG PROTCATLG \n")
    
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
    
def validate_asr_setup(con, orthogroupid):
    """This method validates the results of the method 'setup_asr_analysis'.
    If everything is OK, this method returns with no errors.
    If something is wrong, it exits the program."""
    cur = con.cursor()

    asrdir = "data/" + orthogroupid.__str__()
    if False == os.path.exists(  asrdir  ):
        write_error(con, "I cannot find the ASR working directory " + asrdir)
        exit()

    ntfasta = "data/" + orthogroupid.__str__() + "/" + orthogroupid.__str__() + ".nt.fasta"
    if False == os.path.exists(  ntfasta  ):
        write_error(con, "I cannot find the NT fasta file" + ntfasta)
        exit()
    
    aafasta = "data/" + orthogroupid.__str__() + "/" + orthogroupid.__str__() + ".aa.fasta"
    if False == os.path.exists(  aafasta  ):
        write_error(con, "I cannot find the AA fasta file" + aafasta)
        exit()
    
    configpath = "data/" + orthogroupid.__str__() + "/" + orthogroupid.__str__() + ".config"
    if False == os.path.exists(  configpath  ):
        write_error(con, "I cannot find the configuration file " + configpath)
        exit()
    
    """Now check that all the amino acid and codon sequences match."""
    ntseqs = read_fasta(ntfasta)
    aaseqs = read_fasta(aafasta)
    for taxon in ntseqs:
        if taxon not in aaseqs:
            write_error(con, "I cannot find taxon " + taxon + " in " + aafasta)
            exit()
    for taxon in aaseqs:
        if taxon not in ntseqs:
            write_error(con, "I cannot find taxon " + taxon + " in " + ntfasta)
            exit()
    for taxon in ntseqs:
        check = check_nt_vs_aa(con, aaseqs[taxon], ntseqs[taxon])
        if check == False:
            write_error(con, "The a.a. and n.t. sequences for taxon " + taxon + " do not match. Orthogroup ID " + orthogroupid.__str__())
            exit()   
    
      
def write_asr_scripts(con):
    cur = con.cursor()
    sql = "select groupid from wgd_groups"
    cur.execute(sql)
    x = cur.fetchall()
    commands = []

    #for ii in x:
    #write_log(con, "Restricting the analysis to the first 30 orthogroups only. Checkpoint 663.")    
    for cc in range(0, x.__len__()):
        ii = x[cc] 

        datadir = "data/" + ii[0].__str__()
        scriptpath = datadir + "/runme.sh"

        if os.path.exists(datadir):              
            fout = open(scriptpath, "w")
            fout.write("python ~/EclipseWork/asrpipelinev2/runme.py --configpath " + ii[0].__str__()  + ".config --skip_zorro\n")
            fout.close()
            commands.append("source " + datadir + "/runme.sh")
        else:
            write_error(con, "I cannot find the data directory for orthogroup " + ii[0] + " at " + datadir)
            exit()
            
        if False == os.path.exists(scriptpath):
            write_error(con, "I was unable to write the ASR script " + scriptpath)
            exit()
    
    fout = open("asr_commands.sh", "w")
    for c in commands:
        fout.write(c + "\n")
    fout.close()

def distribute_to_slaves(con):
    cur = con.cursor()
    
    datadir_slave = {}
    
    slaves = ["agassiz", "bago", "darwin", "ericsson"]
    
    sql = "select groupid from wgd_groups"
    cur.execute(sql)
    x = cur.fetchall()
    commands = []
    #write_log(con, "Restricting the analysis to the first 100 orthogroups only. Checkpoint 663.")    
    count_slave = -1
    for cc in range(0, x.__len__()):
        ii = x[cc] 
 
        datadir = "data/" + ii[0].__str__()
        scriptpath = datadir + "/runme.sh"
         
        """Copy the working folder to the slave node's tmp folder"""
        if os.path.exists(scriptpath):
            count_slave += 1
            datadir_slave[ ii[0].__str__() ] = count_slave%(slaves.__len__())
            c = "scp -r " + datadir + " " + slaves[ count_slave%(slaves.__len__()) ] + ":/tmp/"
            print c
            os.system(c)
    
    for ii in range(0, slaves.__len__()):
        fout = open("run." + slaves[ii] + ".sh", "w")
        for d in datadir_slave:
            if datadir_slave[d] == ii:
                fout.write("cd /tmp/" + d + "\n")
                """Run the ASR pipeline for this gene:"""
                fout.write("source runme.sh\n")
                """Each slave sends its data back to the master:"""
                fout.write("scp -r ../" + d + " 10.0.0.100:/Volumes/RAID/victor/fogasr/data\n")
                fout.write("cd -\n")
        fout.close()
        
        os.system("scp " + "run." + slaves[ii] + ".sh " + slaves[ii] + ":/tmp/runme.sh")

def launch_remote_slaves(con):
    slaves = ["agassiz", "bago", "darwin", "ericsson"]
    
    fout = open("hosts.txt", "w")
    fout.write("localhost slots=1\n")
    for s in slaves:
        fout.write(s + " slots=1\n")
    fout.close()
    
    fout = open("asr_commands.sh", "w")
    for s in slaves:
        fout.write("source /tmp/runme.sh\n")
    fout.close()    

    c = "mpirun -np 5 --machinefile hosts.txt /common/bin/mpi_dispatch asr_commands.sh"
    print c
    os.system(c)
 
def collect_from_slaves(con):
    slaves = ["agassiz", "bago", "darwin", "ericsson"]
    
    for s in slaves:
        c = "scp -r " + s + ":/tmp/* ./data/"
        os.system( c )
 
def validate_asr_output(con):
    cur = con.cursor()
    sql = "select groupid from wgd_groups"
    cur.execute(sql)
    x = cur.fetchall()
    
    count_good = 0
    for cc in range(0, x.__len__()):
        ii = x[cc]  
        datadir = "data/" + ii[0].__str__()
        
        """Look for the dN/dS vs. Df comparison"""
        comppath = datadir + "/compare_dnds_Df.txt"
        
        if False == os.path.exists(comppath):
            write_error(con, "I cannot find the dnds-vs.-df comparisons for orthogroup ID " + ii[0].__str__() )
        else:
            count_good += 1
            write_log(con, "OK. I found the dnds-vs.-df comparisons for orthogroup ID " + ii[0].__str__() )
            fin = open(comppath, "r")
            count = 0
            for l in fin.xreadlines():
                if l.__len__() > 1:
                    count += 1
            fin.close()
            print "\n. . . " + count.__str__() + " sites."
    write_log(con, "I found " + count_good.__str__() + " orthogroups with data, and " + (x.__len__()-count_good).__str__() + " without data.")

def read_all_dnds_df_comparisons(con):
    cur = con.cursor()
    sql = "select groupid from wgd_groups"
    cur.execute(sql)
    x = cur.fetchall()
    
    cat1points = []
    negcat1points = []
    cat2points = []
    cat3points = []
    cat23points = []
    cat23points_ancmu = []
    cat23points_pos = []
    cat23points_ancmupos = []
    
    dfpoints = []
    dfranks = []
    dfpoints_ancmu = [] # subset of points with ancestral change
    dfpoints_pos = [] # subset of points with significant positive selection
    dfpoints_ancmupos = [] # subset of points with both previous conditions.
    
    kpoints = []
    ppoints = []
    
    fout = open("outlier_sites.txt", "w")
    
    """ NOTE: we're restricting the analysis to the first 100 orthogroups only."""
    count_good = 0
    for ii in x:
        datadir = "data/" + ii[0].__str__()
        
        """Look for the dN/dS vs. Df comparison"""
        comppath = datadir + "/compare_dnds_Df.txt"
        
        if os.path.exists(comppath):
            count_good += 1
            print ". Reading " + comppath
            
            """FIrst pass - rank df, k, and p"""
            df_sites = []
            p_sites = []
            k_sites = []
            fin = open(comppath, "r")
            for l in fin.xreadlines():
                if l.__len__() < 5:
                    continue
                tokens = l.split()
                if tokens.__len__() < 14:
                    continue
                site = int( tokens[0] )
                df          = abs(  float( tokens[11]) )
                k           = float( tokens[12] )
                p           = float( tokens[13] )
                
                df_sites.append( (df,site)  )
                k_sites.append(  (k,site)   )
                p_sites.append(  (p,site)   )

            df_ranked = sorted(df_sites, reverse=True)
            k_ranked = sorted(k_sites, reverse=True)
            p_ranked = sorted(p_sites, reverse=True)
            
            df_site_rank = {}
            k_site_rank = {}
            p_site_rank = {}
            
            for index, tuple in enumerate(df_ranked):
                df_site_rank[ tuple[1] ] = index
            for index, tuple in enumerate(k_ranked):
                k_site_rank[ tuple[1] ] = index
            for index, tuple in enumerate(p_ranked):
                p_site_rank[ tuple[1] ] = index
            
            fin = open( comppath, "r" )
            for l in fin.xreadlines():
                if l.__len__() > 5:
                    tokens = l.split()
                    if tokens.__len__() >= 14:
                        site = int( tokens[0] )
                        nebcat1     = float( tokens[1] )
                        nebcat2     = float( tokens[2] )
                        nebcat3     = float( tokens[3] )
                        nebancmu    = int(tokens[4]) # did the ancestor mutate at this site?
                        nebsig      = int(tokens[5])
                        
                        bebcat1     = float(tokens[6])
                        bebcat2     = float(tokens[7])
                        bebcat3     = float(tokens[8])
                        bebancmu    = int(tokens[9])
                        if bebancmu < 1:
                            continue
                        
                        bebsig      = int(tokens[10]) # there is significant evidence of positive selection
                        
                        df          = abs(  float( tokens[11]) )
                        k           = float( tokens[12] )
                        p           = float( tokens[13] )
                        
                        """Restrict the analysis to orthogroups in which BEB was used."""
                        if bebcat1 > 0.0 or bebcat2 > 0.0 or bebcat3 > 0.0:
                            #if bebsig > 0 and bebancmu > 0:
                            #    dfpoints_ancmupos.append( df )
                            #    cat23points_ancmupos.append( max(bebcat2,bebcat3) )
                            #elif bebsig > 0:
                            #    dfpoints_pos.append( df )
                            #    cat23points_pos.append( max(bebcat2,bebcat3) )
#                             elif bebancmu > 0:
#                                 dfpoints_ancmu.append( df )
#                                 cat23points_ancmu.append( max(bebcat2,bebcat3) )
                            #else:
                            cat1points.append( bebcat1 )
                            negcat1points.append( 1.0 - bebcat1 )
                            cat2points.append( bebcat2 )
                            cat3points.append( bebcat3 )
                            cat23points.append( max(bebcat2,bebcat3) )
                            dfpoints.append(df)
                            dfrank = df_site_rank[site]
                            dfranks.append( dfrank )
                            kpoints.append(k)
                            ppoints.append(p)
                            
                            
                            is_outlier = 0
                            """Is this site an interesting outlier?"""
                            if df < 20.0 and bebsig > 0:
                                is_outlier = 1

                            if df < 20.0 and bebsig > 0:
                                is_outlier = 1
                                
                            if is_outlier > 0:
                                dbpath = "data/" + ii[0].__str__() + "/asr.db"
                                if False == os.path.exists(dbpath):
                                    write_error(con, "I cannot find the SQL database at " + dbpath)
                                    continue
                            
                                asrcon = lite.connect(dbpath)
                                asrcur = asrcon.cursor()
                                
                                print " * Special Case Type", is_outlier, "in", dbpath
                                
                                fout.write("===============================================================\n")
                                fout.write("comppath: " + comppath + "\n")
                                fout.write("line: " + l + "\n")
                                fout.write("site:" + site.__str__() + "\n")
                                
                                sql = "select id from DNDS_Models where name='Nsites_branch'"
                                asrcur.execute(sql)
                                zz = asrcur.fetchall()
                                if zz.__len__() == 0:
                                    write_log(con, "There are no DNDS_Models in the database, so I'm skipping the comparison of DNDS to Df.")
                                    return
                                nsites_id = zz[0][0]
                                
                                """Get the ancestral PP values for this site."""
                                sql = "select testid from FScore_Sites where site=" + site.__str__() + " and df=" + df.__str__()
                                asrcur.execute(sql)
                                zz = asrcur.fetchall()
                                testid = zz[0][0]
                                #print " \t Fscore testid=" + testid.__str__()
                                sql = "select almethod, phylomodel, ancid1, ancid2 from FScore_Tests where id=" + testid.__str__()
                                asrcur.execute(sql)
                                x = asrcur.fetchall()
                                almethod = x[0][0]
                                phylomodel = x[0][1]
                                ancid1 = x[0][2]
                                ancid2 = x[0][3]
                                sql = "select state, pp from AncestralStates where ancid=" + ancid1.__str__() + " and site=" + site.__str__()
                                asrcur.execute(sql)
                                aaa = asrcur.fetchall()
                                sql = "select state, pp from AncestralStates where ancid=" + ancid2.__str__() + " and site=" + site.__str__()
                                asrcur.execute(sql)
                                bbb = asrcur.fetchall()
                                print "\t Ancestor 1:" + aaa.__str__()
                                print "\t Ancestor 2:" + bbb.__str__()
                                
                                """Get some summary of the codon variation at this site."""
                                sql = "select id from DNDS_Tests where almethod=" + almethod.__str__()
                                sql += " and phylomodel=" + phylomodel.__str__() 
                                sql += " and anc1=" + ancid1.__str__() 
                                sql += " and anc2=" + ancid2.__str__()
                                sql += " and dnds_model=" + nsites_id.__str__()
                                asrcur.execute(sql)
                                testid = asrcur.fetchone()[0]
                                sql = "select * from DNDS_params where testid=" + testid.__str__()
                                asrcur.execute(sql)
                                x = asrcur.fetchall()
                                print "\t DNDS:", x[0][1:]
                                sql = "select * from NEB_scores where testid=" + testid.__str__() + " and site=" + site.__str__()
                                asrcur.execute(sql)
                                x = asrcur.fetchall()
                                nebppcat1 = x[0][2]
                                nebppcat2 = x[0][3]
                                nebppcat3 = x[0][4]
                                nebppcat4 = x[0][5]
                                nebsignificant = x[0][7]
                                print "\t NEB probabilities:", nebppcat1, nebppcat2, nebppcat3, nebppcat4, " -- sig:", nebsignificant
                                sql = "select * from BEB_scores where testid=" + testid.__str__() + " and site=" + site.__str__()
                                asrcur.execute(sql)
                                x = asrcur.fetchall()
                                bebppcat1 = x[0][2]
                                bebppcat2 = x[0][3]
                                bebppcat3 = x[0][4]
                                bebppcat4 = x[0][5]
                                bebsignificant = x[0][7]
                                print "\t BEB probabilities:", bebppcat1, bebppcat2, bebppcat3, bebppcat4, " -- sig:", bebsignificant

                                sql = "select taxonid, alsequence from AlignedSequences where almethod=" + almethod.__str__() + " and datatype=1"
                                asrcur.execute(sql)
                                x = asrcur.fetchall()
                                
                                taxon_aa = {}
                                taxon_codon = {}
                                for ii in x:
                                    aachar = ii[1][site-1]
                                    sql = "select shortname from Taxa where id=" + ii[0].__str__()
                                    asrcur.execute(sql)
                                    taxonname = asrcur.fetchone()[0]
                                    taxon_aa[taxonname] = aachar
                                
                                sql = "select taxonid, alsequence from AlignedSequences where almethod=" + almethod.__str__() + " and datatype=0"
                                asrcur.execute(sql)
                                x = asrcur.fetchall()
                                for ii in x:
                                    codon = ii[1][ (site-1)*3:((site-1)*3)+3 ]
                                    sql = "select shortname from Taxa where id=" + ii[0].__str__()
                                    asrcur.execute(sql)
                                    taxonname = asrcur.fetchone()[0]
                                    asrcur.execute(sql)
                                    taxon_codon[taxonname] = codon                                
                                #print "1079: taxon_aa:", taxon_aa
                                #print "1080: taxon_codon", taxon_codon
                                
                                for taxonname in taxon_aa:
                                    if taxonname not in taxon_aa:
                                        "\t\t", taxonname, " has no a.a."
                                    elif taxonname not in taxon_codon:
                                        "\t\t", taxonname, " has no codon"
                                    else:
                                        print taxon_aa[ taxonname ], taxon_codon[ taxonname ], taxonname                                
            fin.close()
    
    print "\n. I found " + count_good.__str__() + " good orthogroups."
    
    #scatter1(cat23points, dfpoints, xlab="Probability of Positive Selection", ylab="Df")
    
    if dfpoints.__len__() > 0:
        #xsets = [cat23points, cat23points_pos, cat23points_ancmu, cat23points_ancmupos]
        #ysets = [dfpoints,dfpoints_pos, dfpoints_ancmu,dfpoints_ancmupos]
        xsets = [cat23points, cat23points_pos, cat23points_ancmupos]
        ysets = [dfranks,dfpoints_pos,dfpoints_ancmupos]
        scatter1("compare_test", xsets, ysets, xlab="Prob. of Positive Selection", ymin=0, ylab="Functional Score", format="jpeg")
        scatter1("compare_test2_zoom", xsets, ysets, xmin=0.8, logx=False, ymin=0, xlab="Prob. of Positive Selection", ylab="Functional Score", format="jpeg")
    else:
        write_error(con, "I didn't find any data; skipping the scatterplot.")
    
#     scatter_nxm(2, 2, [cat23points, dfpoints], ["P(w>1)", "Df"], "compare_cat23_df", title="dN/dS versus dF", xlab="Probability of Positive Selection", ylab="dF Score", force_square=False, plot_as_rank = [], skip_identity = False, skip_zeros = False)
#     scatter_nxm(2, 2, [cat3points, dfpoints], ["P(cat3)", "Df"], "compare_cat3_df", title="dN/dS versus dF", xlab="Probability of Positive Selection", ylab="dF Score", force_square=False, plot_as_rank = [], skip_identity = False, skip_zeros = False)
#     scatter_nxm(2, 2, [cat2points, dfpoints], ["P(cat2)", "Df"], "compare_cat2_df", title="dN/dS versus dF", xlab="Probability of Positive Selection", ylab="dF Score", force_square=False, plot_as_rank = [], skip_identity = False, skip_zeros = False)
#     scatter_nxm(2, 2, [negcat1points, dfpoints], ["1-P(cat1)", "Df"], "compare_neg1_df", title="dN/dS versus dF", xlab="1.0 - Probability of Positive Selection", ylab="dF Score", force_square=False, plot_as_rank = [], skip_identity = False, skip_zeros = False)
                                  
        
                

def check_again_wgdgroups(con):
    cur = con.cursor()
    groupids = []
    sql = "select groupid from wgd_groups"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        groupids.append( ii[0] )
    
    countgroups = 0
    totalgroups = groupids.__len__()
    for orthogroupid in groupids:
        """Get a list of sequence IDs that are in this orthogroup, AND which can be validated across their aa and nt sequences."""
        sql = "select seqid from group_seq where groupid=" + orthogroupid.__str__() + " and seqid in (select seqid from nt_aa_check where checkval>0)"
        cur.execute(sql)
        x = cur.fetchall()
        seqids = []
        for ii in x:
            seqids.append( ii[0] )
    
        """Now check that all the amino acid and codon sequences match."""
        ntseqs = {}
        aaseqs = {}
        seqname_id = {}
        for seqid in seqids:
            sql = "select sequence from ntseqs where seqid=" + seqid.__str__()
            cur.execute(sql)
            ntsequence = cur.fetchone()[0]
            sql = "select sequence from aaseqs where seqid=" + seqid.__str__()
            cur.execute(sql)
            aasequence = cur.fetchone()[0]
            name = get_seqname(con, seqid)
            seqname_id[name] = seqid
            ntseqs[name] = ntsequence
            aaseqs[name] = aasequence
            
        bad_taxa = []
        for taxon in aaseqs:
            check = check_nt_vs_aa(con, aaseqs[taxon], ntseqs[taxon])
            if check == False:
                bad_taxa.append( taxon )
                write_log(con, "Removing taxon " + taxon + " from wgd_group " + orthogroupid.__str__() + ". Checkping 824.")
        
        for bt in bad_taxa:
            aaseqs.pop(bt)
            ntseqs.pop(bt)
            sql = "delete from nt_aa_check where seqid=" + seqname_id[bt].__str__()
            cur.execute(sql)
        con.commit()
        
        countgroups += 1
        if countgroups%10 == 0:  
            sys.stdout.write("\r    --> %.1f%%"% (100*countgroups/float(totalgroups)) )
            sys.stdout.flush()
            
        

        
        
