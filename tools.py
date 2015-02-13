import os, re
import sqlite3 as lite
from fogasrdb_api import *
from plotutils import *
import time

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


def setup_all_asr(con, new_runmes_only=False):
    cur = con.cursor()
    sql = "select groupid from wgd_groups"
    cur.execute(sql)
    x = cur.fetchall()
    os.system("mkdir data")
    for cc in range(0, x.__len__() ):
        ii = x[cc] 
        print "\n. Building ASR working folder for " + ii[0].__str__()
        setup_asr_analysis(con, ii[0], new_runmes_only=new_runmes_only)
        if new_runmes_only==False:
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


def setup_asr_analysis(con, orthogroupid, new_runmes_only=False):
    """Setup a folder to do the ASR pipeline analysis for one orthogroup.
        new_runmes_only=True will update the *.config and runme.sh files for this orthogroup,
        but it will not validate the sequences. This is a quick hack method for January 26, 2015."""
        
    cur = con.cursor()
    
    """Get a list of sequence IDs that are in this orthogroup, AND which can be validated across their aa and nt sequences."""
    sql = "select seqid from group_seq where groupid=" + orthogroupid.__str__() + " and seqid in (select seqid from nt_aa_check where checkval>0)"
    cur.execute(sql)
    x = cur.fetchall()
    seqids = []
    for ii in x:
        seqids.append( ii[0] )

    if new_runmes_only==False:
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
    fout.write("VERSION = Feb12.2015\n")
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
    fout.write("MAFFT = /usr/local/bin/mafft\n") # important to use this full path for our cluster config. 'mafft' alone causes problems when launched via mpirun.
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
  
    
def write_asr_scripts(con, skip_existing=True, return_ip=None, return_folder=None):
    """This method writes a Python script, named runme.py, in the data directory
        for each orthogroup. The python script will take care of launching the ASR pipeline
        and returning the data to the master node when it's completed."""
    cur = con.cursor()
    sql = "select groupid from wgd_groups"
    cur.execute(sql)
    x = cur.fetchall()
    commands = []

    for cc in range(0, x.__len__() ):
        ii = x[cc] 
        groupid = ii[0]

        datadir = "data/" + groupid.__str__()
        scriptpath = datadir + "/runme.sh"

        if os.path.exists(datadir):              
            fout = open(scriptpath, "w")
            fout.write("#!/bin/bash\n")
            fout.write("source ~/.bash_profile\n")
            fout.write("cd /tmp/data/" + groupid.__str__() + "\n")
            fout.write("python ~/EclipseWork/asrpipelinev2/runme.py --configpath " + groupid.__str__()  + ".config --skip_zorro\n")
            fout.write("scp -r /tmp/data/" + groupid.__str__() + " " + return_ip.__str__() + ":" + return_folder.__str__() + "/\n")
            fout.close()
            commands.append("source /tmp/data/" + groupid.__str__() + "/runme.sh")
        else:
            write_error(con, "I cannot find the data directory for orthogroup " + groupid.__str__() + " at " + datadir)
            exit()
            
        if False == os.path.exists(scriptpath):
            write_error(con, "I was unable to write the ASR script " + scriptpath)
            exit()
    
    fout = open("asr_commands.sh", "w")
    for c in commands:
        fout.write(c + "\n")
    fout.close()


def distribute_and_launch(con, practice_mode=False, skip_tarsend=False):
	"""1. Distribute the configuration file and runme.sh for every orthogroup to each slave /tmp folder
	   2. Use mpi_dispatch to launch the runme.sh scripts."""
	cur = con.cursor()
	datadir_slave = {}    
	slaves = ["agassiz", "bago", "darwin", "ericsson"]

	"""scp_commands is a list of shell commanes, one per orthogroup, that will launch
		the ancestral analysis for that orthogroup."""
	scp_commands = []
	
	if practice_mode == False:
	    if skip_tarsend == False:
	        """Build and send the tar"""
	        os.system("tar -cvf data.tar data")
	        for slave in slaves:
	            os.system("scp data.tar " + slave + ":/tmp/")
	    
	        """Unpack the tar"""     
	        for slave in slaves:
	            os.system("ssh " + slave + " tar xvf /tmp/data.tar -C /tmp/")
	    
		"""Distribute the data"""
		cur = con.cursor()
		sql = "delete from Settings where keyword='last_distributed'"
		cur.execute(sql)
		con.commit()
		sql = "insert into Settings (keyword,value) values('last_distributed'," + time.time().__str__() + ")"
		cur.execute(sql)
		con.commit()
	    	    
	    for slave in slaves:
	        sql = "select groupid from wgd_groups"
	        cur.execute(sql)
	        x = cur.fetchall()
	        for cc in range(0, x.__len__()):
	            groupid = x[cc][0]
	            scp_commands.append("scp data/" + groupid.__str__() + "/runme.sh " + slave + ":/tmp/data/" + groupid.__str__() + "/" )
	            scp_commands.append("scp data/" + groupid.__str__() + "/" + groupid.__str__() + ".config " + slave + ":/tmp/data/" + groupid.__str__() + "/")
	
	"""Write a machinefile for mpirun"""
	fout = open("hosts.txt", "w")
	fout.write("localhost slots=1\n")
	count_nodes = 0
	for s in slaves:
		count_nodes += 1
		fout.write(s + " slots=2\n")
	fout.close()
	
	fout = open("hosts.local.txt", "w")
	fout.write("localhost slots=5\n")
	fout.close()
	
	"""Write a shell script that can launch all the slaves."""
	fout = open("scp_commands.sh", "w")
	for c in scp_commands:
	    fout.write(c + "\n")
	fout.close()
	if False == practice_mode:
	    os.system("source scp_commands.sh")
	
	"""Finally. . . launch the analysis!"""
	if False == practice_mode:
		cur = con.cursor()
		sql = "delete from Settings where keyword='last_launched'"
		cur.execute(sql)
		con.commit()
		sql = "insert into Settings (keyword,value) values('last_launched'," + time.time().__str__() + ")"
		cur.execute(sql)
		con.commit()
		os.system("mpirun -np " + (1+2*count_nodes).__str__() + " --machinefile hosts.txt /common/bin/mpi_dispatch asr_commands.sh")
 
 
def quickcheck_asr_output(con):
    cur = con.cursor()
    sql = "select groupid from wgd_groups"
    cur.execute(sql)
    x = cur.fetchall()
    
    count_good_groups = 0
    for cc in range(0, x.__len__()):
        ii = x[cc]  
        groupid = ii[0]
        datadir = "data/" + groupid.__str__()
        
        """Look for the dN/dS vs. Df comparison"""
        comppath = datadir + "/compare_dnds_Df.txt"
        dbpath = datadir + "/asr.db"
        
        if False == os.path.exists(comppath):
            write_error(con, "I cannot find the dnds-vs.-df comparisons for orthogroup ID " + groupid.__str__() )
        elif False == os.path.exists(dbpath):
            write_error(con, "I cannot find the SQL database for orthogroup ID " + groupid.__str__() )
        else:
            count_good_groups += 1
            write_log(con, "OK. I found the dnds-vs.-df comparisons for orthogroup ID " + groupid.__str__() )
            #fin = open(comppath, "r")
            #count = 0
            #for l in fin.xreadlines():
            #    if l.__len__() > 1:
            #        count += 1
            #fin.close()
            #print "\n. . . " + count.__str__() + " sites."
    write_log(con, "I found " + count_good_groups.__str__() + " orthogroups with data, and " + (x.__len__()-count_good).__str__() + " without data.")


def read_all_dnds_df_comparisons2(con):
    """Note: this method is very messy at the moment."""
    
    cur = con.cursor()
    sql = "select groupid from wgd_groups"
    cur.execute(sql)
    x = cur.fetchall()
    
    """The following lists are just floating-point values that can be used in the scatterplot.
        These lists include points from multiple orthogroups."""
    cat1points = []
    negcat1points = []
    cat2points = []
    cat3points = []
    cat23points = []
    cat23points_sig = []
    
    dfpoints = []
    dfranks = []
    dfranks_sig = []
    
    dfpercentiles = []
    dfpercentiles_sig = []
    
    fromaa_toaa = {} # key = from AA, value = hash; key = to AA, value = count of this mutation
    
    kpoints = []
    ppoints = []
    
    fout = open("outlier_sites.txt", "w")
    
    """ NOTE: we're restricting the analysis to the first 100 orthogroups only."""
    count_good_groups = 0
    count_good_sites = 0
    
    """Iterate over groups"""
    for ii in x:
        groupid = ii[0]
        try:
            datadir = "data/" + groupid.__str__()
            dbpath = datadir + "/asr.db"
            comppath = datadir + "/compare_dnds_Df.txt"
                    
            if False == os.path.exists(comppath):
                write_log(con, "Group " + groupid.__str__() + " has no dnds/df comparison file, skipping: " + comppath)
                continue
            
            last_modified = os.path.getmtime(comppath)
            now = time.time()
            if now - last_modified > 250000:
                print "\n. Group " + groupid.__str__() + " has a stale timestamp. Skipping."
                # skip this orthogroup
                continue
                    
            insert_orthogrup_action_timestamp(con, groupid, timestamp=last_modified, action="last modified")
            
            print ". Reading ", comppath
            
            """FIRST PASS - retrieve and rank df, k, and p"""
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
                bebancmu    = int(tokens[9])
                
                
                """Skip sites in which there was amino acid mutation on the focus branch."""
                if bebancmu < 1:
                    continue
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
                df_site_rank[ tuple[1] ] = index+1
            for index, tuple in enumerate(k_ranked):
                k_site_rank[ tuple[1] ] = index+1
            for index, tuple in enumerate(p_ranked):
                p_site_rank[ tuple[1] ] = index+1
            total_sites = df_sites.__len__()
            
            """site_dfpercentile: key = site, value = df percentile relative to other sites in this orthogroup."""
            site_dfpercentile = {}
            for index, tuple in enumerate(df_ranked):
                site = tuple[1]
                percentile = 1.0 - float(index) / df_ranked.__len__()
                site_dfpercentile[site] = percentile
                        
            """Quality check:"""
            if df_sites.__len__() == 0 or k_sites.__len__() == 0 or p_sites.__len__() == 0:
                write_error(con, "The comparison file seems to lack any data, skipping: " + comppath)
                continue
            else:
                count_good_groups += 1
                
            """Open the orthogroup's SQL database, and retrieve some basic information."""
            if False == os.path.exists(dbpath) and True == os.path.exists(comppath):
                print "\n. What? ", dbpath, comppath
            if False == os.path.exists(dbpath):
                write_error(con, "I cannot find the SQL database at " + dbpath)
                continue
        
            """Connect to the orthogroup's database:"""
            asrcon = lite.connect(dbpath)
            asrcur = asrcon.cursor()
            print "--> SQL OK."
                            
            sql = "select value from Settings where keyword='seedtaxa'"
            asrcur.execute(sql)
            seedname = asrcur.fetchone()[0]
    
            """Second pass -- deeply parse the results."""
            fin = open( comppath, "r" )
            for l in fin.xreadlines():
                if l.__len__() > 5:
                    tokens = l.split()
                    if tokens.__len__() < 14:
                        continue
                    
                    """Parse the line with data for this site"""
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
                    bebsig      = int(tokens[10]) # there is significant evidence of positive selection
                    df          = abs(  float( tokens[11]) )
                    k           = float( tokens[12] )
                    p           = float( tokens[13] )
                    
                    if tokens.__len__() >= 18:
						anc1state = tokens[14]
						anc1pp = tokens[15]
						anc2state = tokens[16]
						anc2pp = tokens[17]
						
						"""Add to the count for this mutation."""
						if anc1state not in fromaa_toaa:
							fromaa_toaa[anc1state] = {}
						if anc2state not in fromaa_toaa[anc1state]:
							fromaa_toaa[anc1state][anc2state] = 0
						fromaa_toaa[anc1state][anc2state] += 1
                    
                    """Restrict the analysis to sites at which an ancestral mutation occurred."""
                    if bebancmu < 1:
                        continue
                    
                    count_good_sites += 1
                    
                    """Restrict the analysis to orthogroups in which BEB was used.
                        i.e., find at least one BEB probability that is greater than 0.0"""
                    if bebcat1 == 0.0 and bebcat2 == 0.0 and bebcat3 == 0.0:
                        continue
    
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
                    
                    dfpercentile = site_dfpercentile[site]
                    dfpercentiles.append( dfpercentile )
                    
                    if bebsig > 0:
                        """If the site is significant, remember the x,y values for a seperate overlay scatterplot."""
                        cat23points_sig.append( max(bebcat2,bebcat3) )
                        dfranks_sig.append( dfrank )
                        dfpercentiles_sig.append( dfpercentile )
                    
                    is_outlier = 0
                    
                    """Is this site an interesting outlier?"""
                    if dfrank == 1 and bebsig > 0:
                        is_outlier = 1
                    elif dfpercentile > 0.95 and bebsig == 0:
                        is_outlier = 2
                    elif dfpercentile < 0.5 and bebsig > 0:
                        is_outlier = 4
                        
                    if is_outlier == 0:
                        """Then we're done processing this site."""
                        continue
                    
                    """Otherwise, analyze this outlier site.
                        Put useful output in the variable 'printline'.  
                        We'll spew its contents at the end of this loop."""
                    printline = ""
                    printline += " * Special Case Type " + is_outlier.__str__() + " in " +  dbpath + "\n"
                    printline += "\t Seed Taxon: " + seedname + "\n"
                    printline += "\t Site: " + site.__str__() + "\n"
                    printline += "\t Df rank: " + dfrank.__str__() + "\tDf percentile: " + dfpercentile.__str__() + "\tBEB significant: " + bebsig.__str__() + "\n"                           

                    printline += "\t" + anc1state + " (" + anc1pp.__str__() + ")"
                    printline += "\t" + anc2state + "(" + anc2pp.__str__() + ")"
                    printline += "\t bebancmu=" + bebancmu.__str__() + ", nebancmu=" + nebancmu.__str__() + "\n"
                    
                    #if bebsig > 0:
                    sql = "select id from DNDS_Models where name='Nsites_branch'"
                    asrcur.execute(sql)
                    zz = asrcur.fetchall()
                    if zz.__len__() == 0:
                        #write_log(con, "There are no DNDS_Models in the database, so I'm skipping the comparison of DNDS to Df.")
                        continue
                    nsites_id = zz[0][0]
                    
                    """Get some summary of the codon variation at this site."""
                    sql = "select id from DNDS_Tests where almethod=" + almethod.__str__()
                    sql += " and phylomodel=" + phylomodel.__str__() 
                    sql += " and anc1=" + ancid1.__str__() 
                    sql += " and anc2=" + ancid2.__str__()
                    sql += " and dnds_model=" + nsites_id.__str__()
                    asrcur.execute(sql)
                    qq = asrcur.fetchall()
                    if qq.__len__() == 0:
                        print "\n. 1088 Something is wrong:"
                        print sql
                        continue
                    testid = qq[0][0]
                    sql = "select * from DNDS_params where testid=" + testid.__str__()
                    asrcur.execute(sql)
                    qq = asrcur.fetchall()
                    if qq.__len__() == 0:
                        print "\n. 1096 Something is wrong:"
                        print sql
                        continue
                    printline += "\t DNDS: " + qq[0][1:].__str__() + "\n"
                    sql = "select * from NEB_scores where testid=" + testid.__str__() + " and site=" + site.__str__()
                    asrcur.execute(sql)
                    qq = asrcur.fetchall()
                    nebppcat1 = qq[0][2]
                    nebppcat2 = qq[0][3]
                    nebppcat3 = qq[0][4]
                    nebppcat4 = qq[0][5]
                    nebsignificant = qq[0][7]
                    printline += "\t NEB probabilities: " + nebppcat1.__str__() + " " + nebppcat2.__str__() + " " +  nebppcat3.__str__() + " " + nebppcat4.__str__() +  " -- sig: " + nebsignificant.__str__() + "\n"
                    sql = "select * from BEB_scores where testid=" + testid.__str__() + " and site=" + site.__str__()
                    asrcur.execute(sql)
                    qq = asrcur.fetchall()
                    if qq.__len__() == 0:
                        print "\n. 1113 Something is wrong:"
                        print sql
                        continue
                    bebppcat1 = qq[0][2]
                    bebppcat2 = qq[0][3]
                    bebppcat3 = qq[0][4]
                    bebppcat4 = qq[0][5]
                    bebsignificant = qq[0][7]
                    printline += "\t BEB probabilities: " + bebppcat1.__str__() + " " + bebppcat2.__str__() + " " + bebppcat3.__str__() + " " + bebppcat4.__str__() + " -- sig: " + bebsignificant.__str__() + "\n"
            
                    """Get the amino acid and codon compositions at this site."""
                    taxon_aa = {}
                    taxon_codon = {}
                    
                    sql = "select taxonid, alsequence from AlignedSequences where almethod=" + almethod.__str__() + " and datatype=1"
                    asrcur.execute(sql)
                    qq = asrcur.fetchall()
                    for ss in qq:
                        aachar = ss[1][site-1:site]
                        sql = "select shortname from Taxa where id=" + ss[0].__str__()
                        asrcur.execute(sql)
                        ww = asrcur.fetchone()
                        if ww == None:
                            print "\n. 1138: An error occurred:"
                            print "Group ID:", groupid.__str__()
                            print sql
                            continue
                        taxonname = ww[0]
                        taxon_aa[taxonname] = aachar
                    
                    sql = "select taxonid, alsequence from AlignedSequences where almethod=" + almethod.__str__() + " and datatype=0"
                    asrcur.execute(sql)
                    qq = asrcur.fetchall()
                    for ss in qq:                                
                        sql = "select shortname from Taxa where id=" + ss[0].__str__()
                        asrcur.execute(sql)
                        taxonname = asrcur.fetchone()[0]
            
                        from Bio.Seq import Seq
                        from Bio.Alphabet import generic_dna
                        
                        codon = ss[1][ ((site-1)*3):(site*3) ]
                        taxon_codon[taxonname] = codon
                        if codon != "---":
                            codon_seq = Seq(codon, generic_dna)
                            translated_codon = codon_seq.transcribe().translate()
            
                            taxon_codon[taxonname] = codon                       
                    
                    taxon_species = {}
                    for taxonname in taxon_aa:
                        sql = "select name from species where id in (select speciesid from seqnames where name='" + taxonname + "')"
                        cur.execute(sql)
                        taxon_species[taxonname] = cur.fetchone()[0].__str__()
                    
                    """Print a mini PHYLIP alignment of the codon"""
                    printline += "\n"
                    printline += " " + taxon_aa.__len__().__str__() + "\t3\n"
                    for taxonname in taxon_aa:
                        if taxonname not in taxon_aa:
                            "\t\t", taxonname, " has no a.a."
                        elif taxonname not in taxon_codon:
                            "\t\t", taxonname, " has no codon"
                        else:
                            printline += taxonname + "\t" + taxon_codon[ taxonname ] + " " + taxon_codon[ taxonname ] + "\n"
                    printline += "\n"
                    
                    sql = "select newick from UnsupportedMlPhylogenies where almethod=" + almethod.__str__() + " and phylomodelid=" + phylomodel.__str__()
                    asrcur.execute(sql)
                    newick = asrcur.fetchone()[0]
                    
                    """Append species info to the newick"""
                    for taxonname in taxon_aa:
                        try:
                            newick = re.sub(taxonname, "'" + taxon_codon[taxonname] + " [" + taxon_aa[taxonname] + "] " + taxon_species[taxonname] + " [" + taxonname + "]'", newick)
                        except KeyError:
                            print "1173:", taxonname, newick
                        #newtaxonname = re.sub("\.", "", taxonname)
                        #newick = re.sub(taxonname, "'" + newtaxonname + "'", newick)
                    #newick = re.sub(" ", "", newick)
                    #newick = re.sub("\.", "", newick)
                    #print newick
                    
                    """Get a string tree."""
                    from cStringIO import StringIO
                    
                    # This code block is complicated way to get the print_plot to a string.
                    # 1. write the newick tree to a file
                    # 2. read that file
                    # 3. print_plot
                    # 4. capture the output from stdout into a string.
                    # 5. return things to normal
                    backup = sys.stdout
                    sys.stdout = StringIO()     # capture output
                    import dendropy
                    from dendropy import Tree
                    this_tree = Tree()
                    fnewickout = open("/tmp/" + groupid.__str__() + ".tre", "w")
                    fnewickout.write( newick + "\n")
                    fnewickout.close()
                    this_tree.read_from_path("/tmp/" + groupid.__str__() + ".tre", "newick")
                    this_tree.print_plot(plot_metric="length", display_width=100)
                    
                    #this_tr = dendropy.Tree(stream=StringIO(newick), schema="newick")
                    #this_tr.print_plot(plot_metric="length", display_width=50)
                    
                    out = sys.stdout.getvalue() # release output
                    sys.stdout.close()  # close the stream 
                    sys.stdout = backup # restore original stdout
                    
                    #print out.upper()   # post processing
        
                    printline += "Newick: " + out + "\n"
                    
                    print printline 
                    fout.write( printline.__str__() )                             
    
            fin.close()
        except:
            print "\n. An error occurred with group", groupid
    fout.close()
    
    print "\n. I found " + count_good_groups.__str__() + " good orthogroups."
    print "\n. I found " + count_good_sites.__str__() + " good sites."
    
    #scatter1(cat23points, dfpoints, xlab="Probability of Positive Selection", ylab="Df")
    
    if dfpoints.__len__() > 0:
        """PLOT 1: probability of positive selection vs. Df rank."""
        xsets = [cat23points, cat23points_sig]
        ysets = [dfranks, dfranks_sig]
        scatter1("scatter.logprobcat23-logdfrank", xsets, ysets, xlab="log(Probability of Positive Selection)", logx=True, logy=True, ylab="log(dF Rank)", format="jpeg")
        scatter1("scatter.probcat23-dfrank", xsets, ysets, xlab="Probability of Positive Selection", logx=False, logy=False, ymin=0,xmin=0.0, xmax=1.0,ylab="dF Rank", format="jpeg")
    
        """PLOT 2: probability of positive selection vs. Df percentile."""
        xsets =[cat23points, cat23points_sig]
        ysets = [dfpercentiles, dfpercentiles_sig]
        scatter1("scatter.logprobcat23-logdfpercentile", xsets, ysets, xlab="log(Probability of Positive Selection)", logx=True, logy=True, ylab="log(dF Percentile)", format="jpeg")
        scatter1("scatter.logprobcat23-dfpercentile", xsets, ysets, xlab="log(Probability of Positive Selection)", logx=True, logy=False, ymin=0.0, ymax=1.0, ylab="dF Percentile", format="jpeg")
        scatter1("scatter.probcat23-dfpercentile", xsets, ysets, xlab="Probability of Positive Selection", logx=False, logy=False, ymin=0, ymax=1.0,xmin=0.0, xmax=1.0, ylab="dF Percentile", format="jpeg")
        
    else:
        write_error(con, "I didn't find any data; skipping the scatterplot.")
      
      
           
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
            
        

        
        
