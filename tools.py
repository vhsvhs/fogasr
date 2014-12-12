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
            #print taxon, aaurl, nturl
            if taxon.__len__() > 1 and aaurl.__len__() > 1 and nturl.__len__() > 1:
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
    
    (taxon_aaurl, taxon_nturl) = get_taxon_urls()
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
        taxa_seq = get_fasta_source( taxon_aaurl[taxon], con)
        for t in taxa_seq:
            sql = "select count(*) from seqnames where speciesid=" + speciesid.__str__() + " and name='" + t + "'"
            cur.execute(sql)
            count = cur.fetchone()[0]
            if count == 0:
                sql = "insert into seqnames (speciesid, name) values(" + speciesid.__str__() + ",'" + t + "')"
                cur.execute(sql)
                con.commit()
            else:
                write_log("Warning: I found a previous copy an AA sequence named " + t)
            
            sql = "select id from seqnames where speciesid=" + speciesid.__str__() + " and name='" + t + "'"
            cur.execute(sql)
            seqid = cur.fetchone()[0]
            
            sql = "insert into aaseqs (speciesid, sequence, seqid) values("
            sql += speciesid.__str__() + ",'" + taxa_seq[t] + "'," + seqid.__str__() + ")"
            cur.execute(sql)
        con.commit()
        
        """nt"""
        taxa_seq = get_fasta_source( taxon_nturl[taxon], con)
        for t in taxa_seq:
            sql = "select count(*) from seqnames where name='" + t + "'"
            cur.execute(sql)
            count = cur.fetchone()[0]
            if count == 0:
                """We're expecting to already have seen this sequence name in the AA data"""
                write_log(con, "I can't find the sequence " + t + ".")
                sql = "delete from seqnames where name='" + t + "'"
                cur.execute(sql)
                con.commit()
                continue
            
            sql = "select id from seqnames where speciesid=" + speciesid.__str__() + " and name='" + t + "'"
            cur.execute(sql)
            seqid = cur.fetchone()[0]
            
            sql = "insert into ntseqs (speciesid, sequence, seqid) values("
            sql += speciesid.__str__() + ",'" + taxa_seq[t] + "'," + seqid.__str__() + ")"
            cur.execute(sql)
        con.commit()
        
        """Ensure all the aa sequences have a companion nt sequence."""
        sql = "delete from aaseqs where seqid not in (select id from ntseqs where speciesid=" + speciesid.__str__() + ")"
        cur.execute(sql)
        con.commit()
        
        sql = "delete from seqnames where id not in (select seqid from ntseqs) or id not in (select seqid from aaseqs)"
        cur.execute(sql)
        con.commit()
        
        sql = "select count(*) from seqnames where speciesid=" + speciesid.__str__()
        cur.execute(sql)
        count = cur.fetchone()[0]
        write_log(con, "Species " + speciesname + " has " + count.__str__() + " verified gene sequences.")

def verify_sequences(con):
    cur = con.cursor()
    sql = "select id from species"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        speciesid = ii[0]
        speciesname = get_species_name(con, speciesid)
        sql = "select count(*) from seqnames where speciesid=" + speciesid.__str__()
        cur.execute(sql)
        count = cur.fetchone()[0]
        write_log(con, "Species " + speciesname + " has " + count.__str__() + " verified gene sequences.")
        

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
    
    lcounter = 0
    for l in lines:
        if l.__contains__(">"):
            lcounter = 0
        else:
            lcounter += 1
        
        
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
                species = item.split("|")[0]
                species_frag = species[0].upper() + ". " + species[1:]
                speciesid = get_species_id(con, species_frag)
                if speciesid == None:
                    write_error("I cannot find the species ID for " + species_frag)
                    exit()
                
                seqname = item.split("|")[1]
                seqid = get_seqid(con, seqname, speciesid)
                
                sql = "insert into group_seq (groupid, seqid) values(" + groupid.__str__() + "," + seqid.__str__() + ")"
                print sql
                
                cur.execute(sql)
            con.commit()
    
    """How many groups did we find?"""
    sql = "select count(*) from orthogroups"
    cur.execute(sql)
    count = cur.fetchone()[0]
    write_log(con, "I found " + count.__str__() + " non-singleton orthogroups.")

    
                
                
                
    
    
    
        