import urllib   
from fogasrdb import *
from tools import *

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


""" MAIN """
con = build_db("fogasr.db")
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
        sql = "insert into seqnames (speciesid, name) values(" + speciesid.__str__() + ",'" + t + "')"
        cur.execute(sql)
        con.commit()
        
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
        
