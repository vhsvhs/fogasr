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
    print htmlSource
    

con = build_db( dbpath = ap.getOptionalArg("--dbpath") )
(taxon_aaurl, taxon_nturl) = get_taxon_urls()
for taxon in taxon_aaurl:
    sql = "select count(*) from species where name='" + taxon + "'"
    cur.execute(sql)
    count = cur.fetchone()[0]
    if count == 0:
        sql = "insert into species (name) values('" + taxon + "')"
        cur.execute(sql)
        con.commit()    
    print taxon, taxon_aaurl[taxon]
    get_fasta_source( taxon_aaurl[taxon], con)
