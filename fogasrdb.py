from fogasrdb_api import *
from version import *
  
def build_db(dbpath = None):
    """Initializes all the tables. Returns the DB connection object.
    If tables already exist, they will NOT be overwritten."""
        
    if dbpath == None or dbpath == False:
        dbpath = "asr.db"
        print "\n. Creating a new database at", dbpath
    else:
        print "\n. Restoring the existing database at", dbpath

    con = lite.connect(dbpath)

    cur = con.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS About(version TEXT)")
    cur.execute("create table if not exists Settings(keyword TEXT, value TEXT)")
    cur.execute("create table if not exists species(id INTEGER primary key autoincrement, name TEXT unique)")
    cur.execute("create table if not exists Log(id INTEGER primary key, time DATETIME DEFAULT CURRENT_TIMESTAMP,  message TEXT, code INT)")
    cur.execute("create table if not exists ErrorLog(id INTEGER primary key, time DATETIME DEFAULT CURRENT_TIMESTAMP,  message TEXT, code INT)")
    cur.execute("create table if not exists seqnames(id integer primary key autoincrement, speciesid INTEGER, name TEXT)")
    cur.execute("create table if not exists ntseqs(id INTEGER primary key autoincrement, speciesid INT, sequence TEXT, seqid INT)")
    cur.execute("create table if not exists aaseqs(id INTEGER primary key autoincrement, speciesid INT, sequence TEXT, seqid INT)")
    cur.execute("create table if not exists nt_aa_check(seqid INTEGER, checkval INT)") # do the nt and aa sequences correspond? If so, then check == 1
    cur.execute("create table if not exists orthogroups(id INTEGER primary key autoincrement, name TEXT unique)")
    cur.execute("create table if not exists orthogroup_action_timestamp(groupid INTEGER, time FLOAT, action TEXT)") # time should be made from Python's time.time() function.
    cur.execute("create table if not exists group_seq(groupid INTEGER, seqid INTEGER)")
    cur.execute("create table if not exists wgd_groups(groupid INTEGER, depth INTEGER)") # which orthogroups are useful for stuyding the whole genome duplication?
    con.commit()
    
    cur.execute("delete from About")
    con.commit()
    cur.execute("insert into About (version) values(" + FOGASR_VERSION.__str__() + ")")
    con.commit()
    cur.execute("delete from Log")
    con.commit()
    cur.execute("delete from ErrorLog")
    con.commit()
    
    return con
    

    