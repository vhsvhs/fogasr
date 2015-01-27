import sqlite3 as lite
import urllib
import os, sys, time

def write_log(con, message, code=None):
    """
    Writes to the log file
    """
    cur = con.cursor()
    sql = "insert into Log (message"
    if code != None:
        sql += ",code"
    sql += ") values(\"" + message
    if code != None:
        sql += "," + code.__str__()
    sql += "\")"
    cur.execute(sql)
    con.commit()
    
    print "\n. " + message
    

def write_error(con, message, code=None):
    cur = con.cursor()
    sql = "insert into ErrorLog (message"
    if code != None:
        sql += ",code"
    sql += ") values(\"" + message
    if code != None:
        sql += "," + code.__str__()
    sql += "\")"
    cur.execute(sql)
    con.commit()
    print "\n. ERROR: " + message
    
def insert_orthogrup_action_timestamp(con, groupid, timestamp=None, action=""):
    cur = con.cursor()
    sql = "insert into orthogroup_action_timestamp (groupid,time,action) "
    if timestamp == None:
        timestamp = time.time()
    sql += "values(" + groupid.__str__() + "," + timestamp.__str__() + ",'" + action + "')"
    cur.execute(sql)
    con.commit()
    
def get_species_name(con, speciesid):
    cur = con.cursor()
    sql = "select name from species where id=" + speciesid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    if x.__len__() == 0:
        return None
    else:
        return x[0][0]

def get_species_id(con, name_contains):
    cur = con.cursor()
    sql = "select id from species where name like '%" + name_contains.__str__() + "%'"
    cur.execute(sql)
    x = cur.fetchall()
    if x.__len__() == 0:
        return None
    else:
        return x[0][0]

def get_seqid(con, seqname, speciesid):
    cur = con.cursor()
    sql = "select id from seqnames where name='" + seqname + "' and speciesid=" + speciesid.__str__()
    #print sql
    cur.execute(sql)
    x = cur.fetchall()
    if x.__len__() == 0:
        return None
    return x[0][0]

def get_seqname(con, seqid):
    cur = con.cursor()
    sql = "select name from seqnames where id=" + seqid.__str__()
    #print sql
    cur.execute(sql)
    x = cur.fetchall()
    if x.__len__() == 0:
        return None
    return x[0][0]