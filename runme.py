from fogasrdb import *
from tools import *

""" MAIN """
con = build_db("fogasr.db")
cur = con.cursor()
import_sequences(con)
verify_sequences(con)
import_orthogroups(con)