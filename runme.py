from fogasrdb import *
from tools import *

""" MAIN """
con = build_db("fogasr.db")
cur = con.cursor()
#import_sequences(con)
#verify_sequences(con)
#import_orthogroups(con)
#filter_orthogroups(con)
build_asr_directories(con)
write_asr_scripts(con)
launch(con, practice_mode=False)

#quickcheck_asr_output(con)
#read_all_dnds_df_comparisons2(con)





