from fogasrdb import *
from tools import *

""" MAIN """
con = build_db("fogasr.db")
cur = con.cursor()
#import_sequences(con)
#verify_sequences(con)
#import_orthogroups(con)
#filter_orthogroups(con)
#check_again_wgdgroups(con)
#setup_all_asr(con)
#write_asr_scripts(con)
#distribute_to_slaves(con)

#launch_remote_slaves(con)
#validate_asr_output(con)
read_all_dnds_df_comparisons(con)