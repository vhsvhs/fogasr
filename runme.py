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

write_asr_scripts(con, skip_existing=False, return_ip="10.0.0.100", return_folder="/Volumes/RAID/victor/fogasr/data")
distribute_to_slaves(con, practice_mode=False, skip_tarsend=True)

#validate_asr_output(con)
#read_all_dnds_df_comparisons2(con)