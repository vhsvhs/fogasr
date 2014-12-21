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

#import os
#os.system("mpirun -np 3 /common/bin/mpi_dispatch asr_commands.sh")

validate_asr_output(con)