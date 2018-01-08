import os
from os import system, listdir
from sys import argv
from tqdm import *


#set working directory
os.chdir('..')


def find_inverted_repeat(in_cons_dir, out_ir_dir, cutoff):

	system('rm -rf %s*'%out_ir_dir)
	system('mkdir %s'%out_ir_dir)
	system('mkdir %sseq'%out_ir_dir)


	for cons_file in tqdm(listdir(in_cons_dir)):

		ir_file = '%s%s' %(out_ir_dir, cons_file)
		seq_ir_file = '%sseq/%s' %(out_ir_dir, cons_file)

		system('einverted -sequence %s%s -gap 12 -threshold %d -match 3 -mismatch -4 -outfile %s -outseq %s -snucleotide 1'%(in_cons_dir,cons_file,cutoff,ir_file,seq_ir_file))
	
		#[delete empty file]
		if os.stat(ir_file).st_size == 0:
			os.remove(ir_file)
			os.remove(seq_ir_file)
		else:
			pass





# ===================using functions===================
# find_inverted_repeat(argv[1],argv[2])
find_inverted_repeat('cons_15_B_softClip/','ir_15_B_softClip/',15)
# find_inverted_repeat('cons_15_fullRead/','ir_15_fullRead/')

