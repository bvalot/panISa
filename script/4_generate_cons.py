import os
from os import system, listdir
from sys import argv
from Bio import	SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import *


#set working directory
os.chdir('..')


def msa_process(in_seq_dir):

	out_msa_dir = 'msa_'+in_seq_dir.split('/')[0].split('seq')[1]
	system('rm -rf %s/*'%out_msa_dir)
	system('mkdir  %s'%out_msa_dir)


	for seq_file in listdir(in_seq_dir):
		system('muscle -in %s%s -out %s/%s.fasta'%(in_seq_dir,seq_file,out_msa_dir,seq_file))

	print('\n---MSA process is completed---\n')


def consensus_porcess(in_msa_dir, out_cons_dir, pc_ident):

	system('rm -rf %s*'%out_cons_dir)
	system('mkdir %s'%out_cons_dir)

	#[call em_cons from emboss]
	# for msa_file in tqdm(listdir(in_msa_dir)):

	# 	name = msa_file.split('.')[0]
	# 	system('em_cons -sequence %s%s -outseq %s%s -identity 80 -name %s -snucleotide 1'%(in_msa_dir,msa_file,out_cons_dir,name,name))

	for msa_file in tqdm(listdir(in_msa_dir)):

		msa_seq_list = []
		len_msa_seq = 0
		consensus_temp = []
		cons_desc = []
		cons_seq = []


		for msa_seq in SeqIO.parse(in_msa_dir+msa_file, "fasta"):
			msa_seq_list.append(msa_seq)


		for pos in range(0,len(msa_seq)):

			allele = set()
			values = []

			for msa_seq in msa_seq_list:
				base = str(msa_seq.seq[pos])
				if base == "-":
					continue

				allele.add(base)
				values.append(base)

			for base in allele:
				if len(values) >= 5 and float(values.count(base))/len(values) >= (float(pc_ident)/100):
					consensus_temp.append(base)

					# if float(values.count(base))/len(values) > 0.9:
					# 	consensus_temp.append(base)
					# else:
					# 	consensus_temp.append(base.lower())


			if len(consensus_temp)-1 != pos:
				consensus_temp.append("N")


		consensus = "".join(consensus_temp)

		#[setup info to write file]
		cons_id = msa_file.split('.')[0]
		cons_desc = 'Start' if 'S' in cons_id.split('_')[-1] else 'End'
		cons_seq = Seq(consensus)

		record_cons_seq = SeqRecord(cons_seq, id = cons_id, description = cons_desc)


		#[write file that not contain only 'N']
		if len(set(cons_seq)) > 1: 
			file_name = '_'.join(msa_file.split('_')[0:-1])

			write_cons_file = open(out_cons_dir+file_name+'|temp','a')
			SeqIO.write(record_cons_seq, write_cons_file, "fasta")
			write_cons_file.close()
		else:
			pass


def concat_cons_file(out_cons_dir):

	for cons_file in listdir(out_cons_dir):

		if 'temp' in cons_file:

			merge_cons_id = ''
			merge_cons_desc = ''
			merge_cons_seq = ''
			file_name = cons_file.split('|temp')[0]


			for cons_data in SeqIO.parse(out_cons_dir+cons_file, 'fasta'):
		
				#[setup info to write file]
				merge_cons_id = file_name
				merge_cons_desc = cons_data.description.split(' ')[1] if len(merge_cons_desc) == 0 else merge_cons_desc+'-'+cons_data.description.split(' ')[1]
				merge_cons_seq = cons_data.seq if len(merge_cons_seq) == 0 else merge_cons_seq+cons_data.seq

				record_seq = SeqRecord(merge_cons_seq, id = merge_cons_id, description = merge_cons_desc)


				write_file = open(out_cons_dir+merge_cons_id,'w')
				SeqIO.write(record_seq, write_file, "fasta")
				write_file.close()

			os.remove(out_cons_dir+cons_file)







def generate_cons(in_seq_dir, out_cons_dir, pc_ident):

	msa_process(in_seq_dir)
	in_msa_dir = 'msa_'+in_seq_dir.split('/')[0].split('seq')[1]+'/'
	consensus_porcess(in_msa_dir, out_cons_dir, pc_ident)
	concat_cons_file(out_cons_dir)



# ===================using functions===================
# generate_cons(argv[1],argv[2],argv[3])
generate_cons('seq15_B_softClip/bwa_8000847_7132693/','cons_15_B_softClip/',80)
# generate_cons('seq15_fullRead/bwa_7139915_7132693/','cons_15_fullRead/')

