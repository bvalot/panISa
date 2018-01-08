from os import listdir, system
import shutil, os
from sys import argv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam
from tqdm import *
from collections import Counter



#[set working directory]
os.chdir('..')


#[read position files]
def read_position(position_directory):


	#[define variables]
	position_dict={'contig':'','seq_start':'','seq_end':'','strand':'', 'genome':''}
	position_list=[]

	#[get all positions in directory]
	for posit_file in listdir(position_directory):
		print 'position file: '+posit_file+'\n'

		open_posit_file = open(position_directory+posit_file, 'r')
		position_readFile = open_posit_file.read().split('\n')
		open_posit_file.close()
		

		for a_line_position in position_readFile:

			if (a_line_position != '') and ('ref' not in a_line_position):

				position_line = a_line_position.split(',')

				position_dict={\
				'contig':position_line[1],\
				'seq_start':int(position_line[2]),\
				'seq_end':int(position_line[3]),\
				'strand':position_line[4],\
				'genome':position_line[5]}

				position_list.append(position_dict)


	return position_list


#[extract sequences by position]
def extract_seq_bam(position_dir, in_bam_dir, out_seq_dir):

	system('mkdir %s'%out_seq_dir)
	system('rm -rf %s*'%out_seq_dir)
	position_info=[]



	for i_read_file in listdir(in_bam_dir):

		genome_dir = out_seq_dir+i_read_file.split('.')[0]+'/'
		system('mkdir %s'%genome_dir)


		#[define variables]
		info_dict={'ref_id':'','read_name':'','read_start':'','read_end':'', 'read_len':'','ref_start':'','ref_end':'','ref_len':'',\
		'cigar':'','cigar_str':'','orient':'','unmapped':'','seq':''}
		info_list=[]
		write_header=True
		clip = set([4])
		count_write_line = 0
		


		#[open bam file to read]
		read_bam = pysam.AlignmentFile(in_bam_dir+i_read_file, 'rb')

		#[count number of reads]
		count_total_read = pysam.AlignmentFile(in_bam_dir+i_read_file, 'rb').count(until_eof=True)

		#[get position info of all genomes from position directory]
		position_info = read_position(position_dir)



		for num_read, each_read in enumerate(tqdm(list(read_bam.fetch(until_eof=True)))):


			#[define value to orientation strand]
			if each_read.is_reverse == True:
				orient_read = "Reverse"
			else:
				orient_read = "Forward"

			#[get reference name from refernece id (tid)]
			if each_read.tid >= 0 :
				ref_name = read_bam.getrname(each_read.tid)
			else:
				ref_name = each_read.tid



			#[save each value into dict]
			info_dict = {'ref_id':ref_name, 'read_name':each_read.qname, 'read_start':each_read.qstart, 'read_end':each_read.qend, 'read_len':each_read.qlen,\
			'ref_start':each_read.reference_start, 'ref_end':each_read.aend, 'ref_len':each_read.alen, \
			'cigar':each_read.cigar, 'cigar_str':each_read.cigarstring, 'orient':orient_read, 'unmapped':each_read.is_unmapped, 'seq':each_read.seq}

			info_list.append(info_dict)



			#[skip for full map and 'None']
			if len(info_dict['cigar']) <= 1:
				continue

			#[looking for seq after end_position]
			value_end = info_dict['cigar'][-1]

			
			#[looking for seq before start_position]
			value_start = info_dict['cigar'][0]


			#[define clip status]
			clip_end = value_end[0]
			if clip_end in clip:
				clip_end = '1'
			else:
				clip_end = '0'

			clip_start = value_start[0]
			if clip_start in clip:
				clip_start = '1'
			else:
				clip_start = '0'



			#[filter for over both start and end postions]
			if (value_start[0] in clip) or (value_end[0] in clip):
				for each_pos in position_info:
					if (info_dict['ref_id'] == each_pos['contig']) and ((info_dict['ref_start'] == each_pos['seq_start']) or (info_dict['ref_end'] == each_pos['seq_end'])):


						#[prepare sequence record to write]
						prep_id = str(each_pos['genome'])+'|'+str(info_dict['ref_id'])
						prep_description = info_dict['cigar_str']+'_pos:'+str(info_dict['ref_start'])+'-'+str(info_dict['ref_end'])+'|'+info_dict['orient']#info_dict['read_name']
						# prep_seq = Seq(info_dict['seq'])


						#[extract only soft clip seq. and write separately]
						if value_start[0] in clip and value_end[0] in clip:
							continue
							# write_info_file = open(genome_dir+each_pos['contig']+'_'+str(each_pos['seq_start'])+'_S','a')
							# prep_seq = Seq(info_dict['seq'][0:int(value_start[1])])
							# write_info_file = open(genome_dir+each_pos['contig']+'_'+str(each_pos['seq_start'])+'_E','a')
							# prep_seq = Seq(info_dict['seq'][150-int(value_end[1]):150])

						elif value_start[0] in clip:
							write_info_file = open(genome_dir+each_pos['contig']+'_'+str(each_pos['seq_start'])+'_S','a')
							prep_seq = Seq(info_dict['seq'][0:int(value_start[1])])

						elif value_end[0] in clip:
							write_info_file = open(genome_dir+each_pos['contig']+'_'+str(each_pos['seq_start'])+'_E','a')
							prep_seq = Seq(info_dict['seq'][150-int(value_end[1]):150])



						prep_record_seq = SeqRecord(prep_seq, id = prep_id, description = prep_description)


						SeqIO.write(prep_record_seq, write_info_file, "fasta")

						count_write_line += 1

						write_info_file.close()


			

		#[print for checking]
		print('======')
		print("Number of reads: %d"%num_read)
		print("Number or reads from count method:%d"%count_total_read)
		print('======\n')
		print("Number of written items: %d"%count_write_line)


		read_bam.close()


#[delete single file]
def check_single_file(target_dir):

	file_name_list = []
	del_list = []

	#[get all files]
	for path, subdirs, files in os.walk(target_dir):

		for each_file in files:
			file_name = '_'.join(each_file.split('_')[0:-1])
			file_name_list.append(file_name)


	count_name_list = Counter(file_name_list)

	
	#[get single file]
	for key, val in count_name_list.iteritems():

		if val < 2:
			S_key = key+'_S'
			del_list.append(S_key)
			E_key = key+'_E'
			del_list.append(E_key)


	#[remove single file]
	for i_file in files:
		if i_file in del_list:
			print(i_file)
			os.remove(os.path.join(path,i_file))
		else:
			pass




# ===================using functions===================
extract_seq_bam('output_B_15/position/','data/B/','seq15_B_softClip/')
# check_single_file('seq15_softClip_sep')


