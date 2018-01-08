import os
from os import listdir, system
from sys import argv
import pysam
from tqdm import *


# [set working directory]
os.chdir('..')
curDir = os.getcwd()


def read_bam_file(in_bam_dir,out_dir):

	system('mkdir %s'%out_dir)
	system('rm %s*'%out_dir)


	for i_read_file in listdir(in_bam_dir):


		# [define variables]
		info_dict={'ref_id':'','read_name':'','read_start':'','read_end':'', 'read_len':'','ref_start':'','ref_end':'','ref_len':'',\
		'cigar':'','cigar_str':'','orient':'','unmapped':'','seq':''}
		info_list=[]
		write_header=True
		clip = set([4,5])
		count_write_line = 0



		# [open bam file to read]
		read_bam = pysam.AlignmentFile(in_bam_dir+i_read_file, 'rb')

		# [count number of reads]
		count_total_read = pysam.AlignmentFile(in_bam_dir+i_read_file, 'rb').count(until_eof=True)

		# [open file to write]
		write_error_file = open(out_dir+'ref_error_'+i_read_file.split('.')[0],'a')
		write_info_file = open(out_dir+'info_'+i_read_file.split('.')[0],'a')

		if write_header == True:
			write_info_file.writelines('\t'.join(['ref_id','read_name','read_start','read_end','read_length','ref_start','ref_end','ref_length',\
				'cigar','clip_start','clip_end','strand','unmapped','seq'])+'\n')
			write_header = False




		for num_read, each_read in enumerate(tqdm(list(read_bam.fetch(until_eof=True)))):


			# [define value to orientation strand]
			if each_read.is_reverse == True:
				orient_read = "Reverse"
			else:
				orient_read = "Forward"

			# [get reference name from refernece id (tid)]
			if each_read.tid >= 0 :
				ref_name = read_bam.getrname(each_read.tid)
			else:
				ref_name = each_read.tid

				# [write error values]
				write_error_file = open(out_dir+'ref_error_'+i_read_file.split('.')[0],'a')
				write_error_file.writelines('\t'.join([str(num_read), str(each_read.tid), str(ref_name), str(each_read.qname), str(each_read.cigarstring)])+'\n')


			# [save each value into dict]
			info_dict = {'ref_id':ref_name, 'read_name':each_read.qname, 'read_start':each_read.qstart, 'read_end':each_read.qend, 'read_len':each_read.qlen,\
			'ref_start':each_read.reference_start, 'ref_end':each_read.aend, 'ref_len':each_read.alen, \
			'cigar':each_read.cigar, 'cigar_str':each_read.cigarstring, 'orient':orient_read, 'unmapped':each_read.is_unmapped, 'seq':each_read.seq}

			info_list.append(info_dict)



			# [skip for full map and 'None']
			if len(info_dict['cigar']) <= 1:
				continue

			# [looking for seq after end_position]
			value_end = info_dict['cigar'][-1]

			
			# [looking for seq before start_position]
			value_start = info_dict['cigar'][0]


			# [define clip status]
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



			# [filter for over both start and end postions]
			if (value_start[0] in clip) or (value_end[0] in clip):
				write_info_file.writelines('\t'.join([str(info_dict['ref_id']), str(info_dict['read_name']), \
					str(info_dict['read_start']), str(info_dict['read_end']), str(info_dict['read_len']), \
					str(info_dict['ref_start']), str(info_dict['ref_end']), str(info_dict['ref_len']), \
					str(info_dict['cigar_str']), clip_start, clip_end, str(info_dict['orient']), \
					str(info_dict['unmapped']), str(info_dict['seq'])])+'\n')
				
				count_write_line += 1


		# [print for checking]
		print('======')
		print("Number of reads: %d"%num_read)
		print("Number or reads from count method:%d"%count_total_read)
		print('======\n')
		print("Number of written items: %d"%count_write_line)


		write_error_file.close()
		write_info_file.close()


		read_bam.close()



# ===================using functions===================
read_bam_file('data/','output/')


