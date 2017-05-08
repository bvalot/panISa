from sys import argv
import subprocess


def simread(refFile, inFastaFile, len_read, cov, out_dir):
    ##[make reads]
    cmdsim = ['dwgsim -1 %s -2 %s -C %s %s %ssim'%(len_read,len_read,cov,inFastaFile,out_dir)]
    subprocess.call(cmdsim, shell=True)

    ##[make index]
    cmdindex = ['bwa index %s'%(refFile)]
    subprocess.call(cmdindex, shell=True)

    ##[make align]
    read1 = out_dir+'sim.bwa.read1.fastq'
    read2 = out_dir+'sim.bwa.read2.fastq'
    alignFile = inFastaFile.split('.fasta')[0]+'_aln.sam'

    ##[alignments with maximal exact matches]
    cmdaln = ['bwa mem %s %s %s > %s'%(refFile,read1,read2,alignFile)]
    subprocess.call(cmdaln, shell=True)

    ##[sam to bam]
    bamFile = inFastaFile.split('.fasta')[0]+'_aln.bam'
    cmdstob = ['samtools view -bS -o %s %s'%(bamFile,alignFile)]
    subprocess.call(cmdstob, shell=True)    

    ##[make sort]
    cmdsort = ['samtools sort %s %s'%(bamFile, inFastaFile.split('.fasta')[0])]
    subprocess.call(cmdsort, shell=True)


# ===================using functions=====================
# simread(argv[1],argv[2],argv[3])
# simread('PA_PAO1.fasta','NC_002516_simIS.fasta', 150, 50) 
