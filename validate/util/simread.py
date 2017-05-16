import sys
import subprocess


def simread(refFile, inFastaFile, len_read, cov, out_dir, thread):
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

    ##[alignments with maximal exact matches to sort bam]
    output = inFastaFile.rstrip('.fasta') + ".bam"
    cmdaln = ['bwa mem -t %s %s %s %s | samtools view -Sb - -o - | samtools sort -o %s -'%(thread,refFile,read1,read2,output)]
    subprocess.call(cmdaln, shell=True)

    # ##[sam to bam]
    # bamFile = inFastaFile.split('.fasta')[0]+'_aln.bam'
    # cmdstob = ['samtools view -bS -o %s %s'%(bamFile,alignFile)]
    # subprocess.call(cmdstob, shell=True)    

    # ##[make sort]
    # cmdsort = ['samtools sort -o %s %s'%(inFastaFile.split('.fasta')[0],bamFile)]
    # subprocess.call(cmdsort, shell=True)
