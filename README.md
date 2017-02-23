# panISa
panISa is a software to search insertion sequence (IS) on resequencing data (bam file) in bacteria

#Idea
The principle of panISa is to search Insertion Sequence on NGS data without knowledge of potential IS present in the bacterial strain.
To achieve that, the software identify a signature of insertion in the alignment with clip read at start and end. 
These reads overlapped on the direct repeat due to IS insertion.
Finally, using a reconstruction of the beginning of both side of the IS, the software valided the IS by searching inverted repeat region

#Development
The program is currently in progress and doesn't work at all.
