# panISa
panISa is a software to search insertion sequence (IS) on resequencing data (bam file) in bacteria

## Idea
The principle of panISa is to search Insertion Sequence on NGS data without knowledge of potential IS present in the bacterial strain.
To achieve that, the software identified a signature of insertion in the alignment by counting clip read at start and end position. 
These reads overlapped on the direct repeat due to IS insertion.
Finally, using a reconstruction of the beginning of both side of the IS, the software valided the IS by searching inverted repeat region

## Requirements
The program used the python library **pysam** (>0.9)

You need to install the **emboss** package:
http://emboss.sourceforge.net

In debian, type:
> sudo apt-get install python-pysam emboss<br />
> sudo apt-get install emboss 

## Command and Options
<pre><code>python panISa.py [options] bam</code></pre>

### Options
<pre>-h           show this help message and exit<br />
-o [<i>str</i>]     Return list of IS insertion by alignment [stdout]<br />
-q [<i>int</i>]     Minimun alignment quality value to conserved a clip read [20]<br />
-m [<i>int</i>]     Minimun number of clip read to look at IS on a position [5]<br />
-s [<i>int</i>]     Maximun size of direct repeat region [15]<br />
-p [<i>float</i>]   Minimun percentage of same base for create consensus [0.8]<br />
-v           show program's version number and exit<br /></pre>

## Recommandation
panISa work well with the alignment from **bwa** software.


## Development
The program is currently in developpemnt, but return some result.
