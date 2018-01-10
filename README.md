# panISa
panISa is a software identifying insertion sequence (IS) on resequencing data (bam file) in bacterial genomes.

## Idea
The aim of panISa is to search for Insertion Sequence on NGS data without knowing the sequence of the potential IS present in the bacterial strain.
To achieve that, the software identify a signature of insertion in the alignment by counting clip reads on the left and right position of the potential IS. 
These reads overlapped on the direct repeat due to IS insertion.
Finally, using a reconstruction of the beginning of both side of the IS, the software valided the IS by searching for inverted repeat region.

## Requirements and Installation
### Requirements
The program used the python library **pysam** (>0.9)

You need to install the **emboss** package:
http://emboss.sourceforge.net

In debian, type:
<pre>sudo apt-get install python-pysam emboss</pre>

### Installation
Download the current tarball and unzip it.

Verify the installation using the test file
<pre>python panISa.py test/test.bam</pre>

## Command and Options
<pre>python panISa.py [options] bam</pre>

### Options
<pre>-h           show this help message and exit<br />
-o [<i>str</i>]     Return list of IS insertion by alignment [stdout]<br />
-q [<i>int</i>]     Minimum alignment quality value to conserve a clip read [20]<br />
-m [<i>int</i>]     Minimum number of clip reads to look at IS on a position [5]<br />
-s [<i>int</i>]     Maximum size of direct repeat region [15bp]<br />
-p [<i>float</i>]   Minimum percentage of same base to create consensus [0.8]<br />
-v           show program's version number and exit<br /></pre>

## Output
PanISa returns result in tabular format with the following columns:
* Chromosome name: location of the potential IS
* Left position: position of the first base of the direct repeat
* Clip reads left: number of clip reads (left position)
* Direct repeat: nucleotidic sequence of the direct repeat
* Right position: position of the last base of the direct repeat
* Clip reads right: number of clip reads (right position)
* Inverted repeats: nucleotidic sequence of inverted repeats and their position
* IS left sequence: reconstruction of the left boundary of the potential IS
* IS right sequence: reconstruction of the right boundary of the potential IS

## Recommandation
panISa works well with the alignment from **bwa** software.
