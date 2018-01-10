panISa
======

panISa is a software identifying insertion sequence (IS) on resequencing
data (bam file) in bacterial genomes.

Idea
----

The aim of panISa is to search for Insertion Sequence on NGS data
without knowing the sequence of the potential IS present in the
bacterial strain. To achieve that, the software identifies a signature of
insertion in the alignment by counting clipped reads on the left and right
position of the potential IS. These reads overlapped on the direct
repeat due to IS insertion. Finally, using a reconstruction of the
beginning of both side of the IS, the software valided the IS by
searching for inverted repeat region.

Requirements and Installation
-----------------------------

Requirements
~~~~~~~~~~~~

The program used the python library **pysam** (>0.9)

You need to install `the emboss package <http://emboss.sourceforge.net>`_

In debian, type:

.. raw:: html

   <pre>sudo apt-get install python-pysam emboss</pre>

Installation
~~~~~~~~~~~~

Download the current tarball and unzip it.

Verify the installation using the test file

.. raw:: html

   <pre>python panISa.py test/test.bam</pre>

Alternativly, you can install from `PyPI repository <https://pypi.python.org/pypi>`_

.. raw:: html

   <pre>pip install panisa</pre>

   
Command and Options
-------------------

.. raw:: html

   <pre>python panISa.py [options] bam</pre>

Options
~~~~~~~

-h     show this help message and exit
-o     Return list of IS insertion by alignment [stdout]
-q     Minimum alignment quality value to conserve a clipped read [20]
-m     Minimum number of clipped reads to look at IS on a position [5]
-s     Maximum size of direct repeat region [15bp]
-p     Minimum percentage of same base to create consensus [0.8]
-v     show program's version number and exit

Output
------

PanISa returns result in tabular format with the following columns: 

Chromosome:
  chromosome id 
Left position:
  position of the last base of the direct repeat and the left bondary of
  the potential IS
Clipped reads left:
  number of clipped reads (left position)
Direct repeat:
  nucleotidic sequence of the direct repeat
Right position:
  position of the first base of the direct repeat and the right
  bondary of the potential IS
Clipped reads right:
  number of clipped reads (right position)
Inverted repeats:
  nucleotidic sequence of inverted repeats and their position
IS left sequence:
  reconstruction of the left boundary of the potential IS
IS right sequence:
  reconstruction of the right boundary of the potential IS

Recommandation
--------------

panISa works well with the alignment from **bwa** software.
