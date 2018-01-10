#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Panisa Treepong
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

def writetabular(output,couples):
    output.write("\t".join(["Chromosome","Left position", "Clipped read left", "Direct repeats","Right position", "Clipped read rigth", \
        "Inverted repeats","Left sequences","Right sequences"]) + "\n")
    for c in couples:
        towrite = []
        towrite.extend([c.chrom, str(c.posend.pos), str(len(c.posend.clipend)), \
                        c.dr, str(c.posstart.pos + 1), str(len(c.posstart.clipstart))])
        if c.ir is None:
            towrite.append("No IR")
        else:
            towrite.append(c.ir)
        towrite.extend([c.cons5prime, c.cons3prime])

        output.write("\t".join(towrite) + "\n")

if __name__=='__main__':
    import doctest
    doctest.testmod()
