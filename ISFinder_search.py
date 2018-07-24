#!/usr/bin/python
# -*- coding: utf-8 -*-

##Copyright (c) 2017 Benoit Valot and Charlotte Couchoud
##benoit.valot@univ-fcomte.fr
##UMR 6249 Chrono-Environnement, Besançon, France
##Licence GPL

"""automate search in ISFinder from panISa output"""

import argparse
import sys
import requests
try:
    from HTMLParser import HTMLParser
except ImportError:
    from html.parser import HTMLParser

desc = "automate search IS homology in ISFinder from panISa output"
command = argparse.ArgumentParser(prog='ISFinder_search.py', \
    description=desc, usage='%(prog)s [options] file')
command.add_argument('-o', '--output', nargs="?", \
    type=argparse.FileType("wt"), default=sys.stdout, \
    help='Return potential IS found in ISFinder (default:stdout)')
command.add_argument('-r', '--remove', nargs="?", \
    type=str, default=".txt", \
    help='Remove suffix at the end of the filename (default:.txt)')
command.add_argument('-l', '--length', nargs="?", \
    type=int, default=30, \
    help='Length of the IRR-IRL search (default:30)')
command.add_argument('-i', '--identity', nargs="?", \
    type=int, default=90, choices=range(0,101), metavar="0-100", \
    help='Percentage of expected identity (default:90)')
command.add_argument('-e', '--evalue', nargs="?", \
    type=float, default=0.001, \
    help='Expected max evalue (default:1E-3)')
command.add_argument('-a', '--alignment', nargs="?", \
    type=float, default=80, choices=range(0,101), metavar="0-100", \
    help='Percentage of expected alignment (default:80)')
command.add_argument('file', nargs="+", type=argparse.FileType("r"), \
    help='PanISa result files to merge')
command.add_argument('-v', '--version', action='version', \
    version='%(prog)s 0.1.4')


class URLParser(HTMLParser):
    urlres = ''
    def handle_starttag(self, tag, attrs):
        if tag == 'meta':
            for i, attr in enumerate(attrs):
                if i == 1:
                    self.urlres=str(attr[1]).lstrip("2;URL=")

class BlastParser(HTMLParser):
    ispre= False
    result = ""
    def handle_starttag(self, tag, attrs):
        if tag == "pre":
            self.ispre=True
    def handle_endtag(self, tag):
        if tag == "pre":
            self.ispre=False
    def handle_data(self, data):
        if self.ispre:
            self.result += data + "\n"

                    
def read_panisa(files, header, remove):
    mergedPanISa = []
    for fi in files: 
        first=fi.readline()
        for line in fi:
            tab = line.rstrip("\n").split("\t") 
            if len(tab) != 9:
                raise Exception("Not a panISa output or incorrect\n"+line)
            result = {k:v for k,v in zip(header,tab)}
            result["sample"] = fi.name.split("/")[-1].rstrip(remove)
            result["chrompos"] = result.get("chrom")+"_"+result.get("lpos")
            mergedPanISa.append(result)
    return mergedPanISa

def get_irl_irr(mergedPanISa, length):
    seqs = {}
    for v in mergedPanISa:
        name = v.get("sample")    
        irl = v.get("irl")[:length]
        irr = v.get("irr")[-length:]
        while len(irr) < length:
            irr = irr+"N"
        while len(irl) < length:
            irl = "N"+irl
        chroPos = v.get("chrompos")
        currentSeq = irl+irr
        if chroPos not in seqs:
            seqs[chroPos] = currentSeq
        elif currentSeq.count("N") < seqs.get(chroPos).count("N"):
            seqs[chroPos] = currentSeq
    return seqs


def search_ISfinder(potential_sequence):
    ##post data
    r = requests.post('https://www-is.biotoul.fr/blast/ncbiIS.php', \
                      data = {'title' : 'test', 'seq' : potential_sequence, 'seqfile' \
                              : '', 'database' : 'ISfindernt', 'prog' : 'blastn', 'blast' \
                              : 'ok', 'alignment' : '7', 'wordsize' : '11','expect' : '10.0' \
                              , 'gapcost' : '5 2'})
    ##get url result
    parser = URLParser()
    parser.feed(r.text)

    ##downoad blast result
    o = requests.get('https://www-is.biotoul.fr/blast/' + parser.urlres)    
    parser = BlastParser()
    parser.feed(o.text)
    return parser.result

def filter_blast(blast_result, identity, evalue, length, alignment):
    blasts = []
    column_name = ["Query", "Subject_ID", "%identity", "alignment_length", "mismatches", \
                   "gap_opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit_score"]
    for line in blast_result.rstrip("\n").split("\n"):
        if line.startswith('#'):
            continue
        blast = line.split("\t")
        if len(blast) != 12:
            sys.stderr.write(line+"\n")
        else:
            resultBlast = {k:v for k,v in zip(column_name, blast)}
            res_evalue = float(resultBlast.get("evalue"))
            res_identity = float(resultBlast.get("%identity"))
            res_alignment = int(resultBlast.get("q.end")) - int(resultBlast.get("q.start"))
            if res_evalue <= evalue and res_identity >= identity and \
               res_alignment >= (alignment*length)/100:
                blasts.append(line)
    return blasts

def write_result(output, mergedPanISa, length, resultISFinder):
    header = 'Sample\tChromosome\tStart_Position\tStop_Position\tPotential_sequence\tPotential_IS\tAlignment\n'
    output.write(header)
    for v in mergedPanISa:
        towrite = [v.get("sample")]
        towrite.append(v.get("chrom"))
        towrite.append(v.get("rpos"))
        towrite.append(v.get("lpos"))
        sequence = seqs.get(v.get("chrompos"))[:length] \
            + "-" + seqs.get(v.get("chrompos"))[-length:]
        towrite.append(sequence)
        isfinder = resultISFinder.get(v.get("chrompos"))
        if len(isfinder) == 0:
            towrite.append("No identity")
            towrite.append("")
        elif len(isfinder) == 1:
            IS = isfinder[0].split("\t")[1]
            towrite.append(IS)
            towrite.append("one side")
        else:
            IS = isfinder[0].split("\t")[1]
            if IS == isfinder[1].split("\t")[1] :
                towrite.append(IS)
                towrite.append("two side")
            else:
                towrite.append(IS)
                towrite.append("one side")
        output.write("\t".join(towrite) + "\n")
    

    
if __name__=='__main__':
    """Performed job on execution script"""
    args = command.parse_args()

    ##read panISa result and associate file name to sample
    header=["chrom", "lpos", "lclipped", 'dr', "rpos", "rclipped", "ir", "irl", "irr"]    
    mergedPanISa = read_panisa(args.file, header, args.remove)

    ##get irr-irl for each chromosome/position independly of the sample
    seqs = get_irl_irr(mergedPanISa, args.length)        

    ##search ISfinder and filter result
    resultISFinder = {}    
    for chrompos, potential_sequence in seqs.items():
        ##search homology in ISFinder
        blast_output = search_ISfinder(potential_sequence)

        ##get filter blast results:
        resultISFinder[chrompos] = filter_blast(blast_output, args.identity, args.evalue, \
                                                args.length, args.alignment)

    ##write result
    write_result(args.output, mergedPanISa, args.length, resultISFinder)

