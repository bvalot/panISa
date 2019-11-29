class Pos:
    """A simple gff class"""
    def __init__(self, bedline):
        self.chro = bedline[0]
        self.feature =  bedline[2]
        self.start = int(bedline[3])
        self.end = int(bedline[4])
        self.gene = bedline[8]
        self.strand = bedline[6]
        self.offset = 0
        if bedline[7] != ".":
            self.offset = int(bedline[7])
        self.annot = {h.split("=")[0]:h.split("=")[1] for h in self.gene.rstrip(";").split(";")}
        
    def is_in(self, chro, pos):
        if chro != self.chro:
            return False
        elif pos>self.start and pos<self.end:
            return True
        else:
            return False

    def is_closer_previous(self, chro, pos, otherpos):
        """Defined if a feature is closer to a position than an other one""" 
        if chro != self.chro:
            return False
        if pos > self.start:
            if otherpos is None:
                return True
            if self.start > otherpos.start:
                return True
        return False

    def is_closer_next(self, chro, pos, otherpos):
        """Defined if a feature is closer to a position than an other one""" 
        if chro != self.chro:
            return False
        if pos < self.end:
            if otherpos is None:
                return True
            if self.end < otherpos.end:
                return True
        return False

    def get_annotation(self, key):
        if key == "Start":
            return str(self.start)
        elif key == "End":
            return str(self.end)
        elif key == "Strand":
            return self.strand
        elif key == "Annotation":
            return self.gene
        else:
            return self.annot.get(key, "")

    def __str__(self):
        return self.chro+"|"+str(self.start)+"_"+str(self.end)
        

def read_gff(fi, feature):
    gffs = []
    for line in fi.readlines():
        if "##FASTA" in line:
            break
        if line[0] == "#":
            continue
        pos = Pos(line.rstrip("\n").split("\t"))
        if pos.feature == feature:
            gffs.append(pos)
    return gffs
