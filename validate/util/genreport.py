from os import listdir, system


class GenReport():
    def __init__(self):
        self.simData = []
        self.panData = []
        self.acc_no = None


    def __readData(self,simulatedFile,panisaoutFile):
        ##[read actual IS file]
        read_simfile = open(simulatedFile,"r")
        simData_read = read_simfile.readlines()
        read_simfile.close()
        ##[read predicted IS file]
        read_panfile = open(panisaoutFile,"r")
        panData_read = read_panfile.readlines()
        read_panfile.close()

        ##[insert into self.simData]
        simkey = ["isname","ispos","ir","dr"]
        if len(str(simData_read[0]).split("\t")) != 10:
            raise Exception("Found some problems in simulated IS file")
        elif len(str(simData_read[0]).split("\t")) == 10:
            simData_read = simData_read[1:]
            for sim in simData_read:
               isname = sim.split("\t")[1]
               dr = int(sim.split("\t")[5])
               ispos = int(sim.split("\t")[8])
               ir = sim.split("\t")[9].split("\n")[0] \
               if "\n" in sim.split("\t")[9] else sim.split("\t")[9]

               simdict = dict(zip(simkey,[isname,ispos,ir,dr]))
               self.simData.append(simdict)

        ##[insert into self.panData]
        pankey = ["isposL","isposR","clipL","clipR","seqL","seqR","ir","dr"]
        if len(str(panData_read[0]).split("\t")) != 9:
            raise Exception("Found some problems in panISa output")
        elif len(str(panData_read[0]).split("\t")) == 9:
            panData_read = panData_read[1:]
            for pan in panData_read:
                self.acc_no = pan.split("\t")[0]
                isposL = int(pan.split("\t")[1])
                clipL = int(pan.split("\t")[2])
                dr = int(len(pan.split("\t")[3])) \
                if "No sequence to build direct repeat" not in pan.split("\t")[3] else None
                isposR = int(pan.split("\t")[4])
                clipR = int(pan.split("\t")[5])
                ir = "Y" if "No IR" not in pan.split("\t")[6] else "N"
                seqL = pan.split("\t")[7]
                seqR = pan.split("\t")[8].split("\n")[0]

                pandict = dict(zip(pankey,[isposL,isposR,clipL,clipR,seqL,seqR,ir,dr]))
                self.panData.append(pandict)


    def __assignCase(self):
        match_panIndex = []
        assignData = []
        pseudo_index = []

        ##[case ISs exist in simulation]
        for s in self.simData:
            for p in self.panData:
                ##[generate interval of (panISa) left&right position]
                minpos = min(p["isposL"],p["isposR"])
                maxpos = max(p["isposL"],p["isposR"])
                if p["dr"] == 0 or p["dr"] is None:
                    range_panpos = range(minpos,maxpos+1)
                elif p["dr"] > 0:
                    range_panpos = range(minpos,minpos+p["dr"]+1)

                ##[case match; ISs exist both simulation and panISa]
                if s["ispos"] in range_panpos:
                    case = "both"
                    assignData.append([case,self.acc_no,s["isname"],str(s["ispos"]),\
                        str(s["dr"]),s["ir"],str(p["isposL"]),str(p["clipL"]),\
                        str(p["dr"]),str(p["isposR"]),str(p["clipR"]),p["seqL"],\
                        p["seqR"],p["ir"],self.simData.index(s)])
                    match_panIndex.append(self.panData.index(p))
                    pseudo_index.append((self.simData.index(s),self.panData.index(p)))
                    break
            ##[case simulation; ISs exist only simulation]
            else:
                case = "simulation"
                assignData.append([case,self.acc_no,s["isname"],str(s["ispos"]),\
                    str(s["dr"]),s["ir"],"","","","","","","","",self.simData.index(s)])

        ##[case panISa; ISs exist only panISa]
        total_panData = set(range(0,len(self.panData)))
        mismatch_panIndex = sorted(total_panData-set(match_panIndex))
        for mm in mismatch_panIndex:
            case = "panISa"
            ##[get new index for sorting]
            mmidx = self.__makeIndex(pseudo_index,mm)
            assignData.append([case,self.acc_no,"","","","",\
                str(self.panData[mm]["isposL"]),str(self.panData[mm]["clipL"]),\
                str(self.panData[mm]["dr"]),str(self.panData[mm]["isposR"]),\
                str(self.panData[mm]["clipR"]),self.panData[mm]["seqL"],\
                self.panData[mm]["seqR"],self.panData[mm]["ir"],mmidx])
        return assignData


    def __makeIndex(self,existList,searchindex):
        ##[make index for sorting from nearest number of existList (default;[(sim,pan),(s,p)])]
        rsindex = min(existList, key=lambda x:abs(x[1]-searchindex)) ##x[1] is pan.index
        newindex = rsindex[0] ##newindex[0] is sim.index
        if searchindex > rsindex[1]:
            newindex = newindex+0.1
        elif searchindex < rsindex[1]:
            newindex = newindex-0.1
        return newindex


    def __formatReport(self,unformdata):
        formData = []
        ##[sort list by index(last column)]
        unformdata = sorted(unformdata, key=lambda x:x[-1])
        ##[delete last column]
        for rm in unformdata:
            formData.append(rm[0:-1])
        return formData


    def processReport(self,simulatedFile,panisaoutFile,outplace):
        ##[read data]
        self.__readData(simulatedFile,panisaoutFile)
        ##[assign case; both or simulation or panISa]
        reportList = self.__assignCase()
        ##[format report; sort data by position]
        formatList = self.__formatReport(reportList)
        ##[write output]
        self.__writeData(outplace,formatList)


    def __writeData(self,outplace,reportdata):
        outplace.write("\t".join(["Case","Accession No.","IS Name","Sim.position",\
            "DR Length","IR","panISa L-position","L-clip","DR Length","panISa R-position",\
            "R-clip","L-Seq.","R-Seq.","IR"]) + "\n")

        for line in reportdata:
            outplace.write("\t".join(line)+"\n")