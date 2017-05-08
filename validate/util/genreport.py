from os import listdir, system

class GenReport():
    def __init__(self):
        self.actGenome = []
        self.actLen_ranDR = []
        self.actDR = []
        self.actISPos = []
        self.actIR = []
        self.actAmount = None
        self.predChrom = []
        self.predDR = []
        self.predLen_DR = []
        self.predISPos_L = []
        self.predISPos_R = []
        self.predIR = []
        self.predClip_L = []
        self.predClip_R = []
        self.predAmount = None
   

    def __mapbyISPos(self,actKey,actVal,predKey,predVal):
        """[map actual IS postion to predicted IS position with identical mapping]"""
        actList = []
        predList = []
        ##[transform list to dict]
        actDict = dict([(actKey[i], actVal[i]) for i in range(len(actKey))])
        predDict = dict([(predKey[j], predVal[j]) for j in range(len(predKey))])
        ##[map value by postion with left and rigth join]
        actList = [(actKey[i], predDict.get(actKey[i])) for i in range(len(actKey))]
        predList = [(predKey[j], actDict.get(predKey[j])) for j in range(len(predKey))]

        yield actList
        yield predList


    def __mapbyNearbyISPos(self,actKey,actVal,predKey,predVal):
        """[map actual IS postion to predicted IS position with nearby mapping; \
        generate range >> (self.predISPos_R,self.predISPos_R+self.predDR+1)]"""
        actList = []
        predList = []
        pKey = []
        pVal = []
        ##[generate predKey interval]
        for p in range(len(predKey)):
            if self.predDR[p] == 0:
                pKey.extend([range(predKey[p],predKey[p]+self.predDR[p]+1)])
                pVal.extend([x for x in [predVal[p]] for repeat in range(1)])
            else:
                pKey.extend([range(predKey[p],predKey[p]+self.predDR[p])])
                pVal.extend([x for x in [predVal[p]] for repeat in range(self.predDR[p])])
        pKey = [y for x in pKey for y in x]
        ##[transform list to dict]
        actDict = dict([(actKey[i], actVal[i]) for i in range(len(actKey))])
        predDict = dict([(pKey[j], pVal[j]) for j in range(len(pKey))])
        ##[map value by postion with left and rigth join]
        actList = [(actKey[i], predDict.get(actKey[i])) for i in range(len(actKey))]
        predList = [(pKey[j], actDict.get(pKey[j])) for j in range(len(pKey))]
        ##[remove some nearest predList data]
        for pd in predList:
            if pd[1] is None:
                predList.remove(pd)

        yield actList
        yield predList


    def __confusionmatrix(self,caseVal,case):
        """[count number of actual&predict cases in confusion matrix]"""
        itemcase = []

        if case == "actual_yes":
            for p in caseVal:
                if p[1] != 0 and p[1] is not None:## p[1] == 1
                    itemcase.append(p[0])
            return sorted(list(set(itemcase)))
        elif case == "actual_no":
            for p in caseVal:
                if p[1] == 0:
                    itemcase.append(p[0])
            return sorted(list(set(itemcase)))
        elif case == "predict_yes":
            for a in caseVal:
                if a[1] != 0 and a[1] is not None:## a[1] == 1
                    itemcase.append(a[0])
            return sorted(list(set(itemcase)))
        elif case == "predict_no":
            for a in caseVal:
                if a[1] == 0:
                    itemcase.append(a[0])
            return sorted(list(set(itemcase)))         


    def __countclassifyterm(self,actVal,predVal,term):
        """[count number of TP/TN/FP/FN items in confusion matrix]"""
        itemterm = []

        if term == "TP":##[count item of True Positive case]
            actual = self.__confusionmatrix(predVal,"actual_yes")
            predict = self.__confusionmatrix(actVal,"predict_yes")
        elif term == "TN":##[count item of True Negative case]
            actual = self.__confusionmatrix(predVal,"actual_no")
            predict = self.__confusionmatrix(actVal,"predict_no")
        elif term == "FP":##[count item of False Positive case]
            actual = self.__confusionmatrix(predVal,"actual_no")
            predict = self.__confusionmatrix(actVal,"predict_yes")
        elif term == "FN":##[count item of False Negative case]
            actual = self.__confusionmatrix(predVal,"actual_yes")
            predict = self.__confusionmatrix(actVal,"predict_no")

        for item in actual:
            counter = predict.count(item)
            if counter >= 1:
                itemterm.append(item)
        return len(itemterm)


    def __readData(self,actualFile,predictedFile):
        ##[read actual IS file]
        read_actfile = open(actualFile,"r")
        actData = read_actfile.readlines()[1:]
        self.actAmount = len(actData)
        read_actfile.close()
        ##[read predicted IS file]
        read_predfile = open(predictedFile,"r")
        predData = read_predfile.readlines()[1:]
        self.predAmount = len(predData)
        read_predfile.close()

        for act in actData:
            self.actGenome.append(act.split("\t")[0])
            self.actDR.append(act.split("\t")[5])
            self.actLen_ranDR.append(1) if int(act.split("\t")[5]) > 0 else self.actLen_ranDR.append(0)
            self.actISPos.append(int(act.split("\t")[8]))
            self.actIR.append(1) if len(act.split("\t")[9]) > 1 else self.actIR.append(0)

        for pred in predData:
            self.predChrom.append(pred.split("\t")[0])
            self.predDR.append(len(pred.split("\t")[3]))
            self.predLen_DR.append(1) if len(pred.split("\t")[3]) > 0 else self.predLen_DR.append(0)
            self.predISPos_L.append(int(pred.split("\t")[1]))
            self.predISPos_R.append(int(pred.split("\t")[4]))
            self.predIR.append(1) if pred.split("\t")[6] != "No IR" else self.predIR.append(0)
            self.predClip_L.append(int(pred.split("\t")[2]))
            self.predClip_R.append(int(pred.split("\t")[5]))


    def processReport(self,actualFile,predictedFile,outfile):
        ##[read data]
        self.__readData(actualFile,predictedFile)

        ##[prepare list of IS position by map with nearest IS position]
        map_actISPos, map_predISPos = self.__mapbyNearbyISPos(self.actISPos,self.actISPos,\
            self.predISPos_R,self.predISPos_R)
        ##[prepare list of IR status by map with nearest IS position]
        map_actIR, map_predIR = self.__mapbyNearbyISPos(self.actISPos,self.actIR,\
            self.predISPos_R,self.predIR)
        ##[prepare list of DR status by map with nearest IS position]
        map_actDR, map_predDR = self.__mapbyNearbyISPos(self.actISPos,self.actLen_ranDR,\
            self.predISPos_R,self.predLen_DR)
    
        ##[write detail to file]
        detailData, FP = self.__saveDetail(map_actISPos)
        self.__writeData(outfile,detailData,"w")

        ##[prepare summary part]
        ##[count correct(T) and wrong(F) cases]
        T_ISpos = self.__countclassifyterm(map_actISPos,map_predISPos,"TP")
        T_IR = self.__countclassifyterm(map_actIR,map_predIR,"TP")+self.__countclassifyterm(map_actIR,map_predIR,"TN")
        T_DR = self.__countclassifyterm(map_actDR,map_predDR,"TP")+self.__countclassifyterm(map_actDR,map_predDR,"TN")
        F_ISpos = FP+(self.actAmount-T_ISpos)
        F_IR = self.__countclassifyterm(map_actIR,map_predIR,"FP")+self.__countclassifyterm(map_actIR,map_predIR,"FN")
        F_DR = self.__countclassifyterm(map_actDR,map_predDR,"FP")+self.__countclassifyterm(map_actDR,map_predDR,"FN")

        ##[find min and max of clipread]
        minPredClip_L = min(self.predClip_L)
        minPredClip_R = min(self.predClip_R)
        maxPredClip_L = max(self.predClip_L)
        maxPredClip_R = max(self.predClip_R)

        ##[write summary to file]       
        summaryData = [("Correct IS: ",T_ISpos),("Correct IR: ",T_IR),("Correct DR: ",T_DR),\
        ("Wrong IS: ",F_ISpos),("Wrong IR: ",F_IR),("Wrong DR: ",F_DR),("Min.L-clip: ",minPredClip_L),\
        ("Min.R-clip: ",minPredClip_R),("Max.L-clip: ",maxPredClip_L),("Max.R-clip: ",maxPredClip_R)]

        summaryReport = self.__saveSummary(summaryData)
        self.__writeData(outfile,summaryReport,"a")


    def __bittochar(self,bit):
        if bit == 1:
            return "Y"
        elif bit == 0:
            return "N"


    def __saveDetail(self,ISdata):
        formDetail = [(["Chromosome","L-clip","R-clip",\
            "L-position","R-position","actual-position","p-IR","actual-IR","p-DR",\
            "actual-DR",-1])]
        pid_found = []

        ##[case report IS according to actual data]
        for i,d in enumerate(ISdata):
            akey = d[0] ##actual IS position; is always in index[0]
            pkey = d[1] ##predicted IS position
            aid = self.actISPos.index(akey)

            if pkey in self.predISPos_R:
                pid = self.predISPos_R.index(pkey)
                pid_found.append(pid)

                formDetail.append([str(self.predChrom[pid]),\
                    str(self.predClip_L[pid]),str(self.predClip_R[pid]),\
                    str(self.predISPos_L[pid]),str(self.predISPos_R[pid]),str(akey),\
                    self.__bittochar(self.predIR[pid]),self.__bittochar(self.actIR[aid]),\
                    str(self.predDR[pid]),str(self.actDR[aid]),(aid)])
            elif pkey not in self.predISPos_R:
                formDetail.append(["Not Found","","","","",\
                    str(akey),"",self.__bittochar(self.actIR[aid]),"",str(self.actDR[aid]),(aid)])
        
        ##[case report IS according to panISa; only False Positive IS]
        all_pid = set(range(0,len(self.predISPos_R )))
        FP_pid = sorted(all_pid-set(pid_found))

        for mispid in FP_pid:
            formDetail.append([str(self.predChrom[mispid]),\
                    str(self.predClip_L[mispid]),str(self.predClip_R[mispid]),\
                    str(self.predISPos_L[mispid]),str(self.predISPos_R[mispid]),"False positive",\
                    self.__bittochar(self.predIR[mispid]),"",\
                    str(self.predDR[mispid]),"",(mispid)])

        ##[sort list by index(last column)]
        sort_formDetail = sorted(formDetail, key=lambda x:x[-1])
        ##[delete last column]
        formDetail = []
        for rm in sort_formDetail:
            formDetail.append("\t".join(rm[0:-1]))
            
        yield "\n".join(formDetail)
        yield len(FP_pid)


    def __saveSummary(self,data):
        formSummary = ['[==========summary===========]']
        for d in data:
            formSummary.append("\t".join([d[0],str(d[1])]))
        return "\n".join(formSummary)


    def __writeData(self,outfile,data,wtype):
        repfile = open(outfile,wtype)
        repfile.write(data+"\n\n")
        repfile.close()