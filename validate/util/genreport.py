from os import listdir, system

class GenReport():
    def __init__(self):
        self.actGenome = []
        self.actLen_ranDR = []
        self.actISPos = []
        self.actIR = []
        self.actAmount = None
        self.predChrom = []
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
        # ##[map value by postion with left and rigth join]
        actList = [(actKey[i], predDict.get(actKey[i])) for i in range(len(actKey))]
        predList = [(predKey[j], actDict.get(predKey[j])) for j in range(len(predKey))]

        yield actList
        yield predList


    def __mapbyNearbyISPos(self,actKey,actVal,predKey,predVal,acceptRange):
        """[map actual IS postion to predicted IS position with nearby mapping; define >> acceptRange]"""
        actList = []
        predList = []
        ##[generate predKey interval]
        predKey[:] = [range(p-acceptRange,p+acceptRange+1) for p in predKey]
        predKey = [y for x in predKey for y in x]
        # predKey = sorted(list(set(predKey)))
        predVal = [x for x in predVal for repeat in range(acceptRange*2+1)]
        ##[transform list to dict]
        actDict = dict([(actKey[i], actVal[i]) for i in range(len(actKey))])
        predDict = dict([(predKey[j], predVal[j]) for j in range(len(predKey))])
        ##[map value by postion with left and rigth join]
        actList = [(actKey[i], predDict.get(actKey[i])) for i in range(len(actKey))]
        predList = [(predKey[j], actDict.get(predKey[j])) for j in range(len(predKey))]
        ##[remove some nearest predList data]
        for pd in predList:
            if pd[0] in range(pd[0]-acceptRange,pd[0]+acceptRange+1) and pd[1] is None:
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
            # self.actLen_ranDR.append(act.split("\t")[5])
            self.actLen_ranDR.append(1) if int(act.split("\t")[5]) > 0 else self.actLen_ranDR.append(0)
            self.actISPos.append(int(act.split("\t")[8]))
            self.actIR.append(1) if len(act.split("\t")[9]) > 1 else self.actIR.append(0)
            # self.actIR.append(1) if "Y" in act.split("\t")[9] else self.actIR.append(0)

        for pred in predData:
            self.predChrom.append(pred.split("\t")[0])
            # self.predLen_DR.append(len(pred.split("\t")[3]))
            self.predLen_DR.append(1) if len(pred.split("\t")[3]) > 0 else self.predLen_DR.append(0)
            self.predISPos_L.append(int(pred.split("\t")[1]))
            self.predISPos_R.append(int(pred.split("\t")[4]))
            self.predIR.append(1) if pred.split("\t")[6] != "No IR" else self.predIR.append(0)
            self.predClip_L.append(int(pred.split("\t")[2]))
            self.predClip_R.append(int(pred.split("\t")[5]))


    def processReport(self,actualFile,predictedFile,acceptRange,outfile):
        ##[read data]
        self.__readData(actualFile,predictedFile)

        ##[prepare list of IS position by map with nearest IS position]
        map_actISPos, map_predISPos = self.__mapbyNearbyISPos(self.actISPos,self.actISPos,\
            self.predISPos_L+self.predISPos_R,self.predISPos_L+self.predISPos_R,acceptRange)
        ##[prepare list of IR status by map with nearest IS position]
        map_actIR, map_predIR = self.__mapbyNearbyISPos(self.actISPos,self.actIR,\
            self.predISPos_L+self.predISPos_R,self.predIR+self.predIR,acceptRange)
        ##[prepare list of IR status by map with nearest IS position]
        map_actDR, map_predDR = self.__mapbyNearbyISPos(self.actISPos,self.actLen_ranDR,\
            self.predISPos_L+self.predISPos_R,self.predIR+self.predLen_DR,acceptRange)
    
        ##[calculate sensitivity: TP/actual_yes]
        sens_ISpos = float(self.__countclassifyterm(map_actISPos,map_predISPos,"TP"))/self.actAmount
        sens_IR = float(self.__countclassifyterm(map_actIR,map_predIR,"TP"))/sum(self.actIR)
        sens_DR = float(self.__countclassifyterm(map_actDR,map_predDR,"TP"))/sum(self.actLen_ranDR)
        
        ##[calculate precision: TP/predicted_yes]
        preci_ISpos = float(self.__countclassifyterm(map_actISPos,map_predISPos,"TP"))/self.predAmount
        preci_IR = float(self.__countclassifyterm(map_actIR,map_predIR,"TP"))/sum(self.predIR)
        preci_DR = float(self.__countclassifyterm(map_actDR,map_predDR,"TP"))/sum(self.predLen_DR)

        ##[calulate false positive rate: FP/actual_no]
        FP_ISpos = self.predAmount-(self.__countclassifyterm(map_actISPos,map_predISPos,"TP"))
        FPrate_IR = float(self.__countclassifyterm(map_actIR,map_predIR,"FP"))/(self.actAmount-sum(self.actIR))
        FPrate_DR = float(self.__countclassifyterm(map_actDR,map_predDR,"FP"))/(self.actAmount-sum(self.actLen_ranDR))

        ##[calculate specificity rate: TN/actual_no]
        speci_IR = float(self.__countclassifyterm(map_actIR,map_predIR,"TN"))/(self.actAmount-sum(self.actIR))
        speci_DR = float(self.__countclassifyterm(map_actDR,map_predDR,"TN"))/(self.actAmount-sum(self.actLen_ranDR))

        ##[find min and max of clipread]
        minPredClip_L = min(self.predClip_L)
        minPredClip_R = min(self.predClip_R)
        maxPredClip_L = max(self.predClip_L)
        maxPredClip_R = max(self.predClip_R)

        ##[write as file]
        fs = self.__FS ##float to string
        report = "\t".join([self.predChrom[0],fs(sens_ISpos),fs(sens_IR),fs(sens_DR),\
            fs(preci_ISpos),fs(preci_IR),fs(preci_DR),str(FP_ISpos),fs(FPrate_IR),fs(FPrate_DR),\
            fs(speci_IR),fs(speci_DR),str(minPredClip_L),str(minPredClip_R),str(maxPredClip_L),str(maxPredClip_R)])
        self.__writeData(outfile,report)


    def __FS(self,value):
        return str(round(value,3))

    def __writeData(self,outfile,data):
        header = "\t".join(["chromosome","sensitivity-pos","sensitivity-IR","sensitivity-DR",\
            "precision-pos","precision-IR","precision-DR","FP-pos","FPR-IR","FPR-DR",\
            "specificity-IR","specificity-DR","min-Lclip","min-Rclip","max-Lclip","max-Rclip"])

        repfile = open(outfile,"w")
        repfile.write("\n".join([header,data]))
        repfile.close()

