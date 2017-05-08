import os
import shutil
from os import listdir

def makedir(newDir):
    if not os.path.exists(newDir):
        os.makedirs(newDir)

def removedir(dirName):
    shutil.rmtree(dirName,ignore_errors=True)

def rmexceptfile(dirName,fileName):
    for f in listdir(dirName):
        if f != fileName:
            os.remove(dirName+f)

def copyfile(oldDir,newDir):
    shutil.copyfile(oldDir,newDir)