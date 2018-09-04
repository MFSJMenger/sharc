#!/usr/bin/env python
from __future__ import print_function

import sys
from general_tools import readFile, writeOutput
from general_tools import NFind, getValues, readFile, scaleCoords
from general_tools import getMoldenAtoms, AtomNumberToName


def getCoordinates(fileName,NAtoms):
    text=readFile(fileName)
    result = getValues(text,"Input orientation",71*NAtoms,322,[1,3,4,5],Vartype=[int,float,float,float])
    if result is not None:
        return result
    else:
        result = getValues(text,"Standard orientation",71*NAtoms,324,[1,3,4,5],Vartype=[int,float,float,float])
        return result

def getNAtoms(fileName):
    with open(fileName,'r') as f:
        text = f.read()
        NAtoms = int( (NFind(text,'NAtoms=',40,0).split("=")[1]
                     ).split()[0])
    return NAtoms


def getGauFreq(text,NAtoms):
    Numbers=map(int,text.splitlines()[0].split())
    shift = 3 
    for i,line in enumerate(text.splitlines()[2:10]): 
        if "Frequencies" in line:
            Freq = map(float,line.split()[shift-1:]) # all others have an additional Space in the keyword!
        elif "Red. masses" in line:
            RedMas = map(float,line.split()[shift:])
        elif "Frc consts" in line:
            FrcConst = map(float,line.split()[shift:])
        elif "IR Inten" in line:
            IRInten = map(float,line.split()[shift:])
        elif "Raman Activ" in line:
            RamanAct = map(float,line.split()[shift:])
        elif "Depolar (P)" in line:
            DepolarP = map(float,line.split()[shift:])
        elif "Depolar (U)" in line:
            DepolarU = map(float,line.split()[shift:])
        elif "Atom  AN" in line:
            lineEnd=i+3
            break
        else:
            raise Exception("Something wrong in getGauFreq")
    Freq1=[]
    Freq2=[]
    Freq3=[]
    for line in text.splitlines()[lineEnd:lineEnd+NAtoms]:
        #if line.strip() == "":
        #    continue
        col = line.split()
        Freq1.append([col[1]]+map(float,col[2:5])) 
        Freq2.append([col[1]]+map(float,col[5:8])) 
        Freq3.append([col[1]]+map(float,col[8:])) 
    
    if len(Freq1) != NAtoms:
            print(len(Freq1))
            raise Exception("Error: getGauFreq, len Freq1 != NAtoms")
    return Numbers,Freq,IRInten,[Freq1,Freq2,Freq3]

def readGaussianFreq(fileName,NAtoms,linear=False,noraman=True, gauversion='g09'):
    text = readFile(fileName)
    if not linear:
        N = (NAtoms-2) 
    else:
        raise Exception("linear Molecules not supported")
    if noraman:
        scale = 4
    else:
        scale = 7
    if gauversion == 'g16':
        val = (2*70+scale*74+78+NAtoms*80)
    elif gauversion == 'g09d01':
        val = (2*70+scale*74+78+NAtoms*80)
    elif gauversion == 'g09':
        val = (2*70+scale*74+78+NAtoms*80)
    elif gauversion == 'g09c01':
        val = (2*69+scale*73+78+NAtoms*80)
    chars = val*N
    Freq = NFind(text,"and normal coordinates:",chars,24)

    Freqs=[ Freq[i*val:(i+1)*val]  for i in range(N) ]
    Numbers = []
    IRInten = []
    GauOut={
        'Numbers' : [],
        'IRInt'   : [],  # IR intensity
        'FreqVal' : [],  # Frequency values
        'FreqDis' : []   # Displacement for Frequencies
        }
    for Freq in Freqs:
        Num,Val,IRIn,FreqDist=getGauFreq(Freq,NAtoms)
        GauOut['Numbers']   += Num
        GauOut['IRInt']     += IRIn
        GauOut['FreqVal']   += Val
        GauOut['FreqDis']   += FreqDist
    return GauOut

def convertGauToMolden(inFile,outFile,Coords,NAtoms,noraman=True,linear=False,unit='Angs', gauversion='g09'):
    formFreq     = "% 14.10F\n"
    formFreqCoord = " %-5s % 14.10F % 14.10F % 14.10F\n"
    formNormCoord= " % 14.10F % 14.10F % 14.10F\n"
    if noraman:
        formInt      = "% 12.4F 0.0000 \n"
    else:
        raise Exception("Raman not included, yet")
   
     
    GauOut=readGaussianFreq(inFile,NAtoms,linear, noraman, gauversion)
    Coords=getCoordinates(inFile,NAtoms)
    if unit.lower() == 'angs':
        Coords=scaleCoords(Coords)

    out = "[Molden Format]\n [Atoms] AU \n" 
    for i,coord in enumerate(Coords):
        out += getMoldenAtoms(i+1,coord)
    
#    out = "[Molden Format]\n [GEOMETRIES] XYZ \n %d\n %s \n" % (NAtoms, inFile)
#    for coord in Coords:
#        out += getXYZ(coord)


    out  += "[FREQ]\n"
    for Freq in GauOut['FreqVal']:
        out += formFreq % Freq

    out += "[INT]\n" 
    for inten in GauOut['IRInt']:
        out += formInt % inten
    
    out += "[FR-COORD]\n"
    for xyz in Coords:
        out += formFreqCoord % (AtomNumberToName[xyz[0]],xyz[1],xyz[2],xyz[3])

    out += "[FR-NORM-COORD]\n" 
    for i,Freqs in enumerate(GauOut['FreqDis']):
        out += "vibration %d \n" % (i+1)
        for xyz in Freqs:
            out += formNormCoord % (xyz[1],xyz[2],xyz[3])

    
    writeOutput(outFile,out)
    return

def convertGauToMoldenFreq(inFile, outFile, cfile, noraman=True, linear=False, unit='Angs', version='g09'):
    """Convert Gaussian log file information to Molden format"""
    # get natoms
    NAtoms = getNAtoms(inFile)
    # get coordinates
    Coords=getCoordinates(inFile,NAtoms)
    # call the main Routine
    convertGauToMolden(inFile,outFile,Coords,NAtoms,noraman=noraman,linear=linear, unit=unit, gauversion=version)
    return

def getCommandoLine():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("fileName", metavar="FileName",type=str,
                        help="Gaussian Frequency file")
    parser.add_argument("-o","--output", metavar="outputfile",type=str,
                        default='freq.molden',
                        help="Output file, if not set it is **freq.molden**!")
    parser.add_argument("-v","--version", metavar="Version",type=str,
                        choices=['g09','g09c01','g09d01','g16'],
                        default='g16',
                        help="Version of the Gaussian code, used!")
    args = parser.parse_args()

    return args.fileName, args.output, args.version


if __name__ == "__main__":
    infile, outfile, version= getCommandoLine()
    convertGauToMoldenFreq(infile,outfile,None,noraman=True,linear=False, version=version)
