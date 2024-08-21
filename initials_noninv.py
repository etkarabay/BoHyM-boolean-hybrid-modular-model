# -*- coding: utf-8 -*-
import sys
sys.path.append('/Users/esratiftik/Documents/')
sys.path.append("/Users/esratiftik/Desktop")
import openpyxl
from pathlib import Path
import boolean2
from boolean2 import util,state,Model,tokenizer
import networkx as nx
import random
from boolean2.plde import helper
import matplotlib.pyplot as plt
import pandas as pd
import xlrd
import math
import numpy as np
from openpyxl import load_workbook
from itertools import permutations
import statistics as st
import csv, StringIO
import string
from itertools import *


#bcatenin=freebcatenin

text="""
1:APC *= IQGAP
1:Cdc42 *= Rac or PIX or Src or cPKC
1:IQGAP *= nPKC or Cdc42 or Rac or Barrestin1
1:CLIP *= MAP1b or not GSK or IQGAP
1:CLASP*= not GSK or IQGAP
1:F-Actin *= Profilin or Actinin or Arp or not Cofilin or not Gelsolin or mDia1 or mDia2 or Vinculin or Talin
1:p190GEF *= FAK
1:Fyn *= ITG
1:LARG *= not PIX or Fyn
1:Paxilin *= FAK or Crk or not GIT1
1:p130Cas *= FAK 
1:Grb2 *= Shc or FAK
1:Src *= FAK or Barrestin2 or Pyk2
1:pdzRhoGEF *= Barrestin1 or Pyk2 or FAK
1:p115RhoGEF *= Barrestin2 or not PIX or not PAK
1:RhoC *= pdzRhoGEF 
1:Cortactin *= Src or Rac or RhoC 
1:p190RhoGap *= Src or Rac
1:PI3K *= RTK or (FAK) or HRas or Barrestin1 or Rap1 or Gbg or Src or CaM
1:p114RhoGEF *= Gbg
1:Rac *=  DOCK180 or Vav or TIAM or PIX or Rac1GEF or Asef  or not FilGAP or not RacGAP1 or p114RhoGEF 
1:Rap1 *= GRF2 or EPAC or PKD
1:RacGAP1 *= Akt
1:SSH *= not GSK or not PKD or Ca or F-Actin 
1:Vav *= PIP3 or Vimentin or Rap1
1:Raf1 *= HRas or Rap1 or not PKA or cPKC
1:cPKC *= Ca and DAG 
1:nPKC  *= DAG
1:Pyk2 *= cPKC or Ca
1:RhoA *=  not PKA or GEFH1 or pdzRhoGEF or p115RhoGEF or p114RhoGEF or not p190RhoGap or p190GEF or LARG or ERK or nPKC
1:Talin *= PKA or cPKC or Rap1 or PIP2 or not Calpain
1:PKD *= nPKC
1:Vimentin *= nPKC or Filamin
1:ppMLC *=  not MLCP or ROCK or PAK or MLCK 
1:PIX *= PAK or not ROCK or PKA
1:RyR *= PKA or Ca 
1:SERCA *= PKA
1:MLCK *= not PAK or not PKA
1:stathmin *= (not PKA or not PAK or not ERK or not CaMK) 
1:PMCA *= CaM or PKA 
1:FAK *= RTK or ((FX or FA ) and PIP2) or PKA or not PTEN  or Src
1:RasGRF *= CaM or PKA or Gbg or FAK or not Barrestin1
1:Calpain *= Ca or not PKA
1:MAP1b *= not GSK 
1:EB *=  not GSK 
1:Vinculin *= Talin 
1:Shc *= RTK or (FAK) or Fyn
1:Gs *= BAR2
1:BAR2 *= actBAR2
1:Gbg *= BAR2 
1:Barrestin1 *= BAR2 
1:Barrestin2 *= BAR2 
1:AC *= Gs or FF
1:cAMP *= AC  
1:EPAC *= cAMP
1:Orai1 *= STIM 
1:STIM *= not SERCA or IP3R or RyR 
1:Ca *= (Orai1 ) or not SERCA or (RyR)  or (IP3R ) or (TRPC ) or (TRPV ) or (TRPM7) or (not PMCA )
1:ROCK *= RhoA or RhoC
1:GSK *= not Akt 
1:Actinin *= not PIP2 
1:ERK *= MEK1  
1:DNA *= ERK or mTOR  or JNK 
1:mTOR*= Akt
1:Crk *= p130Cas
1:DOCK180 *= Crk
1:GRF2 *= Crk
1:PAK *= (Rac or Cdc42)
1:mDia1 *= RhoA 
1:TIAM *= PIP3 or Rap1
1:LIMK *= PAK or ROCK or RhoC 
1:Wave *= Rac
1:Arp *= Wave or Wasp or Cortactin
1:Wasp *= Cdc42 
1:JNK *= Rap1 or Rac
1:SOS *= Grb2
1:HRas *= SOS or RasGRF or Src
1:bRaf *= Rap1
1:MEK1 *= Raf1 or PAK or bRaf 
1:ILK *= PIP3 
1:PDK1 *= PIP3
1:Akt *= PIP3 or PDK1 or ILK
1:IP3 *= PIP2 and (PLCg or PLCe or PLCb) 
1:PLCe *= Rap1
1:DAG *= (PLCg or PLCe or PLCb) and PIP2
1:Rac1GEF *= HRas
1:Gelsolin *= not PIP2 
1:PKA *= cAMP 
1:GEFH1 *= dyMT
1:FilGAP *= ROCK
1:Cofilin *= not LIMK or SSH 
1:PIP5K *= RhoA 
1:MLCP *= not ROCK
1:PLCg *= RTK
1:Profilin *= mDia1 
1:mDia2 *= Cdc42
1:GIT1 *= PIX
1:dyMT *= stathmin or CLIP or CLASP or EB or Gs or Ca or Paxilin
1:stMT *= (EB and mDia1 and APC) or LL5B or Tau or FA or not LIMK 
1:Tau *= Cdc42 or Rac or not GSK
1:LL5B *= PIP3 
1:ERM *= ROCK
1:PIP2 *= PIP5K
1:PIP3 *= (PI3K and PIP2) or not PTEN 
1:Asef *= APC
1:FS *= FS
1:CaMK *= CaM 
1:CaM *= Ca
1:TRPV *= ECM 
1:TRPM7 *= PKA 
1:IP3R *= PKA or IP3 
1:TRPC *= Orai1 
1:PLCb *=Gbg
1:Filamin *= F-Actin
1:FX *= F-Actin and Vinculin and Talin
1:RTK *= actRTK
"""


#Adding sNode



newtext=[]
tokens=tokenizer.tokenize(text)
prem=list(tokenizer.get_nodes(tokens))

nogeneexpression=["Ca","actBAR2","cAMP","FS","FF","actRTK","ECM"]


#nogeneexpression=[]
for token in tokens:
    if token[1].value in nogeneexpression:
        line = tokenizer.tok2line(token[:3]) +"="+ "(" + tokenizer.tok2line(token[4:])+")" 
        newtext.append("\n")
        newtext.append(line)
    
        
        
    else:
        line = tokenizer.tok2line(token[:3]) +"="+ "(" + tokenizer.tok2line(token[4:])+")" + " "+ "and" +" "+"s"+ token[1].value
        newtext.append("\n")
        newtext.append(line)
    

    

A=[]
for token in tokens:
   A.append(token[1].value)


for i in list(tokenizer.get_nodes(tokens)):
    if i in nogeneexpression:
        pass
    
    else:
        if i not in A:
            line = "1:"+" "+ i + " " + "*=" +" "+"s"+i
            newtext.append("\n")
            newtext.append(line)

newtextend=''.join(newtext)


tokens=tokenizer.tokenize(newtextend)
       

m=list(tokenizer.get_nodes(tokens))




sNodes=[]
for i in prem:
    node="s"+ i 
    sNodes.append(node)




#non n
lowexpressed=["sPMCA","sSTIM","sOrai1","sTRPC","sTRPM","sTRPV","sIP3R","sRYR","scPKC"]
#lowexpressed=["sSTIM","sOrai1","sTRPC","sTRPM","sTRPV","sIP3R"]




def main(ECMval):
        d={}
        global m

        for i in range(len(m)):
            d.update({m[i]: {'threshold':0, 'conc':0,'decay':1}})#fp.loc[m[i]]["InitialStates"]
            
            if m[i] in sNodes:
                 d.update({m[i]: {'threshold':0 , 'conc':0.01, 'decay':1}})
            if m[i] in lowexpressed:
                 d.update({m[i]: {'threshold':0 , 'conc':0.01, 'decay':1}})
            if m[i]=="ITG":
                 d.update({m[i]: {'threshold':0 , 'conc':1 , 'decay':1}})
            if m[i]=="FX":
                 d.update({m[i]: {'threshold':0 , 'conc':0 , 'decay':1}})
            if m[i]=="FA":
                 d.update({m[i]: {'threshold':0 , 'conc':0 , 'decay':1}})
            if m[i]=="sFX":
                 d.update({m[i]: {'threshold':0 , 'conc':0 , 'decay':1}})
            if m[i]=="sFA":
                 d.update({m[i]: {'threshold':0 , 'conc':0 , 'decay':1}})
            if m[i]=="sITG":
                 d.update({m[i]: {'threshold':0 , 'conc':0.1 , 'decay':1}})
            if m[i]=="stMT":
                 d.update({m[i]: {'threshold':0 , 'conc':1, 'decay':1}})
            if m[i]=="dyMT":
                 d.update({m[i]: {'threshold':0 , 'conc':0, 'decay':1}})
            if m[i]=="sstMT":
                 d.update({m[i]: {'threshold':0 , 'conc':0.1, 'decay':1}})
            if m[i]=="sdyMT":
                 d.update({m[i]: {'threshold':0 , 'conc':0, 'decay':1}})
            if m[i]=="ECM":
                 d.update({m[i]: {'threshold':0 , 'conc':ECMval, 'decay':1}})
            if m[i]=="sMLCK":
                 d.update({m[i]: {'threshold':0 , 'conc':0.001, 'decay':1}})
            
            
     
          
            
          
            
            
           
           
           
        CONC=d
        

        model = Model( text=newtextend, mode='plde')
        
        model.initialize(missing = helper.initializer(CONC))
        model.iterate( fullt=200, steps=200,module=1,ECMval=ECMval)
    
        data = [model.data]

        #plt.plot(model.data["Ca"],color='orange',label="$Ca^{2+}$")
        #plt.plot(model.data["FA"],color='purple',linestyle='dashed',label="FA")
        #print (model.data["FA"][200])
        #plt.plot(((model.data["ppMLC"]-min(model.data["ppMLC"]))/(max(model.data["ppMLC"])-min(model.data["ppMLC"]))),color='blue',linestyle='dashed',label="ppMLC")
        #plt.plot(((model.data["FA"]-min(model.data["FA"]))/(max(model.data["FA"])-min(model.data["FA"]))),color='purple',linestyle='dashed',label="FA")
        #plt.plot(((model.data["F-Actin"]-min(model.data["F-Actin"]))/(max(model.data["F-Actin"])-min(model.data["F-Actin"]))),color='red',label="F-Actin")

        #print model.data["ECM"]

 


        



        

        return model.data
        
            
           
    

        
        
if __name__ == "__main__":
    main()
  

    



    
    
