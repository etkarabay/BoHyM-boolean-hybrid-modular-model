import sys
import boolean2
from boolean2 import util,state,Model,tokenizer
from boolean2.plde import helper
import matplotlib.pyplot as plt
import pandas as pd


text="""
1:A *= actA
1:B *= A
1:C *= B
1:D *= C
1:E *= D
1:F *= E
1:ITG *= stMT or dyMT or FX or FA or ECM or ppMLC 

"""

global m

#Adding sNodes to main text

newtext=[]
tokens=tokenizer.tokenize(text)

nogeneexpression=["actA"]

"""adding sNodes for nodes have upstream feed"""
for token in tokens:
    if token[1].value in nogeneexpression:
        line = tokenizer.tok2line(token[:3]) +"=" + tokenizer.tok2line(token[4:])
        newtext.append("\n")
        newtext.append(line)       
    else:
        line = tokenizer.tok2line(token[:3]) +"="+ "(" + tokenizer.tok2line(token[4:])+")" + " "+ "and" +" "+"s"+ token[1].value
        newtext.append("\n")
        newtext.append(line)
A=[]
for token in tokens:
   A.append(token[1].value)
"""adding sNodes for nodes do not have upstream feed"""

for i in list(tokenizer.get_nodes(tokens)):
    if i in nogeneexpression:
        pass    
    else:
        if i not in A:
            line = "1:"+" "+ i + " " + "*=" +" "+"s"+i
            newtext.append("\n")
            newtext.append(line)

newtextend=''.join(newtext) #updated text
print(newtextend) 
tokens=tokenizer.tokenize(newtextend)
m=list(tokenizer.get_nodes(tokens)) #allnodesinthelist


#####

import snode_example_initials
init=snode_example_initials.main()



def main():
    for k in range(0,1):           
                d={}
                print ("main module")
                for i in range(len(m)):
                    d.update({m[i]: {'threshold':0, 'conc':init[m[i]][200],'decay':1}})
                    if m[i]=="actA" :
                         d.update({m[i]: {'threshold':0 , 'conc':1, 'decay':1}})
                
                                 
                CONC=d
                
       
                model = Model( text=newtextend, mode='plde')
                model.initialize(missing = helper.initializer(CONC))
                model.iterate( fullt=200, steps=200,ECMval=0.1)
                data = [model.data]

  

               
                plt.plot(model.data["A"],label="A")
                plt.plot(model.data["B"],label="B")
                plt.plot(model.data["C"],label="C")
   

                
                leg = plt.legend(loc='lower right')                                 
                plt.show()      
  
    
    
if __name__ == "__main__":
    main()


 
