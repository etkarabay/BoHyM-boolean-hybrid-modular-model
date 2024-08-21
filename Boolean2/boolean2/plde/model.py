# -*- coding: utf-8 -*-
import sys, os

from itertools import *
from scipy.integrate import odeint
from boolean2.boolmodel import BoolModel
from boolean2 import util, odict, tokenizer
from . import helper 
import math
import pandas as pd
import numpy as np
from numpy import diff
import random
from numpy import exp
import imp



try:
    from matplotlib.pyplot import *
    from scipy.integrate import odeint
    #from pylab import arange, rk4
except ImportError:
    util.error( "matplotlib is missing, install it from: http://matplotlib.sourceforge.net/")


def default_override( node, indexer, tokens ):
    """
    Gets called before the generating each equation.
    If this function returns anything other than None 
    it will override the equation
    """
    return None

def default_equation( tokens, indexer ):
    """
    Default equation generator, override this to generate
    other equations without explicit overrides
    """
    node = tokens[1].value
    text = helper.change(node, indexer) + ' = ' + helper.piecewise(tokens, indexer)
  
    return text

def boolmapper (value):
    if type(value) == tuple:
        return value
    else:
        return util.bool_to_tuple(value)

class PldeModel( BoolModel ):
    """
    This class generates python code that will be executed inside 
    the Runge-Kutta integrator.
    """
    def __init__(self, text, mode='plde'):
        
        # run the regular boolen engine for one step to detect syntax errors
        model = BoolModel(text=text, mode='sync')
        model.initialize( missing=util.randbool )
        model.iterate(steps=1)

        # onto the main initialization
        self.INIT_LINE  = helper.init_line
        self.OVERRIDE   = default_override
        self.DEFAULT_EQUATION = default_equation
        self.EXTRA_INIT = ''

        # setting up this engine
        BoolModel.__init__(self, text=text, mode=mode)
        self.dynamic_code = '*** not yet generated ***'
        self.lazy_data = {}
    
    @property
    def data(self):
        "For compatibility with the async engine"
        return self.lazy_data

    def initialize(self, missing=None, defaults={} ):
        "Custom initializer"
        BoolModel.initialize( self, missing=missing, defaults=defaults )
        
        # will also maintain the order of insertion
        self.mapper  = odict.odict() 
        self.indexer = {}
        
        # this will maintain the order of nodes
        self.nodes = list(self.nodes)
        self.nodes.sort()
        for index, node in enumerate(self.nodes):
            triplet = self.first[node]
            self.mapper [node] = ( index, node, triplet )
            self.indexer[node] = index
        
        # a sanity check
        assert self.nodes == list(self.mapper.keys())

  
    def generate_init( self, localdefs ):
        """
        Generates the initialization lines
        """
        init = [ ]
        
        init.extend( self.EXTRA_INIT.splitlines() )
        init.append( '# dynamically generated code' )
        init.append( '# abbreviations: c=concentration, d=decay, t=threshold, n=newvalue' )
        init.append( '# %s' % list(self.mapper.values()) )

        for index, node, triplet in list(self.mapper.values()):
            conc, decay, tresh = boolmapper(triplet)
            #assert decay > 0, 'Decay for node %s must be larger than 0 -> %s' % (node, str(triplet))   
            store = dict( index=index, conc=conc, decay=decay, tresh=tresh, node=node)
            line = self.INIT_LINE( store )
            init.append( line )
        
        if localdefs:
            init.extend( [ '# custom imports', 'import %s' % localdefs, 'reload(%s)' % localdefs, 'from %s import *' % localdefs ]   )

        init_text = '\n'.join( init )
        return init_text

    
       
 
    
    def create_equation( self, tokens ):
        """
        Creates a python equation from a list of tokens.
        """
        original = '#' + tokenizer.tok2line(tokens)
        node  = tokens[1].value
        lines = [ '', original ]
        
        line  = self.OVERRIDE(node, indexer=self.indexer, tokens=tokens)
        line
        if line is None:
            line = self.DEFAULT_EQUATION( tokens=tokens, indexer=self.indexer )
        
        if isinstance(line, str):
            line = [ line ]
        


        lines.extend( [ x.strip() for x in line ] )
        return lines
        
    


    def generate_function(self ):

        """
        Generates the function that will be used to integrate
        """
        sep = ' ' * 4

        indices = [ x[0] for x in self.mapper.values() ]
        
        assign  = [ 'c%d' % i for i in indices ]
       
    
        retvals = [ 'n%d' % i for i in indices ]
        
        zeros   = map(lambda x: '0.0', indices ) 
        assign  = ', '.join(assign)
        retvals = ', '.join(retvals)
        zeros   = ', '.join( zeros )

        body = []
        body.append( 'x0 = %s' % assign )
        body.append( 'def derivs( x, t):' )
         
        body.append( '    %s = x' % assign )
        body.append( '    %s = %s' % (retvals, zeros) )
        for tokens in self.update_tokens:
            equation = self.create_equation( tokens )
            equation = [ sep + e for e in equation ]
            body.append( '\n'.join( equation)  )
        body.append( '' )
        body.append( "    return ( %s ) " % retvals )
        

        text = '\n'.join( body )
    
        
        return text




    #def iterate( self, fullt,steps,kbeta_cons, kgamma_cons1,
            #kdel_cons1,kdel_cons2,kdel_cons3,
        #kepsi_cons1,kepsi_cons2,kepsi_cons3,autogen_fname=None, localdefs=None, autogen='autogen'):
     
    def iterate( self, fullt,steps,ECMval,autogen_fname=None, localdefs=None, autogen='autogen'):

        esra = []
       

     
       
        """
        Iterates over the system of equations 
        """
        
        if autogen_fname is not None:
            autogen = autogen_fname
            del autogen_fname
            util.warn("parameter 'autogen_fname' is deprecated. Use 'autogen' instead." )
        
        # setting up the timesteps
        dt = fullt/float(steps)
        self.t  = [ dt * i for i in range(steps+1) ]


        # generates the initializator and adds the timestep
        
        self.init_text  = self.generate_init( localdefs=localdefs )
        self.init_text += '\ndt = %s' % dt
     

        
        # generates the derivatives
        self.func_text = self.generate_function()
        
       
        self.dynamic_code = self.init_text + '\n' + self.func_text
       
        
        try:
            fp = open( '%s.py' % autogen, 'wt')
            fp.write( '%s\n' % self.init_text )
            fp.write( '%s\n' % self.func_text )
            fp.close()
            autogen_mod = __import__( autogen )
        

            try:
                os.remove( '%s.pyc' % autogen )
            except OSError:
                imp.reload( autogen_mod )
        except Exception as exc:
            msg = "'%s' in:\n%s\n*** dynamic code error ***\n%s" % ( exc, self.dynamic_code, exc )
            util.error(msg)
    
        """BoHyM Extension"""

         
        from boolean2.plde.helper import line2,line3,line4,linew,line5,line6,line7,line8,line9,line10,line11,line12,line13,line14
        from boolean2.plde.helper import seq1,linew,seq2
        diffw=[]
        diffn=[]
        weightchange=[]

        #Node concentrations for rate constants
        concx=[]

        #FAx=[]

        #loop time and steps
        interloopsteps=2000
        loop_time=dt
        timex=np.linspace(0,loop_time,interloopsteps)

        #Initials List
        alldata1a=[] 
        alldata1=[]
        alldata1a.append(list(np.array(autogen_mod.x0)))
        alldata1.append(list(np.array(autogen_mod.x0)))
        alldatay=[0]*len(self.nodes) #initial value for odeint

        
        #initial and total values of FA 
        totalFA=((alldata1a[0][self.nodes.index("sITG")]+alldata1a[0][ self.nodes.index("sFX")]+alldata1a[0][self.nodes.index("sFA")]))
        totalactin_cyt =((alldata1a[0][self.nodes.index("sg_act_cyt")]+alldata1a[0][ self.nodes.index("sg_ox_cyt")]+alldata1a[0][self.nodes.index("sf_act_cyt")]))
        totalactin_nuc = ((alldata1a[0][self.nodes.index("sg_act_nuc")]+alldata1a[0][ self.nodes.index("sg_ox_nuc")]+alldata1a[0][self.nodes.index("sf_act_nuc")]))
        totalMT=(alldata1a[0][self.nodes.index("sdyMT")]+alldata1a[0][ self.nodes.index("sstMT")])
        #print totalFA


        alldata1[0][self.nodes.index("ITG")]=alldata1a[0][self.nodes.index("ITG")]
        alldata1[0][self.nodes.index("FX")]=alldata1a[0][self.nodes.index("FX")]
        alldata1[0][self.nodes.index("FA")]=alldata1a[0][self.nodes.index("FA")]
        alldata1[0][self.nodes.index("g_act_cyt")]=alldata1a[0][self.nodes.index("g_act_cyt")]
        alldata1[0][self.nodes.index("f_act_cyt")]=alldata1a[0][self.nodes.index("f_act_cyt")]
        alldata1[0][self.nodes.index("g_act_nuc")]=alldata1a[0][self.nodes.index("g_act_nuc")]
        alldata1[0][self.nodes.index("f_act_nuc")]=alldata1a[0][self.nodes.index("f_act_nuc")]
        alldata1[0][self.nodes.index("g_ox_nuc")]=alldata1a[0][self.nodes.index("g_ox_nuc")]
        alldata1[0][self.nodes.index("g_ox_cyt")]=alldata1a[0][self.nodes.index("g_ox_cyt")]

           


        #initial and total values of MT 
        
        
        

     
        alldata1[0][self.nodes.index("dyMT")]=alldata1a[0][self.nodes.index("dyMT")]
        alldata1[0][self.nodes.index("stMT")]=alldata1a[0][self.nodes.index("stMT")]



        
            

        #alldata1[0][self.nodes.index("Ca")]=0.5
        alldata1[0][self.nodes.index("ECM")]=ECMval
        
        #Internal loops 
        for i in range(1,len(self.t)):
          

    
            #Extracting node values from helper
            for m in range(len(self.nodes)):
               exec("%s = autogen_mod.%s"%("d%d"%m,"d%d"%m)) in locals ()

            
              
                
            for m in range(len(self.nodes)):

                exec("%s = %.6f" %("c%d"%m,alldata1[i-1][m])) in locals ()
            
            #print(alldata1[i-1][22])

            #if i==2000:
               # print("ppMLC",alldata1[i-1][self.nodes.index("ppMLC")])

            Tmin=0.5
            Tmax=500
            #print alldata1[i-1][self.nodes.index("ECM")]
            ECM=alldata1[i-1][self.nodes.index("ECM")]
            T=(alldata1[i-1][self.nodes.index("ppMLC")]*ECM)*Tmax+Tmin
            esra.append(T)

            if T>10:
                T=10

      
    

            wFX=eval(line4)
            #print eval(linew)
            #print eval(line4)
            #constants, factors and weight factors
            

            global wPolarity #calling weight function of polarity value for time 200
            wPolarity=[]
            if i==200:
                wPolarity=eval(line13)
            

            global wKcell #calling weight function of kcell value for time 200
            wKcell=[]
            if i==200:
                wKcell=eval(line14)
            


            

           

            kbeta_cons=1

            kbeta=kbeta_cons*wFX

            #if i==2000:print("kbeta",kbeta)
            #if i==2000:print("wFX",wFX)

            #print(kbeta)



            kgamma=(-exp(0.5*25) + exp(-25*(T/5-0.5)))/((1-exp(0.5*25))*(1+exp(-25*(T/5-0.5))))+0.15

            #if i==2000:
                #print("kgamma",kgamma)

            kdel_cons=100
            kepsi_cons=100
           

            kdelta=((3309*exp(-T)+(3.942*10**-7*exp(T)+5.819*10**-2*exp(-T))**-1)**-1)*kdel_cons
            #if i==2000:
                #print("kdelta",kdelta)

            
            kepsilon=((3309*exp(-T)+(3.942*10**-7*exp(T)+5.819*10**-2*exp(-T))**-1)**-1)*kepsi_cons



            k1=eval(line2) #dyMT
            
            k2=eval(line3) #stMT
           
            ###########
            #Focal Adhesion Loop
            def rxnb(F,timex):
                ITG=F[0]
                FX=F[1]
                FA=F[2]

                dITGdt=-kbeta*ITG\
                +kdelta*FX\
                -kgamma*FX*ITG\
                +kepsilon*FA
                
                        
                dFXdt=kbeta*ITG\
                -kdelta*FX\
                -kgamma*FX*ITG\
                +kepsilon*FA

                dFAdt=kgamma*FX*ITG\
                -kepsilon*FA
                                                            
                return [ dITGdt,dFXdt,dFAdt]

            #FA internal loop initials
            F0=[0,0,0]

        

            F0[0]=alldata1[i-1][ self.nodes.index("ITG")]
            F0[1]=alldata1[i-1][ self.nodes.index("FX")]
            F0[2]=alldata1[i-1][ self.nodes.index("FA")]
         
            F=odeint(rxnb,F0,timex)

            if (F[interloopsteps-1][1]+F[interloopsteps-1][0]+F[interloopsteps-1][2])==0:
                alldatay[self.nodes.index("ITG")]=0
                alldatay[ self.nodes.index("FX")]=0
                alldatay[ self.nodes.index("FA")]=0
            
            else:
            
                alldatay[self.nodes.index("ITG")]=totalFA*F[interloopsteps-1][0]/((F[interloopsteps-1][1]+F[interloopsteps-1][0]+F[interloopsteps-1][2]))
                alldatay[ self.nodes.index("FX")]=totalFA*F[interloopsteps-1][1]/((F[interloopsteps-1][1]+F[interloopsteps-1][0]+F[interloopsteps-1][2]))
                alldatay[ self.nodes.index("FA")]=totalFA*F[interloopsteps-1][2]/((F[interloopsteps-1][1]+F[interloopsteps-1][0]+F[interloopsteps-1][2]))

            ##########
            #factin to gactin oxidation loop
            m=10
            k5=eval(line5)*m #RhoGTPase cyt plus upstream
            k6=eval(line6)*m #RhoGTPase nuc
            k7=eval(line7)*m #cofilin cyt
            k8=eval(line8)*m #cofilin nuc
            k9=eval(line9)*m #mical2 cyt
            k10=eval(line10)*m #mical2 nuc
            k11=eval(line11)*m #MSRB_cyt
            k12=eval(line12)*m #MSRB_nuc

            

            #WPolarity=eval(line13) #weight function polarity
            #print(WPolarity)
   
            
            #Put in mass action equations
            def rxnb(G,timex):
                f_act_cyt=G[0]
                f_act_nuc=G[1]
                g_act_cyt=G[2]
                g_act_nuc=G[3]
                g_ox_cyt=G[4]
                g_ox_nuc=G[5]

               
                
                d_f_act_cyt_dt = k5*g_act_cyt-f_act_cyt*k9-f_act_cyt*k7

                d_f_act_nuc_dt = -f_act_nuc*k10-f_act_nuc*k8+k6*g_act_nuc
                
                d_g_act_cyt_dt = f_act_cyt*k7-g_act_cyt*k5+k11*g_ox_cyt

                d_g_act_nuc_dt = f_act_nuc*k8-g_act_nuc*k6+k12*g_ox_nuc

                d_g_ox_cyt_dt = f_act_cyt*k9-k11*g_ox_cyt

                d_g_ox_nuc_dt = f_act_nuc*k10-k12*g_ox_nuc

                '''print(k5,'k5')
                print(k6,'k6')
                print(k7,'k7')
                print(k8,'k8')
                print(k9,'k9')
                print(k10,'k10')
                print(k11,'k11')
                print(k12,'k12')
                print(g_act_cyt,'g_act_cyt')
                print(g_act_nuc,'g_act_nuc')
                print(f_act_cyt,'f_act_cyt')
                print(f_act_nuc,'f_act_nuc')
                print(g_ox_cyt,'g_ox_cyt')
                print(g_ox_nuc,'g_ox_nuc')'''


                return [d_f_act_cyt_dt,d_f_act_nuc_dt,d_g_act_cyt_dt,d_g_act_nuc_dt,d_g_ox_cyt_dt,d_g_ox_nuc_dt]
            
            
            G0= [alldata1[i-1][self.nodes.index("f_act_cyt")], alldata1[i-1][self.nodes.index("f_act_nuc")], alldata1[i-1][self.nodes.index("g_act_cyt")], alldata1[i-1][self.nodes.index("g_act_nuc")], alldata1[i-1][self.nodes.index("g_ox_cyt")], alldata1[i-1][self.nodes.index("g_ox_nuc")]]


            G=odeint(rxnb,G0,timex)

            
            alldatay[self.nodes.index("f_act_cyt")]=totalactin_cyt*G[steps-1][0]/((G[steps-1][0]+G[steps-1][2]+G[steps-1][4]))
            alldatay[ self.nodes.index("f_act_nuc")]= totalactin_nuc*G[steps-1][1]/((G[steps-1][1]+G[steps-1][3]+G[steps-1][5]))
            alldatay[ self.nodes.index("g_act_cyt")]=totalactin_cyt*G[steps-1][2]/((G[steps-1][0]+G[steps-1][2]+G[steps-1][4]))
            alldatay[ self.nodes.index("g_act_nuc")]= totalactin_nuc*G[steps-1][3]/((G[steps-1][1]+G[steps-1][3]+G[steps-1][5]))
            alldatay[ self.nodes.index("g_ox_cyt")]=totalactin_cyt*G[steps-1][4]/((G[steps-1][0]+G[steps-1][2]+G[steps-1][4]))
            alldatay[ self.nodes.index("g_ox_nuc")]= totalactin_nuc*G[steps-1][5]/((G[steps-1][1]+G[steps-1][3]+G[steps-1][5]))
            #explain alldatay restores all data values for every step

            #print(alldatay[self.nodes.index("g_ox_nuc")])            
            
        



            #Microtubules loop        
            def rxna(MT,timex):
                dyMT=MT[0]
                stMT=MT[1]
                               
                dstMTdt=k2*dyMT-k1*stMT
                ddyMTdt=-k2*dyMT+k1*stMT
                
                return [ddyMTdt,dstMTdt]

            D0=[0,0]
          
  
            D0[0]=alldata1[i-1][ self.nodes.index("dyMT")]
            D0[1]=alldata1[i-1][ self.nodes.index("stMT")]

            ##########
           

            #ODE solver
            MT=odeint(rxna,D0,timex)

            if (MT[interloopsteps-1][1]+MT[interloopsteps-1][0])==0:
                alldatay[self.nodes.index("dyMT")]=0
                alldatay[self.nodes.index("stMT")]=0
                
            else:
                alldatay[self.nodes.index("dyMT")]=totalMT*MT[interloopsteps-1][0]/((MT[interloopsteps-1][1]+MT[interloopsteps-1][0]))
                alldatay[self.nodes.index("stMT")]=totalMT*MT[interloopsteps-1][1]/((MT[interloopsteps-1][1]+MT[interloopsteps-1][0]))
            
               
            #Keeps the nodes in the internal loop
            loop_comp=[self.nodes.index("dyMT"),self.nodes.index("stMT"),self.nodes.index("FA"),self.nodes.index("FX"),self.nodes.index("ITG"),

            self.nodes.index("f_act_cyt"),self.nodes.index("f_act_nuc"),self.nodes.index("g_act_cyt"),self.nodes.index("g_act_nuc"),self.nodes.index("g_ox_cyt"),self.nodes.index("g_ox_nuc")]


            #Temporary list and extracting value from previous step,updating nodes other than loop_comp
            for index, node in enumerate( self.nodes ):
                if index not in loop_comp:                 
                    alldatay[index]=alldata1[i-1][index]


            alldata=odeint(autogen_mod.derivs,alldatay,np.linspace(0,1,10)) #until here you are going to update
         


            fin=[0]*len(self.nodes)                    

            #alldata[1][self.nodes.index("Ca")]=1
            #Calcium Oscillation test
            #alldata[1][self.nodes.index("actBAR2")]=(math.sin(i*math.pi/500)+1)/2
            
          
            #alldata[1][self.nodes.index("Ca")]=0.5
            alldata[1][self.nodes.index("ECM")]=ECMval
            #alldata[1][self.nodes.index("ppMLC")]=0
            #alldata[1][self.nodes.index("Ca")]=random.uniform(0.45,0.55)
            #alldata[1][self.nodes.index("actBAR2")]=random.uniform(0.45,0.55)
            
            #Adding Transcription Factors, 0.01 until step 100, then 0.1
            if i<=99:
                alldata[1][self.nodes.index("sActinin_cyt")]=0.01
                alldata[1][self.nodes.index("sSEMA")]=0.01
                alldata[1][self.nodes.index("sArp_cyt")]=0.01
                alldata[1][self.nodes.index("sERM")]=0.01
                alldata[1][self.nodes.index("TGLN")]=0.01
                alldata[1][self.nodes.index("sVinculin_cyt")]=0.01
                alldata[1][self.nodes.index("sWDR1")]=0.01
            if i>=100:
                alldata[1][self.nodes.index("sActinin_cyt")]=1
                alldata[1][self.nodes.index("sSEMA")]=1
                alldata[1][self.nodes.index("sArp_cyt")]=1
                alldata[1][self.nodes.index("sERM")]=1
                alldata[1][self.nodes.index("TGLN")]=1
                alldata[1][self.nodes.index("sVinculin_cyt")]=1
                alldata[1][self.nodes.index("sWDR1")]=1
            
            #small example TFs
            #sNodeVals = [0.0001,0.001,0.01,0.1,1]
            
            '''if i<=99:
                alldata[1][self.nodes.index("sG")]=sTF

            if i>=100:
                alldata[1][self.nodes.index("sG")]=TF'''
            
            

                
            
            
            #print(type(alldata[1][self.nodes.index("Actinin_cyt")]))

            '''conc_val = []

            for alldata in range(0,200):
                alldata[1][self.nodes.index("Actinin_cyt")]=conc_val
    
                if alldata in range(0,100):
                    conc_val = [0.01]
                if alldata in range(100,200):
                    conc_val = [0.1]'''

                
                
            
           
         
                    
            
         


            
            for index, node in enumerate( self.nodes ):
                fin[index]=alldata[1][index]

            #updating only MT and FA with previous loop values
            for z in loop_comp:
                fin[z]=alldatay[z]
           
            #all nodes adding by time   
            alldata1.append(fin)


    

            

          
        FAx=[]
        for row in alldata1:
            FAx.append(row[self.nodes.index("FA")])
        dFAxdt=diff(FAx)/diff(self.t)
            
        
        '''for m in range(len(self.t)-1):
            alldata1[m][self.nodes.index("FS")]=dFAxdt[m]'''
        for index, node in enumerate( self.nodes ):
            self.lazy_data[node] = [ row[index] for row in alldata1]
      


        
            
        
if __name__ == '__main__':
    text = """
    #
    # this is a comment
    #
    # conc, decay, treshold
    # 100%
    A = (1, 1, 0.5)
    B = (1, 1, 0.5)
    C = (1, 1, 0.5 )
    1: A* = not A 
    2: B* = A and B
    3: C* = C
    """
    
    model = PldeModel( text=text )
    model.initialize()
    model.iterate( fullt=5, steps=100 )

    from pylab import *
    plot( model.data['A'] ,'o-' )
    plot( model.data['B'] ,'o-' )
    plot( model.data['C'] ,'o-' )
    show()
