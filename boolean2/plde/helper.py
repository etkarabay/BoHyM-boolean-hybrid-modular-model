# -*- coding: utf-8 -*-
"""
Helper functions
"""

import csv, io
from itertools import *
import math
import pandas as pd
import numpy as np
# these function get injected into the generated code

helper_modules = """
try:
    from %s import *
except:
    pass
"""

from .defs import *

def change(node, indexer):
    "Returns the change for a node"
    index = indexer[node]
    return ' n%d ' % index 

def newval(node, indexer):
    """
    Returns the change for a node,
    
    chose a bad name originally still here for compatibility reasons, 
    will be deprecated
    """
    return change( node, indexer)
 
def conc( node, indexer):
    "Returns the concentration for a node"
    index = indexer[node]
    return ' c%d ' % index 

def decay( node, indexer):
    "Returns the decay for a node"
    index = indexer[node]
    return ' d%d ' % index 

def threshold( node, indexer):
    "Returns the threshold for a node"
    index = indexer[node]
    return ' t%d ' % index 

def default( node, indexer, tokens):
    "Default equation builder"
    newval = change( node, indexer )
    piece  = piecewise( tokens, indexer)
    return '%s = %s' % ( newval, piece )

def hill_func( node, indexer, par):
    """
    Generates a hill function call based on the parameters 
    """
    index = indexer[node]
    try:
        text = ' hill( c%d, h=%s, n=%s ) ' % ( index, par[node]['h'], par[node]['n'] )
    except Exception as exc:
        msg = "error creating hill function for node %s -> %s" % (node, exc)
        raise Exception(msg)
    return text

def prop_func( node, indexer, par):
    """
    Generates a proportion function call based on the parameters        
    """
    index = indexer[node]
    try:
        nconc = conc(node, indexer)
        text = ' prop( r=%s, rc=%s ) - %s ' % ( par[node].r, par[node].rc, nconc )
    except Exception as exc:
        msg = "error creating proportion function for node %s -> %s" % (node, exc)
        raise Exception(msg)
    return text


"""BoHyM Extension""" "kat start from here"

line2=[]
line3=[]
line4=[]
seq1=[]
seq2=[]
linew=[]
line5=[] #f-actin_cyt
line6=[] #F-actin_nuc
line7=[] #g-actin_cyt
line8=[] #g-actin_nuc
line9=[] #g-ox_cyt
line10=[] #g-ox_nuc
line11=[] #k11
line12=[] #k12
line13=[] #Polarity
line14=[] #Kcell



def piecewise(tokens,indexer):

    global seq1
    global seq2

    zag=["dyMT","stMT","ITG","FX","FA","f_act_cyt","f_act_nuc","g_act_cyt","g_act_nuc","g_ox_cyt","g_ox_nuc","K11","K12","Polarity","Kcell"]

    
    base_node  = tokens[1].value  
    base_index = indexer[base_node]

   
    suma=0  
    sumb=0
    line=[]

   
    
    nodes = [ t.value for t in tokens[4:] ]
    wnode="BAR2"
    if base_node==wnode:
        for node in nodes:
            if node in indexer:
                seq1.append(node)
                seq2.append(indexer[node])
   
        
    #if indexer[base_node]==285: print (base_node)

    
   
    a=1
    b=1
    num=0

    for node in nodes:
        if node in indexer:
            if nodes[nodes.index(node)-1]=='not':
                sumb=1+sumb
            if nodes[nodes.index(node)-1]!='not' :
                suma=1+suma


    #sensitivity analysis

    #num=0

    #from sens_initials_inv import a,b

    #print (a,b)
    

    '''
    for node in nodes:
        if node in indexer:
            if nodes[nodes.index(node)-1]=='not':
                sumb=b+sumb
            if nodes[nodes.index(node)-1]!='not' :
                suma=a+suma'''





    #activation part of weight function
    line1=[]
    active=0    
    for node in nodes:
            # replace each node with the comparison
            if node in indexer:
                 if nodes[nodes.index(node)-1]!='not' :

                     """Biased signaling test"""
                     
                     '''if node=="BAR2" :
                         if base_node=="Barrestin1": a=0.1
                         if base_node=="Barrestin2": a=0.1
                         if base_node=="Gs": a=0.1
                         if base_node=="Gbg": a=0.1
                         
                     else: a=1'''
                     active=1
                     num=1+num                  
                     if num <=1:
                         line.append ( 'float(sum([')               
                     index = indexer[node]             
                     value = " %f*c%d, " % (a,index)                                   
                     line.append ( value )
                     #suma=a+suma
                     #print(suma,node)
    
    
    if active==1:
        line.append('])/%f)'%(suma+sumb))


     
    #inhibition part of weight function
    
    num=0
    inhib=0
    for node in nodes:
        if nodes[nodes.index(node)-1]=='not':
            num=1+num
            inhib=1
            if num <=1:
                if len(line)==0 and nodes[0]=='not':
                    line.append ( 'float(1-sum([')
                else :
                    line.append ( '*float(1-sum([')                    
            index = indexer[node]   
            value = "  %f*c%d, " % (b,index)
            line.append ( value )     
    if inhib==1:
        line.append('])/%f)'%(suma+sumb))

    if base_node not in zag:
        if active==1:        
            line.append ("/")
            line.append ("(%f/((%f+%f)*d%d))" %(suma,suma,sumb,base_index))

    
    if base_node  in zag:
        if active==1:
            line.append ("/(%f/(%f+%f))" % (suma,suma,sumb))


    #weight function

    w="".join(line)
    #w="1"
    line=[]
    
    
    
    
    line1.append(w)
    
    line.append(w)
    

    line.append ( '*float(' )
    nodes = [ t.value for t in tokens[4:] ]
    for node in nodes:
        if node in indexer:        
            index = indexer[node]
            value = " ( c%d>t%d ) " % ( index,index)
        else:
            value = node
        line.append ( value )
    line.append ( ')' )
    
    #Removing the decay function from internal loop components' equations
    
    if base_node in zag:
        line.append("")        
    else:line.append ("- d%d *c%d" % ( base_index,base_index))

   


   
    A="".join(line1)
    
 
    global line2
    global line3
    global line4
    global linew
    global line5
    global line6
    global line7
    global line8 
    global line9 
    global line10 
    global line11
    global line12
    global line13
    global line14

    #internal loop weight functions

    if base_node=="dyMT": line2=A
    if base_node=="stMT": line3=A
    if base_node==wnode: linew=A
    if base_node=="FX": line4=A
    if base_node=="f_act_cyt": line5=A 
    if base_node=="f_act_nuc": line6=A
    if base_node=="g_act_cyt": line7=A
    if base_node=="g_act_nuc": line8=A
    if base_node=="g_ox_cyt": line9=A
    if base_node=="g_ox_nuc": line10=A
    if base_node=="K11": line11=A
    if base_node=="K12": line12=A
    if base_node=="Polarity": line13=A
    if base_node=="Kcell": line14=A
  


    return ' '.join(line)

def init_line( store ):

    #print (store)
    """
    Store is an incoming dictionary prefilled with parameters
    """
    patt = 'c%(index)d, d%(index)d, t%(index)d = %(conc).10f, %(decay)f, %(tresh)f #%(node)s' 
    return patt % store

def init_from_conc_max_threshold( param ):
    """
    Store is an incoming dictionary prefilled with parameters
    """
    patt = 'c%(index)d, d%(index)d, t%(index)d = %(conc)f, 1.0/%(decay)f, %(tresh)f # %(node)s' 
    return patt % param

def initializer(data,labels=None,**kwds):
    """
    Function factory that returns an initializer 
    that can initialize based on a parameter row (and naming convention)

    If a node is missing the function will raise an error. If a
    default parameter is passed to the function factory the
    function will return this value upon errors
    """

    # allow to override the parameter labels, the order is important
    labels = labels or 'conc decay threshold'.split()
    
    def func( node):
        # the extra work is to produce helpful error messages
        try:
            values = [ data[node][label] for label in labels ]
            return tuple(values) 
        except KeyError as exc:
            
            if 'default' in kwds:
                return kwds['default']
            else:
                if node not in data:
                    raise KeyError( 'could not find parameter %s' % node )
                for label in labels:
                    if label not in data[node]:
                        raise KeyError( "could not find parameter %s['%s']" % (node, label) )   

        
    return func 

class Parameter(object):
    """
    Allows attribute access to the parameters (Bunch)
    """
    def __init___(self, **kwds):
        self.__dict__.update( kwds )
    
    def __getattr__(self, attr):
        return self[attr]

    def __getitem__(self, key):
        if key not in self.__dict__:
            raise KeyError("parameter field '%s' not present" % key)
        else:
            return self.__dict__[key]
    
    def __setitem__(self, key, value):
        self.__dict__[key] = value
    
    def __repr__(self):
        return str(self.__dict__)
    
    def __contains__(self, key):
        return key in self.__dict__

    def setdefault(self, key, default):
        return self.__dict__.setdefault(key, default)

def read_parameters( fname ):
    """
    Reads parameters from a comma separated file and 
    returns a bunch object with attributes corresponding
    to the second line in the file
    """

    #
    # needs extra error checking because files created with 
    # Excel may contain artifacts => empty columns, empty lines, invisible values (spaces)
    #

    
    def something( row ):
        # skips rows with empty elements
        return filter(lambda x:x, map(string.strip, row ))
    
    # load the file, skipping commented or empty rows
    lines = filter( something, csv.reader( CommentedFile(fname))) 

    # check file size 
    assert len(lines) > 2, "file '%s' needs to have more than two lines" % fname
    
    # same number of columns in each line
    colnum = len( lines[0] )
    def coltest( elems ) :
        size = len(elems)
        if size != colnum:
            raise Exception( "column number mismatch expected %d, found %s, at line '%s'" % (size, colnum, ', '.join(elems)))
        return True
    
    lines = filter( coltest, lines )
    
    # nodes and attributes
    nodes, attrs = lines[0:2]
    
    # tries to coerce the value into a datastructure, float tuples, float or string
    def tuple_cast( word ):
        try:
            values = map( float, word.split(',') )
            if len( values ) > 1 :
                return tuple(values)
            else:
                return values[0]
        except ValueError:
            return word

    # generate the datastructure
    store  = []
    for elems in lines[2:]:
        param = Parameter()
        for index, attr, node in zip(count(), attrs, nodes):
            value = elems[index]
            param.setdefault( node, Parameter() )[attr] = tuple_cast( value )
        store.append( param )
    
    return store

class CommentedFile:
    """
    A file reader that skips comments in files
    """
    def __init__(self, fp):
        if isinstance(fp, str):
            fp = file(fp, 'rU')
        self.fp = fp

    def next(self):
        line = self.fp.next()
        while line.startswith('#'):
            line = self.fp.next()
        return line

    def __iter__(self):
        return self


