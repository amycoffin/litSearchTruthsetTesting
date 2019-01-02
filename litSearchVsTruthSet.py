# Amy Coffin
# Updated January 2, 2019

"""
Tests the output of any version of LitSearch/pubMunch (provided it is in the appropriate format) to produce a dictionary 
data structure, stored as an object "pubs" in file outPubs.py. 
Found variant aliases are stored together as a single tuple objec, which represents a single variant.  
The pubs dictionary takes on the following structure:
pubs = {'PMID': { ('variantalias1', 'variantalias2', ...) : {'inMunchOutput': boolean, 'inTruthSet': boolean} }, 
'PMID':{():{}}, ...  }
This structure can be used to research performance based on specific PMIDs, and is also used to provide performance 
statistics.  
stats is a function run on pubs just prior to the file being written to outPubs.py, which prints performance information 
about pubMunch compared to the truth sets provided. 

pubs is passed to all functions, obtaining more information with each function run. disabling functions or changing 
the order in which they are run may create false data. 

variable names such as pcol and vcol are used to refer to PMID column, and variant columns (where any aliases of all 
types are found) respectively.
pcol does NOT refer to protein identifier columns 
variant columns are searched using cPatternList and pPatternList.
cPatternList corresponds with regexes associated with cDNA nomenclature, while pPatternList uses regexes that 
correspond with unusual protein nomenclature
patternLists were compiled to interrogate specific regexes that were used to improve PubMunch, so they do NOT
reflect all regexes used in PubMunch/LitSearch

patternLists should be updated according to which regexes you are trying to ensure are working properly, per the truth set

To alter the columns interrogated in truthset or pubMunchOutput files, refer to load functions 
(loadLovdBRCA1, loadLovdBRCA2, loadMunch) 
and change arguments associated with truthSet (designed for LOVD truth sets) 
and munchOut (designed for a specific format of LitSearch/PubMunch output). 
"""

import sys
import re
import argparse 

pmidPat = re.compile(r".*PMID(\d{8}).*")  
nucPat0 = re.compile(r"(c\.\d+[ACTG]>[ACTG])")
nucPat1 = re.compile(r"(c\.\d+\+\d+[ACTG]>[ACTG]+)")
nucPat2 = re.compile(r"(c\.\d+\-\d+[ACTG]>[ACTG]+)")
nucPat3 = re.compile(r"(c\.\d+[_\-\+]\d+del[ACTG]+)")
nucPat4 = re.compile(r"(IVS\d+[\-\+]\d+[ACTG]>[ACTG]+)") 
nucPat5 = re.compile(r"(IVS\d+[\-\+]\d+del[ACTG]+)")
protPat = re.compile(r"([A-Z]\d+[A-Z])")
nucPat6 = re.compile(r"(c\.\d+[ACTG]+>None)")
nucPat7 = re.compile(r"(IVS\d+[\-\+]\d+[ACTG]+>None)")

cPatternList = [nucPat0, nucPat1, nucPat2, nucPat3, nucPat4, nucPat5, nucPat6, nucPat7]
pPatternList = [nucPat4, nucPat5, nucPat6, nucPat7,  protPat]
"""organize the patternLists according to the variant formats found in columns. In this case,
 cPatternList is used for variants of somewhat predictable c. formats, 
 while pPatternList searches for variants with IVS, deletion, protein, etc. formats"""


def parseArgs(): 
    parser = argparse.ArgumentParser(description='processes your files')
    parser.add_argument('ts1', help='your first truth set file')
    parser.add_argument('ts2', help='your second truth set file')
    parser.add_argument('pOut', help='your pubMunch output file') 
    return parser.parse_args()

def loadLovdBrca1(f, pubs):
    """opens then passes truth set file to truthSet function with integer arguments that indicate what column the pubmed ids (pcol) and variants (vcol, vcol2) will be found in"""
    truth1F = open(f)
    truthSet(truth1F, 8, 6, 5, pubs)
    
 
def loadLovdBrca2(f, pubs):
    """opens then passes a second truth set file to truthSet function with integer arguments that indicate what column the pubmed ids and variants will be found in"""
    truth2F = open(f) 
    truthSet(truth2F, 11, 6, 5, pubs) #args tell truthSet which columns to look in for PMID (arg1, and variant aliases (arg2, arg3)

def loadMunch(f, pubs): 
    """opens then passes pubMunch output file to munchOut function with integer arguments that indicate what column the pubmed ids and variants will be found in"""
    munchF = open(f)
    munchOut(munchF, 0, 10, 15, pubs) #args tell munchOut


def truthSet(f, pcol, vcol, vcol2, pubs):  
    """processes truth set data and conditionally enters it into the pubs data structure."""
    header = f.readline() 
    cols = f.readline()
    for line in f:
        line = line.strip().split('\t') #each line gets made into a list, and the arguments are used as indecies to retrieve pubs data in each line
        try:
           pmid = pmidPat.match(line[pcol]).group(1) #truthset used required a regex to obtain pmid 
        except: #handles lines without pmids/variant data 
           continue
        p = []  
        for pat in pPatternList: #iterates through the regexes to test variant against each one  
            try:
                p += pat.findall(line[vcol2]) #creates a list based on pattern matches 
            except:
                p = None #ensures that an empty list won't be added to pubs
        c = []
        for pat in cPatternList: #iterates through a pattern list to test for different variant formats in a different column
            try: 
                c += pat.findall(line[vcol])
            except: 
                c = None
        variant = [] #'final version' of the variant, with its differnt formats, to be added to pubs as one variant object
        if p is not None: 
            for pvar in p:  
                if pvar not in variant:#ensures no redundancies 
                    if pvar:  
                        variant.append(pvar) #they are appended to the master list
        if c is not None:  
            for cvar in c: 
                if cvar not in variant:
                    if cvar:
                        variant.append(cvar)
        if pmid in pubs:
            for var in variant: 
                varFound = False
                for vtup in pubs[pmid]: 
                    if var in vtup: #checks for current var formats in existing pubs   
                        pubs[pmid][vtup]["truth"] = True
                        varFound = True 
                if varFound == False:
                    pubs[pmid][tuple(variant)] = {"truth": True, "munch": False} #catches new variants
        else:
            pubs[pmid] = {}
            if variant: 
                pubs[pmid][tuple(variant)] = {"truth":True,"munch":False} #ensures that the first variant in the first iteration of a pmid is caught

def munchOut(f, pcol, vcol, vcol2, pubs): 
    """processes pubMunch output data and conditionally enters it into the pubs data structure. Works similarly to truthSet, but with different boolean assignments"""
    header = f.readline()
    for line in f: 
        line = line.strip().split('\t') 
        if len(line) <  2:
            continue
        try:
           pmid = line[pcol] 
        except:  
            continue
           
        c = []   
        for pat in cPatternList: 
            try:
                c += pat.findall(line[vcol]) 
            except: 
                c = None 
        p = []
        for pat in pPatternList:
            try:
                p += pat.findall(line[vcol2]) 
            except:
                p = None
        variant = []
        if p is not None: 
            for pvar in p:
                if pvar not in variant:
                    if pvar:
                        variant.append(pvar)
        if c is not None:  
            for cvar in c: 
                if cvar not in variant:
                    if cvar:
                        variant.append(cvar)
        if pmid in pubs:            
            for var in variant: 
                varFound = False
                for vtup in pubs[pmid]:
                    if var in vtup:              
                        pubs[pmid][vtup]["munch"] = True
                        varFound = True
                if varFound == False:
                    pubs[pmid][tuple(variant)] = {"truth": False, "munch": True} 
     
        else:
            pubs[pmid] = {}
            if variant:
                pubs[pmid][tuple(variant)] = {"truth":False,"munch":True}
             

def stats(p, mv):
    """This function will calculate performance of pubMunch output against the truthset. """   
    truthMunch = 0
    truthNotMunch = 0
    notTruthMunch = 0
    notTruthNotMunch = 0     
    neg = {}
    pos = {}
    totVar = 0
    #missedVars = []
    
    for pmid, rest in  pubs.iteritems(): 
        for variant in rest: 
            totVar += 1 
            if rest[variant]["munch"]: 
                if rest[variant]["truth"]:
                   truthMunch += 1
                else:
                   notTruthMunch += 1
            else:
                 if rest[variant]["truth"]:
                     truthNotMunch += 1
                     tempdict = {pmid: variant}
                     missedVars.append(tempdict)
                 else:
                     notTruthNotMunch += 1
    
    #print totVar  
    print "Truthset Variants Found by PubMunch:", truthMunch
    print "Truthset Variants Not Found by PubMunch:", truthNotMunch
    print "PubMunch Variants Not Present in TruthSet:", notTruthMunch
    print "Variants in Neither Set:", notTruthNotMunch
    #print missedVars

#def main(args):
pubs = {}
missedVars = []
arguments = parseArgs()
truth1 = arguments.ts1
truth2 = arguments.ts2
pmOutput = arguments.pOut
loadLovdBrca1(truth1, pubs)
loadLovdBrca2(truth2, pubs) 
loadMunch(pmOutput, pubs) 
stats(pubs, missedVars)
f = open("outPubs.py", 'w')
f.write('pubs = ' + str(pubs) + '\n') #pubs is now available through the module outPubs.py
g = open('missedVars.py', 'w')
g.write('missedVars = ' + str(missedVars) + '\n')
