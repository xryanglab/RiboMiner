#!/usr/bin/env python2.7
'''
NAME:        Seq2Logo Module
AUTHOR:      Martin Thomsen
DESCRIPTION:
   This module contains all relevant classes an functions needed to retrieve
   sequence data from an inputfile (see SeqParse for a complete list of
   supported formats), calculate a PSSM and create a sequence logo of the data.
REFERENCE:
   Martin C F Thomsen and Morten Nielsen,
   'Seq2Logo: a method for construction and visualization of amino acid binding
   motifs and sequence profiles including sequence weighting, pseudo counts and
   two-sided representation of amino acid enrichment and depletion.'
   Nucleic Acids Research (2012) 40 (W1): W281-W287
DEPENDENCIES:
   gs2.07 - Ghostscript version 2.07+ (www.ghostscript.com/download/gsdnld.html)
CONTENT:
   CLASES:
      Seqparse     - Sequence Parser
      SequenceData - Sequence Data Processing
      MakeLogo     - EPS Logo Generation
      Reg          - Extended Regular Expression Handler
      Error        - Systemic Error Handler
   FUNCTIONS:
       Hobohm      - Hobohm clustering algorithm 1 (Hobohm et al., 1992)
       Heuristic   - Heuristic clustering algorithm (S. HeniKoff and J. HeniKoff, 1994)
       EPS2Formats - Prepares the conversion from EPS to the requested formats (JPEG, PNG or PDF)
       ExecuteGS   - Employ Ghost Script to convert the EPS
'''
import sys, os, re, numpy as np, datetime as dt
from argparse import ArgumentParser, RawDescriptionHelpFormatter, SUPPRESS
from math import ceil, sqrt
import textwrap

################################################################################
##########                       Class Library                        ##########
################################################################################
class SeqParse:
   '''
   NAME:        SeqParse - Sequence Parser
   AUTHOR:      Martin Thomsen
   DESCRIPTION:
      This class is used to parse different kinds of sequence data, such as a
      list of raw peptides, fasta sequences or PSSM, to a sequence data class.
   SUPPORTED FORMATS:
      Raw peptides    (prealigned)
      Fasta sequences (prealigned)
      ClustalW
      Position Specific Scoring Matrices (PSSM):
         Frequency table
         BLAST PSSM
         Weight matrix
   DEPENDENCIES:
      class SequenceData
      class Reg
      class Error
   USAGE:
      >>> seqpar = SeqParse()
      >>> seqdat = SequenceData() # Object where the sequence data is stored
      >>> alphabet, bgfreq = seqpar.parseBackgroundFrequencies('bg_freqs.txt')
      >>> seqpar.SetAlphabet(seqdat)
      >>> inputType, debug, showcalc = seqpar.parsefile(inputfile)
   '''
   def __init__(self):
      # INITIALIZE VARS
      self.debug        = False
      self.showcalc     = False
   def SetAlphabet(self, seqdat, alphabet='FAMILVWQSTYNCGPHKRDE', gaps="_-.*"):
      # INITIALIZE VARS
      self.gaps         = gaps
      self.aaRegEx      = alphabet
      self.gapRegEx     = re.escape(gaps)
      self.unknownRegEx = 'X' if not 'X' in alphabet else ''
      self.allowed      = np.array(list(alphabet + gaps + self.unknownRegEx)) # A list of all allowed symbols
      # Set Regular Expressions
      self.seqRegEx     = Reg('^[%s]+'%(self.aaRegEx+self.gapRegEx+self.unknownRegEx), re.I) #re.compile('^[%s]+'%(aaRegEx+gapRegEx+unknownRegEx), re.I) 
      self.reASCII      = Reg('[^a-z0-9\_\-]', re.I)
      # SPECIAL case STAR and DOT can't be used as gaps in CLUSTAL, only dash is gaps
      self.reCLUSTAL    = Reg('(\S+)\s+([%s]+)(\s+[0-9]+)?\s*$'%(self.aaRegEx+re.escape('-')+self.unknownRegEx), re.I) #re.compile('\S+\s+[%s]+(\s+[0-9]+)?\s*'%(aaRegEx+gapRegEx+unknownRegEx), re.I)
      self.reWHEAD      = Reg('^\s*(\S\s+){2,50}(\S)\s*$') #Reg('^\s*([\w%s]\s+){2,50}([\w%s])\s*$'%(self.gapRegEx, self.gapRegEx)) #re.compile('^\s*([\w%s]\s+){2,50}\S'%(gapRegEx), re.I)
      # Set Sequence Data Object
      self.seqdat = seqdat
   def parsefile(self, inputfile):
      ''' This function determines the type of the input file and parses the sequence data to sequnce data class. '''
      # DETERMINE INPUT TYPE
      inputType = ""
      with open(inputfile, 'r') if inputfile != '' else sys.stdin as fobj:
         f = readstream(fobj)
         for lc, l in enumerate(f,1): #lc = line count
            l=l.strip()
            if not l: continue
            if l[0]!='#':
               if l[0] == ">": # FASTA NAMESTRING
                  inputType = "FASTA"
                  f.seekline(0) # reset file cursor
                  self.__parseFASTA(f) # PARSE SEQUENCE DATA
                  break
               elif l[:22].upper() == "LAST POSITION-SPECIFIC" or self.reWHEAD.match(l): # WEIGHT MATRIX HEADERS
                  inputType = "WEIGHT"
                  if l[:22].upper() != "LAST POSITION-SPECIFIC": f.seekline(0) # reset file cursor
                  self.__parseWEIGHT(f) # PARSE SEQUENCE DATA
                  break
               elif l[:7].upper() == "CLUSTAL" or self.reCLUSTAL.match(l): # CLUSTAL HEADER OR DATA LINE
                  inputType = "CLUSTAL"
                  if l[:7].upper() != "CLUSTAL": f.seekline(0) # reset file cursor
                  self.__parseCLUSTAL(f) # PARSE SEQUENCE DATA
                  break
               elif self.seqRegEx.match(l) and len(l.split()[0]) > 1: # RAW PEPTIDE SEQUENCES
                  inputType = "RAW"
                  f.seekline(0) # reset file cursor
                  self.__parseRAW(f) # PARSE SEQUENCE DATA
                  break
               else: ParsingError('unknown format', (lc, l.strip())) # Could not determin filetype
            elif l[:6].upper() == "#DEBUG":    self.debug    = True # DEBUG PARAMETER SPECIFIED IN INPUT DATA
            elif l[:9].upper() == "#SHOWCALC": self.showcalc = True # SHOWCALC PARAMETER SPECIFIED IN INPUT DATA
      return inputType, self.debug, self.showcalc
   def parseBlosum(self, fp):
      ''' Retrieve and return Blosum Matrix, used as prior probability for
      pseudo count calculations. Needed to adjust data sets with few datapoint.
      '''
      alphabet = []
      mat = []
      reVal = Reg('^\s*(\s*\-?[0-9]+(\.?\,?[0-9]+)?(e-?[0-9]+)?)+', re.I)
      with open(fp, 'r') as f:
         for l in f:
            if not l.strip(): continue
            if l[0]=='#': continue
            tmp = l.strip().split()
            if reVal.match(l): mat.append(tmp) # values found
            else:
               alphabet = tmp # Alphabet found
      return alphabet, mat
   def parseBackgroundFrequencies(self, fp):
      ''' Retrieve and return Background Frequencies. Used to calculate the
      Kullback-Leibler distance, which compares the probability of the element
      occuring in the current environment to the probability of it occurs by
      random for a non-uniform distribution of elements. '''
      gapRegEx = re.escape("_-.*")
      reWHEAD      = Reg('^\s*(\S\s+){2,50}(\S)\s*$') #Reg('^\s*([\w%s]\s+){2,50}([\w%s])\s*$'%(gapRegEx, gapRegEx))
      foundalphabet = False
      alen = 0
      alphabet = []
      weights = []
      with open(fp, 'r') as f:
         for lc, l in enumerate(f,1): #lc = line count
            l=l.strip()
            if not l: continue
            if l[0]=='#': continue
            if foundalphabet:
               # Find weights
               if reWeights.match(l):
                  g = reWeights.getgroup(0).split()
                  # Remove position number and concensus sequence if found
                  for i in xrange(len(g)-alen): g.pop(0)
                  # Replacing , with . if any is found
                  if ',' in l: g=[x.replace(',', '.') for x in g]
                  # Adding weights to the weight matrix
                  weights.append(g)
               else: ParsingError('ref_weight_vary', lc) #ERROR
            else:
               # Find alphabet
               if reWHEAD.match(l):
                  alphabet = [a.strip() for a in reWHEAD.getgroup(0).split()]
                  alen = len(alphabet)
                  reWeights = Reg('^\s*(\-?[0-9]+\s+)?([%s%s]+\s+)?(\s*\-?[0-9]+(\.?\,?[0-9]+)?(e-?[0-9]+)?){%d}'%(''.join(alphabet), gapRegEx, alen), re.I)
                  foundalphabet = True
      return alphabet, weights
   def __parseCLUSTAL(self, f):
      ''' This function parses a clustalw alignment file to the alignment. '''
      alignment = {}
      firstname = ''
      for lc, l in enumerate(f,1): #lc = line count
         l=l.strip()
         if not l: continue
         if l[0]=='#': continue
         if self.reCLUSTAL.match(l):
            k = self.reASCII.sub('', self.reCLUSTAL.getgroup(1)) # name
            v = self.reCLUSTAL.getgroup(2) # sequence
            if k in alignment:
               if firstname == k:
                  if curseq != numseq: ParsingError('missing sequence', lc) #ERROR incurrect number of sequences in line lc
                  else: curseq = 0
               alignment[k] += v
               curseq +=1
            elif firstname == '':
               alignment[k] = v
               curseq =1
               numseq =1
               firstname = k
            else:
               alignment[k] = v
               numseq +=1
               curseq +=1
         elif curseq == numseq: continue # Annotation line
         else: ParsingError('clustal_err') # Raising Error due to unrecognised line
      # Adding the found sequences to the alignment
      for s in alignment.values():
         e = self.seqdat.add_sequence(s)
         if e==1: ParsingError('seq_empty') #ERROR
         elif e==2: ParsingError('seq_vary') #ERROR
   def __parseWEIGHT(self, f):
      ''' This function parses a PSSM input file. '''
      foundalphabet = False
      alen = 0
      blast = False
      for lc, l in enumerate(f,1): #lc = line count
         l=l.strip()
         if not l: continue
         if l[0]=='#': continue
         if foundalphabet:
            # Find weights
            if reWeights.match(l):
               g = reWeights.getgroup(0).split()
               diff = len(g)-alen
               if diff > 0: self.seqdat.append_position(int(g.pop(0))) # Save original position number #TODO fix originalPosition
               if diff == 2:
                  # Store consensus sequence
                  self.seqdat.append_consensus(g.pop(0))
               if ',' in l: g=[x.replace(',', '.') for x in g] # replacing , with . if any is found
               self.seqdat.append_weights(g) # Adding weights to the sequence data class
            elif blast and 'lambda' in l.lower(): break # Blast end details found...
            else: ParsingError('weight_vary', lc) #ERROR
         else:
            # Find alphabet
            if self.reWHEAD.match(l):
               # CHECK FOR DUPLICATES -> IF FOUND IT IS BLAST
               alphabet = "".join([a.strip() for a in self.reWHEAD.getgroup(0).split()])
               alen = len(alphabet)
               cutoff = -1
               for i in xrange(alen):
                  a = alphabet[i]
                  for j in xrange(i+1,alen):
                     if alphabet[j] == a:
                        cutoff=j
                        break
                  if cutoff > 0: break
               if cutoff > 0 and cutoff < alen:
                  # Alphabet contained doublicates -> cutting away the excess
                  alphabet = alphabet[:cutoff]
                  # blast matrix
                  blast = True
               self.seqdat.set_alphabet(alphabet)
               foundalphabet = True
               alen = self.seqdat['alen']
               if blast:
                  reWeights = Reg('^\s*(\-?[0-9]+\s+[a-zA-Z%s]+\s+)?(\s*\-?[0-9]+(\.?\,?[0-9]+)?(e-?[0-9]+)?){%d}'%(self.gapRegEx,alen), re.I)
               else:
                  reWeights = Reg('^\s*(\-?[0-9]+\s+)?([%s%s]+\s+)?(\s*\-?[0-9]+(\.?\,?[0-9]+)?(e-?[0-9]+)?){%d}'%(alphabet,self.gapRegEx,alen), re.I)
   def __parseFASTA(self, f):
      ''' This function parses a fasta sequence file input to the alignment. '''
      s='' # string containing the sequence
      for lc, l in enumerate(f,1): #lc = line count
         l=l.strip()
         if not l: continue
         if l[0]=='#': continue
         if l[0] == ">": # Fasta namestring
            if s != "":
               e = self.seqdat.add_sequence(s) # Appending alignment
               if e==1: ParsingError('seq_empty', line=(lc, l)) #ERROR
               elif e==2: ParsingError('seq_vary', line=(lc, l)) #ERROR
               s = "" # resetting sequence string
         elif self.seqRegEx.match(l):
            l=l.split()[0].strip().upper() #Removing anything but the sequence in the first column
            s += l # appending string
         else: ParsingError('unrecognised line', lc, self.__ShowError(l, l.split()[0].strip())) # Raising Error due to unrecognised symbol
      if s != "":
         e = self.seqdat.add_sequence(s) # Adding the last sequence to the alignment
         if e==1: ParsingError('seq_empty', line=(lc, l)) #ERROR
         elif e==2: ParsingError('seq_vary', line=(lc, l)) #ERROR
   def __parseRAW(self, f):
      ''' This function parses a raw sequence file input to the alignment. '''
      for lc, l in enumerate(f,1): #lc = line count
         l=l.strip()
         if not l: continue
         if l[0]=='#': continue
         if self.seqRegEx.match(l):
            l=l.split()[0].strip().upper() #Removing anything but the sequence in the first column
            e = self.seqdat.add_sequence(l)
            if e==1: ParsingError('seq_empty', line=(lc, l)) #ERROR
            elif e==2: ParsingError('seq_vary', line=(lc, l)) #ERROR
         else: ParsingError('unrecognised line', lc, self.__ShowError(l, l.split()[0].strip())) #Raising Error due to unrecognised symbol
   def __ShowError(self, line, seq):
      c = (line.index(seq) + 
          np.argmin(np.any(self.allowed==np.array(list(seq))[np.newaxis].T, 1)))
      return "%s\n%s^\n"%(line, '-'*(c))


class SequenceData:
   '''
   NAME:        SequenceData - Sequence Data Manager
   AUTHOR:      Martin Thomsen
   DESCRIPTION:
      This class defines a structure for storing and calculating all relevant
      sequence data needed to create a sequence logo.
   DEPENDENCIES:
      class Error
   USAGE:
      >>> seqdat = SequenceData()      # Object where the sequence data is stored
      >>> seqdat.set_alphabet(alphabet, blosum, gaps, unknown)
      >>> seqdat.add_sequence(sequence)         # Adding a sequence
      >>> seqdat.set_sw(weights)                # Set Sequence Weighting
      >>> seqdat.calc_freq()                    # Calculate Observed Frequencies
      >>> seqdat.calc_prob(alpha, beta)         # Calculate Probabilities
      >>> seqdat.calc_weights(yLabel, logoType) # Calculate PSSM
      >>> seqdat.set_heights(logoType)          # Calculate Symbol Heights
      >>> seqdat.set_positions(startPos)        # Add Position Numbers
      >>> seqdat.set_width(inputType, minwidth) # Calculate and set stack width
   '''
   def __init__(self):
      self.alphabet = []      # Symbols contained in the alignment
      self.alen = 0           # Alphabet length
      self.a2i = {}           # Dictionary converting symbol to index
      self.bg = {}            # Background frequencies {aa: bgfreq}
      self.bl = None          # Blosum matrix (nested dictionary)
      self.gaps = None        # contains a list of gaps, used to determin the width of a stack
      self.aln = None         # alignment / sequences
      self.seqlen = -1        # Length of sequences
      self.numseq = 0         # Number of sequences
      self.f = []             # frequencies    (N*M  N=alen (column), M=seqlen (row))
      self.p = None           # probabilities  (N*M  N=alen (column), M=seqlen (row))
      self.q = None           # background     (N*M  N=alen (column), M=seqlen (row))
      self.w = None           # weights        (N*M  N=alen (column), M=seqlen (row))
      self.h = None           # heights        (N*M  N=alen (column), M=seqlen (row))
      self.sw = None          # array of sequence weights
      self.pos = np.array([], dtype='i4') # array of position numbers
      self.consensus = None
      self.width = None        
      self.maxHeight = 0
      self.minHeight = 0
   def set_alphabet(self, alphabet, blosum=None, gaps='', unknown=''):
      ''' This sub-routine appends the provided arrays to the alphabet array, bg array and the a2i dictionary. '''
      self.alphabet = list(alphabet)
      self.alen = len(self.alphabet)
      if blosum: self.bl = np.array(blosum, dtype='f8')
      # Adding all characters in the alphabet to the dictionary a2i pointing to their respective indexes in the alpahbet
      self.a2i =  dict([(k,v) for v,k in list(enumerate(self.alphabet))])
      # Adding all gap characters to the dictionary a2i pointing to 'NaN'
      if gaps+unknown: map(lambda g,d=self.a2i.__setitem__: d(g, 'NaN'), list(gaps+unknown))
   def set_sw(self, weights):
      ''' This sub-routine appends the provided arrays to the alphabet array, bg array and the a2i dictionary. '''
      try: self.sw = np.array([weights], dtype='f8').T
      except Exception, e: DataError('set_sw', e=e.message) #ERROR
   def add_sequence(self, seq):
      ''' This sub-routine appends the provided string to the alignment matrix. '''
      if len(seq) == self.seqlen:
         self.__add_sequence(seq)
      elif len(seq) == 0: return 1 # Return error
      elif self.seqlen == -1:
         self.seqlen = len(seq)
         self.__add_sequence(seq)
      else: return 2
      return 0
   def __add_sequence(self, seq):
      ''' This is a private function used to limit the amount of repeated script. '''
      # Converting the symbols in seq to their corresponding indexes in the alphabet
      try: self.aln = np.concatenate((self.aln, np.array([[self.a2i[a] for a in list(seq)]],dtype='f8')), axis=0)
      except Exception, e:
         if self.aln is None:
            try: self.aln = np.array([[self.a2i[a] for a in list(seq)]],dtype='f8') #Init aln
            except KeyError, e: ParsingError('seq_unknown', e=e.message) #ERROR
         else: ParsingError('seq_unknown', e=e.message) #ERROR
   def calc_freq(self):
      ''' calculate frequencies: Count occurences of each symbol in the alphabet
          per position. Works by making a reference array (a) and then counting
          the occurences of all references in the alignment per position (pos).
          This returns an array with symbol counts according to the alphabet.
          Sequence weighting adjust the count. After the counting, the
          frequencies are calculated by normalising according to the row sum.'''
      a = np.array(range(self.alen))
      pos = np.hsplit(self.aln,self.seqlen)
      # Counting (weighted) observations
      for i in xrange(self.seqlen):
         # a==pos[i] ==> 2d count array, rows=sequences, cols=alphabet
         # *sw applies the sequence weights for the counts
         # np.sum adds up all the counts column-wise
         try: self.f.append(np.sum((a==pos[i])*self.sw, 0))
         except Exception, e: DataError('freq_err', e=e.message) #ERROR
      # Frequencies are calculated
      self.f /= np.sum(self.f,1)[np.newaxis].T # [np.newaxis].T ==> convert 1d array to transposed 2d array
   def append_frequencies(self, freqs): #USED?
      ''' This sub-routine appends the provided array to the frequency matrix. '''
      if len(freqs) == self.alen: self.f = np.append(self.f, np.array(freqs, dtype='f8'), 0)
   def calc_prob(self, a=0, b=0):
      ''' Calculate probabilities, using pseudo count, a is alpha which is the
          weight on the observed frequencies. b is beta which is the weight on
          prior knowledge acquired from the blosum substitution matrix. '''
      if not isinstance(a, np.ndarray):
         if isinstance(a, (float, int)):
            a = np.ones(self.seqlen) * a
         else:
            print("Unhandled alpha format: %s"%type(a))
      if (a >= 0).all and b > 0 and self.bl is not None: # Calculate pseudo counts
         if isinstance(a, np.ndarray) and len(a.shape) == 1:
            # transpose a
            a = a[:,np.newaxis]
         f =  self.f  # observed frequencies
         bl = self.bl # blosum substitution matrix
         g = np.dot(f,bl) # Pseudo counts
         self.p = np.array((np.multiply(a,f) + b*g) / (a + b))
      else: self.p = self.f # transfer frequencies to probabilities
   def calc_weights(self, unitType='Bits', logoType=None, bgfreqs=None):
      ''' Calculate weights. w=log2(p/q) '''
      # Handle input arguments
      unitType = unitType.lower() # lowering capital letters
      if not bgfreqs is None:
         # Verify correct usage of reference distribution and Adjust
         # distribution accordingly if not
         bgfreqs = np.array(bgfreqs, dtype='f8')
         # Check shape has more than 1 dimension
         bshape = bgfreqs.shape
         pshape = self.p.shape
         if bshape[0] > 1:
            # Check if it fits p.shape
            if bshape[0] != pshape[0]:
               # Fix bad shape
               print(bshape, pshape)
               if bshape[1] != pshape[1]:
                  # Set uniform reference distribution
                  q = np.ones(pshape, dtype='f8')/pshape[1]
                  # Warn User
                  RaiseWarning('ref_dist_uniform', '') #WARN
               elif bshape[0] > pshape[0]:
                  # Cut off the excess
                  q = bgfreqs[:pshape[0],:]
                  # Warn User
                  RaiseWarning('ref_dist_too_many_pos', '') #WARN
               else:
                  # Extend with uniform distribution
                  q = np.vstack([bgfreqs, np.ones((pshape[0]-bshape[0],pshape[1]))/pshape[1]])
                  # Warn User
                  RaiseWarning('ref_dist_too_few_pos', '') #WARN
            else:
               q = bgfreqs
         else:
            # Global reference distribution
            q = np.ones(pshape)*bgfreqs[0]
         # Check for negative or zero sum positions
         if np.any(q==0):
            # Warn User
            RaiseWarning('ref_dist_zero', '0 found at [row, col] = %s'%(', '.join([str([y+1 for y in x]) for x in zip(*np.nonzero(q==0))]))) #WARN
            # Fix sum=zero rows
            for x in np.nonzero(q==0)[0]:
               q[x,:] = 1./pshape[1]
         bgvalues = q[0]
         bgsum = np.sum(bgvalues)
         bgzeros = np.sum(bgvalues<=0)
         # Check background frequencies are valid
         if (bgsum < 0.95 or bgsum > 1.05) or bgzeros > 0 or len(bgvalues) != self.alen:
            DataError('bgfreq_invalid') #ERROR
      else:
         # Set uniform reference distribution
         q = np.ones(pshape, dtype='f8')/pshape[1]
      # INITIALIZE w
      self.w = np.zeros(pshape)
      pnot0 = self.p>0
      # HANDLE WEIGHTS FOR p==0 CASES
      self.w[self.p==0]=-99.999
      # CALCULATE WEIGHTS FOR THOSE p>0
      if unitType == 'bits':
         self.w[pnot0]=np.log2(self.p[pnot0]/q[pnot0])
         if np.any(self.w<-50): # upholding cut-off on -50
            self.w[self.w<-50]= -50
            RaiseWarning('weight_neg_inf', '-50 bits') #WARN
      elif unitType == 'halfbits':
         self.w[pnot0]=np.log2(self.p[pnot0]/q[pnot0])*2
         if np.any(self.w<-99.999): # upholding cut-off on -99.999
            self.w[self.w<-99.999]= -99.999
            RaiseWarning('weight_neg_inf', '-99.999 halfbits') #WARN
      else: # Calculate as Bits
         self.w[pnot0]=np.log2(self.p[pnot0]/q[pnot0])
         if np.any(self.w<-50): # upholding cut-off on -50
            self.w[self.w<-50]= -50
            RaiseWarning('weight_neg_inf', '-50 bits') #WARN
      if self.consensus is not None and self.consensus == []: # Find Consensus
         self.consensus = [self.alphabet[i] for i in np.argmax(self.p, axis=1)]
   def append_weights(self, weights):
      ''' This sub-routine appends the provided array to the weight matrix. '''
      if len(weights) == self.alen:
         try: self.w = np.concatenate((self.w, np.array([weights], dtype='f8')), axis=0)
         except Exception, e:
            if self.w is None: self.w = np.array([weights], dtype='f8')
            else: DataError('weight_err', msg=(e, self.w, weights)) #ERROR
      else: DataError('weight_vary') #ERROR
   def set_heights(self, heighttype):
      ''' This sub-routine calculates the height of the symbols in the alphabet
          from the sequence data depending on the height type and stores them
          in the heights matrix. I=p*w, bh = sum(I) '''
      if heighttype == 5: # PSSM-LOGO
         # CHECK WEIGHTS EXISTS
         if self.w.std() == 0: DataError('weight_zero') #ERROR
         # SET HEIGHT
         self.h = self.w
      else:
         # CHECK WEIGHTS AND PROBABILITIES EXISTS
         if self.w.std() == 0: DataError('weight_zero') #ERROR
         if self.p.std() == 0: DataError('prob_zero')   #ERROR
         # DEFINING HEIGHT FUNCTIONS
         fct_sha  = lambda p, w: np.sum(w*p, 1)[np.newaxis].T * p
         fct_kl   = lambda p, w: np.sum(w*p, 1)[np.newaxis].T * p * np.sign(w)
         fct_wkl  = lambda p, w: np.sum(w*p, 1)[np.newaxis].T * w / np.sum(np.absolute(w), 1)[np.newaxis].T
         fct_pwkl = lambda p, w: np.sum(w*p, 1)[np.newaxis].T * p*w / np.sum(np.absolute(w)*p, 1)[np.newaxis].T
         h={ 1 : fct_sha,   # bh * p,           # SHANON
             2 : fct_kl,    # bh * p * sign(w), # KULLBACK-LEIBLER
             3 : fct_wkl,   # bh * w / ws,      # WEIGHTED KULLBACK-LEIBLER
             4 : fct_pwkl } # bh * p * w / pws  # P-WEIGHTED KULLBACK-LEIBLER
         # CALCULATE HEIGHTS
         self.h = h[heighttype](self.p, self.w)
         # CHECK HEIGHTS WERE CALCULATED
         if self.h.std() == 0: DataError('height_zero')   #ERROR
   def set_positions(self, startPos=1, exclude_zero=True):
      ''' This method does several things:
          * Sets seqlen if unset.
          * Checks if positions have been assigned properly.
            - Assigning linear position numbers with a user set offset.
            - Circumvents 0, as this is nonsense in respect to positions
      '''
      if startPos > 0: startPos -=1
      if self.seqlen == -1: self.seqlen = len(self.h)
      if len(self.pos) < self.seqlen:
         self.pos = np.arange(self.seqlen) + startPos
         if exclude_zero:
            self.pos[self.pos>=0]+=1
   def append_position(self, newPos):
      ''' Appending a position number to the position array.
      (needed to parse positions from weight matrices.) '''
      self.pos = np.append(self.pos, newPos)
   def set_width(self, inputType, minwidth):
      ''' This sub-routine calculates the width of the stacks in the sequence
      logo. Counts sum of collective gaps per position and deduct the gap fraction
      compared to the sequence length. '''
      if inputType == "WEIGHT": # PSSM
         self.width = np.ones(self.seqlen)
      else: # MSA or RAW
         if self.sw is not None:
            self.width = np.around(np.sum((self.aln>=0)*self.sw, 0) / np.sum(self.sw, dtype='|f8'), decimals=2)
         else:
            self.width = np.around(np.sum(self.aln>=0,0) / float(self.numseq), decimals=2)
         self.width[self.width < minwidth] = minwidth
   def set_consensus(self):
      ''' Enable consensus retrieval and presentation. '''
      self.consensus = []
   def append_consensus(self, cs):
      ''' Add symbols to the consensus sequence. '''
      if self.consensus is not None: self.consensus.append(cs)
   def __getitem__(self, key):
      ''' This function returns the internal variable 'key' '''
      return getattr(self,key)
   def __setitem__(self, key, value):
      ''' This function sets the internal variable 'key' to 'value'. '''
      return setattr(self,key, value)


class MakeLogo:
   '''
   NAME:        MakeLogo - Logo Generation Manager
   AUTHOR:      Martin Thomsen
   DESCRIPTION:
      This class provides the methods to create a sequence logo in Encapsulated
      PostScript (EPS) format from sequence data inputs.
   DEPENDENCIES:
      class Error
      EPS template file (See he bin directory)
      math (ceil, sqrt)
      datetime as dt
   USAGE:
      >>> logo = MakeLogo(outputfile)               # Create a new logo file
      >>> logo.SetParameters(title, pageSize, spl, lpp, userYaxis, colors, xTic,
                             heights, showX, showY, rotNum, xLabel, yLabel,
                             eLabel)                # Set logo environment
      >>> logo.PrepareHeader(eps_template)          # Copy and edit the template
      >>> logo.AddData(widths, positions, alphabet) # Add the logo data
      >>> logo.Finish()                             # End and close the file
   '''
   def __init__(self, filename, pos, title, pageSize, spl, lpp, userYaxis, colors, xTic, heights, showFine, showX, showY, rotNum, xLabel, yLabel, eLabel, eps_template, widths, alphabet, consensus):
      # Initialising the Class Variables
      self.fn = filename+".eps"
      self.pH = 0         # Page Height of the EPS file
      self.pW = 0         # Page Width of the EPS file 
      self.title = ''     # Shown title on the EPS pages
      self.fine  = ''     # Shown fineprint
      self.pages = 0      # Number of pages in the logo
      self.lpp = 0        # Number of lines per page
      self.spl = 0        # Number of stacks per line
      self.minHeight = 0  # Y-axis lowest value
      self.maxHeight = 0  # Y-axis highest value
      self.yTicLarge = 0  # Y-axis Large Tic Interval
      self.yTicSmall = 0  # Y-axis Small Tic Interval
      self.yLabel = ''    # Y-axis Unit Label
      self.xLabel = ''    # X-axis Unit Label
      self.eLabel = ''    # Sequence End Labels (eg. 5') - d: DNA, p: PROTEIN (make automatic if needed!)
      self.xTic = 0       # X-axis Numbered Tic interval
      self.rotNum = None  # Rotate Numbers
      self.showY = True   # Show Y-axis on Logo
      self.showX = True   # Show Y-axis on Logo
      self.cDict = None   # Color Dictionary
      self.date = dt.date.today() # EPS Creation Date
      self.heights = None # Heights from Sequence Data
      # Opening the EPS File for Writing
      self.OUT = open(self.fn, 'w')
      # MAKE LOGO
      self.SetParameters(pos, title, pageSize, spl, lpp, userYaxis, colors, xTic, heights, showFine, showX, showY, rotNum, xLabel, yLabel, eLabel)
      self.PrepareHeader(eps_template)
      self.AddData(widths, pos, alphabet, consensus)
      self.Finish()
   def SetParameters(self, pos, title='', pageSize='640x480', spl=40, lpp=4, userYaxis=[0.0, 0.0], colors=None, xTic=0, heights=None, showFine=True, showX=True, showY=True, rotNum=False, xLabel='', yLabel='Bits', eLabel=''):
      ''' This class method does several things:
          *  Sets Several User-defined Variables.
          *  Sets Page Size.
          *  Correction of Stacks Per Line and Lines Per Page.
          *  Sets Best Matching Numbered Tic Interval for the X-axis
          *  Determines and Sets Boundries on the Y-axis.
          *  Calculating Best Matching Y-axis Tic Intervals.
          *  Sets the Color Dictionary Used in the Logo.
      '''
      # Set User Defined Variables
      self.title  = title
      self.lpp    = lpp
      self.spl    = spl
      self.showX  = str(showX).lower()
      self.showY  = str(showY).lower()
      self.rotNum = str(rotNum).lower()
      self.yLabel = yLabel
      self.xLabel = xLabel
      self.eLabel = eLabel
      if showFine: self.fine = 'Created by Seq2Logo'
      # Set Sequence Defined Variables
      self.heights = heights
      self.seqlen  = len(heights)
      # Call internal methods
      self.__SetPageBoundries(pageSize)
      self.__SetYaxis(userYaxis)
      self.__SetXaxis(pos, xTic)
      self.__SetColor(colors)
   def PrepareHeader(self, eps_template):
      ''' Retrieve, modify and write EPS template to output file. '''
      tmpPH = Reg('\{\%\=(\w+)\%\}', re.I) # RE Object which matches place-holders in the EPS template
      with open(eps_template, 'r') as f:
         for l in f:
            while tmpPH.match(l): l = tmpPH.sub(str(getattr(self,tmpPH.getgroup(1))), l, 1)
            self.OUT.write(l)
   def AddData(self, widths, positions, symbols, consensus):
      ''' Add Logo data to the output file. '''
      ls = len(positions)-1 # last stack
      done = False
      if consensus == []: consensus = None
      # Find valid positions - where there is found height information (not nan's)
      heights = self.heights
      if np.any(np.isnan(heights)):
         # Warn User
         RaiseWarning('no_data_pos', 'nan found at [row, col] = %s'%(', '.join([str([y+1 for y in x]) for x in zip(*np.nonzero(heights==np.nan))]))) #WARN
         # Remove rows without data
         nans = np.any(np.isnan(heights),1)
         valid = np.logical_not(nans)
         heights = heights[valid]
         positions = positions[valid]
         widths = widths[valid]
         print consensus
         if consensus != None: consensus = consensus[valid]
         # Readjust last position
         ls -= nans.sum()
      # Add the Data
      for p in xrange(self.pages): # l # Page Number in Logo
         self.OUT.write(("%%Page: %d %d\n\n")%(p+1, self.pages))
         self.OUT.write("StartLogo\n")
         for l in xrange(self.lpp): # i # Line Number in Page
            self.OUT.write(" StartLine\n")
            for s in xrange(self.spl): # j # Stack Number in Line
               cs = (p*self.lpp+l)*self.spl+s # Current Stack
               w, pos = widths[cs], positions[cs] # Set Width and Position Number
               hs = sorted(zip(heights[cs], symbols)) # Heights and Symbols zipped and sorted
               if pos%self.xTic != 0: pos = ''
               elif consensus is not None: pos = str(pos) + consensus[cs]
               self.OUT.write("   (%s) StartStack\n"%(pos))
               neg = [] # array of negative symbols
               for h, s in hs:
                  if h >= 0: self.OUT.write("      %f %f (%s) ShowSymbol\n"%(w, h, s))
                  else:          neg.append("      %f %f (%s) ShowSymbol\n"%(w, h, s))
               if len(neg) > 0: #  and showDepletion  #removed since Shannon can't have negative heights
                  self.OUT.write("   DepletedStack\n")
                  self.OUT.write(''.join(reversed(neg)))
               self.OUT.write("   EndStack\n")
               if cs >= ls: # last stack reached
                  done = True
                  break
            self.OUT.write(" EndLine\n")
            if done: break
         self.OUT.write("EndLogo\n\n")
         if done: break
   def Finish(self):
      ''' Closing Logo File and Returns the Path. '''
      self.OUT.write("%%EOF")
      self.OUT.close()
      return self.fn
   def __SetPageBoundries(self, pageSize='640x480'):
      ''' Set Page Size and correction of stacks per line and lines per page, to
      fit the actually observed, if less than specified (private method). '''
      stdPageSizes = { 'A4' : (617.14, 872.82) }
      if pageSize in stdPageSizes: self.pW, self.pH = stdPageSizes[pageSize]
      else:
         try:    self.pW, self.pH = [int(x) for x in pageSize.split('x')]
         except: self.pW, self.pH = 640, 480 # Default
      linesInLogo = ceil(self.seqlen/float(self.spl))
      if self.spl > self.seqlen:
         self.spl = self.seqlen
         self.lpp = 1
      elif self.lpp > linesInLogo:
         self.lpp = int(linesInLogo)
      self.pages = int(ceil(linesInLogo / self.lpp)) # number of pages
   def __SetYaxis(self, userYaxis=[0.0, 0.0]):
      ''' Set Y-axis (private method). '''
      # Determine and set boundries on the Y-AXIS
      ph = np.copy(self.heights); ph[ph<0]=0 # collection of positive heights
      nh = np.copy(self.heights); nh[nh>0]=0 # collection of negative heights
      self.minHeight, self.maxHeight = abs(np.nanmin(np.sum(nh,1))), np.nanmax(np.sum(ph,1))
      if self.minHeight < userYaxis[0]: self.minHeight = userYaxis[0]
      if self.maxHeight < userYaxis[1]: self.maxHeight = userYaxis[1]
      # Calculating best matching Y-axis tics.
      th = self.minHeight + self.maxHeight
      self.yTicLarge, self.yTicSmall = th/5, th/25 # Default
      d={
         1  : (0.2, 0.04),
         1.5: (0.5, 0.1),
         2  : (0.5, 0.1),
         3  : (0.5, 0.25),
         4  : (1.0, 0.2),
         5  : (1.0, 0.2),
         6  : (1.0, 0.2),
         7  : (1.0, 0.25),
         8  : (2.0, 0.4),
         9  : (2.0, 0.4) # 3.0, 0.5
      }
      s = str(th).strip('0.').replace('.','')
      if len(s)>1: x=float(s[:2])/10
      else: x = float(s)
      # Matching best fitting tic interval
      if x in d: self.yTicLarge, self.yTicSmall = d[x]
      else:
         x=int(x*2)/2.0
         if x in d: self.yTicLarge, self.yTicSmall = d[x]
         elif int(x) in d: self.yTicLarge, self.yTicSmall = d[int(x)]
      # Calculating the proportion p
      c = str(th/x).split('.')
      if int(c[0]) > 0: p = 10.0**(len(c[0])-1)
      else: p = 1/10.0**(len(c[1])-len(c[1].lstrip('0'))+1)
      if x in d or int(x) in d:
         self.yTicLarge *= p
         self.yTicSmall *= p
   def __SetXaxis(self, pos, xTic=0):
      ''' Set Numbered Tic Interval for the X-axis (private method). '''
      if xTic < 1:
         if self.seqlen <= 10:   self.xTic = 1
         elif self.seqlen <= 20: self.xTic = 2
         else:                   self.xTic = 5
      else: self.xTic = xTic
      # ENABLE X-AXIS NUMBER ROTATION IF THERE IS NOT ENOUGH SPACE
      if self.xTic*2 < len(str(np.nanmax(np.absolute(pos)))):
         self.rotNum = 'true'
         RaiseWarning('X_num_rot') #WARN
   def __SetColor(self, colors=None):
      ''' Set Color dictionary (private method). '''
      color={
         'black' : "[0.0 0.0 0.0]",
         'red'   : "[0.9 0.0 0.0]",
         'blue'  : "[0.0 0.0 1]",
         'green' : "[0.0 0.85 0.0]",
         'yellow': "[0.9 0.9 0.1]",
         'purple': "[0.6 0.1 0.9]",
         'orange': "[1 0.6 0]"
      }      
      if colors == None: colors = [('DE', color['red']), ('NQSGTY', color['green']), ('HKR', color['blue'])] # Default
      else: colors = dict([[self.hex2pscolor(c) if not c in color else color(c), ss] for c, ss in colors.items()])
      self.cDict = ''.join(["(%s) %s\n"%(s, c) for c, ss in colors.items() for s in ss])
   def hex2pscolor(self, s):
      return "[%s]"%(' '.join(str(x) for x in self.hex2percent(s)))
   def hex2percent(self, s, p=2):
      return [round((int(s[i:i+p],16))/255.0,2) for i in range(0,len(s),p)]
   def __getitem__(self, key):
      ''' This function returns the internal variable 'key' (private method) '''
      return getattr(self,key)


class Reg:
   '''
   NAME:        Reg - Extended Regular Expression Handler
   AUTHOR:      Martin Thomsen
   DESCRIPTION:
      This class enables a simplistic usage of regular expression to get
      contained groups in a match statement. But it also allows to do some of
      the normal re call, such as findall and sub.
   DEPENDENCIES:
      re (regular expression module)
   USAGE:
      >>> RegEx = Reg(pattern, flag)
      >>> if RegEx.match(string):
      >>>    RegEx.getgroup(index)
   EXAMPLE:
      >>> RegEx = Reg('[^a]*(a)[^b]*(b)[^Y]*(Y)(a)?', re.I)
      >>> if RegEx.match('aBcdefgHIJKLmnOpqrstuvwxyz'):
      ...    print(RegEx.getgroup(0), # index=0 -> full match
      ...          RegEx.getgroup(1),
      ...          RegEx.getgroup(2),
      ...          RegEx.getgroup(3),
      ...          RegEx.getgroup(4))
      ... 
      ('aBcdefgHIJKLmnOpqrstuvwxy', 'a', 'B', 'y', None)
      # NIFTY SUBSTITUTION LOOP
      >>> string = 'There are {%=count%} {%=animal%} on the {%=location%}!'
      >>> # Dictionary containing place-holders and values
      ... # (make sure all placeholders in the string is included!)
      ... d = { 'count': 5, 'animal': 'cows', 'location': 'Battle Field' }
      >>> # RE Object which matches place-holders in the string
      ... tmpPH = Reg('\{\%\=(\w+)\%\}', re.I)
      >>> # substitute all placeholders
      ... while tmpPH.match(string): string = tmpPH.sub(str(d[tmpPH.getgroup(1)]), string, 1)
      ...
      >>> print(string)
      There are 5 cows on the Battle Field!
      '''
   def __init__(self, pattern, flag=''):
      if flag: self.re = re.compile(pattern, flag)
      else: self.re = re.compile(pattern)
      self.matches = None
   def sub(self, replace, string, count=0):
      ''' returns new string where the matching cases (limited by the count) in
      the string is replaced. '''
      return self.re.sub(replace, string, count)
   def findall(self, s):
      ''' Finds all matches in the string and returns them in a tuple. '''
      return self.re.findall(s)
   def match(self, s):
      ''' Matches the string to the stored regular expression, and stores all
      groups in mathches. Returns False on negative match. '''
      self.matches = self.re.search(s)
      return self.matches
   def getgroup(self, x):
      ''' Returns requested subgroup. '''
      return self.matches.group(x)
   def getgroups(self):
      ''' Returns all subgroups. '''
      return self.matches.groups()

class Logging:
   '''
   NAME:        Logging - Global Logging Manager
   AUTHOR:      Martin Thomsen
   DESCRIPTION:
      This class is used to manage all Seq2Logo errors, warning and debug
      loggings, and present them to the user in a graceful manner.
   DEPENDENCIES:
      sys (system module)
   USAGE:
      >>> RaiseError = Error(webmode=False, webbin='/path/to/web/bin/').RaiseError
      >>> RaiseError("ErrorType", LineCount, die=True)
   '''
   def __init__(self, ltype, *args):
      if ltype == 'error':
         if args: self.alphabet, self.gaps, self.unknown = args
         else: self.alphabet, self.gaps, self.unknown = None, None, None
         self.type = None
         self.cause = None
         self.help = None
      if ltype == 'warn':
         self.warnings = []
   def SetAlphabet(self, alphabet, gaps='', unknown=''):
      self.alphabet = alphabet
      self.gaps = gaps
      self.unknown = unknown
   def ParsingError(self, eid, line='', e=''):
      '''  '''
      msg = {
         'missing sequence'    : {'cause': 'The error was caused by a missing or unrecognized sequence.', 'help': 'To fix this problem, try checking your input file for symbols which are not in the provided alphabet, and if it is Clustal data then check if all the sequences are repeated correctly.', 'alphabet': True},
         'unknown format'      : {'cause': '%sSeq2Logo could not determine the data type.'%e, 'help': 'Please check that you provided a correct alphabet. Refer to the instructions and examples for additional help.', 'alphabet': True},
         '*too many sequences' : {'cause': 'The error was caused by an inconsistent number of sequences or unrecognised symbols.', 'help': 'To fix this problem, try checking your input file for symbols which are not in the provided alphabet, and if it is Clustal data then check if all the sequences are repeated correctly.', 'alphabet': True},
         'unrecognised line'   : {'cause': 'The error was caused by the use of an unknown symbol.\n%s'%e, 'help': 'To fix this problem, try checking your input file for symbols which are not in the provided alphabet.', 'alphabet': True},
         'weight_err'          : {'cause': 'The error was caused by the use of an unknown symbol or incorrect structure.', 'help': 'To fix this problem, compare your input file with the examples to find the error.', 'alphabet': False},
         'weight_vary'         : {'cause': 'The error was caused by the number of weights not matching the alphabet length.', 'help': 'Please make sure there is the same number of weights as there is symbols in the header line.', 'alphabet': False},
         '*nan'                : {'cause': 'The error was caused by a non numerical weight.', 'help': 'Please make sure there are no unaccepted symbols in the format, and that the correct structure is followed.', 'alphabet': False},
         'seq_vary'            : {'cause': 'The error was caused by one or more sequences which vary in length.', 'help': 'Please make sure all sequences have equal length, and that only valid symbols are used.', 'alphabet': True},
         '*empty alignment'    : {'cause': 'No alignment was found.', 'help': 'Please refer to the examples for help.', 'alphabet': False},
         'seq_empty'           : {'cause': 'Empty sequence found', 'help': 'To fix this problem, try checking your input file for missing sequences or symbols which are not in the provided alphabet.', 'alphabet': False},
         'seq_unknown'         : {'cause': 'Unknown symbol (%s) encountered in a sequence.'%e, 'help': 'You could either remove the symbol, or provide a background frequencies for all symbols.', 'alphabet': False},
         'clustal_err'         : {'cause': 'The error was caused by a missing or unrecognized sequence.', 'help': 'To fix this problem, try checking your input file for symbols which are not in the provided alphabet.', 'alphabet': True},
         'ref_weight_vary'     : {'cause': 'The error was caused by an incorrect number of weights in the reference distribution compared to the provided alphabet.', 'help': 'Please make sure there is the same number of weights as there is symbols in the header line.', 'alphabet': False},
         '' : {'cause': '', 'help': '', 'alphabet': False}
      }
      msg = msg[eid]
      self.type    = 'Error occurred during parsing of your input file.'
      self.cause   = msg['cause']
      self.help    = msg['help']
      if isinstance(line, tuple): self.cause += '\n   Line %d ==> "%s"'%(line)
      elif isinstance(line, int): self.cause += ' (Line %d)'%(line)
      if msg['alphabet'] and self.alphabet + self.gaps + self.unknown != '':
         self.help += '\nAlphabet: %s\nGaps:     %s\nUnknown:  %s'%(self.alphabet, self.gaps, self.unknown)
      self.__RaiseError()
   def DataError(self, eid, line='', e=''):
      '''  '''
      msg = {
         'set_sw'         : {'cause': 'Seq2Logo could not apply the sequence weighting.', 'help': e},
         'freq_err'       : {'cause': 'Calculation of frequencies was not possible due to:\n%s'%e, 'help': 'Please contact the support officer by email [mcft@cbs.dtu.dk].'},
         'bgfreq_invalid' : {'cause': 'Invalid Background frequencies.', 'help': 'Make sure the background frequencies sum up to 1 and there are no values equal to or below 0. Also check that there is a background frequency for each symbol in the alphabet used. If you do not have any background frequencies for your data, then switch to another logo type. Eg. Shannon or PSSM.'},
         'blosum_invalid' : {'cause': 'The Blosum Matrix did not contain all Symbols found in the background frequency file.', 'help': 'Please fix this!\nA quick fix would be to set beta to 0.'},
         'gs_invalid'     : {'cause': 'Call to executeGS failed.', 'help': 'Invalid inputs.'},
         'gs_err'         : {'cause': 'Ghostscript Error.', 'help': 'Could not convert to %s'%e},
         'weight_zero'    : {'cause': 'No PSSM data was found/calculated.', 'help': 'If your alignment does infact contain data, please contact the support officer by email [mcft@cbs.dtu.dk].'},
         'prob_zero'      : {'cause': 'No probalities data was found/calculated.', 'help': 'If your alignment does infact contain data, please contact the support by email [mcft@cbs.dtu.dk].'},
         'height_zero'    : {'cause': 'No height data was calculated due to internal error.', 'help': 'Please contact the support officer by email [mcft@cbs.dtu.dk].'},
         'aln_invalid'    : {'cause': 'The alignment is invalid.', 'help': 'Please see instructions for help.'},
         '' : {'cause': '', 'help': ''}
      }
      msg = msg[eid]
      self.type    = 'Error occurred during data processing.'
      self.cause   = msg['cause']
      if line: self.cause += '(Line %d)'%(line)
      self.help    = msg['help']
      self.__RaiseError()
   def RuntimeError(self, eid, e=''):
      '''  '''
      msg = {
         'logo_invalid' : {'cause': 'Invalid logo options was found.', 'help': 'Please contact the support officer by email [mcft@cbs.dtu.dk].'},
         '' : {'cause': '', 'help': ''}
      }
      msg = msg[eid]
      self.type    = 'Error occurred during parsing of your input file.'
      self.cause   = msg['cause']
      self.help    = msg['help']
      if e: self.cause += '\n%s'%(e)
      self.__RaiseError()
   def ArgError(self, eid, e=''):
      '''  '''
      msg = {
         'segment_invalid' : {'cause': 'Selected segment is invalid.', 'help': 'Use the format start-end eg. 5-10'},
         '' : {'cause': '', 'help': ''}
      }
      msg = msg[eid]
      self.type    = 'Error occurred during parsing of the arguments.'
      self.cause   = msg['cause']
      self.help    = msg['help']
      self.__RaiseError()
   def __RaiseError(self):
      '''  '''
      sys.exit('\n'.join(['', '='*80, 'ERROR', self.__Wrap(self.type), self.__Wrap(self.cause), '', self.__Wrap(self.help), '='*80, '']))
   def RaiseWarning(self, wid, other=''):
      '''  '''
      msg = {
         'weight_neg_inf'        : {'cause': 'One or more weights was less than %s.'%(other), 'action': 'To avoid such amino acid dominating the logo position, these was reduced to %s.'%(other), 'solution': 'To remove this warning you could try choosing the option of pseudo counts.'},
         'X_num_rot'             : {'cause': 'There was not enough space to fit x-axis labels without overlapping.', 'action': 'The x-axis numbers was rotated 90 degrees to fix this.', 'solution': 'To avoid rotation of the x-axis numbers, please change the frequency by which the x-axis number is shown.'},
         'bgs_nomatch'           : {'cause': 'The alphabet from the background frequencies did not match that found in the input file.', 'action': 'A flat background distribution was used as default.', 'solution': 'To avoid this error, please provide a background frequency pssm which matches the the alphabet in the input file.'},
         'flat_bgs_found'        : {'cause': 'Conversion of the weights back to frequencies had the best result using a flat background distribution.', 'action': 'A flat background distribution was used as conversion to probabilities, and the weight were then recalculated using the chosen background distribution.', 'solution': 'If this action is unwanted then please upload a fitting background distribution which were used in the generation of the log-odds scores.'},
         'formats_invalid'       : {'cause': 'The chosen formats were invalid.', 'action': 'The EPS was not converted to other formats.', 'solution': 'Please fix the format argument. Check the instructions for additional help.'},
         'ref_dist_zero'         : {'cause': 'One or more positions had an invalid reference distribution containing zero or negative value. (%s)'%(other), 'action': 'The reference distribution for these position was appointed a uniform distribution.', 'solution': 'Please fix your reference distribution. Values of zero or below is impossible, all distributions should add up to 1 for each position.'},
         'ref_dist_too_few_pos'  : {'cause': 'Provided reference distribution contained insufficient positions.', 'action': 'the missing positions was appointed a uniform distribution.', 'solution': 'To avoid this error, please provide a reference distribution which matches your input data.'},
         'ref_dist_too_many_pos' : {'cause': 'Provided reference distribution contained too many positions.', 'action': 'The excess was removed.', 'solution': 'Please fix your reference distribution to match your input file.'},
         'ref_dist_uniform'      : {'cause': 'Provided reference distribution was invalid, wrong number of columns.', 'action': 'An uniform distribution was used for all positions.', 'solution': 'Please fix your reference distribution to match your input file.'},
         'no_data_pos'           : {'cause': 'One or more positions had no data associated. (%s)'%(other), 'action': 'These positions was removed from the logo.', 'solution': 'Please fix your input file. Every position should contain data.'},
         '' : {'cause': '', 'action': '', 'solution': ''}
      }
      msg = msg[wid]
      self.warnings.append([msg['cause'], msg['action'], msg['solution']])
   def GetWarnings(self):
      '''  '''
      if self.warnings: return '\n'.join(['='*80, 'WARNINGS\n', '\n\n'.join(['\n'.join([self.__Wrap(x) for x in xs]) for xs in self.warnings]),'='*80, ''])
      else: return ''
   def __Wrap(self, s):
      '''  '''
      return '\n'.join([textwrap.fill(x, 80) for x in s.split('\n')])


# INITIALIZE ERROR HANDLING
Error = Logging('error')
DataError = Error.DataError
ParsingError = Error.ParsingError
ArgError = Error.ArgError
RuntimeError = Error.RuntimeError
Warning = Logging('warn')
RaiseWarning = Warning.RaiseWarning


class readstream():
   '''
   DESCRIPTION:
      This class can wrap txt files and stdin streams in a line seekable
      structure. Which makes it esier to work across file input and standard
      input.
   EXAMPLE:
      with open('test_data/fasta.txt') as fobj:
      f = readstream(fobj)
      for l in f:
         print l.strip()
         if l.strip() == '>sequencename 5': 
            print 'hello'
            break
      for l in f: print l.strip()
   '''
   def __init__(self, fileobj):
      self.f = fileobj
      self.lines = []
      self.ln = 0
   def readline(self):
      l = self.f.readline()
      self.lines.append(l)
      self.ln += 1
      return l
   def __iter__(self):
      ll = len(self.lines)
      if ll > self.ln:
         for i in xrange(self.ln, ll):
            yield self.lines[i]
      eof = 0
      while not eof:
         l = self.readline()
         if not l: eof = 1
         yield l
   def seekline(self, ln):
      if ln <= len(self.lines):
         self.ln = ln
      else:
         sys.exit('ERROR: seeked line (%s) does not exist!\n'%ln)


################################################################################
##########                      Function Library                      ##########
################################################################################
def Hobohm(alignment, alphabet='FAMILVWQSTYNCGPHKRDE', gaps='-', threshold=0.63):
   '''
   NAME:        Hobohm - Hobohm Clustering Algorithm 1 (Sequence Weighting)
   AUTHOR:      Martin Thomsen
   DESCRIPTION:
      This algorithm finds similar sequences and clusters them together. After
      which it calculates the sequence weights as 1 divided by the number of
      sequences in the corresponding cluster.
      The method returns position specific alphas (number of clusters - gaps -1)
      and a list of weights corresponding to the sequences in the alignment.
   REFERENCE:
      Uwe Hobohm, Michael Scharf, Reinhard Schneider and Chris Sander
      'Selection of representative protein data sets'
      Protein Science, 1992, I: p409-417
   DEPENDENCIES:
      Numpy as np (Numeric array module)
   EXAMPLE:
      >>> import numpy as np
      >>> Hobohm(['ABCDE','ABCDE','ABCDF','-BDEF','A--EF','A--EF'], 'ABCDEF',
      ...        '-', 0.6)
      (array([ 1.,  1.,  1.,  2.,  2.]), (0.3333333333333333,
      0.3333333333333333, 0.3333333333333333, 1.0, 0.5, 0.5))
      >>> Hobohm(['ALAKAAAAM', 'ALAKAAAAN', 'ALAKAAAAR', 'ALAKAAAAT',
      ... 'ALAKAAAAV', 'GMNERPILT', 'GILGFVFTM', 'TLNAWVKVV', 'KLNEPVLLL',
      ... 'AVVPFIVSV'])
      (array([ 5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.]),
      (0.2, 0.2, 0.2, 0.2, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0))
      >>> Hobohm([[0,1,0,2,0,0,0,0,3], [0,1,0,2,0,0,0,0,4], [0,1,0,2,0,0,0,0,5],
      ... [0,1,0,2,0,0,0,0,6], [0,1,0,2,0,0,0,0,7], [8,3,4,9,5,10,11,1,6],
      ... [8,11,1,8,12,7,12,6,3], [6,1,4,0,13,7,2,7,7], [2,1,4,9,10,7,1,1,1],
      ... [0,7,7,10,12,11,7,14,7]])
      (array([ 5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.]),
      (0.2, 0.2, 0.2, 0.2, 0.2, 1.0, 1.0, 1.0, 1.0, 1.0)
   '''
   a = np.array(alignment)
   if a.dtype.type == np.str_: # Convert alignment to indexes
      s2i = dict(zip(list(alphabet+gaps),range(len(alphabet))+[float('NaN'),]*len(gaps)))
      a = np.array([[s2i[s] for s in seq] for seq in a])
   if len(a.shape)==1: a = np.array([list(seq) for seq in a])
   if len(a.shape)==2:
      numseq, seqlen = a.shape
      # SORT SEQUENCES according to their gap percentage ascending.
      sort=list(zip(*sorted(zip(np.sum(a!=a,1),np.arange(a.shape[0]))))[1])
      a=a[sort]
      # CALCULATE DISTANCE MATRIX, based on differences not including gaps
      d=np.zeros([numseq]*2) # Init. distance matrix
      for i in xrange(numseq): d[i+1:,i] = d[i,i+1:] = 1-np.sum(a[i]==a[i+1:],1)/float(np.sum(a[i]==a[i]))
      # APPOINT SEQUENCES TO CLUSTERS DEPENDING ON THEIR DISTANCES
      c = {} # Cluster
      s = np.ones(numseq) # Array showing which sequences aren't assigned to a cluster yet
      for i in xrange(numseq):
         if s[i]: # Unassigned => Find matches which are unassigned and below threshold
            m = np.nonzero(np.logical_and(d[i,i+1:]<1-threshold, s[i+1:]))[0]+i+1
            c[i]= list(m) if m.size else [] # Add matching sequences to cluster
            if m.size: s[m]=0 # Change status of matched sequences to assigned
      # COMPUTING CLUSTER WEIGHTS
      w = zip(*sorted([(x,1.0/len([k]+v)) for k,v in c.items() for x in [k]+v]))[1]
      # REVERSE SORTING *Note not reversing c since we only use the length
      w = zip(*sorted(zip(sort, w)))[1]
      return len(c)-1-np.sum((a!=a).T*w,1), w
   else: DataError('aln_invalid') #ERROR


def Heuristic(alignment, alphabet='FAMILVWQSTYNCGPHKRDE', gaps='-'):
   '''
   NAME:        Heuristic - Heuristic Sequence Weighting Algorithm
   AUTHOR:      Martin Thomsen
   DESCRIPTION:
      This is a sequence weighting algorithm which builds on a heuristic
      approach developed by Henikoff and Henikoff (see Ref).
      It returns alpha (average number of unique symbols at each position -1),
      and a list of sequence weights corresponding to the sequences in the
      alignment.
   REFERENCE:
      Steven Henikoff and Jorja G. Henikoff,
      'Position-based Sequence Weights'
      Journal of Molecular Biology, 1994, 243: p574-578
   DEPENDENCIES:
      Numpy as np (Numeric array module)
   EXAMPLE:
      >>> import numpy as np
      >>> Heuristic(['ALAKAAAAM', 'ALAKAAAAN', 'ALAKAAAAR', 'ALAKAAAAT',
      ...            'ALAKAAAAV', 'GMNERPILT', 'GILGFVFTM', 'TLNAWVKVV',
      ...            'KLNEPVLLL', 'AVVPFIVSV'])
      (3.777777777777778, [0.414047619047619, 0.49738095238095237,
         0.49738095238095237, 0.414047619047619, 0.38626984126984121,
         1.3583333333333334, 1.4583333333333333, 1.2746031746031745,
         1.1857142857142859, 1.5138888888888891])
      >>> Heuristic(alignment=['AGCTC','AGCTT','AGCTC','-GTCT','A--TC'],
      ...           alphabet='AGTC', gaps='-')
      (0.6000000000000001, [0.95833333333333326, 1.0416666666666665,
      0.95833333333333326, 1.5, 0.54166666666666663])
      >>> # NUMERIC SEQUENCE ALIGNMENT -----------------------------------------
      >>> Heuristic([[0,1,0,2,0,0,0,0,3], [0,1,0,2,0,0,0,0,4],
      ... [0,1,0,2,0,0,0,0,5], [0,1,0,2,0,0,0,0,6], [0,1,0,2,0,0,0,0,7],
      ... [8,3,4,9,5,10,11,1,6], [8,11,1,8,12,7,12,6,3], [6,1,4,0,13,7,2,7,7],
      ... [2,1,4,9,10,7,1,1,1], [0,7,7,10,12,11,7,14,7]])
      (3.7777777777777777, [0.414047619047619, 0.49738095238095237,
      0.49738095238095237, 0.414047619047619, 0.38626984126984121,
      1.3583333333333334, 1.4583333333333333, 1.2746031746031745,
      1.1857142857142859, 1.5138888888888891])
   '''
   a = np.array(alignment)
   if a.dtype.type == np.str_: # Convert alignment to indexes
      s2i = dict(zip(list(alphabet+gaps),range(len(alphabet))+[float('NaN'),]*len(gaps)))
      a = np.array([[s2i[s] for s in seq] for seq in a])
   if len(a.shape)==1: a = np.array([list(seq) for seq in a])
   if len(a.shape)==2:
      seqlen = a.shape[1]
      ai = np.arange(len(alphabet)) # alphabet index
      # Calculate Wkps
      Wkp = []
      rsum = 0
      for p in range(seqlen):
         pos = np.hsplit(a,seqlen)[p] # Alignment segment at position p 
         s = np.sum(ai==pos,0)          # Symbols count at position p
         r = np.sum(s!=0)                     # number of unique symbols at position p
         with np.errstate(divide='ignore'):
            irs = 1.0/(s*r)                   # inverted r*s
         rsum += r
         Wkp.append([irs[s] if s==s else 0 for s in pos.ravel()])
      # Sum Wks
      Wk = np.sum(np.array(Wkp).T, 1)         # Sequence weights
      return float(rsum)/seqlen-1, list(Wk)
   else: DataError('aln_invalid') #ERROR


# EPS CONVERSION MODULE
def EPS2Formats(outputfile, formats, pageWidth, pageHeight, gsPath='gs'):
   '''
   NAME:        EPS2Formats - EPS Conversion Manager
   AUTHOR:      Martin Thomsen
   DESCRIPTION:
      This method manages the conversion of the provided EPS file to the
      requested formats.
   SUPPORTED FORMATS:
      EPS -> JPEG
      EPS -> PNG
      EPS -> PDF
      EPS -> SVG (only works for Ghostscript v9.05 or later)
   DEPENDENCIES:
      sys (system module)
      ExecuteGS
      Reg
   USAGE:
      >>> EPS2Formats(outputfile, formats, pageWidth, pageHeight)
   '''
   # CHECKING INPUT
   if not isinstance(outputfile, str) or not isinstance(formats, str) or not isinstance(pageWidth, (str,int,float)) or not isinstance(pageHeight, (str,int,float)):
      RuntimeError('logo_invalid', 'Call to EPS2Formats failed. Invalid inputs!')
   else:
      pageWidth = str(pageWidth)
      pageHeight = str(pageHeight)
   # Execute ghostscript conversions
   if "JPEG" in formats: #JPEG - always show at least jpeg image if not png requested
      #  -dEPSCrop Crop an EPS file to the bounding box. This is useful when converting an EPS file to a bitmap.
      specialArgs =  " -dGraphicsAlphaBits=4 -dTextAlphaBits=4 -dAlignToPixels=0"
      ExecuteGS(outputfile, pageWidth, pageHeight, "jpeg", "-%03d.jpg", specialArgs, "JPEG", gsPath=gsPath)
   if "PNG" in formats:
      specialArgs =  " -dGraphicsAlphaBits=4 -dTextAlphaBits=4 -dAlignToPixels=0"
      ExecuteGS(outputfile, pageWidth, pageHeight, "png16m", "-%03d.png", specialArgs, "PNG", gsPath=gsPath)
   if "PDF" in formats: 
      specialArgs =  " -dAutoRotatePages=/None -c '<< /PageSize ["+pageWidth+" "+pageHeight+"] /Orientation 0>> setpagedevice'"
      ExecuteGS(outputfile, pageWidth, pageHeight, "pdfwrite", ".pdf", specialArgs, "PDF", gsPath=gsPath)
   if "SVG" in formats: 
      specialArgs =  ''#" -dAutoRotatePages=/None -c '<< /PageSize ["+pageWidth+" "+pageHeight+"] /Orientation 0>> setpagedevice'"
      if ExecuteGS(outputfile, pageWidth, pageHeight, "svg", ".svg", specialArgs, "SVG", gsPath=gsPath):
         SVGerr = Reg('([^=]\'\#[0-9a-fA-F]{6})', re.I)
         # Fixing GS Error, where ' fill=' is missing in some instances of a 'g' element
         with open(outputfile+'.svg', 'r') as f, open('out_fixed.svg', 'w') as o:
            for l in f:
               if SVGerr.match(l): # SVG Error found
                  g = SVGerr.getgroup(1)
                  l = l.replace(g, "%s fill='%s'"%(g[:2],g[2:])) # Correct Error
               o.write(l)
         os.rename('out_fixed.svg', outputfile+'.svg')
   if not 'specialArgs' in locals() and formats != '': RaiseWarning('formats_invalid')


def ExecuteGS(outputfile, pageWidth, pageHeight, device, suffix, specialArgs='', name='', gsPath = 'gs'):
   '''
   NAME:        ExecuteGS - Ghostscript Execution Wrapper
   AUTHOR:      Martin Thomsen
   DESCRIPTION:
      This function executes a ghostscript command.
   DEPENDENCIES:
      sys (system module)
      os  (operation system module)
   USAGE:
      >>> executeGS(outputfile, formats, pageWidth, pageHeight, "jpeg",
      ...           "-%03d.jpg", specialArgs, "JPEG")
   '''
   # CHECKING INPUTS
   if not isinstance(outputfile, str) or not isinstance(pageWidth, (str,int,float)) or not isinstance(pageHeight, (str,int,float)) or not isinstance(device, str) or not isinstance(suffix, str):
      DataError('gs_invalid') #ERROR
   else:
      pageWidth = str(pageWidth)
      pageHeight = str(pageHeight)
   # SETTING GHOSTSCRIPT COMMAND
   cmd = (gsPath+" -dNOPAUSE -dBATCH -dSAFER -dQUIET -sDEVICE="+device+" -sOutputFile="+outputfile+suffix+ specialArgs+
          " -dDEVICEWIDTHPOINTS="+pageWidth+" -dDEVICEHEIGHTPOINTS="+pageHeight+
          " -g"+pageWidth+"x"+pageHeight+" "+outputfile+".eps")
   try: os.system(cmd)
   except OSError:
      DataError('gs_err', msg=name) #ERROR
      return False
   else: return True


def savePSSM(filename, pos, alphabet, weights, probs, args, yLabel, logoType, header="Last position-specific scoring matrix computed, values are in halfbits\n"):
   '''
   NAME:        savePSSM - PSSM File Creator
   AUTHOR:      Martin Thomsen
   DESCRIPTION:
      This function writes a PSSM file.
   DEPENDENCIES:
      ComputeWeightMatrix
   USAGE:
      >>> savePSSM("PSSM.txt", positions, alphabet, weights, CMDlineArgs)
   '''
   with open(filename,'w') as f:
      # WRITE HEADERS
      f.write("#Seq2Logo.py %s\n"%(' '.join(args))) # CMDLINE
      f.write(header)
      # WRITE WEIGHT MATRIX
      f.write(ComputeWeightMatrix(pos, alphabet, weights, probs, yLabel, logoType))


def ComputeWeightMatrix(pos, alphabet, weights, probs, yLabel, logoType):
   '''
   NAME:        ComputeWeightMatrix - Weight Matrix Renderer
   AUTHOR:      Martin Thomsen
   DESCRIPTION:
      This function renders and returns a weight matrix which consist of an
      alphabet, position numbers, concensus sequence and weights.
   DEPENDENCIES:
      Numpy as np (Numeric array module)
   EXAMPLE:
      >>> print( ComputeWeightMatrix(np.array([1,2,3]), np.array(list('ABC')),
      ...                            np.array([[3,1,2],[4,6,5],[8,7,9]])) )
      A B C
      1 A 6 2 4
      2 B 8 12 10
      3 C 16 14 18
   '''
   # COMPUTE WEIGHT MATRIX - weights converted to halfbits and rounded to 3 decimals
   p = len(pos)
   if not probs is None: consensus = np.array([alphabet[i] for i in np.argmax(probs, axis=1)], dtype='|S1')
   else: consensus = np.array(['']*p, dtype='|S1')
   if yLabel.lower()!='halfbits' and logoType!=5: # skipping PSSM logos and weights already in halfbits
      weights*=2
      weights[weights<-99.999]= -99.999
   w = weights.round(decimals=3).astype('|S10')
   wm = np.concatenate((pos.reshape(p,1), consensus.reshape(p,1), w), axis=1)
   return '\n'.join([' '.join(alphabet)]+[' '.join([str(c) for c in r]) for r in wm])

def ParseArguments(dirbin):
   # INITIALIZE THE PARSER
   parser = ArgumentParser(
      formatter_class=RawDescriptionHelpFormatter,
      description=('''
********************************************************************************

NAME:        Seq2Logo
VERSION:     2.1 (Last modified: 20. Jun 2014)
AUTHOR:      Martin Christen Frolund Thomsen

DESCRIPTION:
Seq2Logo takes a multiple sequence alignment (MSA) or a weight matrix (PSSM)
and creates a sequence logo.

INPUT:
Following formats are accepted:
   MSAs:   Fasta and ClustalW
   PSSMs:  BLASTmatrix, weightmatrix, frequency table
   OTHER:  Aligned raw peptide sequences
To see a full list, and examples of these go to:
http://www.cbs.dtu.dk/biotools/Seq2Logo-2.1/instructions.php#formats

OUTPUT:
The program produces at least two files:
   1. A sequence logo in eps format. (eg. output.eps)
   2. A Weightmatrix as txt file.    (eg. output.txt)
In addition to that, the user can also specify to get the eps, converted to:
PDF, JPEG, SVG and/or PNG.

INSTRUCTIONS:
1. Select a logo type. ('-I', 5 different options)
   - If you are submitting a PSSM and thus just want to show the data 'as is',
      it is recommended to use the PSSM logo option 5
2. Select clustering type. ('-C', 3 options)
   - This step is recommended if you have a lot of unfiltered data which
      includes dublicates etc. and you want to reduce the impact of the
      redundancy.
   - Otherwise set this to 0 to skip clustering.
3. Set Threshold for the Hobohm algorithm. ('-T', percentage)
   - default is 0.63, which means that two sequences are clustered togather if
      more than 63% of the longest sequence is identical with the other
      sequence.
   - This can be ignored if the Hobohm clustering algorithm was not selected.
4. Set the weight of prior. ('-b', number)
   - This option is recommended if there are fewer than 50 sequences, a beta on
      200 is reccomended for Hobohm, and 50 is recommended for heuristics
      clustering.
   - If you are submitting a PSSM this step can be ignored.

This was a short guide to the usage of Seq2Logo, to see a full guide, please
visit:
http://www.cbs.dtu.dk/biotools/Seq2Logo-2.1/instructions.php

********************************************************************************
'''), epilog=('''
********************************************************************************
   
'''))
   
   # ADD ARGUMENTS
   parser.add_argument("-f", type=str, dest="inputfile", default='',
               help="Input file path")
   parser.add_argument("-o", type=str, dest="outputfile", default="output",
               help="Output file name eg. 'mylogo' would result in mylogo.eps, mylogo.txt and mylogo-001.jpg etc. [%(default)s]")
   parser.add_argument("-I", type=int, dest="logoType", default=2,
               help="Logo type, defines how the data is represented on the motif. 1:Shanon, 2:Kullback-Leibler, 3:Weighted Kullback-Leibler, 4:P-Weighted Kullback-Leibler, 5:PSSM-Logo. [%(default)s]")
   parser.add_argument("-b", type=float, dest="beta", default=200.0,
               help="Weight on prior, defines how much the prior knowledge of amino acid substitution probabilities (the Blosum matrix) impacts the estimated proabilities. This is great for small datasets as it reduces the impact of outliers and corrects for missing values. [%(default)s]")
   parser.add_argument("-C", type=int, dest="sequenceWeighting", default=2,
               help="Clustering: 0:None, 1:Heuristics and 2:Hobohm. [%(default)s]")
   parser.add_argument("-T", type=float, dest="threshold", default=0.63,
               help="Threshold for hobohm in fractions. [%(default)s]")
   parser.add_argument("-u", type=str, dest="yLabel", default="Bits",
               help="The unit on the Y-axis. Eg. Bits, Halfbits, Z-score etc. NOTE: Only PSSM-logos can show unit types different than Bits and Halfbits. [%(default)s]")
   parser.add_argument("-s", type=int, dest="spl", default=40,
               help="Number of stacks per line. [%(default)s]")
   parser.add_argument("-l", type=int, dest="lpp", default=3,
               help="Number of lines per page. [%(default)s]")
   parser.add_argument("-p", type=str, dest="pageSize", default="640x480",
               help="PageSize Width x Height eg. 640x480. [%(default)s]")
   parser.add_argument("--colors", type=str, dest="colors", default='E60000:DE,00D900:QSTYNG,0000FF:HKR',
               help="String of colors and symbols. Eg. 'red:DE' or 'E60000:DE' to color D and E red. [%(default)s]")
   parser.add_argument("-t", type=str, dest="title", default='',
               help="Title of the logo. [%(default)s]")
   parser.add_argument("-m", type=float, dest="minwidth", default=0.5,
               help="Minimum stackwidth for positions with many gaps in fraction. [%(default)s]")
   parser.add_argument("-S", type=int, dest="startPos", default=1,
               help="Starting position number. [%(default)s]")
   parser.add_argument("-Z", type=str, dest="exclude_zero", default='',
               help="Defines whether to included the zeroes position in the x-axis values, set to 'on' to enable it. [%(default)s]")
   parser.add_argument("-i", type=int, dest="xTic", default=0,
               help="The tic interval on the X-axis, eg. to show the number of every position set this to 1. [%(default)s]")
   parser.add_argument("-c", type=str, dest="segment", default=None,
               help="Cut out a segment of the alignment by setting a start and end position of the segment, eg. 10-35. [%(default)s]")
   parser.add_argument("-y", type=floatrange, dest="userYaxis", default="0:0", action = 'store',
               help="Define the boundries of the y-axis, eg. 4.32:4.32 for an yaxis ranging from -4.32 to 4.32 bits. [%(default)s]")
   parser.add_argument("--consensus", type=str, dest="consensus", default='',
               help="Defines whether to included the consensus in the x-axis values, set to 'on' to enable it. [%(default)s]")
   parser.add_argument("-H", type=str, dest="hide", default='',
               help="Defining which elementents to hide list as 'yaxis,xaxis,yaxis_label,fineprint,rotate_numbers,ends'. [%(default)s]")
   parser.add_argument("--format", type=str, dest="formats", default="JPEG",
               help="Output formats: JPEG, PNG, SVG and/or PDF. To get multiple formats use comma separation, eg. 'JPEG,PNG,SVG,PDF'. [%(default)s]")
   parser.add_argument("--bg", type=str, dest="bgPath", default=dirbin+'B62_Reference_Distribution.pssm',
               help="Path to the file listing the background frequencies. [%(default)s]")
   parser.add_argument("--blosum", type=str, dest="blosumPath", default=dirbin+'Blosum62.txt',
               help="Path to the blosum frequency matrix. [%(default)s]")
   parser.add_argument("-w", action="store_true", dest="webmode", default=False, help=SUPPRESS)
   
   # VALIDATING ARGUMENTS
   args = parser.parse_args()
   if args.segment:
      try: args.segment = [int(x)-1 for x in args.segment.split('-',1)]
      except: ArgError('segment_invalid')
   
   if args.exclude_zero != '': args.exclude_zero = True
   else: args.exclude_zero = False
   
   if args.consensus: args.consensus = []
   else: args.consensus = None
   
   # MANAGE ESCAPED COLOR SYMBOLS
   string = args.colors
   # Dictionary containing place-holders and values
   d = dict([(str(x),chr(x)) for x in range(256)])
   # RE Object which matches place-holders in the string
   tmpPH = Reg('\%(\d+)\%', re.I)
   # substitute all placeholders
   while tmpPH.match(string): string = tmpPH.sub(str(d[tmpPH.getgroup(1)]), string, 1)
   args.colors = string
   
   # MANAGING ARGUMENTS
   args.title = args.title.replace("_", " ")
   args.colors = dict([x.split(':') for x in args.colors.split(',')])
   args.hide = args.hide.split(',')
   args.showX    = not "xaxis" in args.hide           # Show X-axis on Logo
   args.showY    = not "yaxis" in args.hide           # Show Y-axis on Logo
   args.eLabel   = '' if "ends" in args.hide else 'p' # Sequence End Labels (eg. 5') - d: DNA, p: PROTEIN
   args.rotNum   = "rotate_numbers" in args.hide      # Rotate Numbers on X-axis
   args.showFine = not "fineprint" in args.hide
   args.yLabel = args.yLabel.replace("_", " ") if not "yaxis_label" in args.hide else ''
   del(args.hide)
   
   # RETURN ARGUMENT LIST
   return [(arg, getattr(args, arg)) for arg in dir(args) if arg[:1]!='_']

def floatrange(value):
    values = value.split(':')
    if len(values) != 2: raise argparse.ArgumentError
    values = map(float, values)
    return values