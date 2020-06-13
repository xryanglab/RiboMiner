#!/usr/bin/env python2.7
'''
NAME:        Seq2Logo - Sequence Logo Generator
AUTHOR:      Martin Thomsen
DESCRIPTION:
   This program uses the classes and functions in the Seq2Logo module to
   retrieve sequence data from an inputfile (see Seq2Logo_module.SeqParse for a
   complete list of supported formats), calculate a PSSM and create a sequence
   logo of the data.
PLEASE CITE:
   Martin C F Thomsen and Morten Nielsen,
   'Seq2Logo: a method for construction and visualization of amino acid binding
   motifs and sequence profiles including sequence weighting, pseudo counts and
   two-sided representation of amino acid enrichment and depletion.'
   Nucleic Acids Research (2012) 40 (W1): W281-W287
ACKNOWLEDGEMENT:
   WebLogo3.0 played a big inspiration in the creation of Seq2Logo.
      [Crooks et al.,2004]
   Also the original creators of sequence logos should be credited.
      [Schneider & Stephens,1990].
'''
import sys

# CHANGES in version 2.1
# * Now using position specific alpha instead of global alpha, which gives a
#   more correct influence of pseudo counts at positions with many gaps.
# * The alphabet order can now be mixed between input matrices (input pssm,
#   Blosum matrix, background frequencies)
# * Background frequencies can now be provided in any of the supported
#   alignment formats. The Background frequencies will be computed from the
#   alignemt as position specific frequencies
# * Fixed a bug where a flat background distribution was used for PSSM
#   calculations of weight inputs, when it was observed that the weight input
#   was created using a flat background. Now the weights are recalculated using
#   the provided background frequencies.
# * Fix Display of fineprint, so it doesn't overlap the x-axis

# Fixed / TODO
# * Hobohm
# Heuristic
# no clustering
# use of alpha

# alphabet rearrangement function
# alphabet check INPUT
# alphabet check BLOSUM
# alphabet check BG

# BG Parser using a seqpar object

# Display of fineprint



# EXTERNAL APPLICATION PATHS AND ENVIRONMENT VARIABLES
gsPath='gs'
dirbin = "%s/etc/"%(sys.path[0])
sys.path.append(dirbin)

import  numpy as np
from Seq2Logo_module import SeqParse, SequenceData, MakeLogo, EPS2Formats, Heuristic, Hobohm, Logging, savePSSM, ParseArguments, Error, Warning, RaiseWarning

# -------------------- MAP COMMANDLINE ARGUMENTS TO GLOBALS --------------------
map(lambda x, d=locals().__setitem__: d(*x), ParseArguments(dirbin))

################################################################################
#                                     MAIN                                     #
################################################################################
def main():
   global startPos, eLabel, xTic
   global seqpar, seqdat
   # ----------------------- INITIALIZE AND SET PARAMETERS ------------------------
   # SEQ2LOGO VARS
   gaps = '_-.*'
   unknown = 'X'
   eps_template = "%s/logo_template.eps"%(dirbin)
   xLabel = '' # X-axis Unit Label (Not Used)
   
   # INITIALIZE SEQUENCE PARSE AND DATA MANAGER
   seqpar = SeqParse()
   seqdat = SequenceData()
   if consensus is not None: seqdat.set_consensus()
   
   # PARSE BACKGROUND FREQUENCIES
   alphabet, bgfreq = seqpar.parseBackgroundFrequencies(bgPath)
   
   # PARSE BLOSUM MATRIX
   blosum = None
   if beta > 0:
      alphabet2, blosum = seqpar.parseBlosum(blosumPath)
      # Check that there is a BLOSUM entry for each Symbol in the alphabet (the alphabet is defined from the background frequencies)
      if len(np.setdiff1d(alphabet, alphabet2)) > 0: Error.DataError('blosum_invalid')
   
   # INITIALIZE SEQUENCE DATA ENVIRONMENT
   unknown = '' if 'X' in alphabet else 'X'
   seqpar.SetAlphabet(seqdat, ''.join(alphabet), gaps)
   seqdat.set_alphabet(alphabet, blosum, gaps, unknown)
   Error.SetAlphabet(''.join(alphabet), gaps, unknown)
   
   # ------------------------------ PARSE INPUTFILE -------------------------------
   inputType, debug, showcalc = seqpar.parsefile(inputfile)
   sys.stderr.write('inputType: %s\n'%inputType)
   
   # -------------------------------- COMPUTE PSSM --------------------------------
   if inputType == 'WEIGHT':
      xTic = 1 # Show position number for all positions
      if logoType != 5: # not PSSM logo
         # NOTE: it is assumed the weights represents log-odds values either in bits
         # or halfbit, or that they are frequencies/probabilities. This imply that
         # a sparse weight matrix, could prove incompatible. Use PSSM logos for these.
         # MAKE SEGMENT
         if segment:
            seqdat.w = seqdat.w[segment[0]:segment[1]+1]
            seqdat.pos = seqdat.pos[segment[0]:segment[1]+1]
         # Set Number of Sequences
         seqdat.numseq = len(seqdat.w)
         # Get weights of first position andbackground frequencies, to determine weight type.
         weights = seqdat.w[0] # Weights of first position 
         # Calculate sums for testing
         alphabet2 = alphabet
         alphabet = seqdat.alphabet
         if np.in1d(alphabet2, alphabet).all() and np.in1d(alphabet, alphabet2).all():
            alen = len(alphabet)
            bgfreq = np.array(bgfreq, dtype='f8')
            bshape = bgfreq.shape
            seqlen = len(seqdat.pos)
            if bshape[0] == 1:
               bgfreq = bgfreq[0]
            elif bshape[1] < alen:
               # Set uniform reference distribution
               bgfreq = np.ones(alen, dtype='f8')/alen
               RaiseWarning('ref_dist_uniform', '') #WARN
            elif bshape[0] > seqlen:
               # Cut off the excess
               bgfreq = bgfreq[:alen,:]
               RaiseWarning('ref_dist_too_many_pos', '') #WARN
            elif bshape[0] < seqlen:
               # Extend with uniform distribution
               bgfreq = np.vstack([bgfreq, np.ones((seqdat.pos.shape[0]-bshape[0], alen))/alen])
               RaiseWarning('ref_dist_too_few_pos', '') #WARN
         else:
            bgfreq = np.ones(len(alphabet), dtype='f8')/len(alphabet)
            Warning.RaiseWarning('bgs_nomatch')
         bgfreq_flat = np.ones(len(alphabet), dtype='f8')/len(alphabet)
         # Estimate PSSM type
         sumf = np.sum(weights)
         # CHECK IF WEIGHTS ARE FREQUENCIES - Max of 2% deviation from sum==1, and frequencies can't be negative
         if not np.any(weights<0) and sumf >= 0.90 and sumf <= 1.10: # Frequencies/probabilities
            sys.stderr.write('pssmType: Frequencies\n')
            # SET PROBABILITIES - The frequencies are asumed to be probabilities, since pseudo counts are not applicable
            seqdat.p = seqdat.w
            # CHANGE bgfreq FOR SHANNON LOGOS TO A FLAT DISTRIBUTION
            if logoType == 1: bgfreq = bgfreq_flat
            # Re-dimension 1-dim bgfreq
            if len(bgfreq.shape) == 1: bgfreq=np.array([bgfreq])
            # CALCULATE PSSM (WEIGHTS)
            seqdat.calc_weights(yLabel, logoType, bgfreq)
         else:
            # Decide on PSSM type
            p = np.power(2,np.array([weights.copy() for i in range(4)])/np.array([[2,1,2,1]]).T)*np.array([bgfreq,bgfreq,bgfreq_flat,bgfreq_flat])
            sump =  np.sum(p,1)
            pssmType = np.argmin((np.array(sump)-1)**2) # Finding the best sump == 1 match
            sys.stderr.write('pssmType: %s\n'%['Halfbits','Bits','Halfbits_flat','Bits_flat'][pssmType])
            if pssmType == 0 or pssmType == 2: seqdat.w/=2           # Halfbits
            if pssmType <= 1: seqdat.p = np.power(2,seqdat.w)*bgfreq
            else:
               # Calculate Probabilities using flat bg distribution
               seqdat.p = np.power(2,seqdat.w)*bgfreq_flat
               Warning.RaiseWarning('flat_bgs_found')
               # Re-dimension 1-dim bgfreq
               if len(bgfreq.shape) == 1: bgfreq=np.array([bgfreq])
               # Recalculate weights using user-defined bgfreq
               seqdat.calc_weights(yLabel, logoType, bgfreq)
      else: # PSSM logo
         pass # nothing is needed done
   else:
      # MAKE SEGMENT
      if segment:
         seqdat.aln = seqdat.aln.T[segment[0]:segment[1]+1].T
         seqdat.seqlen = segment[1]-segment[0]+1
         startPos+=segment[0]
      # Set Number of Sequences
      seqdat.numseq = len(seqdat.aln)
      # CALCULATE SEQUENCE WEIGHTING
      if sequenceWeighting == 1: # Heuristics
         alpha, weights = Heuristic(seqdat.aln, alphabet, gaps)
      elif sequenceWeighting == 2: # Hobohm
         alpha, weights = Hobohm(seqdat.aln, alphabet, gaps, threshold)
      else: alpha, weights = seqdat.numseq - 1, [1]*seqdat.numseq # Set weights to 1
      # APPLY SEQUENCE WEIGHTING
      seqdat.set_sw(weights)
      # CALCULATE FREQUENCIES
      seqdat.calc_freq()
      # CALCULATE PROBABILITIES
      seqdat.calc_prob(alpha, beta)
      # CHANGE bgfreq FOR SHANNON LOGOS TO A FLAT DISTRIBUTION
      if logoType == 1: bgfreq = [np.ones(len(alphabet), dtype='f8')/len(alphabet)]
      # CALCULATE PSSM (WEIGHTS)
      seqdat.calc_weights(yLabel, logoType, bgfreq)
   
   # ----------------------------- PREPARE LOGO DATA ------------------------------
   # CALCULATE HEIGHTS
   seqdat.set_heights(logoType)
   # SET POSITIONS
   seqdat.set_positions(startPos, exclude_zero)
   # CALCULATE AND SET STACK WIDTH
   seqdat.set_width(inputType, minwidth)
   # CHECK IF ENDS MARK-UP IS REQUESTED AND SEQUENCE-TYPE IS DNA
   if eLabel and len(seqdat.alphabet)==4 and len(np.setdiff1d(list('ATGC'), seqdat.alphabet))==0:
      eLabel = 'd' # ENDS ARE SHOW AS DNA 3' and 5'
   
   # -------------------------------- CREATE LOGO ---------------------------------
   logo = MakeLogo(outputfile, seqdat.pos, title, pageSize, spl, lpp, userYaxis,
                   colors, xTic, seqdat.h, showFine, showX, showY, rotNum, xLabel,
                   yLabel, eLabel, eps_template, seqdat.width, seqdat.alphabet,
                   seqdat.consensus)
   
   # ----------------------------- SAVE WEIGHT MATRIX -----------------------------
   if not logoType==5: savePSSM(outputfile+".txt", seqdat.pos, seqdat.alphabet, seqdat.w, seqdat.p, sys.argv[1:], yLabel, logoType)
   # ----------------------------- SAVE FREQ MATRIX -------------------------------
   if not logoType==5: savePSSM(outputfile+"_freq.mat", seqdat.pos, seqdat.alphabet, seqdat.p, seqdat.p, sys.argv[1:], yLabel, 5, "#Position Specific Frequency Matrix\n")
   
   # ---------------------- CONVERT EPS TO REQUESTED FORMATS ----------------------
   EPS2Formats(outputfile, formats, logo.pW, logo.pH, gsPath)
   
   # ------------------------------ Present warnings ------------------------------
   return Warning.GetWarnings()

def GetArg(arg):
   ''' return a local variable. '''
   if arg in globals(): return globals()[arg]
   else: return None

if __name__ == '__main__':
   try: warnings = main()
   except Exception, e: raise #sys.exit(e)
   else: sys.stderr.write(warnings+'\n')
