﻿* Encoding: UTF-8.

RELIABILITY
  /VARIABLES=NOS1 NOS2 NOS3
  /SCALE('ALL VARIABLES') ALL
  /MODEL=ALPHA.

RELIABILITY
  /VARIABLES=SKEP1 SKEP2 SKEP3 SKEP4
  /SCALE('ALL VARIABLES') ALL
  /MODEL=ALPHA.

RELIABILITY
  /VARIABLES=SCN1 SCN2 SCN3 SCN4
  /SCALE('ALL VARIABLES') ALL
  /MODEL=ALPHA.

RELIABILITY
  /VARIABLES=SUPPORT1 SUPPORT2 SUPPORT3
  /SCALE('ALL VARIABLES') ALL
  /MODEL=ALPHA.

RELIABILITY
  /VARIABLES=PA1 PA2
  /SCALE('ALL VARIABLES') ALL
  /MODEL=ALPHA.

RELIABILITY
  /VARIABLES=NA1 NA2
  /SCALE('ALL VARIABLES') ALL
  /MODEL=ALPHA.

COMPUTE Nostalgia=(NOS1 + NOS2 + NOS3) / 3.
EXECUTE.

COMPUTE Skepticism=(SKEP1 + SKEP2 + SKEP3 + SKEP4) / 4.
EXECUTE.

COMPUTE Connectedness=(SCN1 + SCN2 + SCN3 + SCN4) / 4.
EXECUTE.

COMPUTE Support=(SUPPORT1 + SUPPORT2 + SUPPORT3) / 3.
EXECUTE.

COMPUTE Paffect=(PA1 + PA2) / 2.
EXECUTE.

COMPUTE Naffect=(NA1 + NA2) / 2.
EXECUTE.

COMPUTE Adoption=ADOPT1 + ADOPT2 + ADOPT3 + ADOPT4 + ADOPT5 + ADOPT6.
EXECUTE.

ONEWAY Nostalgia Skepticism Connectedness Support Paffect Naffect BY Mnos
  /STATISTICS DESCRIPTIVES 
  /MISSING ANALYSIS.

NPAR TESTS
  /K-S(NORMAL)=Adoption
  /K-S(POISSON 1)=Adoption
  /MISSING ANALYSIS
  /KS_SIM CIN(99) SAMPLES(10000).

GENLIN Adoption BY Mnos (ORDER=ASCENDING)
  /MODEL Mnos INTERCEPT=YES
 DISTRIBUTION=POISSON LINK=LOG
  /CRITERIA METHOD=FISHER(1) SCALE=1 COVB=MODEL MAXITERATIONS=100 MAXSTEPHALVING=5 
    PCONVERGE=1E-006(ABSOLUTE) SINGULAR=1E-012 ANALYSISTYPE=3(WALD) CILEVEL=95 CITYPE=WALD 
    LIKELIHOOD=FULL
  /MISSING CLASSMISSING=EXCLUDE
  /PRINT CPS DESCRIPTIVES MODELINFO FIT SUMMARY SOLUTION (EXPONENTIATED).


