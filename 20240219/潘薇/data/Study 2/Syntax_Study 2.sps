* Encoding: UTF-8.

DATASET ACTIVATE DataSet1.
RELIABILITY
  /VARIABLES=Nostalgia1 Nostalgia2 Nostalgia3 Nostalgia4 Nostalgia5 Nostalgia6 Nostalgia7
  /SCALE('ALL VARIABLES') ALL
  /MODEL=ALPHA.

RELIABILITY
  /VARIABLES=skepticism1 skepticism2 skepticism3 skepticism4
  /SCALE('ALL VARIABLES') ALL
  /MODEL=ALPHA.

RELIABILITY
  /VARIABLES=Connect1 Connect2 Connect3 Connect4
  /SCALE('ALL VARIABLES') ALL
  /MODEL=ALPHA.

RELIABILITY
  /VARIABLES=Support_AI1 Support_AI2 Support_AI3
  /SCALE('ALL VARIABLES') ALL
  /MODEL=ALPHA.

RELIABILITY
  /VARIABLES=Support_5G1 Support_5G2 Support_5G3
  /SCALE('ALL VARIABLES') ALL
  /MODEL=ALPHA.

COMPUTE Nostalgia=(Nostalgia1 + Nostalgia2 + Nostalgia3 + Nostalgia4 + Nostalgia5 + Nostalgia6 + 
    Nostalgia7) / 7.
EXECUTE.

COMPUTE Skepticism=(skepticism1 + skepticism2 + skepticism3 + skepticism4) / 4.
EXECUTE.

COMPUTE Connectedness=(Connect1 + Connect2 + Connect3 + Connect4) / 4.
EXECUTE.

COMPUTE Sup_AI=(Support_AI1 + Support_AI2 + Support_AI3) / 3.
EXECUTE.

COMPUTE Sup_5G=(Support_5G1 + Support_5G2 + Support_5G3) / 3.
EXECUTE.

CORRELATIONS
  /VARIABLES=Nostalgia Skepticism Connectedness Sup_AI Sup_5G
  /PRINT=TWOTAIL NOSIG
  /MISSING=PAIRWISE.


