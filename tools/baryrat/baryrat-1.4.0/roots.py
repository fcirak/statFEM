import sys
import numpy as np
from interfaceBrasil import InterfaceBrasil


m = int(sys.argv[1])                     # degree
exponent = float(sys.argv[2])            # exponent
interpolation = InterfaceBrasil(m)       # create and initalise object
interpolation.set_exponent(exponent)     # set the exponent of power function
p,q = interpolation.rational()           # top and bottom factors

np.savetxt('top.out', p)
np.savetxt('bot.out', q)

