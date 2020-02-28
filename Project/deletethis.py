import math
import numpy as np 
import matplotlib.pyplot as plt
from ProtonBunch import ProtonBunch
from scipy.stats import norm

def paritcleEnergies():
    mean, var = 0.047, 0.1
    x = norm.pdf()



myBunch = ProtonBunch(0.047, 5)
for i in myBunch.bunch:
    print(i)
