from scipy.stats import norm
from scitools.easyviz.matplotlib_ import plot
from numpy import linspace

y = 0.12*norm.ppf(linspace(0.01, 0.99, 20)) + 2.269

plot(y)

raw_input("press enter")
