import matplotlib.pyplot as plt
import numpy
import os
from paraview.simple import *
from subprocess import check_output

cwd_base = os.getcwd().rsplit(sep = os.sep, maxsplit = 2)
cwd = os.path.join(cwd_base[0], 'results','')
stdout = check_output('ls' +' {}'.format(cwd), shell=True)
files = stdout.decode().split(sep='\n')
reader = OpenDataFiles([cwd + a for a in files])
