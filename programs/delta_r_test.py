import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from my_functions import*

"""
Let's test how the delta_r function works. For this purpose we're going to invent a set of data. This data needs to have a certain
structure. The data must be arranged in pairs (etapairs and phipairs) like: [2,3] or like [[2,3],[2,1],[4,2]]. That means the 
function supports data containing just one pair or data containing arrays of pairs.
"""
test_eta=[[2,3],[2,1],[4,2]]
test_phi=[[1,2],[2,3],[-1,-2]]
result=delta_r(test_eta,test_phi)

print(r"The computed delta R values are:", result)