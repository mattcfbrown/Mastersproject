import unittest
import numpy as np

from my_functions import max_prior
from my_functions import search
from my_functions import index
from my_functions import posterior

#Testing it actuall gets the maximum value
class TestSum(unittest.TestCase):
    def test_max_prior(self):
        #Testing if we are getting the maximum prior
        genes = ['T2','T1']
        prior_info = np.array([0,0.1,0.2,0.3,0,0.4,0.6,0.5,0]).reshape(3,3)
        result = max_prior(genes,prior_info)
        self.assertEqual(result,0.5*2.2)

#Testing it will search properly
class TestSum(unittest.TestCase):
    def test_search(self):
        #Testing if the search algorithm works (It should)
        list1 = ['T1----T2']
        list2 = ['T0----T1']
        gene_list = np.array(['T0', 'T1', 'T1', 'T2', 'T0', 'T2']).reshape(3,2)
        result = search(gene_list,list1,list2)
        self.assertEqual(result,[2,1,3])

#Testing if we get the agreed index
class TestSum(unittest.TestCase):
    def test_index(self):
        connection_type = [1,2,3,2,1,3,1,1]
        x = 1
        result = index(connection_type,x)
        self.assertEqual(result,[0,4,6,7])

#Now testing that the posterior calculation works
class TestSum(unittest.TestCase):
    def test_posterior(self):
        priors = [0.1,0.5]
        data = np.array([0.2,0.3,0.6,0.4,1.2,0.15]).reshape(2,3)
        index = [0,2]
        result = posterior(priors,data,index)
        print(result)
        self.assertEqual(result,[0.95,-1])        

if __name__ == '__main__':
    unittest.main()