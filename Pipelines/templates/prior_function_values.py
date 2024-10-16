#Here are some plotting tools to help my thesis

import math
import matplotlib.pyplot as plt
import numpy as np

#Write out the prior function:
def prior_fn(x,w0,multi):
    if x == 0.0:
        return math.exp(w0) / ( math.exp(w0) + math.exp(x) )
    elif 0.0*2.2 < x <0.1*2.2:
        return math.exp(w0) / ( math.exp(w0) + 0*math.exp(x) )
    elif 0.1*2.2 <= x < 0.5*2.2:
        return math.exp(w0) / ( math.exp(w0) + multi*0.5*math.exp(x) )        
    else:
       return  math.exp(w0) / ( math.exp(w0) + multi*math.exp(x) )  
    
x = np.linspace(0, 2.2, 300)
function_1 = []
function_10 = []
function_100 = []
for x_val in x:
    function_1.append(prior_fn(x_val,2.2,1.0))
    function_10.append(prior_fn(x_val,2.2,10.0))
    function_100.append(prior_fn(x_val,2.2,100.0))

#Plotting it all
plt.plot(x,function_1, label = "multi = 1")
plt.plot(x,function_10, label = "multi = 10")
plt.plot(x,function_100, label = "multi = 100")
plt.xlim([-0.001, 2.2])
plt.title('Prior function')
plt.ylabel('Alpha value')
plt.xlabel('x')
plt.legend()
plt.show()