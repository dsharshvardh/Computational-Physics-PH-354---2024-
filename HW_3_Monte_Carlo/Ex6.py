'''
Question 6(a)

The range of θ is 0 to pi, and ɸ is 0 to 2pi
Given p(θ)=sin(θ)/2. 
          integral of (p(θ)dθ) from the range of 0 to pi  
          = -cos(θ)/2 from 0 to pi
          = -cos(pi/2) + cos(0)
          = 1
Now p(ɸ)=1/2pi, 
          integral of (p(ɸ)dɸ) from range of 0 to 2pi
          = ɸ/2pi from range of 0 to 2pi
          = (2pi-0)/2pi
          = 1   
Therefore, both the individual distributions are normalised.

'''




# Question 6(b)

#   Consider a transformation from some uniform random variable x to θ. Now considering continuity, if x(θ)=t, then x(θ+dθ) is close to x and we can consider it as x+dx
#   So, the Jacobian of this transformation will be given as :-
#                           p(x)dx = p(θ)dθ for "probability conservation"
#                           p(x) = 1 for a uniform distribution from 0 to 1.
#                           dx = sin(θ)/2 dθ
#   Let x(0)=0. Then integration gives us x=(1-cos(θ))/2 which can be written as θ=arccos(1-2x)   [As math.acos gives value from 0 to pi]
#   Therefore, if we have a uniform random variable in 0 to 1, we can get θ and ɸ according to densities θ and ɸ by taking ɸ=2pi*x and θ=math.acos(1-2x)






# Question (c)

import math
import random

phi=2*math.pi*random.random()
theta=math.acos((1-2*random.random()))
print("phi=",phi)
print("theta=",theta)
  
