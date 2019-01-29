
# coding: utf-8

# In[ ]:


import numpy as np


# In[ ]:


#using the potential described in https://arxiv.org/pdf/1611.00222.pdf section 3.2, identical to that used in their
#recent Nature paper

#defining the constants in pc
a = 6.5 * 10 **(3)
b = .26 * 10 **(3)
c = 0.7 * 10 **(3)
d = 12 * 10 **(3)

#potential for each portion of the galaxy, where r and z are from your data
halo = (173.2)**2 * np.log(1 + (r**2 / d**2) + (z**2 / d**2))
disk = - (4.302  *  10**(-3) * 6.3 * 10**(10)) / (np.sqrt(r**2 + (a + np.sqrt(z**2 + b**2))**2))
bulge = - (4.302 * 10**(-3) * 2.1 * 10**(10))/ (np.sqrt(r**2 + z**2) + c)

KE = .5 * (V**2 + U**2 + W**2)

Energy = -(halo + disk + bulge + KE)

L_z = r * V_phi

