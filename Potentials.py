import numpy as np


def Spin_V(Z_1, Z_2, mu_12, V_0, R, a, l, r):
  """A function which calculates a nuclear potential for a system of clusters with charges Z_1 and Z_" and reduced mass of mu_12, bound with orbital angular momentum of l.

  Args:
      Z_1 (_type_): Charge of cluster 1.
      Z_2 (_type_): Charge of cluster 2.
      mu_12 (_type_): Reduced mass of the two clusters.
      V_0 (_type_): Depth of the Wood-Saxon potential.
      R (_type_): Nuclear distance between the two clusters.
      a (_type_): Diffusness parameter of the Wood-Saxon potential.
      l (_type_): Orbital angular momentum quantum number
      r (_type_): Radial distance value, can be an array, in which case the output will return an np.array-s in which the index correspond to an r value with the same indes.

  Returns:
      _type_: Tuple consisting of values for [Wood-Saxon (0), Coulomb (1), Centrifugal barrier (2), Overall potential (3)]
  """
  #Constants
  h_bar = 197.326980 #MeV fm from NIST
  alpha = 7.2973525643*10**(-3) #fine structure constant from NIST
  e = (h_bar*alpha)**(0.5) #elementary charge

  #if r is a number:
  if isinstance(r, np.ndarray) != True:
    #Wood-Saxon potential:
    V_WS = -1.0*V_0/(1+np.exp((r-R)/a))
    #Centrifugal barrier:
    V_l = ((h_bar**2)/(2*mu_12))*l*(l+1)/(r**2)
    #Coulomb potential:
    if r < R:
      V_C = Z_1*Z_2*(e**2)*(1.5-(((r/R)**2)/2))*(1/R)
    elif r >= R:
      V_C = Z_1*Z_2*(e**2)*(1/r)
    else:
      print("Something went wrong :/...")

  elif isinstance(r, np.ndarray):
    V_WS = np.empty(0)
    V_l = np.empty(0)
    V_C = np.empty(0)
    for i in r:
      #Wood-Saxon potential:
      V_WS = np.append(V_WS,-1.0*V_0/(1+np.exp((i-R)/a)))
      #Centrifugal barrier:
      V_l = np.append(V_l,((h_bar**2)/(2*mu_12))*l*(l+1)/(i**2))
      #Coulomb potential:
      if i < R:
       V_C = np.append(V_C,Z_1*Z_2*(e**2)*(1.5-(((i/R)**2)/2))*(1/R))
      elif i >= R:
       V_C = np.append(V_C,Z_1*Z_2*(e**2)*(1/i))
      else:
       print("Something went wrong :/...")


  V = V_WS + V_C + V_l

  return V_WS, V_C, V_l, V


def square_well(V_0, L, r, V_h = 0):
  """Function giving a square well potential of depth V_0 and ranging from -L/2 to L/2. The height can be adjusted using V_h parameter.

  Args:
      V_0 (_type_): Depth of the potential.
      L (_type_): Width of the potential, cetered on r = 0.
      r (_type_): Distance variable, should be ranging from <-L/2 to >L/2.
      V_h (int, optional): Height of the potential. Defaults to 0.

  Returns:
      V _type_: Either an np.array containing the potential values corresponding to distance values or the potential value corresponding to a given distance value.
  """
  if isinstance(r, np.ndarray) != True:
    if r < -1*L/2:
     V = V_h
    elif -1*L/2 < r < L/2:
     V = V_0
    elif r>L/2:
     V=V_h
  elif isinstance(r, np.ndarray):
    V = np.empty(0)
    for i in r:
      if i < -1*L/2:
       V = np.append(V,V_h)
      elif -1*L/2 < i < L/2:
       V = np.append(V,V_0)
      elif i>L/2:
       V=np.append(V,V_h)
  return V

def WoodSaxon( V_0, R, a, r):
  """Calculates a Wood-Saxon potential with given parameters.

  Args:
      V_0 (_type_): Potential depth in MeV.
      R (_type_): The nuclear radius given by r_0(A1^(1/3)+A2^(1/3)), with A1 and A2 being the masses of the nuclei in amu.
      a (_type_): surface thickness of the nucleous
      r (_type_): radial distancem can be either an array or a single value.

  Returns:
     V_WS _type_: Values of the potentatial at the provided r, if r is an array returns an array of corresponding values.
  """
  #Constants
  h_bar = 197.326980 #MeV fm from NIST
  alpha = 7.2973525643*10**(-3) #fine structure constant from NIST
  e = (h_bar*alpha)**(0.5) #elementary charge

  #if r is a number:
  if isinstance(r, np.ndarray) != True:
    #Wood-Saxon potential:
    V_WS = -1.0*V_0/(1+np.exp((r-R)/a))

  elif isinstance(r, np.ndarray):
    V_WS = np.empty(0)
    for i in r:
      #Wood-Saxon potential:
      V_WS = np.append(V_WS,-1.0*V_0/(1+np.exp((i-R)/a)))
  return V_WS

def CentrifugalBarrier(mu_12, l, r):
  """Calculates the orbital angular momentum barrier.

  Args:
      mu_12 (_type_): Reduced mass of the two clusters in the composite nucleus in MeV/c.
      l (_type_): orbital angular momentum number between two clusters.
      r (_type_): radial distancem can be either an array or a single value.
  Returns:
      _type_: Values of the "potentatial" at the provided r, if r is an array returns an array of corresponding values.
  """
  #Constants
  h_bar = 197.326980 #MeV fm from NIST
  a = 7.2973525643*10**(-3) #fine structure constant from NIST
  e = (h_bar*a)**(0.5) #elementary charge

  #if r is a number:
  if isinstance(r, np.ndarray) != True:
    #Centrifugal barrier:
    V_l = ((h_bar**2)/(2*mu_12))*l*(l+1)/(r**2)

  elif isinstance(r, np.ndarray):
    V_l = np.empty(0)
    for i in r:
      #Centrifugal barrier:
      V_l = np.append(V_l,((h_bar**2)/(2*mu_12))*l*(l+1)/(i**2))
  return V_l

def Coulomb(Z_1, Z_2, R, r):
  """Calculates the Coulomb potential.

  Args:
      Z_1 (_type_): Charge of cluster 1.
      Z_2 (_type_): Charge of cluster 2.
      R (_type_): The nuclear radius given by r_0(A1^(1/3)+A2^(1/3)), with A1 and A2 being the masses of the nuclei in amu.
      r (_type_): radial distancem can be either an array or a single value.

  Returns:
      V_c _type_: Values of the potentatial at the provided r, if r is an array returns an array of corresponding values.
  """
  #Constants
  h_bar = 197.326980 #MeV fm from NIST
  alpha = 7.2973525643*10**(-3) #fine structure constant from NIST
  e = (h_bar*alpha)**(0.5) #elementary charge

  #if r is a number:
  if isinstance(r, np.ndarray) != True:
    #Coulomb potential:
    if r < R:
      V_C = Z_1*Z_2*(e**2)*(1.5-(((r/R)**2)/2))*(1/R)
    elif r >= R:
      V_C = Z_1*Z_2*(e**2)*(1/r)
    else:
      print("Something went wrong :/...")

  elif isinstance(r, np.ndarray):
    V_C = np.empty(0)
    for i in r:
      #Coulomb potential:
      if i < R:
       V_C = np.append(V_C,Z_1*Z_2*(e**2)*(1.5-(((i/R)**2)/2))*(1/R))
      elif i >= R:
       V_C = np.append(V_C,Z_1*Z_2*(e**2)*(1/i))
      else:
       print("Something went wrong :/...")
  return V_C

def SpinOrbit(V_sl0, j, l, s, R, a, r):
  """_summary_

  Args:
      V_sl0 (_type_): Potential depth in MeV.
      j (_type_): J value of the compound nucleus.
      l (_type_): Orbital angular momentum of the two clusters.
      s (_type_): Spin of one of the two clusters.
      R (_type_): The nuclear radius given by r_0(A1^(1/3)+A2^(1/3)), with A1 and A2 being the masses of the nuclei in amu.
      a (_type_): surface thickness of the nucleous
      r (_type_): radial distancem can be either an array or a single value.
  Returns:
      _type_: _description_
  """
  #Constants
  h_bar = 197.326980 #MeV fm from NIST
  alpha = 7.2973525643*10**(-3) #fine structure constant from NIST
  e = (h_bar*alpha)**(0.5) #elementary charge
  m_pi = 139.570 #pion mass in MeV/c

  ls = j*(j+1)-l*(l+1)-s*(s+1)

  #if r is a number:
  if isinstance(r, np.ndarray) != True:
    #Spin-orbit potential:
    V_sl = -1.0*ls*V_sl0*np.exp((r-R)/a) * ((h_bar)**2) / (2*((m_pi)**2)*r*a*((1+np.exp((r-R)/a))**2))

  elif isinstance(r, np.ndarray):
    V_sl = np.empty(0)
    for i in r:
      #Spin-orbit potential:
      V_sl = np.append(V_sl,-1.0*ls*V_sl0*np.exp((i-R)/a) * ((h_bar)**2) / (2*((m_pi)**2)*i*a*((1+np.exp((i-R)/a))**2)))

  return V_sl
  