def gamma2LK(gamma):
  """Converts Luttinger (Gamma) parameters to Luttinger-Kohn parameters.""" 

  L = -gamma[0] - 4 * gamma[1]
  M = 2 * gamma[1] - gamma[0]
  Nm = M - 1
  N = -6 * gamma[2]
  Np = N - Nm

  # Alternative, but equivalent formulation (by Foreman)
  #gdelta = 1 / 9. * (1 + gamma[0] + gamma[1] - 3 * gamma[2])
  #gmu = 1 / 2. * (gamma[2] - gamma[1])
  #gbar = 1 / 2. * (gamma[2] + gamma[1])
  #gsigma = gbar - 1 / 2. * gdelta
  #gpi = gmu + 3 / 2. * gdelta

  #L = -(gamma[0] + 4 * gamma[1])
  #M = -(gamma[0] - 2 * gamma[1])
  #Np = -6 * (gsigma - gdelta)
  #Nm = -6 * gpi
  #N = Nm + Np
  
  return [L, M, N, Np, Nm]

def LK2gamma(lk_lmn):
  """Converts Luttinger (Gamma) parameters to Luttinger-Kohn parameters."""

  g1 = -1. / 3. * (lk_lmn[0] + 2 * lk_lmn[1])
  g2 = -1. / 6. * (lk_lmn[0] - lk_lmn[1])
  g3 = -1. / 6. * lk_lmn[2] 

  return (g1, g2, g3)

def gamma2dresselhaus(gamma):
  """Converts Luttinger (Gamma) parameters to Dresselhaus parameters.""" 

  L = -gamma[0] - 4 * gamma[1] - 1 
  M = 2 * gamma[1] - gamma[0] - 1
  Nm = M
  N = -6 * gamma[2]
  Np = N - Nm

  return [L, M, N, Np, Nm]

