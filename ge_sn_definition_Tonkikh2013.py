# Definition of germanium material parameters (@ 300K)
Ge_basics = {'EG':0.8927, 'EL':0.7517, 'EX': 0.9467, 'Eav':0.00, 'delta_so':0.29, \
          'a0':0.56573, 'C11':128.53, 'C12':48.28, 'C44':0.0, \
	  'av':1.24, 'bv':-2.90, 'aG':-10.06, 'aL':-1.54, 'aX':0.14, 'xiX':9.42}

# From gamma parameters: 13.38, 4.24, 5.69 
Ge_ema = {'meffG': 0.038, 'meffL': (1.57,0.0807), 'meffX': (1.0, 1.0), 'meffhh': 0.2041, 'mefflh': 0.0457, 'meffso': 1.0}
Ge_kp6 = {'L':-30.34, 'M':-4.90, 'N':-34.14, 'Np':-28.24, 'Nm':-5.90}


# Definition of silicon material parameters (@ 300K)
Sn_basics = {'EG':0.5437, 'EL':0.9627, 'EX': 1.8667, 'Eav':0.69, 'delta_so':0.80, \
                 'a0':0.64892, 'C11':69.0, 'C12':29.3, 'C44':0.0, \
                 'av':1.58, 'bv':-2.70, 'aG':-6.00, 'aL':-2.14, 'aX':-0.46, 'xiX':9.42}

# From gamma parameters: -14.97, -10.61, -8.52
Sn_ema = {'meffG': 0.058, 'meffL': (1.478, 0.075), 'meffX': (1.0, 1.0), 'meffhh': 0.1600, 'mefflh': -0.0276, 'meffso': 1.0}
Sn_kp6 = {'L':57.41, 'M':-6.25, 'N': 51.12, 'Np':58.37, 'Nm':-7.25}

# Definition of bowing parameters for SiGe alloy
GeSn_bowing = {'EG': 2.61, 'EL': 0.80, 'EX': 0.10}

