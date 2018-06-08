import math

def alloy_gradients_1d(profile_type, param):
  
  if profile_type == 'linear':
    class LinearGradient:
      def __init__(self, param):
        startpos = param[0]
        thickness = param[1]
        x_init = param[2]
        x_final = param[3]

        self.offset = x_init
        self.slope = (x_final - x_init) / thickness
        self.center = startpos

      def __call__(self, coord):
        return self.offset + self.slope * (coord[0] - self.center)

    return LinearGradient(param)

  elif profile_type == 'fermi':
    class FermiGradient:
      def __init__(self, param):
        startpos = param[0]
        thickness = param[1]
        x_init = param[2]
        x_final = param[3]
        slope = abs(param[4])

        self.scaling = x_final - x_init
        self.center = startpos + 0.5 * thickness
        self.offset = x_init
	self.slope = slope

      def __call__(self, coord):
        # Check if exponential cannot cause an overflow
	if -self.slope * (coord[0] - self.center) > 709:
	  fermi_value = self.offset
	else:
          fermi_value = self.offset + self.scaling / (math.exp(-self.slope * (coord[0] - self.center)) + 1)
	
	return fermi_value

    return FermiGradient(param)

  elif profile_type == 'expmodgauss':
    class ExpModGaussGradient:
      def __init__(self, param):
        startpos = param[0]
	y0 = param[1]
	A = param[2]
	xc = param[3]
	w = param[4]
	t0 = param[5]
        
	self.xs = startpos
	self.y0 = y0
	self.A = A
	self.xc = xc
	self.w = w
	self.t0 = t0

      def __call__(self, coord):
        x = coord[0] - self.xs
	z = (x - self.xc) / self.w - self.w / self.t0
	erfw = 0.5 * (1 + math.erf(z / math.sqrt(2)))

        expmodgauss_value = self.y0 + self.A / self.t0 * \
	  math.exp(0.5 * (self.w / self.t0) ** 2 - (x - self.xc) / self.t0) * erfw

	return expmodgauss_value / 100.0

    return ExpModGaussGradient(param)


  else:
    raise Exception("Alloy profile '" + profile_type + "' unknown.")

