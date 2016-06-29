import math
import numpy as np

from SEN import SEN, so2_rotation


def hat(twist_coords, dt=1):
	if np.shape(twist_coords) == (3, 1):
		# rot =
		return np.matrix([[0, -twist_coords[2], twist_coords[0]], []])
	elif np.shape(twist_coords) == (6, 1):
		return np.matrix([[]])


# def expSE2(twist_coords, dt):
# 	J = np.matrix([[0, 1], [-1, 0]])
#
# 	if parameterized_twist[2] != 0:
# 		d = - 1/xi(3)*(eye(2)- so2_rotation(xi(3)))*J*xi(1:2);
#
# 	return 0
#
# def expSE3(parameterized_twist, dt):
# 	return 0

class ExponentialSEN(object):
	def __init__(self, twist_coords, dt=1):
		if np.shape(twist_coords) == (3, 1):
			print "asdf"
		elif np.shape(twist_coords) == (6, 1):
			print "asdf"
		else:
			print "invalid parameterized twist dimensions"
