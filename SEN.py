import math
import numpy as np

class SEN(object):
	translation = np.matrix('0; 0')     # or np.matrix('0; 0; 0')           x, y, z
	rotation = np.matrix('0; 0')        # or np.matrix('0; 0; 0')   rads in x, y, z
	type = 'SE2'                        # or type = 'SE3'

	def __init__(self, translation, rotation):
		if np.shape(translation) == (2, 1):
			self.translation = translation
		elif np.shape(translation) == (3, 1):
			self.translation = translation
			self.type = 'SE3'
		else:
			print "translation matrix doesn't have valid dimensions"

		if np.shape(rotation) == (1, 1) and type == 'SE2':
			self.rotation = rotation
		elif np.shape(rotation) == (3, 1) and type == 'SE3':
			self.rotation = rotation
		else:
			print "rotation matrix doesn't have valid dimensions"

		# convert translations and rotations to resulting G matrix

	# returns corresponding rotation matrix for SE2 and SE3
	def so2_rotation(self, theta):
		return np.matrix([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])

	def so3_rotation(self, x_theta, y_theta, z_theta):
		x = np.matrix([[1, 0, 0], [0, math.cos(x_theta), -math.sin(x_theta)], [0, math.sin(x_theta), math.cos(x_theta)]])
		y = np.matrix([[math.cos(y_theta), 0, -math.sin(y_theta)], [0, 1, 0], [math.sin(y_theta), 0, math.cos(y_theta)]])
		z = np.matrix([[math.cos(z_theta), -math.sin(z_theta), 0], [math.sin(z_theta), math.cos(z_theta, 0)], [0, 0, 1]])
		return x * y * z

