import math
import numpy as np


# returns corresponding rotation matrix for SO2 and SO3
def so2_rotation(theta):
	return np.matrix([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])


def so3_rotation(x_theta, y_theta, z_theta):
	x = np.matrix([[1, 0, 0], [0, math.cos(x_theta), -math.sin(x_theta)], [0, math.sin(x_theta), math.cos(x_theta)]])
	y = np.matrix([[math.cos(y_theta), 0, -math.sin(y_theta)], [0, 1, 0], [math.sin(y_theta), 0, math.cos(y_theta)]])
	z = np.matrix([[math.cos(z_theta), -math.sin(z_theta), 0], [math.sin(z_theta), math.cos(z_theta), 0], [0, 0, 1]])
	return x * y * z


class SEN(object):
	translation = np.matrix('0; 0')                 # or np.matrix('0; 0; 0')               pos in x, y, z
	rotation = np.matrix('0 0; 0 0')                # or np.matrix('0 0 0; 0 0 0; 0 0 0')   mat in x, y, z
	se_type = 'SE2'                                 # or se_type = 'SE3'
	g_matrix = np.matrix('0 0 0; 0 0 0; 0 0 1')

	def __init__(self, translation, rotation, is_rotation_matrix=False):
		if np.shape(translation) == (2, 1):
			self.translation = translation
		elif np.shape(translation) == (3, 1):
			self.translation = translation
			self.se_type = 'SE3'
		else:
			print "translation matrix doesn't have valid dimensions"

		if not is_rotation_matrix:
			if np.shape(rotation) == (1, 1) and self.se_type == 'SE2':
				self.rotation = so2_rotation(rotation[0])
			elif np.shape(rotation) == (3, 1) and self.se_type == 'SE3':
				self.rotation = so3_rotation(rotation[0], rotation[1], rotation[2])
			else:
				print "inputted rotation values don't have valid dimensions"
		else:
			if np.shape(rotation) == (2, 2) and self.se_type == 'SE2':
				self.rotation = rotation
			elif np.shape(rotation) == (3, 3) and self.se_type == 'SE3':
				self.rotation = rotation
			else:
				print "rotation matrix don't have valid dimensions"

		# convert translations and rotations to resulting G matrix
		if self.se_type == 'SE2':
			rot = np.concatenate((self.rotation, np.matrix('0 0')))
			tra = np.concatenate((self.translation, np.matrix('1')))
			self.g_matrix = np.concatenate((rot, tra), axis=1)
		elif self.se_type == 'SE3':
			rot = np.concatenate((self.rotation, np.matrix('0 0 0')))
			tra = np.concatenate((self.translation, np.matrix('1')))
			self.g_matrix = np.concatenate((rot, tra), axis=1)

		# converting rotation matrices to quaternions
		# how do we take a Jacobian and keep it as a quaternion


def inverse(input_sen):
	return SEN(-np.transpose(input_sen.rotation) * input_sen.translation, np.transpose(input_sen.rotation),
	           is_rotation_matrix=True)

