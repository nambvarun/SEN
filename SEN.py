import math
import numpy as np


"""
Set of classes to handle operations with Special Euclidean 2 and 3 groups.
"""


class SE2(object):
	"""Class for a Special Euclidean Group 2 (SE2) matrix."""
	def __init__(self, translation=np.matrix('0; 0'), rotation=np.matrix('1 0; 0 1')):
		"""
		Initializer for a SE2 object.

		Parameters
		----------
		:param translation: np.matrix with shape (2, 1)
			The translation used in the SE2 object.
		:param rotation: np.matrix with shape (2, 2) or a float that is a rotation angle.
			The rotation used in the SE2 object.
		"""
		if np.prod(np.shape(rotation)) == 4:        # if i give a rotation matrix
			self._mat = np.vstack((rotation, np.matrix('0 0')))
			self._mat = np.hstack((self._mat, np.vstack((translation, np.matrix('1')))))
		else:                                       # if i just give a rotation angle
			rot_mat = np.matrix([[math.cos(rotation), -math.sin(rotation)], [math.sin(rotation), math.cos(rotation)]])
			self._mat = np.vstack((rot_mat, np.matrix('0 0')))
			self._mat = np.hstack((self._mat, np.vstack((translation, np.matrix('1')))))

	def __mul__(self, other):
		"""
		Allows you to use the (*) operator with SE2 objects.

		Parameter
		---------
		:param other: SE2 object or np.matrix with shape (2, n) or (3, n)
			If the np.matrix object has a shape of (2, n) a matrix of (1, n) ones are added to the bottom of the matrix.
			For a (3, n) input matrix, there isn't any modification.

		Return
		------
		:return: SE2 object or np.matrix with shape (2, n) or (3, n).
			For a SE2 object input, the function returns an SE2 object. If a np.matrix is passed as a parameter the
			function returns a np.matrix of the same size.
		"""
		if isinstance(other, SE2):
			output_mat = self._mat * other._mat
			return SE2(output_mat[0:2, 2], output_mat[0:2, 0:2])
		elif isinstance(other, np.matrix):
			if np.size(other, 0) == 2:
				other = np.vstack((other, np.ones((1, np.size(other, 1)))))
				out = self._mat * other
				return out[0:3, :]
			elif np.size(other, 0) == 3:
				return self._mat * other
			else:
				return None
		else:
			return None

	def adjoint(self, x):
		"""
		Allows you to perform the adjoint operation with the current SE2 object on vectors and SE2 objects.

		Parameter
		---------
		:param x: SE2 object or np.matrix with shape (2, n) or (3, n)
			If the np.matrix object has a shape of (2, n) a matrix of (1, n) ones are added to the bottom of the matrix.
			For a (3, n) input matrix, there isn't any modification.

		Return
		------
		:return: SE2 object or np.matrix with shape (2, n) or (3, n).
			For a SE2 object input, the function returns an SE2 object. If a np.matrix is passed as a parameter the
			function returns a np.matrix of the same size after applying the adjoint operation.
		"""
		if isinstance(x, SE2):
			return self * x * inv(self)
		elif isinstance(x, np.matrix):
			if np.size(x, 0) == 3:
				z = np.vstack((np.hstack((self.rotation(), np.matrix('0 1; -1 0') * self.translation())), np.matrix('0 0 1')))
				return z * x
			elif np.size(x, 0) == 2:
				x = np.vstack((x, np.zeros((1, x.shape[1]))))
				z = np.vstack((np.hstack((self.rotation(), np.matrix('0 1; -1 0') * self.translation())), np.matrix('0 0 1')))
				z *= x
				return z[0:2, :]
			else:
				return None
		else:
			return None

	def rotation(self):
		"""
		Get the SE2 object's rotation matrix.

		Return
		------
		:return: np.matrix with shape (2, 2)
			Returns the rotation component of the SE2 matrix.
		"""
		return self._mat[0:2, 0:2]

	def translation(self):
		"""
		Get the SE2 object's translation matrix.

		Return
		------
		:return: np.matrix with shape (2, 1)
			Returns the translation component of the SE2 matrix.
		"""
		return self._mat[0:2, 2]

	def matrix(self):
		"""
		Get the matrix representation of the SE2 matrix.

		Return
		------
		:return: np.matrix with shape (3, 3)
			Returns the internal representation of the matrix.
		"""
		return self._mat

	def __repr__(self):
		return np.round(self._mat, 4).__str__()


class SE3(object):
	"""Class for a Special Euclidean Group 3 (SE3) matrix."""
	def __init__(self, translation=np.matrix('0; 0; 0'), rotation=np.matrix('1 0 0; 0 1 0; 0 0 1')):
		"""
		Initializer for a SE3 object.

		Parameters
		----------
		:param translation: np.matrix with shape (3, 1)
			The translation used in the SE3 matrix.
		:param rotation: np.matrix with shape (3, 3) or with shape (1, 3).
			The rotation used in the SE3 matrix. If a matrix with shape (1, 3) is passed in the matrix should be encoded
			as follows [theta_x in radians, theta_y in radians, theta_z in radians].
		"""
		if np.prod(np.shape(rotation)) == 9:
			self._mat = np.vstack((rotation, np.matrix('0 0 0')))
			self._mat = np.hstack((self._mat, np.vstack((translation, np.matrix('1')))))
		else:
			self._mat = np.vstack((EulerXYZtoR(rotation[0, 0], rotation[0, 1], rotation[0, 2]), np.matrix('0 0 0')))
			self._mat = np.hstack((self._mat, np.vstack((translation, np.matrix('1')))))

	def __mul__(self, other):
		"""
		Allows you to use the (*) operator with SE3 objects.

		Parameter
		---------
		:param other: SE3 object or np.matrix with shape (3, n) or (4, n)
			If the np.matrix object has a shape of (3, n) a matrix of (1, n) ones are added to the bottom of the matrix.
			For a (4, n) input matrix, there isn't any modification.

		Return
		------
		:return: SE3 object or np.matrix with shape (3, n) or (4, n).
			For a SE3 object input, the function returns an SE3 object. If a np.matrix is passed as a parameter the
			function returns a np.matrix of the same size.
		"""
		if isinstance(other, SE3):
			output_mat = self._mat * other._mat
			return SE3(output_mat[0:3, 3], output_mat[0:3, 0:3])
		elif isinstance(other, np.matrix):
			if np.size(other, 0) == 3:
				other = np.vstack((other, np.ones((1, np.size(other, 1)))))
				out = self._mat * other
				return out[0:3, :]
			elif np.size(other, 0) == 4:
				return self._mat * other

	def adjoint(self, x):
		"""
		Allows you to perform the adjoint operation with the current SE3 object on vectors and SE3 objects.

		Parameter
		---------
		:param x: SE3 object or np.matrix with shape (3, n) or (4, n)
			If the np.matrix object has a shape of (3, n) a matrix of (1, n) ones are added to the bottom of the matrix.
			For a (4, n) input matrix, there isn't any modification.

		Return
		------
		:return: SE3 object or np.matrix with shape (3, n) or (4, n).
			For a SE3 object input, the function returns an SE3 object. If a np.matrix is passed as a parameter the
			function returns a np.matrix of the same size after applying the adjoint operation.
		"""
		if isinstance(x, SE3):
			return self * x * inv(self)
		elif isinstance(x, np.matrix):
			if np.size(x, 0) == 6:
				z = vec_hat(x) * self.rotation()
				z = np.hstack((self.rotation(), z))
				z = np.vstack((z, np.hstack((np.zeros((3, 3)), self.rotation()))))
				return z * x
			elif np.size(x, 0) == 3:
				x = np.vstack((x, np.zeros((3, x.shape[1]))))
				t = vec_hat(self.translation())
				z = np.hstack((self.rotation(), t * self.rotation()))
				z = np.vstack((z, np.hstack((np.zeros((3, 3)), self.rotation()))))
				z = z * x
				return z[0:3, :]
			elif np.shape(x) == (4, 4):
				return np.matrix('0')

	def rotation(self):
		"""
		Get the SE3 object's rotation matrix.

		Return
		------
		:return: np.matrix with shape (3, 3)
			Returns the rotation component of the SE3 matrix.
		"""
		return self._mat[0:3, 0:3]

	def translation(self):
		"""
		Get the SE3 object's translation matrix.

		Return
		------
		:return: np.matrix with shape (3, 1)
			Returns the translation component of the SE3 matrix.
		"""
		return self._mat[0:3, 3]

	def matrix(self):
		"""
		Get the matrix representation of the SE3 matrix.

		Return
		------
		:return: np.matrix with shape (4, 4)
			Returns the internal representation of the matrix.
		"""
		return self._mat

	def __repr__(self):
		return np.round(self._mat, 4).__str__()


def inv(sen):
	"""
	Get the inverse of a SE2 or SE3 object.

	Parameter
	---------
	:param sen: SE2 or SE3 object.

	Return
	------
	:return: SE2 or SE3 object.
		Returns the corresponding Special Euclidean inverse for a given input.
	"""
	if isinstance(sen, SE2):
		return SE2(-np.transpose(sen.rotation()) * sen.translation(), np.transpose(sen.rotation()))
	elif isinstance(sen, SE3):
		return SE3(-np.transpose(sen.rotation()) * sen.translation(), np.transpose(sen.rotation()))


def prod_range(sen_list, start_index, end_index):
	"""
	Returns the product of a list of special euclidean groups for a given range. The start and end indices are
	inclusive.

	Parameters
	----------
	:param sen_list: list of SE2s or SE3s
	:param start_index: integer
	:param end_index: integer

	Return
	------
	:return: A SE2 or SE3 object. Depends on the Special Euclidean objects in the list.
		The product of the elements of the list from start to end indices inclusive.
	"""
	if isinstance(sen_list[start_index], SE2):
		product = SE2()
	elif isinstance(sen_list[start_index], SE3):
		product = SE3()
	else:
		product = None

	for i in range(start_index, end_index + 1):
		product = product * sen_list[i]

	return product


def vec_hat(xi):
	"""
	Performs the hat operation on an input vector.

	Parameter
	---------
	:param xi: np.matrix with minimum dimension of (3, 1)

	Return
	------
	:return: np.matrix with dimension (3, 3)
	"""
	r1 = np.matrix([0, -xi[2, 0], xi[1, 0]])
	r2 = np.matrix([xi[2, 0], 0, -xi[0, 0]])
	r3 = np.matrix([-xi[1, 0], xi[0, 0], 0])
	return np.vstack((r1, r2, r3))


def rotX(theta):
	"""
	Returns a matrix representing a 3D x-axis rotation.

	Parameter
	---------
	:param theta: float
		A value in radians for the x-axis rotation.

	Return
	------
	:return: np.matrix with dimension (3, 3)
	"""
	return np.matrix([[1, 0, 0],
					 [0, math.cos(theta), -math.sin(theta)],
					 [0, math.sin(theta), math.cos(theta)]])


def rotY(theta):
	"""
	Returns a matrix representing a 3D y-axis rotation.

	Parameter
	---------
	:param theta: float
		A value in radians for the y-axis rotation.

	Return
	------
	:return: np.matrix with dimension (3, 3)
	"""
	return np.matrix([[math.cos(theta), 0, math.sin(theta)],
					 [0, 1, 0],
					 [-math.sin(theta), 0, math.cos(theta)]])


def rotZ(theta):
	"""
	Returns a matrix representing a 3D z-axis rotation.

	Parameter
	---------
	:param theta: float
		A value in radians for the z-axis rotation.

	Return
	------
	:return: np.matrix with dimension (3, 3)
	"""
	return np.matrix([[math.cos(theta), -math.sin(theta), 0],
					 [math.sin(theta), math.cos(theta), 0],
					 [0, 0, 1]])


def EulerXYZtoR(theta_x, theta_y, theta_z):
	"""
	Returns a 3D rotation matrix.

	Parameters
	----------
	:param theta_x: float
		A value in radians for the x-axis rotation.
	:param theta_y: float
		A value in radians for the y-axis rotation.
	:param theta_z: float
		A value in radians for the z-axis rotation.

	Return
	------
	:return: np.matrix with dimension (3, 3)
		The 3D rotation matrix as a np.matrix object.
	"""
	return rotX(theta_x) * rotY(theta_y) * rotZ(theta_z)

