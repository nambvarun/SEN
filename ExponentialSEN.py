import math
import numpy as np

from SEN import SEN, so2_rotation, inverse


# change returns so that they are SENs rather than matrices
def hat(twist_coords, dt=1):
	if np.shape(twist_coords) == (3, 1):
		return np.matrix([[0, -twist_coords[2], twist_coords[0]],
		                  [twist_coords[2], 0, twist_coords[1]],
		                  [0, 0, 0]]) * dt
	elif np.shape(twist_coords) == (6, 1):
		return np.matrix([[0, -twist_coords[5], twist_coords[4], twist_coords[0]],
		                  [twist_coords[5], 0, -twist_coords[3], twist_coords[1]],
		                  [-twist_coords[4], twist_coords[3], 0, twist_coords[2]],
		                  [0, 0, 0, 0]]) * dt
	else:
		print "invalid twist coordinates"


def rotation_hat(twist_omegas, dt=1):
	if np.shape(twist_omegas) == (1, 1):
		return twist_omegas * np.matrix('0 -1; 1 0') * dt
	elif np.shape(twist_omegas) == (3, 1):
		return np.matrix([[0, -twist_omegas[2], twist_omegas[1]],
		                  [twist_omegas[2], 0, -twist_omegas[0]],
		                  [-twist_omegas[1], twist_omegas[0], 0]]) * dt
	else:
		print "invalid omega dimensions"


def unhat(mat, dt=1):
	if np.shape(mat) == (3, 3):
		return np.row_stack((mat[0:2, 2], mat[1, 0])) / dt
	elif np.shape(mat) == (4, 4):
		return np.row_stack((mat[0:3, 3], mat[2, 1], mat[0, 2], mat[1, 0])) / dt
	else:
		print "invalid input matrix for unhatting"


def rotation_unhat(r_mat, dt=1):
	if np.shape(r_mat) == (2, 2):
		return np.row_stack((r_mat[1, 0])) / dt
	elif np.shape(r_mat) == (3, 3):
		return np.row_stack((r_mat[2, 1], r_mat[0, 2], r_mat[1, 0])) / dt
	else:
		print "invalid input rotation matrix for unhatting"


class ExponentialSEN(object):
	twist_coords = np.matrix('0; 0; 0')
	sen = SEN(np.matrix('0; 0'), np.matrix('0 0; 0 0'), is_rotation_matrix=True)

	def __init__(self, twist_coords, dt=1):
		if np.shape(twist_coords) == (3, 1):
			self.twist_coords = twist_coords

			J = np.matrix('0 1; -1 0')

			# calculating the translation vector
			if twist_coords(2) != 0:
				T = -1 / twist_coords[1] * (np.identity(2) - so2_rotation(twist_coords[2])) * J * twist_coords[:1]
			else:
				T = twist_coords[:1] * dt

			# calculating the rotation matrix
			R = so2_rotation(twist_coords(2))
		elif np.shape(twist_coords) == (6, 1):
			self.twist_coords = twist_coords

			w = twist_coords[3:]
			v = twist_coords[:3]

			if np.linalg.norm(w) == 0:
				R = np.identity(3)
				T = v * dt
			else:
				R = np.identity(3) + rotation_hat(w) / np.linalg.norm(w) * math.sin(np.linalg.norm(w) * dt) + \
				         rotation_hat(w) ** 2 / np.linalg.norm(w) ** 2 * (1 - math.cos(np.linalg.norm(w) * dt))
				T = (np.identity(3) - R) * rotation_hat(w) * np.transpose(v) / np.linalg.norm(w) ** 2 + \
				         w * np.transpose(w) / np.linalg.norm(w) ** 2 * np.transpose(v) * dt
		else:
			print "invalid twist coordinates"

		self.sen = SEN(T, R, is_rotation_matrix=True)


def adjoint(sen, twist):
	if twist.isclass(SEN):
		return sen.g_matrix * twist * inverse(sen.g_matrix)
	elif sen.isclass(SEN) and sen.se_type == 'SE2':
		J = np.matrix('0, 1; -1, 0')

		if np.shape(twist)[0] == 3:
			return np.row_stack((np.column_stack((sen.rotation, sen.translation)), [0, 0, 1]))
		elif np.shape(twist)[0] == 2:
			return sen.rotation * twist
		else:
			print "twist has incorrect dimensions for SE2 adjoint operation"
	elif sen.isclass(SEN) and sen.se_type == 'SE3' and np.shape(twist) == (6, 1):
		quad1 = sen.rotation
		quad2 = rotation_hat(sen.translation) * sen.rotation
		quad3 = np.zeros((3, 3))
		final1 = np.concatenate((quad1, quad2), axis=1)
		final2 = np.concatenate((quad3, quad1), axis=1)
		return np.row_stack((final1, final2)) * twist
	else:
		print "incorrect adjoint param dimensions inputted"


def log(sen, dt=1):
	if sen.se_type == 'SE2':
		J = np.matrix('0 1; -1 0')
		g_mat = sen.g_matrix
		omega = np.arctan2(g_mat[1, 0], g_mat[0, 0]) / dt
		if omega != 0:
			v = omega * J * np.linalg.inv(np.identity(2) - sen.rotation) * sen.translation
		else:
			v = sen.translation / dt

		return np.row_stack((v, omega))
	elif sen.se_type == 'SE3':
		wn = (1 / dt) * np.arccos((np.trace(sen.rotation) - 1) / 2)
		g_mat = sen.g_matrix

		if wn == 0:
			omega = np.matrix('0; 0; 0')
		else:
			omega = wn / (2 * np.sin(wn * dt)) * np.matrix([[g_mat[2, 1] - g_mat[1, 2]],
			                                                [g_mat[0, 2] - g_mat[2, 0]],
			                                                [g_mat[1, 0] - g_mat[0, 1]]])

		if np.count_nonzero(omega):
			v = ((wn ** 2) * np.linalg.inv((np.identity(3) - sen.rotation) * rotation_hat(omega) + dt *
			                               (omega * np.transpose(omega)))) * sen.translation
		else:
			v = sen.translation / dt

		return np.row_stack(v, omega)
	else:
		print "invalid SEN type inputted into log function"
