import math
import numpy as np


class Quaternion(object):
	quaternion_components = [0, 0, 0, 0]
	rotation_matrix = np.identity(3)

	# input in the form q = w + xi + yj + zk
	# matrix should look like [w, x, y, z]
	def __init__(self, components):
		if len(components) == 4:
			self.quaternion_components = components
			w = self.quaternion_components(0)
			x = self.quaternion_components(1)
			y = self.quaternion_components(2)
			z = self.quaternion_components(3)

			n = w * w + x * x + y * y + z * z
			s = 0 if n == 0 else 2 / n

			wx = s * w * x
			wy = s * w * y
			wz = s * w * z

			xx = s * x * x
			xy = s * x * y
			xz = s * x * z

			yy = s * y * y
			yz = s * y * z
			zz = s * z * z

			self.rotation_matrix = np.matrix([1 - (yy + zz), xy - wz, xz + wy], [xy + wz, 1 - (xx + zz), yz - wx], [xz - wy, yz + wx, 1 - (xx + yy)])
		else:
			print "Inputted quaternion length is incorrect."

	# still not done implementing....
	def rotation_to_quaternion(self, components):
		w = components(0)
		x = components(1)
		y = components(2)
		z = components(3)

		n = 1.0 / math.sqrt(x * x + y * y + z * z + w * w)

		x *= n
		y *= n
		z *= n
		w *= n