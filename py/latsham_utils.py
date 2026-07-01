import numpy as np
from math import log, exp, sqrt, cos, sin, pi

def spherical_to_linear_parameters(hyperparams):
    if len(hyperparams) == 6:
        theta1a, theta1b, theta1c, theta2a, theta2b, theta2c = hyperparams
        return[cos(theta1a), sin(theta1a)*cos(theta1b), sin(theta1a)*sin(theta1b)*cos(theta1c), sin(theta1a)*sin(theta1b)*sin(theta1c),
                    cos(theta2a), sin(theta2a)*cos(theta2b), sin(theta2a)*sin(theta2b)*cos(theta2c), sin(theta2a)*sin(theta2b)*sin(theta2c)]
    if len(hyperparams) == 2:
        theta1, theta2 = hyperparams
        return [cos(theta1), sin(theta1), cos(theta2), sin(theta2)]
    else:        
        raise ValueError("Hyperparams must be either 2 or 6 elements long.")
    

def project_to_hypersphere(vec, radius=1.0):
    norm = sqrt(sum(v*v for v in vec))
    if norm == 0:
        raise ValueError("Cannot project the zero vector.")
    return [radius * v / norm for v in vec]


def linear_to_spherical_parameters(linear_params):
    if len(linear_params) == 8:
        x1, y1, z1, w1, x2, y2, z2, w2 = linear_params

        vec1 = project_to_hypersphere([x1, y1, z1, w1])
        vec2 = project_to_hypersphere([x2, y2, z2, w2])

        x1, y1, z1, w1 = vec1
        x2, y2, z2, w2 = vec2

        theta1a = np.arccos(x1)
        theta1b = np.arctan2(np.sqrt(z1*z1 + w1*w1), y1)
        theta1c = np.arctan2(w1, z1)

        theta2a = np.arccos(x2)
        theta2b = np.arctan2(np.sqrt(z2*z2 + w2*w2), y2)
        theta2c = np.arctan2(w2, z2)

        # Make them all positive and in the range [0, 2*pi) (some are really in [0, pi] and some are in [0, 2*pi])
        theta1a = theta1a % (2 * np.pi)
        theta1b = theta1b % (2 * np.pi)
        theta1c = theta1c % (2 * np.pi)
        theta2a = theta2a % (2 * np.pi)
        theta2b = theta2b % (2 * np.pi)
        theta2c = theta2c % (2 * np.pi)

        return [
            theta1a, theta1b, theta1c,
            theta2a, theta2b, theta2c
        ]
    if len(linear_params) == 4:
        x1, y1, x2, y2 = linear_params

        # Project the 2D vectors onto the unit circle
        vec1 = project_to_hypersphere([x1, y1])
        vec2 = project_to_hypersphere([x2, y2])

        x1, y1 = vec1
        x2, y2 = vec2

        theta1 = np.arccos(x1)
        theta2 = np.arccos(x2)

        return [theta1, theta2]
    else:
        raise ValueError("Linear params must be either 4 or 8 elements long.")
    

assert spherical_to_linear_parameters(linear_to_spherical_parameters([1, 0, 0, 0, 1, 0, 0, 0])) == [1, 0, 0, 0, 1, 0, 0, 0]


angles = np.random.rand(6)
angles[0] *= np.pi
angles[1] *= np.pi
angles[2] *= 2*np.pi
angles[3] *= np.pi
angles[4] *= np.pi
angles[5] *= 2*np.pi

lin = spherical_to_linear_parameters(angles)
angles2 = linear_to_spherical_parameters(lin)
lin2 = spherical_to_linear_parameters(angles2)

assert(np.max(np.abs(np.array(lin) - np.array(lin2))) < 1e-10)