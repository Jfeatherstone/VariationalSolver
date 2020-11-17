import sympy as sp

# Radial exponential form
def RadialExponential():
    r, b = sp.symbols('r b', real=True, positive=True)

    var = r
    bounds = [0, sp.oo]
    volumeElement = r**2
    params = b
    form = sp.exp(-b * r)
    return form, var, bounds, volumeElement, params

# Spherical Exponential (radial but normalized over theta and phi)
def SphericalExponential():
    r, b, theta, phi = sp.symbols(r'r b \theta \phi', real=True, positive=True)

    var = [r, theta, phi]
    bounds = [[0, sp.oo], [0, sp.pi], [0, 2*sp.pi]]
    volumeElement = r**2 * sp.sin(theta)
    params = b
    form = sp.exp(-b * r)
    return form, var, bounds, volumeElement, params

# Radial gaussian form
def RadialGaussian():
    r, b = sp.symbols('r b', real=True, positive=True)

    var = r
    bounds = [0, sp.oo]
    volumeElement = r**2
    params = b
    form = sp.exp(-b * r**2)
    return form, var, bounds, volumeElement, params

# Spherical Exponential (radial but normalized over theta and phi)
def SphericalGaussian():
    r, b, theta, phi = sp.symbols(r'r b \theta \phi', real=True, positive=True)

    var = [r, theta, phi]
    bounds = [[0, sp.oo], [0, sp.pi], [0, 2*sp.pi]]
    volumeElement = r**2 * sp.sin(theta)
    params = b
    form = sp.exp(-b * r**2)
    return form, var, bounds, volumeElement, params

# Spherical Exponential (radial but normalized over theta and phi)
def CubedExponential():
    r, b, theta, phi = sp.symbols(r'r b \theta \phi', real=True, positive=True)

    var = [r, theta, phi]
    bounds = [[0, sp.oo], [0, sp.pi], [0, 2*sp.pi]]
    volumeElement = r**2 * sp.sin(theta)
    params = b
    form = sp.exp(- b * r**3)
    return form, var, bounds, volumeElement, params
# Comination of spherical gaussian and exponential
#def SphericalGaussExponential():
#    r, a, b, theta, phi = sp.symbols(r'r a b \theta \phi', real=True, positive=True)
#
#    var = [r, theta, phi]
#    bounds = [[0, sp.oo], [0, sp.pi], [0, 2*sp.pi]]
#    volumeElement = r**2 * sp.sin(theta)
#    params = [a, b]
#    form = sp.exp(-b * r**2 - a * r)
#    return form, var, bounds, volumeElement, params
