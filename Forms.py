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
def SphericalGaussianProduct():
    r, b, theta, phi = sp.symbols(r'r b \theta \phi', real=True, positive=True)

    var = [r, theta, phi]
    bounds = [[0, sp.oo], [0, sp.pi], [0, 2*sp.pi]]
    volumeElement = r**2 * sp.sin(theta)
    params = b
    form = sp.exp(-b * r**2) * r
    return form, var, bounds, volumeElement, params

# Spherical Exponential (radial but normalized over theta and phi)
def SphericalExponentialTrigThetaSin():
    r, b, theta, phi = sp.symbols(r'r b \theta \phi', real=True, positive=True)

    var = [r, theta, phi]
    bounds = [[0, sp.oo], [0, sp.pi], [0, 2*sp.pi]]
    volumeElement = r**2 * sp.sin(theta)
    params = b
    form = sp.exp(-b * r) * sp.sin(theta)
    return form, var, bounds, volumeElement, params

# Spherical Exponential (radial but normalized over theta and phi)
def SphericalExponentialTrigThetaCos():
    r, b, theta, phi = sp.symbols(r'r b \theta \phi', real=True, positive=True)

    var = [r, theta, phi]
    bounds = [[0, sp.oo], [0, sp.pi], [0, 2*sp.pi]]
    volumeElement = r**2 * sp.sin(theta)
    params = b
    form = sp.exp(-b * r) * sp.cos(theta)
    return form, var, bounds, volumeElement, params

# Spherical Exponential (radial but normalized over theta and phi)
def SphericalExponentialTrigPhi():
    r, b, theta, phi = sp.symbols(r'r b \theta \phi', real=True, positive=True)

    var = [r, theta, phi]
    bounds = [[0, sp.oo], [0, sp.pi], [0, 2*sp.pi]]
    volumeElement = r**2 * sp.sin(theta)
    params = b
    form = sp.exp(-b * r) * sp.exp(-sp.I * phi)
    return form, var, bounds, volumeElement, params

# Spherical Exponential (radial but normalized over theta and phi)
def SphericalLorentzian():
    r, b, theta, phi = sp.symbols(r'r b \theta \phi', real=True, positive=True)

    var = [r, theta, phi]
    bounds = [[0, sp.oo], [0, sp.pi], [0, 2*sp.pi]]
    volumeElement = r**2 * sp.sin(theta)
    params = b
    form = b**2 / (r**2 + 1/4 * b**2)
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
def SphericalGaussExponential():
    r, a, b, alpha, beta, theta, phi = sp.symbols(r'r a b \alpha \beta \theta \phi', real=True, positive=True)

    var = [r, theta, phi]
    bounds = [[0, sp.oo], [0, sp.pi], [0, 2*sp.pi]]
    volumeElement = r**2 * sp.sin(theta)
    params = [a, b, alpha, beta]
    form = beta * sp.exp(-b * r**2) + alpha * sp.exp(-a * r)
    return form, var, bounds, volumeElement, params

def CartesianGaussian():
    x, y, z = sp.symbols('x y z', real=True)
    b = sp.symbols('b', real=True, positive=True)

    var = [x, y, z]
    bounds = [[-sp.oo, sp.oo], [-sp.oo, sp.oo], [-sp.oo, sp.oo]]
    volumeElement = 1
    params = b
    form = sp.exp(-b * (x**2 + y**2 + z**2))
    return form, var, bounds, volumeElement, params


def CartesianLorentzian():
    x, y, z = sp.symbols('x y z', real=True)
    b = sp.symbols('b', real=True, positive=True)

    var = [x, y, z]
    bounds = [[-sp.oo, sp.oo], [-sp.oo, sp.oo], [-sp.oo, sp.oo]]
    volumeElement = 1
    params = b
    form = b**2/(x**2 + .25*b**2) * b**2/(y**2 + .25*b**2) * b**2/(z**2 + .25*b**2)
    return form, var, bounds, volumeElement, params
