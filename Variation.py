import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

# Normalize a given symbolic form by integrating over the relevant dimensions
def normalize(form, var, bounds, volumeElement=1):
    """
    Examples:

    Integrate an exponential form over the radial coordinate r
    >>> r, b = sp.symbols(r'r b', real=True, positive=True)
    >>> expForm = sp.exp(-b * r)
    >>> normExpForm = normalize(expForm, r, [0, sp.oo], volumeElement=r**2)

    >> 2*b**(3/2)*exp(-b*r)

    Integrate an exponential form over all three spherical coordinates
    >>> r, b, theta, phi = sp.symbols(r'r b \theta \phi', real=True, positive=True)
    >>> expForm = sp.exp(-b * r)
    >>> normExpForm = normalize(expForm, [r, theta, phi], [[0, sp.oo],
                                [0, sp.pi], [0, 2*sp.pi]], volumeElement=r**2 * sp.sin(theta))

    >> b**(3/2)*exp(-b*r)/sqrt(pi)
    """
    # Begin with the integrad as the usual psi psi* with volume element if provided
    integrand = form * sp.conjugate(form) * volumeElement

    # If we are passed a list of variables, we integrate over each one successively
    if isinstance(var, list):
        for i in range(len(var)):
            integrand = sp.integrate(integrand, (var[i], bounds[i][0], bounds[i][1]))
    # Otherwise we just integrate over the one
    else:
        integrand = sp.integrate(integrand, (var, bounds[0], bounds[1]))

    # Make sure our function is normalizable
    if integrand == 0 or abs(integrand) == sp.oo:
        raise Exception(f'Unnormalizable form passed to normalize: {form}')

    # Return the original form with normalization constant
    normalization = sp.sqrt(1/integrand)
    return normalization * form

# Parameterize a given Hamilontian in terms of a given form
# Similar argument signature to minimize, and should work very nicely
# with any of the return signatures of Forms.*
def parameterizeHamiltonian(hamiltonian, form, var, bounds, volumeElement):

    # Begin with the integrad as the usual psi psi* with volume element if provided
    integrand = form * hamiltonian(sp.conjugate(form), var) * volumeElement

    # If we are passed a list of variables, we integrate over each one successively
    if isinstance(var, list):
        for i in range(len(var)):
            integrand = sp.integrate(integrand, (var[i], bounds[i][0], bounds[i][1]))
    # Otherwise we just integrate over the one
    else:
        integrand = sp.integrate(integrand, (var, bounds[0], bounds[1]))

    return integrand

def minimizeHamiltonian(parameterizedHamiltonian, params, plot=False, initialGuess=None):
    if isinstance(params, list):
        # Define a function that we can evaluate to minimize
        def eval_H(x):
            subForm = dict(zip(params, x))
            return parameterizedHamiltonian.evalf(subs=subForm)

        # If we aren't provided a list of guess, just take 1
        if initialGuess == None:
            initialGuess = np.ones(len(params))

        # Now optimize numerically with scipy
        res = minimize(eval_H, initialGuess)

        finalSol = res.x

        return finalSol, eval_H(finalSol)
    else:
        # If it's just a single variable, we should be able to solve
        deriv = sp.diff(parameterizedHamiltonian, params)
        solutions = sp.solve(deriv, params)

        # Take the first positive value, since negative values most likely don't make sense
        finalSol = None
        if len(solutions) > 1:
            for s in solutions:
                if s > 0:
                    finalSol = s
                    break
        else:
            finalSol = solutions[0]

        if finalSol == None:
            raise Exception(f'No minimization of Hamiltonian {parameterizedHamiltonian} available')

        if plot:
            paramArr = np.linspace(float(finalSol*.5), float(finalSol*1.5), 100)
            energyArr = [parameterizedHamiltonian.evalf(subs={params:pi}) for pi in paramArr]
            optimalEnergy = parameterizedHamiltonian.evalf(subs={params:finalSol})

            plt.plot(paramArr, energyArr)
            plt.axvline(finalSol, linestyle='--', label=f'{params} = {finalSol:.5}', color='tab:orange')
            plt.axhline(y=optimalEnergy, linestyle='--', label=f'E = {optimalEnergy:.5}', color='tab:gray')
            plt.legend()
            plt.ylabel(r'$\langle \hat H \rangle$')

        return finalSol, parameterizedHamiltonian.evalf(subs={params:finalSol})

