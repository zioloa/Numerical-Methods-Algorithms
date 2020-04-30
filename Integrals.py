""" Adam Ziółkowski

	Simple program containing three methods for computing integrals (analytical, Romberg's method and Gaussian quadrature).
	Here it is computed for f(x) = 0.02721 x^4 + 0.5415 x^3 + 0.1211 x^2 + 4.21 x -2.21
"""
import numpy as np

x0 = -20
x1 = 5

coefficients = np.array([0.02721, 0.5415, 0.1211, 4.27, -2.21])
powers = np.arange(4, -1, -1)


def f(x):
    return np.sum(coefficients * pow(x, powers))


"""
    Analytical method
"""


def integral_f(x):
    return np.sum((coefficients / (powers + 1)) * pow(x, (powers + 1)))
# [ 0.005442    0.135375    0.04036667  2.135      -2.21      ]


def analytical():
    return integral_f(x1) - integral_f(x0)


analytical = analytical()

"""
    Romberg's method
"""


def trapeze(xp, xk):
    return (xk - xp) * ((f(xp) + f(xk)) / 2)


i = 1
oh = [[trapeze(x0, x1)]]
done = False
while not done:
    temp = []
    x = np.linspace(x0, x1, 2**i + 1)

    for j in range(2 ** i):
        temp.append(trapeze(x[j], x[j+1]))

    oh[0].append(np.sum(temp))

    for k, el in enumerate(oh):
        k += 1
        if len(el) >= 2:
            if len(oh) > k:
                oh[k].append((4**k * el[-1] - el[-2]) / (4**k - 1))
            else:
                oh.append([(4 ** k * el[-1] - el[-2]) / (4 ** k - 1)])

    if np.abs((analytical - oh[-1][0])/analytical) * 100 <= 2:
        done = True

    i += 1

"""
    Gaussian quadrature
"""


def gaussian_quadrature(t, c):
    temp = []
    for i in range(3):
        temp.append(((x1-x0)/2)*(c[i] * f(((x1 + x0)+(x1 - x0)*t[i])/2)))

    return np.sum(temp)


t = [-np.sqrt(3/5), 0, np.sqrt(3/5)]
c = [5/9, 8/9, 5/9]
quadrature = gaussian_quadrature(t, c)


"""
    Results
"""

print("-" * 70)
print("Analytical method: ", analytical)
print("Romberg's method: ",oh[-1][0])
print("Gaussian quadrature: ", quadrature)
print("-" * 70)
print("Result polynomial coefficients in analytical method:")
print(coefficients/(powers + 1))
print("-" * 70)
print("Romberg's intermediate matrix :\n")
for el in np.transpose(oh):
    print(el)
print("-" * 70)
