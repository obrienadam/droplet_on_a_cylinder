from math import sin, cos, acos, sqrt, pi
from scipy.optimize import newton
import matplotlib.pyplot as plt

def compute_droplet_area(rs, rl, lc):
    return pi*rl**2 \
           - rs**2*acos((lc**2 + rs**2 - rl**2)/(2*lc*rs)) \
            - rl**2*acos((lc**2 + rl**2 - rs**2)/(2*lc*rl)) \
            + sqrt((-lc + rs + rl)*(lc + rs - rl)*(lc - rs + rl)*(lc + rs + rl))/2.

rs = 5e-3
rl0 = 5e-3
l0 = 5e-3
theta = 150.

## Begin
theta *= pi/180.

# Find the starting area. This is preserved
area = compute_droplet_area(rs, rl0, l0)

print "Original area:", pi*rl0**2, area

d = 2*rs

# The root of this function is the solution

dl = lambda lc: d*cos(theta) + sqrt(d**2*cos(theta)**2 - d**2 + 4*lc**2)
alpha = lambda lc: acos((d**2 + 4*lc**2 - dl(lc)**2)/(4*d*lc))
f = lambda lc: area - (alpha(lc) + theta)*(dl(lc)/2)**2 + alpha(lc)*(d/2)**2 - (dl(lc)*d/4)*sin(theta)

# Finde the spacing
lc = newton(f, l0)

print "Droplet center location: {}, radius = {}".format(lc, dl(lc)/2)
print "Check error:", f(lc)
print "Check error:", compute_droplet_area(rs, dl(lc)/2, lc)

cylinder = plt.Circle((0, 0), rs, color='r')
droplet = plt.Circle((0, lc), dl(lc)/2, color='b')

fig, ax = plt.subplots()

ax.add_artist(droplet)
ax.add_artist(cylinder)

plt.axis('equal')
plt.show()
