#!/usr/bin/env python
"""Routines for running an n-body simulation.

This implementation uses the standard Runge-Kutta method for integrating
orbits. It's not state of the art, but it's a big improvement over the simpler
Euler method.

The original version of this file was created by Thanassis Tsiodras. Blog
post at:

    http://users.softlab.ntua.gr/~ttsiod/gravity.html
    http://ttsiodras.github.com/gravity.html
"""

import math
import random
import sys

from PIL import Image

# The window size
WIDTH, HEIGHT = 900, 600
WIDTHD2, HEIGHTD2 = WIDTH/2., HEIGHT/2.

# The gravity coefficient - it's my universe, I can pick whatever I want :-)
GRAVITYSTRENGTH = 1.e4

# The global list of planets
planets = []


class State:
    """Class representing position and velocity."""
    def __init__(self, x, y, vx, vy):
        self._x, self._y, self._vx, self._vy = x, y, vx, vy

    def __repr__(self):
        return 'x:{x} y:{y} vx:{vx} vy:{vy}'.format(
            x=self._x, y=self._y, vx=self._vx, vy=self._vy)


class Derivative:
    """Class representing velocity and acceleration."""
    def __init__(self, dx, dy, dvx, dvy):
        self._dx, self._dy, self._dvx, self._dvy = dx, dy, dvx, dvy

    def __repr__(self):
        return 'dx:{dx} dy:{dy} dvx:{dvx} dvy:{dvy}'.format(
            dx=self._dx, dy=self._dy, dvx=self._dvx, dvy=self._dvy)


class Planet:
    """Class representing a planet."""
    def __init__(self, x, y, vx, vy, mass):
        self._st = State(x, y, vx, vy)
        self._m = mass
        # Generate a random bright color.
        color = [0, 255, random.randint(0, 255)]
        random.shuffle(color)
        self._color = tuple(color)

    def __repr__(self):
        return repr(self._st)

    def acceleration(self, state, unused_t):
        """Calculate acceleration caused by other planets on this one."""
        ax = 0.0
        ay = 0.0
        for p in planets:
            if p is self:
                continue
            dx = p._st._x - state._x
            dy = p._st._y - state._y
            dsq = dx*dx + dy*dy
            dr = math.sqrt(dsq)
            force = GRAVITYSTRENGTH*self._m*p._m/dsq if dsq>1e-10 else 0.
            ax += force*dx/dr
            ay += force*dy/dr
        return (ax, ay)

    def initialDerivative(self, state, t):
        """Part of Runge-Kutta method."""
        ax, ay = self.acceleration(state, t)
        return Derivative(state._vx, state._vy, ax, ay)

    def nextDerivative(self, initialState, derivative, t, dt):
        """Part of Runge-Kutta method."""
        state = State(0., 0., 0., 0.)
        state._x = initialState._x + derivative._dx*dt
        state._y = initialState._y + derivative._dy*dt
        state._vx = initialState._vx + derivative._dvx*dt
        state._vy = initialState._vy + derivative._dvy*dt
        ax, ay = self.acceleration(state, t+dt)
        return Derivative(state._vx, state._vy, ax, ay)

    def updatePlanet(self, t, dt):
        """Runge-Kutta 4th order solution to update planet's pos/vel."""
        a = self.initialDerivative(self._st, t)
        b = self.nextDerivative(self._st, a, t, dt*0.5)
        c = self.nextDerivative(self._st, b, t, dt*0.5)
        d = self.nextDerivative(self._st, c, t, dt)
        dxdt = 1.0/6.0 * (a._dx + 2.0*(b._dx + c._dx) + d._dx)
        dydt = 1.0/6.0 * (a._dy + 2.0*(b._dy + c._dy) + d._dy)
        dvxdt = 1.0/6.0 * (a._dvx + 2.0*(b._dvx + c._dvx) + d._dvx)
        dvydt = 1.0/6.0 * (a._dvy + 2.0*(b._dvy + c._dvy) + d._dvy)
        self._st._x += dxdt*dt
        self._st._y += dydt*dt
        self._st._vx += dvxdt*dt
        self._st._vy += dvydt*dt


def main():
    image = Image.new('RGB', (WIDTH, HEIGHT), (0, 0, 0))

    global planets

    # And God said: Let there be lights in the firmament of the heavens...
    planets = []
    for i in xrange(0, 10):
        sr = 1.5
        planets.append(Planet(random.uniform(0, WIDTH),
                                      random.uniform(0, HEIGHT),
                                      random.uniform(-sr, sr),
                                      random.uniform(-sr, sr), 0.014137))

    sun = Planet(WIDTHD2, HEIGHTD2, 0, 0, 5)
    planets.append(sun)
    for p in planets:
        if p is sun:
            continue

    # t and dt are unused in this simulation, but are in general, 
    # parameters of engine (acceleration may depend on them)
    t, dt = 0., 1.

    for frame_number in range(1000):
        t += dt
        for p in planets:
            x = int(WIDTHD2+WIDTHD2*(p._st._x-WIDTHD2)/WIDTHD2)
            y = int(HEIGHTD2+HEIGHTD2*(p._st._y-HEIGHTD2)/HEIGHTD2)
            if x < 0 or y < 0 or x >= WIDTH or y >= HEIGHT:
                continue
            image.putpixel((x, y), p._color)

        # Update all planets' positions and speeds (should normally double
        # buffer the list of planet data, but turns out this is good enough :-)
        for p in planets:
            if p is sun:
                continue
            p.updatePlanet(t, dt)

    image.save('orbits.png')


if __name__ == "__main__":
    main()
