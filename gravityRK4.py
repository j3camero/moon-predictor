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

    def acceleration(self, state, unused_t, other_planets, G):
        """Calculate acceleration caused by other planets on this one."""
        ax = 0.0
        ay = 0.0
        for p in other_planets:
            if p is self:
                continue
            dx = p._st._x - state._x
            dy = p._st._y - state._y
            dsq = dx*dx + dy*dy
            dr = math.sqrt(dsq)
            force = G*self._m*p._m/dsq if dsq>1e-10 else 0.
            ax += force*dx/dr
            ay += force*dy/dr
        return (ax, ay)

    def initialDerivative(self, state, t, other_planets, G):
        """Part of Runge-Kutta method."""
        ax, ay = self.acceleration(state, t, other_planets, G)
        return Derivative(state._vx, state._vy, ax, ay)

    def nextDerivative(self, initialState, derivative, t, dt,
                       other_planets, G):
        """Part of Runge-Kutta method."""
        state = State(0., 0., 0., 0.)
        state._x = initialState._x + derivative._dx*dt
        state._y = initialState._y + derivative._dy*dt
        state._vx = initialState._vx + derivative._dvx*dt
        state._vy = initialState._vy + derivative._dvy*dt
        ax, ay = self.acceleration(state, t+dt, other_planets, G)
        return Derivative(state._vx, state._vy, ax, ay)

    def update(self, t, dt, other_planets, G):
        """Runge-Kutta 4th order solution to update planet's pos/vel."""
        a = self.initialDerivative(self._st, t, other_planets, G)
        b = self.nextDerivative(self._st, a, t, dt*0.5, other_planets, G)
        c = self.nextDerivative(self._st, b, t, dt*0.5, other_planets, G)
        d = self.nextDerivative(self._st, c, t, dt, other_planets, G)
        dxdt = 1.0/6.0 * (a._dx + 2.0*(b._dx + c._dx) + d._dx)
        dydt = 1.0/6.0 * (a._dy + 2.0*(b._dy + c._dy) + d._dy)
        dvxdt = 1.0/6.0 * (a._dvx + 2.0*(b._dvx + c._dvx) + d._dvx)
        dvydt = 1.0/6.0 * (a._dvy + 2.0*(b._dvy + c._dvy) + d._dvy)
        self._st._x += dxdt*dt
        self._st._y += dydt*dt
        self._st._vx += dvxdt*dt
        self._st._vy += dvydt*dt


class NBodySimulation:
    """Represents an n-body simulation consisting of several Planets."""
    def __init__(self, G):
        self.G = G
        self.planets = []
        self.t = 0

    def tick(self, dt):
        self.t += dt
        for p in self.planets:
            p.update(self.t, dt, self.planets, self.G)

    def run(self, max_t, dt, image_filename, image_size, plot_radius):
        if image_filename:
            image = Image.new('RGB', (image_size, image_size), (0, 0, 0))
        while self.t < max_t:
            self.tick(dt)
            if image_filename:
                for p in self.planets:
                    x = int(0.5 * image_size * (p._st._x / plot_radius + 1))
                    y = int(0.5 * image_size * (p._st._y / plot_radius + 1))
                    if x < 0 or y < 0 or x >= image_size or y >= image_size:
                        continue
                    image.putpixel((x, y), p._color)
        if image_filename:
            image.save(image_filename)


def main():
    simulation = NBodySimulation(G=1.e4)
    simulation.planets.append(Planet(300, -40, 2, 0, 0.5))
    simulation.planets.append(Planet(300, 40, -2, 0, 0.5))
    simulation.planets.append(Planet(-300, -40, 2, 0, 0.5))
    simulation.planets.append(Planet(-300, 40, -2, 0, 0.5))
    simulation.run(1000, 0.1, 'orbits.png', 600, 400)


if __name__ == "__main__":
    main()
