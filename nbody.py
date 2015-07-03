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

from PIL import Image

class Ray:
    """Class representing position and velocity."""
    def __init__(self, x, y, vx, vy):
        self.x, self.y, self.vx, self.vy = x, y, vx, vy

class Planet:
    """Class representing a planet in an n-body simulation."""
    def __init__(self, x, y, vx, vy, mass):
        self.state = Ray(x, y, vx, vy)
        self.m = mass
        # Generate a random bright color.
        color = [0, 255, random.randint(0, 255)]
        random.shuffle(color)
        self.color = tuple(color)

    def Acceleration(self, state, sim):
        """Calculate acceleration caused by other planets on this one."""
        ax = 0.0
        ay = 0.0
        for p in sim.planets:
            if p is self:
                continue
            dx = p.state.x - state.x
            dy = p.state.y - state.y
            dsq = dx*dx + dy*dy
            dr = math.sqrt(dsq)
            acceleration = sim.G*p.m/dsq if dsq>1e-10 else 0.
            ax += acceleration * dx / dr
            ay += acceleration * dy / dr
        return (ax, ay)

    def InitialDerivative(self, sim):
        """Part of Runge-Kutta method."""
        ax, ay = self.Acceleration(self.state, sim)
        return Ray(self.state.vx, self.state.vy, ax, ay)

    def NextDerivative(self, derivative, dt, sim):
        """Part of Runge-Kutta method."""
        state = Ray(self.state.x + derivative.x*dt,
                    self.state.y + derivative.y*dt,
                    self.state.vx + derivative.vx*dt,
                    self.state.vy + derivative.vy*dt)
        ax, ay = self.Acceleration(state, sim)
        return Ray(state.vx, state.vy, ax, ay)

    def Update(self, dt, sim):
        """Runge-Kutta 4th order solution to update planet's pos/vel."""
        a = self.InitialDerivative(sim)
        b = self.NextDerivative(a, dt*0.5, sim)
        c = self.NextDerivative(b, dt*0.5, sim)
        d = self.NextDerivative(c, dt, sim)
        dxdt = 1.0/6.0 * (a.x + 2.0*(b.x + c.x) + d.x)
        dydt = 1.0/6.0 * (a.y + 2.0*(b.y + c.y) + d.y)
        dvxdt = 1.0/6.0 * (a.vx + 2.0*(b.vx + c.vx) + d.vx)
        dvydt = 1.0/6.0 * (a.vy + 2.0*(b.vy + c.vy) + d.vy)
        self.state.x += dxdt*dt
        self.state.y += dydt*dt
        self.state.vx += dvxdt*dt
        self.state.vy += dvydt*dt


class NBodySimulation:
    """Represents an n-body simulation consisting of several Planets."""
    def __init__(self, G):
        self.G = G
        self.planets = []
        self.t = 0

    def AddPlanet(self, x, y, vx, vy, mass):
        """"Shortcut method for adding a planet to the simulation."""
        self.planets.append(Planet(x, y, vx, vy, mass))

    def Tick(self, dt):
        """Advance the simulation by one tick."""
        self.t += dt
        for p in self.planets:
            p.Update(dt, self)

    def Run(self, max_t, dt, image_filename=None, image_size=0, plot_radius=0):
        """Advance the simulation to the specified time max_t."""
        if image_filename:
            image = Image.new('RGB', (image_size, image_size), (0, 0, 0))
        while self.t < max_t:
            self.Tick(dt)
            if image_filename:
                for p in self.planets:
                    x = int(0.5 * image_size * (p.state.x / plot_radius + 1))
                    y = int(0.5 * image_size * (p.state.y / plot_radius + 1))
                    if x < 0 or y < 0 or x >= image_size or y >= image_size:
                        continue
                    image.putpixel((x, y), p.color)
        if image_filename:
            image.save(image_filename)

    def Normalize(self):
        """Centers the system at (0,0) and zeroes the average momentum.

        Only run this once at the beginning of the simulation. Running it
        periodically might cause really strange artifacts in the simulation.
        """
        barycenter = Ray(0, 0, 0, 0)
        total_mass = 0
        for p in self.planets:
            barycenter.x += p.state.x * p.m
            barycenter.y += p.state.y * p.m
            barycenter.vx += p.state.vx * p.m
            barycenter.vy += p.state.vy * p.m
            total_mass += p.m
        barycenter.x /= total_mass
        barycenter.y /= total_mass
        barycenter.vx /= total_mass
        barycenter.vy /= total_mass
        for p in self.planets:
            p.state.x -= barycenter.x
            p.state.y -= barycenter.y
            p.state.vx -= barycenter.vx
            p.state.vy -= barycenter.vy
