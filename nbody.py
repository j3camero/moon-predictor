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
import os
import random
import string

from PIL import Image, ImageDraw

from progress import ProgressBar

pi = math.pi

class Ray(object):
    """Represents position and velocity."""
    def __init__(self, x, y, vx, vy):
        self.x, self.y, self.vx, self.vy = x, y, vx, vy

class Planet(Ray):
    """Represents a planet in an n-body simulation."""
    def __init__(self, x, y, vx, vy, mass, name=''):
        super(Planet, self).__init__(x, y, vx, vy)
        self.m = mass
        # Generate a random bright color.
        color = [0, 255, random.randint(0, 255)]
        random.shuffle(color)
        self.color = tuple(color)
        self.name = name

    def Acceleration(self, state, sim):
        """Calculate acceleration caused by other planets on this one."""
        ax = 0.0
        ay = 0.0
        for p in sim.planets:
            if p is self:
                continue
            dx = p.x - state.x
            dy = p.y - state.y
            dsq = dx*dx + dy*dy
            dr = math.sqrt(dsq)
            acceleration = sim.G*p.m/dsq if dsq>1e-10 else 0.
            ax += acceleration * dx / dr
            ay += acceleration * dy / dr
        return (ax, ay)

    def InitialDerivative(self, sim):
        """Part of Runge-Kutta method."""
        ax, ay = self.Acceleration(self, sim)
        return Ray(self.vx, self.vy, ax, ay)

    def NextDerivative(self, derivative, dt, sim):
        """Part of Runge-Kutta method."""
        state = Ray(self.x + derivative.x*dt,
                    self.y + derivative.y*dt,
                    self.vx + derivative.vx*dt,
                    self.vy + derivative.vy*dt)
        ax, ay = self.Acceleration(state, sim)
        return Ray(state.vx, state.vy, ax, ay)

    def CalculateUpdate(self, sim, dt):
        """Runge-Kutta 4th order solution to update planet's pos/vel."""
        a = self.InitialDerivative(sim)
        b = self.NextDerivative(a, dt*0.5, sim)
        c = self.NextDerivative(b, dt*0.5, sim)
        d = self.NextDerivative(c, dt, sim)
        dx = 1.0/6.0 * (a.x + 2.0*(b.x + c.x) + d.x)
        dy = 1.0/6.0 * (a.y + 2.0*(b.y + c.y) + d.y)
        dvx = 1.0/6.0 * (a.vx + 2.0*(b.vx + c.vx) + d.vx)
        dvy = 1.0/6.0 * (a.vy + 2.0*(b.vy + c.vy) + d.vy)
        return Ray(dx, dy, dvx, dvy)

    def ApplyUpdate(self, update, dt):
        self.x += update.x * dt
        self.y += update.y * dt
        self.vx += update.vx * dt
        self.vy += update.vy * dt

    def Speed(self):
        """Calculates the instantaneous speed of the planet."""
        return math.sqrt(self.vx**2 + self.vy**2)


def RandomWord(length):
   return ''.join(random.choice(string.lowercase) for i in range(length))


class NBodySimulation:
    """Represents an n-body simulation consisting of several Planets."""
    def __init__(self, G):
        self.G = G
        self.planets = []
        self.t = 0

    def AddPlanet(self, x, y, vx, vy, mass, name=''):
        """"Shortcut method for adding a planet to the simulation."""
        self.planets.append(Planet(x, y, vx, vy, mass, name))

    def Tick(self, dt):
        """Advance the simulation by one tick."""
        self.t += dt
        updates = []
        for p in self.planets:
            u = p.CalculateUpdate(self, dt)
            updates.append(u)
        for p, u in zip(self.planets, updates):
            p.ApplyUpdate(u, dt)

    def PlotPlanets(self, image, image_size, plot_radius):
        """Plots the positions of the planets onto an existing PIL image.

        Each planet is plotted by a single pixel. Call this function for
        many ticks re-using the same image and the planets leave traces,
        plotting their paths over time.
        """
        for p in self.planets:
            x = int(0.5 * image_size * (p.x / plot_radius + 1))
            y = int(0.5 * image_size * (p.y / plot_radius + 1))
            if x < 0 or y < 0 or x >= image_size or y >= image_size:
                continue
            image.putpixel((x, y), p.color)

    def DrawCirclePlanets(self, image, image_size, plot_radius, planet_radius):
        """Plots the positions of the planets onto an existing PIL image.

        The planets are drawn as large colored circles.
        """
        for p in self.planets:
            x = int(0.5 * image_size * (p.x / plot_radius + 1))
            y = int(0.5 * image_size * (p.y / plot_radius + 1))
            if x < 0 or y < 0 or x >= image_size or y >= image_size:
                continue
            draw = ImageDraw.Draw(image)
            r = planet_radius
            draw.ellipse((x-r, y-r, x+r, y+r), fill=p.color)

    def LogOppositionsAndConjunctions(self, log_dir, index_file, traces_image,
                                      image_size, plot_radius, dt):
        """Detect conjunctions and oppositions and write them to an HTML log.

        Also output handy images to help visualize each event.
        """
        b = self.Barycenter()
        angular = []
        # Calculate angle and angular speed for each planet w.r.t. barycenter.
        for p in self.planets:
            dx = p.x - b.x
            dy = p.y - b.y
            r2 = dx * dx + dy * dy
            dot = p.vx * (-dy) + p.vy * dx
            angle = math.atan2(dy, dx)
            angular_speed = dot / r2
            pair = (angle, angular_speed)
            angular.append(pair)
        # Compare each pair of planets to detect conjunction and/or opposition.
        n = len(self.planets)
        for i in range(n):
            p1 = self.planets[i]
            a1, v1 = angular[i]
            for j in range(i + 1, n):
                p2 = self.planets[j]
                a2, v2 = angular[j]
                angle_rate = v2 - v1
                angle_diff = a2 - a1
                if math.fabs(angle_rate) < 1.e-20:
                    continue
                # Detect conjunctions. Transform angle to the range [-pi,+pi].
                angle_diff -= 2 * pi * int((angle_diff + pi) / (2 * pi))
                conjunction = -angle_diff / angle_rate
                if conjunction <= 0 and conjunction >= -dt:
                    con_angle = a2 + conjunction * v2
                    con_time = self.t + conjunction
                    debug_image = traces_image.copy()
                    self.DrawCirclePlanets(debug_image, image_size,
                                           plot_radius, 10)
                    image_filename = RandomWord(6) + '.png'
                    debug_image.save(os.path.join(log_dir, image_filename))
                    html = '<tr><td>CON</td>'
                    html += '<td>%s</td><td>%s</td>' % (p1.name, p2.name)
                    html += '<td>%f</td><td>%f</td>' % (con_angle, con_time)
                    html += '<td><img src=%s></td>' % image_filename
                    html += '</tr>\n'
                    index_file.write(html)
                # Detect oppositions. Transform angle to the range [0,2*pi].
                angle_diff -= 2 * pi * int(angle_diff / (2 * pi))
                opposition = -(angle_diff - pi) / angle_rate
                if opposition <= 0 and opposition >= -dt:
                    opp_angle = a2 + opposition * v2
                    opp_time = self.t + opposition
                    debug_image = traces_image.copy()
                    self.DrawCirclePlanets(debug_image, image_size,
                                           plot_radius, 10)
                    image_filename = RandomWord(6) + '.png'
                    debug_image.save(os.path.join(log_dir, image_filename))
                    html = '<tr><td>OPP</td>'
                    html += '<td>%s</td><td>%s</td>' % (p1.name, p2.name)
                    html += '<td>%f</td><td>%f</td>' % (opp_angle, opp_time)
                    html += '<td><img src=%s></td>' % image_filename
                    html += '</tr>\n'
                    index_file.write(html)

    def Run(self, max_t, dt, image_filename=None, image_size=0, plot_radius=0,
            eta_report_frequency=0, oppconj_dir=None):
        """Advance the simulation to the specified time max_t."""
        if image_filename:
            image = Image.new('RGB', (image_size, image_size), (0, 0, 0))
        if oppconj_dir:
            try:
                os.mkdir(oppconj_dir)
            except:
                pass
            index_filename = os.path.join(oppconj_dir, 'index.html')
            index_file = open(index_filename, 'w')
            index_file.write('<table>\n')
        progress = ProgressBar(eta_report_frequency)
        while self.t < max_t:
            self.Tick(dt)
            if image_filename:
                self.PlotPlanets(image, image_size, plot_radius)
            if oppconj_dir:
                self.LogOppositionsAndConjunctions(oppconj_dir, index_file,
                                                   image, image_size,
                                                   plot_radius, dt)
            progress.MaybeReport(float(self.t) / max_t)
        if image_filename:
            image.save(image_filename)
        if oppconj_dir:
            index_file.write('</table>\n')
            index_file.close()

    def Barycenter(self):
        """Calculate the barycenter of the system, and also total momentum."""
        x, y, vx, vy, m = 0, 0, 0, 0, 0
        for p in self.planets:
            x += p.x * p.m
            y += p.y * p.m
            vx += p.vx * p.m
            vy += p.vy * p.m
            m += p.m
        s = 1.0 / m
        return Ray(x * s, y * s, vx * s, vy * s)

    def Normalize(self):
        """Centers the system at (0,0) and zeroes the average momentum.

        Only run this once at the beginning of the simulation. Running it
        periodically might cause really strange artifacts in the simulation.
        """
        b = self.Barycenter()
        for p in self.planets:
            p.x -= b.x
            p.y -= b.y
            p.vx -= b.vx
            p.vy -= b.vy
