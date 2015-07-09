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
from math import pi
import os
import random
import string

from PIL import Image, ImageDraw

from progress import ProgressBar

class Ray(object):
    """A position and velocity."""
    def __init__(self, x, y, vx, vy):
        self.x, self.y, self.vx, self.vy = x, y, vx, vy

class Planet(Ray):
    """A planet in an n-body simulation."""
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


def PolarToCartesian(r, angle):
    """Converts a radius and angle in radians to (x,y) coords."""
    return r * math.cos(angle), r * math.sin(angle)


class RailPlanet(Planet):
    """A special planet that follows a circular closed-form orbit."""
    def __init__(self, orbital_radius, orbital_period, phase, mass, name):
        self.orbital_radius = orbital_radius
        self.orbital_period = orbital_period
        self.orbital_speed = 2 * pi * orbital_radius / orbital_period
        self.phase = phase
        s = self.PositionAndVelocityAtTime(0)
        super(RailPlanet, self).__init__(s.x, s.y, s.vx, s.vy, mass, name)

    def PositionAndVelocityAtTime(self, t):
        orbit_count = float(t) / self.orbital_period + self.phase
        angle = orbit_count * 2 * pi
        x, y = PolarToCartesian(self.orbital_radius, angle)
        vx, vy = PolarToCartesian(self.orbital_speed, angle + 0.5 * pi)
        return Ray(x, y, vx, vy)

    def CalculateUpdate(self, sim, unused_dt):
        return self.PositionAndVelocityAtTime(sim.t)

    def ApplyUpdate(self, update, unused_dt):
        self.x = update.x
        self.y = update.y
        self.vx = update.vx
        self.vy = update.vy


class EllipticalOrbit(object):
    """A parameterization of an idealized elliptical two-body orbit."""
    def __init__(self, semi_major_axis, eccentricity):
        self.semi_major_axis = semi_major_axis
        self.eccentricity = eccentricity

    @staticmethod
    def CalculateFromMotion(ray, gm):
        r = math.sqrt(ray.x**2 + ray.y**2)
        v2 = ray.vx**2 + ray.vy**2
        # Calculate the semi-major axis using vis-visa equation.
        a = gm * r / (2 * gm - r * v2)
        # Angular momentum (not sure)
        h = ray.x * ray.vy - ray.y * ray.vx
        ecc = math.sqrt(1 - h**2 / (gm * a))
        return EllipticalOrbit(a, ecc)


def RandomWord(length):
    """Generate a random string of lowercase letters."""
    return ''.join(random.choice(string.lowercase) for i in range(length))


class NBodySimulation:
    """Represents an n-body simulation consisting of several Planets."""
    def __init__(self, G):
        self.G = G
        self.planets = []
        self.particles = []
        self.t = 0

    def AddPlanet(self, x, y, vx, vy, mass, name=''):
        """Shortcut method for adding a planet to the simulation."""
        self.planets.append(Planet(x, y, vx, vy, mass, name))

    def AddRailPlanet(self, orbital_radius, orbital_period, phase,
                      mass, name=''):
        """Shortcut method for adding a rail planet to the simulation."""
        self.planets.append(RailPlanet(orbital_radius, orbital_period, phase,
                                       mass, name))

    def AddParticle(self, orbital_radius, phase, speed_multiplier=1):
        """Add a low-mass particle to the simulation."""
        x, y = PolarToCartesian(orbital_radius, phase)
        orbital_speed = math.sqrt(self.G * self.TotalMass() / orbital_radius)
        orbital_speed *= speed_multiplier
        vx, vy = PolarToCartesian(orbital_speed, phase + 0.5 * pi)
        p = Planet(x, y, vx, vy, 1)
        p.color = (255, 255, 255)
        self.particles.append(p)

    def AddRandomParticle(self, lo, hi, speed_range=0):
        """Add a particle in a random circular orbit."""
        # Generate random radius using sqrt to get a uniform density.
        #orbital_radius = max_orbital_radius * math.sqrt(random.random())
        orbital_radius = random.random() * (hi - lo) + lo
        phase = random.random() * 2 * pi
        kick = random.uniform(1 - speed_range, 1 + speed_range)
        self.AddParticle(orbital_radius, phase, kick)

    def AddRandomParticles(self, num_particles, lo, hi, speed_range=0):
        """Add particles in random circular orbits."""
        for i in range(num_particles):
            self.AddRandomParticle(lo, hi, speed_range)

    def TotalMass(self):
        """The total mass of the planets in this simulation."""
        total = 0
        for p in self.planets:
            total += p.m
        return total

    def Tick(self, dt, deletion_distance):
        """Advance the simulation by one time step."""
        self.t += dt
        # Calculate and update particles in one pass since they have
        # negligible effect on larger planets.
        remaining_particles = []
        for p in self.particles:
            u = p.CalculateUpdate(self, dt)
            p.ApplyUpdate(u, dt)
            if p.x**2 + p.y**2 < deletion_distance**2:
                remaining_particles.append(p)
        self.particles = remaining_particles
        # Calculate and update planets in two passes so they're
        # order-independent.
        updates = []
        for p in self.planets:
            u = p.CalculateUpdate(self, dt)
            updates.append(u)
        for p, u in zip(self.planets, updates):
            p.ApplyUpdate(u, dt)

    def PlotParticles(self, image, image_size, plot_radius):
        """Plots the positions of the particles onto an existing PIL image."""
        for p in self.particles:
            x = int(0.5 * image_size * (p.x / plot_radius + 1))
            y = int(0.5 * image_size * (p.y / plot_radius + 1))
            if x < 0 or y < 0 or x >= image_size or y >= image_size:
                continue
            image.putpixel((x, y), p.color)

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
        non_resonators = ['Pluto', 'Charon', 'Jupiter']
        n = len(self.planets)
        for i in range(n):
            p1 = self.planets[i]
            if p1.name in non_resonators:
                continue
            a1, v1 = angular[i]
            for j in range(i + 1, n):
                p2 = self.planets[j]
                if p2.name in non_resonators:
                    continue
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

    def Run(self, max_t, dt, deletion_distance, image_filename=None,
            image_size=0, plot_radius=0, eta_report_frequency=0,
            oppconj_dir=None, snapshot_dir=None, snapshot_period=86400,
            particle_csv_filename=None):
        """Advance the simulation to the specified time max_t."""
        if image_filename or snapshot_dir:
            image = Image.new('RGB', (image_size, image_size), (0, 0, 0))
        if oppconj_dir:
            try:
                os.mkdir(oppconj_dir)
            except:
                pass
            index_filename = os.path.join(oppconj_dir, 'index.html')
            index_file = open(index_filename, 'w')
            index_file.write('<table>\n')
        if snapshot_dir:
            try:
                os.mkdir(snapshot_dir)
            except:
                pass
        next_snapshot = 0
        snapshot_count = 0
        progress = ProgressBar(eta_report_frequency)
        while self.t < max_t:
            self.Tick(dt, deletion_distance)
            if image_filename:
                self.PlotPlanets(image, image_size, plot_radius)
            if oppconj_dir:
                self.LogOppositionsAndConjunctions(oppconj_dir, index_file,
                                                   image, image_size,
                                                   plot_radius, dt)
            if snapshot_dir and self.t >= next_snapshot:
                snapshot_image = image.copy()
                self.DrawCirclePlanets(snapshot_image, image_size,
                                       plot_radius, 5)
                self.PlotParticles(snapshot_image, image_size, plot_radius)
                filename = os.path.join(snapshot_dir,
                                        '%06d.png' % snapshot_count)
                snapshot_image.save(filename)
                snapshot_count += 1
                next_snapshot += snapshot_period
            progress.MaybeReport(float(self.t) / max_t)
        if image_filename:
            self.PlotParticles(image, image_size, plot_radius)
            image.save(image_filename)
        if oppconj_dir:
            index_file.write('</table>\n')
            index_file.close()
        if particle_csv_filename:
            self.ParticleOrbitsToCsv(particle_csv_filename)

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

    def MoonSpread(self, names):
        """Maximum angle separating two moons from conjunction or opposition.
        """
        filtered_planets = [p for p in self.planets if p.name in names]
        n = len(filtered_planets)
        moon_spread = pi
        b = self.Barycenter()
        for i in range(n):
            p1 = filtered_planets[i]
            x1 = p1.x - b.x
            y1 = p1.y - b.y
            r1 = math.sqrt(x1*x1 + y1*y1)
            for j in range(i + 1, n):
                p2 = filtered_planets[j]
                x2 = p2.x - b.x
                y2 = p2.y - b.y
                r2 = math.sqrt(x2*x2 + y2*y2)
                cos_angle = (x1 * x2 + y1 * y2) / (r1 * r2)
                angle = math.acos(abs(cos_angle))
                moon_spread = min(moon_spread, angle)
        return moon_spread

    def CalculateParticleOrbits(self):
        """Calculates the idealized elliptical orbital parameters.

        This is not actually used for physics purposes. It's just for stats
        and analysis. Mainly, looking for clusters of objects with similar
        orbital characteristics.
        """
        gm = self.TotalMass() * self.G
        return [EllipticalOrbit.CalculateFromMotion(p, gm)
                for p in self.particles]

    def ParticleOrbitsToCsv(self, csv_filename):
        with open(csv_filename, 'w') as csv_file:
            csv_file.write('semi_major_axis,eccentricity\n')
            orbits = self.CalculateParticleOrbits()
            for orbit in orbits:
                if orbit.eccentricity >= 1:
                    continue
                line = '%f,%f\n' % (orbit.semi_major_axis, orbit.eccentricity)
                csv_file.write(line)
