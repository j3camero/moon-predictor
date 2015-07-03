from nbody import NBodySimulation

sim = NBodySimulation(G=6.674e-11)
sim.AddPlanet(0, 0, 0, 0, 5.972e24) # Earth
sim.AddPlanet(3.844e8, 0, 0, 1023, 7.348e22) # Moon
sim.AddPlanet(0, 1.496e11, -30000, 0, 1.989e30) # Sun
sim.Normalize()
sim_time = 365 * 24 * 3600
sim.Run(sim_time, 100, 'orbits.png', 2000, 1.6e11)
