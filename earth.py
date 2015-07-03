from nbody import NBodySimulation

sim = NBodySimulation(G=6.674e-11)
sim.AddPlanet(0, 0, 0, 0, 5.972e24)
sim.AddPlanet(3.844e8, 0, 0, 1023, 7.348e22)
sim.Normalize()
sim_time = 30 * 24 * 3600
sim.Run(sim_time, 300, 'orbits.png', 600, 4.e8)
