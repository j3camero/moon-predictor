from gravityRK4 import NBodySimulation

sim = NBodySimulation(G=1.e4)
sim.AddPlanet(300, -40, 2, 0, 0.5)
sim.AddPlanet(300, 40, -2, 0, 0.5)
sim.AddPlanet(-300, -40, 2, 0, 0.5)
sim.AddPlanet(-300, 40, -2, 0, 0.5)
sim.Run(1000, 0.1, 'orbits.png', 600, 400)
