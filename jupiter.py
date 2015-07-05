from nbody import NBodySimulation

sim = NBodySimulation(G=6.67384e-11)
sim.AddPlanet(0, 0, 0, 0, 1.89813e27, 'Jupiter')
sim.AddPlanet(4.218e8, 0, 0, 17330, 8.9298e22, 'Io')
sim.AddPlanet(-6.711e8, 0, 0, -13739, 4.7987e22, 'Europa')
sim.AddPlanet(1.0704e9, 0, 0, 10879, 1.4815e23, 'Ganymede')
sim.AddPlanet(0, 1.8827e9, -8203, 0, 1.0757e23, 'Callisto')
sim.Normalize()
for p in sim.planets:
    print p.name, p.Speed()
sim.Run(16 * 24 * 3600, 60, 'jupiter.png', 600, 2.e9, 5)
for p in sim.planets:
    print p.name, p.Speed()