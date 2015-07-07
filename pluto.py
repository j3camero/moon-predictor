from nbody import NBodySimulation

sim = NBodySimulation(G=6.67384e-11)
sim.AddRailPlanet(2.035e6, 551856.7, 0, 1.305e22, 'Pluto')
sim.AddRailPlanet(1.7536e7, 551856.7, 0.5, 1.587e21, 'Charon')
sim.AddPlanet(0, 4.2656e7, -153.8588, 0, 1.e16, 'Styx')
sim.AddPlanet(0, -4.8694e7, 142.4736, 0, 5.e16, 'Nix')
sim.AddPlanet(5.7783e7, 0, 0, 130.6316, 2.e16, 'Kerberos')
sim.AddPlanet(0, -6.4738e7, 123.2372, 0, 5.e16, 'Hydra')
sim.Normalize()
sim.Run(10 * 365.25 * 24 * 3600, 60, 'pluto.png', 400, 7.e7, 1)
