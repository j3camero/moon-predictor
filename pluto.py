from nbody import NBodySimulation

sim = NBodySimulation(G=6.67384e-11)
sim.AddRailPlanet(2.035e6, 551856.7, 0, 1.305e22, 'Pluto')
sim.AddRailPlanet(1.7536e7, 551856.7, 0.5, 1.587e21, 'Charon')
sim.AddRailPlanet(4.2656e7, 1741957.92, 0.25, 1.e16, 'Styx')
sim.AddRailPlanet(4.8694e7, 2147440.032, 0.75, 5.e16, 'Nix')
sim.AddRailPlanet(5.7783e7, 2779277.184, 0.13, 2.e16, 'Kerberos')
sim.AddRailPlanet(6.4738e7, 3300632.928, 0.75, 5.e16, 'Hydra')
sim.AddRandomParticles(1000, 3.e7, 4.e7, 0.01)
sim.Normalize()
sim.Run(38 * 24 * 3600, 4 * 3600, 7.e7, 800, 7.e7, 'pluto_snapshots', 24*3600,
        'particles.csv')
