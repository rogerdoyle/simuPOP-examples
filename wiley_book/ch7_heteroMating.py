import simuPOP as sim
pop = sim.Population(size=[10000, 10000], loci=1)
pop.setVirtualSplitter(sim.ProportionSplitter([0.8, 0.2]))
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(freq=[0.5, 0.5])
    ],
    preOps=[
        sim.Stat(homoFreq=0, subPops=[0,1], vars='homoFreq_sp'),
        sim.PyEval(r"'(%.2f, %.2f)\n' % (subPop[0]['homoFreq'][0], "
            "subPop[1]['homoFreq'][0])"),
    ],
    matingScheme=sim.HeteroMating(matingSchemes=[
        sim.RandomMating(subPops=[(0, 0), 1]),
        sim.SelfMating(subPops=[(0, 1)]),
    ]),
    gen = 3
)
