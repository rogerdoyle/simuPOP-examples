import simuOpt
simuOpt.setOptions(quiet=True, alleleType='long')
import simuPOP as sim
pop = sim.Population(size=[2500]*10, loci=1)
simu = sim.Simulator(pop, rep=2)
simu.evolve(
    initOps=[
        sim.InitSex(),
        sim.InitGenotype(genotype=20),
    ],
    preOps=[
        sim.StepwiseMutator(rates=0.0001, reps=0),
        sim.KAlleleMutator(k=10000, rates=0.0001, reps=1),
    ],
    matingScheme=sim.RandomMating(),
    postOps=[
        # Use vars=['alleleFreq_sp'] to calculate allele frequency for
        # each subpopulation
        sim.Stat(alleleFreq=0, vars=['alleleFreq_sp'], step=200),
        sim.PyEval('gen', step=200, reps=0),
        sim.PyEval(r"'\t%.2f' % (sum([len(subPop[x]['alleleFreq'][0]) "
            "for x in range(10)])/10.)", step=200),
        sim.PyOutput('\n', reps=-1, step=200)
    ],
    gen=601
)
