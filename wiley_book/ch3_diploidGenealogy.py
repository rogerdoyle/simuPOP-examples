import simuPOP as sim
pop = sim.Population(1000, ancGen=-1, 
    infoFields=['ind_id', 'father_id', 'mother_id'])
pop.evolve(
    initOps=[
        sim.InitSex(),
        sim.IdTagger(),
    ],
    matingScheme=sim.RandomMating(
        ops=[
            sim.IdTagger(),
            sim.PedigreeTagger()
        ],
    ),
    gen = 1000
)
# a pedigree with only paternal information
pop.asPedigree()
IDs = pop.identifyAncestors()
allIDs = [ind.ind_id for ind in pop.allIndividuals()]
removedIDs = list(set(allIDs) - set(IDs))
pop.removeIndividuals(IDs=removedIDs)
# number of ancestors...
sizes = [pop.popSize(ancGen=x) for x in range(pop.ancestralGens())]
print(sizes[0], sizes[100], sizes[500], sizes[999])
