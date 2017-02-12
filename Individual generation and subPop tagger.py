"""
Roger W. Doyle
Two small subpopulations go through an acute bottleneck and other changes
in census number. Pedigrees of all individuals are saved along with their generations, subpopulation
membership and genotypes. (gen_id and sp_id are generation and subpopulation, respectively.)
"""

from simuPOP import *
import simuPOP as sim
gen_evolve=20 #gens of evolution and retention

#sim.setRNG(seed=12481632) #useful for repeating a simulation run
seed_used = sim.getRNG().seed()
print '\nThe random seed used in the simulation was', seed_used, '\n'

def censuscontrol(gen):
    if gen>=10 and gen<13:
        return [10,20]
    elif gen>=13 and gen<15:
        return [20,40]
    elif gen==15:
        return [22,44]
    else:
        return [30,60]
pop = sim.Population(size=[8,16], ploidy=2, loci=[0,1,2], infoFields=['ind_id', 'father_id',
    'mother_id','gen_id','sp_id'],
    lociPos = (1, 2, 3), 
    ancGen=gen_evolve)  
pop.evolve(
    initOps=[
        sim.IdTagger(begin=0, end=-1),
	    sim.InitSex(maleProp=0.5),
        sim.InitGenotype(freq=[0.2, 0.2, 0.2, 0.2, 0.2], loci=[0,1,2]), 
        sim.PedigreeTagger(output='>>simp_Pedigree.ped', outputLoci=[0,1,2],
            outputFields=['gen_id', 'sp_id']) 
    ], #end of initOps
    preOps=[PyOperator(lambda pop: [pop.setIndInfo(x, "sp_id", x) for x in range(pop.numSubPop())] is not None),
    ],
    matingScheme=sim.MonogamousMating(subPopSize=censuscontrol,numOffspring=8,
        sexMode=(sim.NUM_OF_MALES, 2),                                      
    ops=[
        sim.InfoExec('gen_id = gen'),  
        sim.MendelianGenoTransmitter(),
        sim.IdTagger(),
        sim.InheritTagger(infoFields='sp_id'),
        sim.PedigreeTagger(output='>>simp_Pedigree.ped', outputLoci=[0,1,2],
            outputFields=['gen_id', 'sp_id']),
    ], #end of Ops
    ), #end of matingScheme   
    gen = gen_evolve,
) #end of pop.evolve
#Pedigree plus info sent to the monitor
print "ind_id, Sire, Dam, Sex, Affection, gen_id, sp_id, Loc1a,Loc1b, Loc2a,Loc2b, Loc3a, Loc3b,"  
print(open('simp_Pedigree.ped').read()) 


