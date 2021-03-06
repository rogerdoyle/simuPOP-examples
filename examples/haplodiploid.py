#!/usr/bin/env python

'''
Author: Bo Peng (bpeng@mdanderson.org)

Purpose:
  This script demonstrates how to handle recombination in haplodiploid populations.

Date: 2008-12-14
'''

from simuPOP import *

#
# This script demonstrate how to construct a genotype transmitter (during
# mating operator) that can be used to handle recombination in a
# haplodiploid population.
#
# Recall that, in a haplodiploid population
# - Female individuals have two sets of homologous chromosomes.
# - Male individuals have one set of homologous chromosomes.
# - Female offspring get one set of chromosome from her mother (possibly
#   with recombination), and one set of chromosome from her father (no
#   recombination).
# - Male offspring et on set of chromosome from his mother (possibly
#   with recombination).
#
# This is implemented in the following Python during-mating operator:
#
class HaplodiploidRecombinator(PyOperator):
    def __init__(self, rates, intensity, loci, convMode, *args, **kwargs):
        '''
        Create an instance of a Python operator, which will call
        ``self.transmitGenotype`` to create offspring. For performance
        considerations, this example uses two existing operators, namely
        ``Recombinator`` and ``GenoTransmitter`` to recombine and copy
        genotype. It is of course possible to use functions such as
        ``Individual::setGenotype()`` directly if no existing operator
        fits your need.
        '''
        # This operator is used to recombine maternal chromosomes
        self.Recombinator = Recombinator(rates, intensity, loci, convMode)
        # this operator is used to copy paternal chromosomes
        self.copier = GenoTransmitter()
        # With no *param*, this operator expects a function
        # in the form of ``(pop, off, dad, mom)``. If *param* is given, the
        # function should have the form ``(pop, off, dad, mom, param)``. If
        # *offspringOnly* is set to ``True``, the function can be simplied
        # to ``(off)`` or ``(off, param)``.
        PyOperator.__init__(self, func=self.transmitGenotype, *args, **kwargs)

    def transmitGenotype(self, pop, off, dad, mom):
        # Form the first homologous copy of offspring.
        self.Recombinator.transmitGenotype(mom, off, 0)
        # If the offspring is male, copy the second homologous copy from
        # her father. Male individuals only have one homologous set.
        if off.sex() == FEMALE:
            self.copier.copyChromosomes(dad, 0, off, 1)
        return True


def HaplodiploidRecMating(replacement=True, rates=[], intensity=-1, loci=[],
        convMode=NO_CONVERSION, numOffspring = 1., sexMode = RANDOM_SEX,
        ops = [], subPopSize = [], subPops = [], weight = 0):
    '''
    Return a mating scheme that uses random parents chooser and a customized
    during mating operator. A large number of parameters are provided to support
    number of offspring, sex-specification and population size changes. If none
    of these is needed, and you do not need a fancy Recombinator, this function
    can be simplied to:

    def haplodiploidRecMating(rate):
        return HomoMating(
            chooser = RandomParentsChooser(False),
            generator = OffspringGenerator(
                ops = [haplodiploidRecombinator(rate=rate)]
            )
        )
    '''
    return HomoMating(
        chooser = RandomParentsChooser(replacement),
        generator = OffspringGenerator(
            [HaplodiploidRecombinator(rates, intensity, loci, convMode)] + ops,
            numOffspring, sexMode),
        subPopSize = subPopSize,
        subPops = subPops,
        weight = weight)


def simuHaplodiploid(N, numMito=3, gen=10):
    '''
    The default genotype transmitter (HaplodiploidGenoTransmitter used in
    mating cheme haplodiploidMaiting) does not support recombination.
    We cannot use a Recombinator directly because it will also recombine
    maternal chromosomes. This example defines a Python during mating
    operator that actually uses a Recombinator to recombine maternal
    chromosomes, and then a GenoTransmitter to copy paternal chromosomes.
    '''
    pop = Population(N, ploidy=HAPLODIPLOID, loci=[20]*2,
        # record indexes of parents for verification purpose
        ancGen=1, infoFields=['father_idx', 'mother_idx'])

    pop.evolve(
        initOps=[
            InitSex(),
            # initialize alleles 0, 1, 2, 3 with different frequencies
            InitGenotype(freq=[0.4] + [0.2]*3),
        ],
        matingScheme = HaplodiploidRecMating(rates=0.1,
            # this operator is appended, instead of replace the haplodiploid
            # Recombinator
            ops = [ParentsTagger()]),
        postOps = Dumper(structure=False),
        gen = gen
    )
    return pop


if __name__ == '__main__':
    simuHaplodiploid(10, 3, 2)

# A possible output:
#
#SubPopulation 0 (unnamed):
#    0: MU 00120212010013211123 33003123030303201020 | ____________________ ____________________ |  9 3
#    1: FU 22233212010013211033 01003123030303201002 | 20002010320030301002 30030021133211001113 |  5 3
#    2: FU 00310232201103033000 00032103001203103020 | 00003100331001033002 30230200113331320030 |  4 1
#    3: MU 00310232201103032000 01231103010103303021 | ____________________ ____________________ |  8 1
#    4: FU 22233011331013211123 01003123030303201020 | 22010033023001203010 32132203031200003313 |  2 3
#    5: FU 01033012330020032010 00032120010103132120 | 22010033023001203010 32132203031200003313 |  2 1
#    6: MU 02233212010013211123 01003123030303201002 | ____________________ ____________________ |  6 3
#    7: FU 02220212031220111030 33120123030303201020 | 00003100331001033002 30230200113331320030 |  4 3
#    8: FU 00310212201120033000 01231320010103132120 | 00311103020110213030 30001033001002233012 |  7 1
#    9: MU 22233011310010111033 33120102130323320002 | ____________________ ____________________ |  0 3
# End of individual genotype.
# 
# Genotype of individuals in the present generation:
# SubPopulation 0 (unnamed):
#    0: MU 20033210320030301002 30033123030303001113 | ____________________ ____________________ |  0 1
#    1: FU 02220210331001033000 33120123030303200030 | 00310232201103032000 01231103010103303021 |  3 7
#    2: MU 00020212031001033002 30120123030303200030 | ____________________ ____________________ |  6 7
#    3: FU 00311103020110033000 01231323001002233122 | 02233212010013211123 01003123030303201002 |  6 8
#    4: FU 00310232201103033000 30030200113203120020 | 00120212010013211123 33003123030303201020 |  0 2
#    5: MU 00003100201103033000 30032103001203103020 | ____________________ ____________________ |  3 2
#    6: FU 00010232201103033002 00230200113333103020 | 22233011310010111033 33120102130323320002 |  9 2
#    7: MU 22010033023013203010 01003123030303201020 | ____________________ ____________________ |  0 4
#    8: MU 22233010320030301002 00033121133211001113 | ____________________ ____________________ |  3 1
#    9: FU 00310232201103033002 30232100113303103020 | 00310232201103032000 01231103010103303021 |  3 2
# End of individual genotype.

# 
# - Genotype of all male individuals are unused.
# - Individual 0 (male),
#   20033210320030301002 is a recombined copy of maternal chromosomes 
#       (22233212010013211033 20002010320030301002)
# - Individual 1 (Female),
#   The second homologous sets (00310232201103032000 01231103010103303021) are copied
#   from her father, without recombination.
#   
