import osprey
from ..settings import MUTATIONS, RAM, GPU, CPU, EPSILON, RESIDUES_SIDE1, RESIDUES_SIDE2


def run(pdb, interest_aa, neighbours_side1, neighbours_side2):
    osprey.start(heapSizeMiB=RAM)

    # choose a forcefield
    ffparams = osprey.ForcefieldParams()

    # read a PDB file for molecular info
    mol = osprey.readPdb(pdb)
    # make sure all strands share the same template library
    templateLib = osprey.TemplateLibrary(ffparams.forcefld)

    # define the protein strand
    protein = osprey.Strand(
        mol,
        templateLib=templateLib,
        residues=[RESIDUES_SIDE1.get('FIRST'), RESIDUES_SIDE1.get('LAST')]
    )
    protein.flexibility[interest_aa].setLibraryRotamers(osprey.WILD_TYPE,
                                                        *MUTATIONS).addWildTypeRotamers().setContinuous()

    map(lambda x: protein.flexibility[x].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous(),
        neighbours_side1)


    # define the ligand strand
    ligand = osprey.Strand(
        mol,
        templateLib=templateLib,
        residues=[RESIDUES_SIDE2.get('FIRST'), RESIDUES_SIDE2.get('LAST')])

    map(lambda x: ligand.flexibility[x].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous(),
        neighbours_side2)


    # make the conf space for the protein
    proteinConfSpace = osprey.ConfSpace(protein)

    # make the conf space for the ligand
    ligandConfSpace = osprey.ConfSpace(ligand)

    # make the conf space for the protein+ligand complex
    complexConfSpace = osprey.ConfSpace([protein, ligand])

    # how should we compute energies of molecules?
    # (give the complex conf space to the ecalc since it knows about all the templates and degrees of freedom)
    parallelism = osprey.Parallelism(cpuCores=CPU, gpus=GPU)
    minimizingEcalc = osprey.EnergyCalculator(complexConfSpace, ffparams, parallelism=parallelism, isMinimizing=True)

    # BBK* needs a rigid energy calculator too, for multi-sequence bounds on K*
    rigidEcalc = osprey.SharedEnergyCalculator(minimizingEcalc, isMinimizing=False)

    # how should we define energies of conformations?
    def confEcalcFactory(confSpace, ecalc):
        eref = osprey.ReferenceEnergies(confSpace, ecalc)
        return osprey.ConfEnergyCalculator(confSpace, ecalc, referenceEnergies=eref)

    # configure BBK*
    bbkstar = osprey.BBKStar(
        proteinConfSpace,
        ligandConfSpace,
        complexConfSpace,
        numBestSequences=2,
        epsilon=EPSILON,  # you proabably want something more precise in your real designs
        writeSequencesToConsole=True,
        writeSequencesToFile=f'bbkstar.results.{interest_aa}.tsv'
    )

    # configure BBK* inputs for each conf space
    for info in bbkstar.confSpaceInfos():
        # how should we define energies of conformations?
        eref = osprey.ReferenceEnergies(info.confSpace, minimizingEcalc)
        info.confEcalcMinimized = osprey.ConfEnergyCalculator(info.confSpace, minimizingEcalc, referenceEnergies=eref)

        # compute the energy matrix
        ematMinimized = osprey.EnergyMatrix(info.confEcalcMinimized, cacheFile='emat.%s.dat' % info.id)

        # how should confs be ordered and searched?
        # (since we're in a loop, need capture variables above by using defaulted arguments)
        def makeAStarMinimized(rcs, emat=ematMinimized):
            return osprey.AStarTraditional(emat, rcs, showProgress=False)

        info.confSearchFactoryMinimized = osprey.BBKStar.ConfSearchFactory(makeAStarMinimized)

        # BBK* needs rigid energies too
        confEcalcRigid = osprey.ConfEnergyCalculatorCopy(info.confEcalcMinimized, rigidEcalc)
        ematRigid = osprey.EnergyMatrix(confEcalcRigid, cacheFile='emat.%s.rigid.dat' % info.id)

        def makeAStarRigid(rcs, emat=ematRigid):
            return osprey.AStarTraditional(emat, rcs, showProgress=False)

        info.confSearchFactoryRigid = osprey.BBKStar.ConfSearchFactory(makeAStarRigid)

        # how should we score each sequence?
        # (since we're in a loop, need capture variables above by using defaulted arguments)
        def makePfunc(rcs, confEcalc=info.confEcalcMinimized, emat=ematMinimized):
            return osprey.PartitionFunction(
                confEcalc,
                osprey.AStarTraditional(emat, rcs, showProgress=False),
                osprey.AStarTraditional(emat, rcs, showProgress=False),
                rcs
            )

        info.pfuncFactory = osprey.KStar.PfuncFactory(makePfunc)

    # run BBK*
    scoredSequences = bbkstar.run(minimizingEcalc.tasks)

    # make a sequence analyzer to look at the results
    analyzer = osprey.SequenceAnalyzer(bbkstar)

    # use results
    for scoredSequence in scoredSequences:
        print("result:")
        print("\tsequence: %s" % scoredSequence.sequence)
        print("\tK* score: %s" % scoredSequence.score)

        # write the sequence ensemble, with up to 10 of the lowest-energy conformations
        numConfs = 10
        analysis = analyzer.analyze(scoredSequence.sequence, numConfs)
        print(analysis)
        analysis.writePdb(
            'seq.%s.pdb' % scoredSequence.sequence,
            'Top %d conformations for sequence %s' % (numConfs, scoredSequence.sequence)
        )
