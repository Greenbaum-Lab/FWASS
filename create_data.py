import msprime


def run_simulation(num_indv, num_snps):
    demography = msprime.Demography()
    demography.add_population(name='pop', initial_size=10000)

    ts = msprime.sim_ancestry(samples={'a': num_indv}, demography=demography, recombination_rate=1e-8, random_seed=1,
                              sequence_length=num_snps)
    mts = msprime.sim_mutations(ts, model=msprime.BinaryMutationModel(), rate=8e-7, random_seed=1)
    return mts


if __name__ == '__main__':
    mutations = run_simulation(20, 100)
    print(type(mutations))
    print(mutations)