#!/usr/bin/env python
# coding: utf-8

# In[1]:
#Import packages
import msprime
import numpy as np
import matplotlib.pyplot as plt


# In[2]:
#Set up demographic history
demography = msprime.Demography()
taxa = {  # popsize, age-in-generations-ago
    "root" : (1000, 500000),
    "ABCD" : (1000, 400000),
    "BCD" : (1000, 300000),
    "par_ari" : (1000, 200000),
    "CD" : (1000, 100000),
    "D" : (1000, 90000),
    "cal_lon" : (1000, 80000),
    "pun" : (1000, 70000),
    "mainland" : (1000, 70000),
    "SD" : (1000, 60000),
    "cle" : (1000, 0),
    "gra" : (1000, 0),
    "ari" : (1000, 0),
    "par" : (1000, 0),
    "aur" : (1000, 0),
    "island_lon" : (1000, 0),
    "mainland_lon" : (1000, 0),
    "cal" : (1000, 0),
    "OC" : (1000, 0),
    "red" : (1000, 0),
    "yellow" : (1000, 0),
}
assert taxa["par_ari"][1] > taxa["CD"][1], "par-ari split must be older than everything-else split"

for n in taxa:
    demography.add_population(name=n, initial_size=taxa[n][0])
    
demography.add_population_split(time=taxa["root"][1], derived=["cle", "ABCD"], ancestral="root")
demography.add_population_split(time=taxa["ABCD"][1], derived=["gra", "BCD"], ancestral="ABCD")
demography.add_population_split(time=taxa["BCD"][1], derived=["CD", "par_ari"], ancestral="BCD")   
demography.add_population_split(time=taxa["par_ari"][1], derived=["par", "ari"], ancestral="par_ari")    
demography.add_population_split(time=taxa["CD"][1], derived=["aur", "D"], ancestral="CD") 
demography.add_population_split(time=taxa["D"][1], derived=["cal_lon", "pun"], ancestral="D") 
demography.add_population_split(time=taxa["cal_lon"][1], derived=["island_lon", "mainland"], ancestral="cal_lon")
demography.add_population_split(time=taxa["mainland"][1], derived=["cal", "mainland_lon"], ancestral="mainland")
demography.add_population_split(time=taxa["pun"][1], derived=["SD", "OC"], ancestral="pun")
demography.add_population_split(time=taxa["SD"][1], derived=["red", "yellow"], ancestral="SD")

# gene flow between island longiflorus starting at half the time back to the island-mainland split
# note that times are "time-ago", so first one applies to now; and second  one applies to longer-ago-than the given time
demography.set_symmetric_migration_rate(rate=0, populations=["island_lon", "par"])
#demography.add_symmetric_migration_rate_change(time=taxa["cal_lon"][1] / 2, rate=1e-2, populations=["island_lon", "par"])
#demography.add_symmetric_migration_rate_change(time=taxa["cal_lon"][1] / 2 + 1000, rate=0, populations=["island_lon", "par"])


# gene flow between parviflorus-ancestor and the branch above the MRCA of aur-cal-lon as long as both lineages exist
demography.add_symmetric_migration_rate_change(time=taxa["CD"][1], rate=0, populations=["par", "CD"])
demography.add_symmetric_migration_rate_change(time=150000, rate=0.1, populations=["par", "CD"])
demography.add_symmetric_migration_rate_change(time=150000+1, rate=0, populations=["par", "CD"])


demography.sort_events()

samples = [
    msprime.SampleSet(num_samples=6, population=n) for n in taxa if taxa[n][1] == 0
]

# In[3]:


dd = demography.debug()
# dd.print_history()


# In[4]:
#Run simulation
ts = msprime.sim_ancestry(samples, 
                          demography=demography,
                          sequence_length=5e8,
                          ploidy=2,
                          recombination_rate=1e-8
                         )
ts = msprime.sim_mutations(ts, rate=1e-8)
ts


# In[5]:
#Assign pop ids and sample names
pop_ids = { p.metadata['name'] : p.id for p in ts.populations() }
pop_samples = { n : ts.samples(population = pop_ids[n]) for n in pop_ids }


# TODO: check ordering of arguments to f4, at https://tskit.dev/tskit/docs/latest!

# In[6]:
#Calculate f4
windows = np.linspace(0, ts.sequence_length, 101)
f4_old = ts.f4(sample_sets=[
    pop_samples["par"],
    pop_samples["ari"],
    pop_samples["island_lon"],
    pop_samples["cle"]
],
     windows = windows)
f4_recent = ts.f4(sample_sets=[
    pop_samples["island_lon"],
    pop_samples["mainland_lon"],
    pop_samples["par"],
    pop_samples["cle"]
],
     windows = windows)


# In[7]:
#plot f4
fig, ax = plt.subplots()
ax.plot(windows[1:], f4_old, label="old")
ax.plot(windows[1:], f4_recent, label="recent")
ax.set_xlabel("genome")
ax.set_ylabel("f4")
ax.legend();


# In[8]:
#Print mean f4
print(f"mean recent f4: {np.mean(f4_recent)}")
print(f"mean old + recent f4: {np.mean(f4_old)}")
print(f"mean old f4: {np.mean(f4_old)-np.mean(f4_recent)}")


# In[9]:
#Output VCF
indiv_names = [f"{ts.population(ts.node(ind.nodes[0]).population).metadata['name']}_{ind.id}" for ind in ts.individuals()]
with open("sim_1_gens_gene_flow_0.1_only_ancient_x10_divergence_time.vcf", "w") as f:
    ts.write_vcf(f, individual_names=indiv_names)


# In[ ]:




