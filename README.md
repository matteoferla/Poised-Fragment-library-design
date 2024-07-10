# Sociably poised fragment library design
The design of a fragment library that is sociably poised for the XChem in house synthesis robot

## Aims

* New fragment library
* Updated focused in house database for elaboration

Requirements:

* Sociable (=elaboratable) compounds, with focus limited reaction repertoire
* Functionally diversity
* The building blocks for the fragment do _not_ need to be synthetically accessible

## Workflow

### Groundwork
What has come before

1. Analysis of historic data, expanding off [Carbery et al. 2022](https://pubs.acs.org/doi/10.1021/acs.jmedchem.2c01004).

### Validation
Given a library assess how much it fits our criteria

* Sociability
* Diversity

I experimented with Gold docking in the hope that there would be a modicum of validity,
but it opens more issues that it resolves.

### Design

1. Starting from a collection of amidation/Suzuki–decomposed Murcko scaffolds
2. Count superstructures at each vector in the greater building blocking and screening collection
3. Count superstructures with reactive groups we care for
4. Assign a bioactive frequency akin to [Bueler and Raymond 2022](https://pubs.acs.org/doi/10.1021/acs.jcim.3c01096) and some other metrics
5. Perform boostrapped Pareto sampling with all the metrics
6. Ranking, Clustering by Tanimoto followed by subclustering by Graph Edit Distance
7. For each subcluster, enumerate compounds with a minimum moiety from the reaction, say a N-methylamide as opposed to the huge two-part fragments we have in the current DSiPosed
8. Cluster and pick best from each cluster
9. Test 
10. Iterate

### User reference guide
Make a set of pages (html or pdf) that describes each compound, its sociability, its analogues and its presence in Chembl.

## DSiPoised analysis

> See [DSiPoised analysis](DSiPoised_analysis.md)

## Misc

Enamine downloads on cluster see https://github.com/matteoferla/Fragment-hit-follow-up-chemistry
Enamine subsampled with [enamine_subsample.py](enamine_subsample.py)

Enamine building block scaffold vectors
Note that scaffolds were not fragmented by common reaction moieties, just Murcko decomposition.
No effort was made to address isomorphic compounds, so index order will affect the values.
Protection groups were not removed, so there are many Fmocs.

![top_BB_scaffold.png](images%2Ftop_BB_scaffold.png)
![second_top_BB_scaffold.png](images%2Fsecond_top_BB_scaffold.png)