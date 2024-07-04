# Poised-Fragment-library-design
The design of a fragment library that is sociably poised for the XChem in house synthesis robot

## Aims

* New fragment library
* Updated focused in house database for elaboration

Requirements:

* Sociable (=elaboratable) compounds, with focus limited reaction repertoire
* Functionally diversity
* The building blocks for the fragment do _not_ need to be synthetically accessible

## Workflow

## Groundwork
What has come before

1. Analysis of historic data, expanding off [Carbery et al. 2022](https://pubs.acs.org/doi/10.1021/acs.jmedchem.2c01004).

## Validation
Given a library assess how much it fits our criteria

* Sociability
* Diversity

I experimented with Gold docking in the hope that there would be a modicum of validity,
but it opens more issues that it resolves.

## Design

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

## Limitations of DSiPoised

The DSiPoised is great for elaborations, but does not cover as much interaction diversity as
the 3D libraries for example.

The DSiPoised was made by shortlisting compounds by scoring for reaction moieties among various things. Four issues arise:
1.	These were not the smallest expansions, for example, a N-methylamide. A reaction moiety joining two big scaffolds means that it is not an expansion vector and instead introduces the problem that the two need to be assesses for which is the largest contributor.
2.	How the building blocks were made is irrelevant, instead it kills diversity
3.	Broad reaction set. Most in-house robotic synthesis combinatorial approaches rely on a very small subset of reactions (Suzuki, Schotten–Baumann and Buchwald–Hartwig), only for certain series are specific reactions employed. Cycloaddition products could certainly be helpful (e.g. Hantzsch, Diels–Alder, Büchi–Paterno cycloadditions), but generally the ring may need changing.
4.	Poor sociability. This is because the scaffold could be made via the reaction used in scoring, but this does not reflect the number of close analogues actually in catalogue space. The filtering criteria also could be said to have caused issues.
5.	Overly large. Large compounds are problematic as they open up the box of what is most relevant
6.	Tanimoto. The similarity filtering was done by Tanimoto similarity as a result there are many smaller compounds differ by one atom. Graph edit distance is by far a better metric (but very expensive).

## Steps

Enamine downloads on cluster see https://github.com/matteoferla/Fragment-hit-follow-up-chemistry
Enamine subsampled.
