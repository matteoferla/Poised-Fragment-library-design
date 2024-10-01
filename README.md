# Sociably poised fragment library design
The design of a fragment library that is sociably poised for the XChem in house synthesis robot

## Aims

* New fragment screening library
* Updated focused in-house database for elaboration

Requirements:

* Sociable (=elaboratable) compounds, with focus limited reaction repertoire
* Functionally diversity
* The building blocks for the fragment do _not_ need to be synthetically accessible

## Workflow

### Groundwork

> See [DSiPoised analysis](DSiPoised_analysis/DSiPoised_analysis.md)

What has come before: analysis of historic data, expanding off [Carbery et al. 2022](https://pubs.acs.org/doi/10.1021/acs.jmedchem.2c01004).

### Validation
Given a library assess how much it fits our criteria

* Sociability
* Diversity

I experimented with Gold docking in the hope that there would be a modicum of validity,
but it opens more issues that it resolves.

### Design (version 1)
This is no longer applicable as the updated focused in-house database is more urgent.


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

Example of minimum moiety (although this would be as super strict)

![minimum-2nd-amide.png](images/minimum-2nd-amide.png)

### User reference guide
Make a set of pages (html or pdf) that describes each compound, its sociability, its analogues and its presence in Chembl.

## DSiPoised analysis

> See [DSiPoised analysis](DSiPoised_analysis/DSiPoised_analysis.md)

## ChEMBL bioactive

> See [ChEMBL bioactive](ChEMBL_actives/ChEMBL_bioactive.md)

## Subsetting screening library

For Steph's fragment expansion algorithm, https://github.com/stephwills/FragmentKnitwork,
an Astex-inspired fragment network is needed.
Due to the size of catalogue space, a subset of the library is needed.

Initial notes are in [initial-details-on-subsetting.md](initial-details-on-subsetting.md).
These may be outdated, but may contain more info.

This was done in 3 steps.

The module `library_subsetting_module` contains the class `CompoundSieve`
to effectively give a verdict on a compound based on thresholds from a randomly drawn subset.
The cluster scripts are in `library_subsetting_cluster_scripts`.

For the sake of speed if a compound fails a check further checks are not performed.
To do a diagnostic, set the class property `cutoffs` to an empty dictionary.
The cutoffs are in the format min_ / max_ + value calculated, allowing easy subclassing
(as the method `assess` is runs on the key in `cutoffs` and `verdict` dictionary).


The first step (property `mode` set to `SieveMode.basic`) 
is a cutoff based on the properties of the CXSMILES row.
Removes by default the lower quartile of number of HBonds/HAC and Rotatable bonds/HAC.
This is because Enamine REAL is rich in alkanes and greasy compounds with no HBond donors.

The second step (property `mode` set to `SieveMode.substructure`)
* N rings ≥ 1
* HAC ≤ 35
* (HBD+HBA)/HAC ≥ 1 / 5 (~75% quantile)
* Rotatable-bonds/HAC ≤ 1 / 5 (~75% quantile)
* N methylene ≤ 6
* largest ring size ≤ 8
* no unwanted substructures:
  * protection groups,
  * PAINS and 
  * some exocyclic groups disliked by Medicinal Chemists (exocyclic carbamate, ester, imine, and hydrazines)
* 'boringness' ≤ 0 (~80% quantile)
* 'synthon_score_per_HAC' ≥ 0.138 (~75% quantile)

The third step (property `mode` set to `SieveMode.synthon`) generates a 3D conformer
and generates a weighted Zscore combining the above metrics, a 'get-out-of-jail' value for compounds
that are superstructures of library compounds, and a pharmacophore score based 
on uniqueness of pharmacophore distance trios (see pharmacophore).
This unified Zscore allows ranking.

## Other

* Enamine downloads on cluster see https://github.com/matteoferla/Fragment-hit-follow-up-chemistry
* Enamine subsampled with [enamine_subsample.py](library_subsetting_cluster_scripts/enamine_random_subsample.py)