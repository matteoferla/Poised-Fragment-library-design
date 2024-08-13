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
What has come before

1. Analysis of historic data, expanding off [Carbery et al. 2022](https://pubs.acs.org/doi/10.1021/acs.jmedchem.2c01004).

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

Vendor catalogue space is too big for the in-house database.
There are two ways to subset the vendor catalogue space:

* **Projection**. Decompose and enumerate analogues of the compounds in the libraries, 
    and subset all catalogue compounds that include them
* **Trim**. Subset the catalogues by a set of rules, including
    a metric that approximates how 'amicable' the compounds are to each other and accessible
    to a restricted reaction repertoire.

Due to the urgency, the latter has to be done.

This requires the following steps:

* **Decomposition**. Decompose the compounds into synthons that can be combined via the robot's restricted reaction repertoire,
* **Amicability**. Assign a metric that approximates how 'amicable' the compounds are

The first is done by the class `RoboDecomposer`.
The reactions had to be written from scratch to avoid problems with lactams and sulfams.
This reverses the following bonds:

* _Amides_: fully exocyclic, or exocyclic carbonyls with alicylic nitrogen, secondary or tertiary amides,
    avoids ureido and lactam.
    This results in carboxylic acid and amine.
    For simplicity this lumps together amidation and Schotten-Baumann reaction (cf. tertiary amides).
    The Buckwald-Hartwig amination is not included, hence the lack of aryl nitrogens.
* _Sulfonamides_: fully exocyclic, or exocyclic carbonyls, secondary or tertiary sulfonamides.
    Avoids sulfams.
    This results in sulfonyl chloride and amine. This is sulfo Schotten-Baumann reaction product.
* _Biaryls_: two arenes connected by a single bond.
    This results in two aryl halides (for simplicity: one would be a boronic acid/ester).
    This is the Suzuki–Miyaura reaction product.

Other reactions, such as Buckwald-Hartwig amination, Sonogashira coupling, Buch reductive amination, Williamson ether synthesis, Ugi reaction, Chan-Lam coupling etc. are not included.
As a results the distribution of synthons heavy atom count is positively skewed.

![](images/HAC-synthon.png)

However, the synthon amicability score is higher for smaller synthons.

![amiHAC](images/amicability-HAC.png)

Then the 'synthon amicability' is calculated.
Due to time constraints and because I was instructed to use USRCAT,
I used USRCAT as a pharmacophoric similarity metric.
There is an issue in that USRCAT is fast as positions the molecules relative to its moments of itertia,
and uses the positions of the pharmacophores to make a vector, that can be used to calcuate the USRCATScore,
without having to superimpose the molecules (`Open3D`) which would be a combinatorial problem.
This is unsuitable for fragments/synthons, as the PMIs are easily shifted,
whereas methods that use the internal distances between pharmacophores, may be more suitable  (see section).

The cutoff of 0.7 was chosen because it is traditional in USRCAT, but also is in the tail range of the distribution of USRCAT scores.

This synthon illuestrates
![img.png](images/example_synthon.png)

![img.png](images/most_amicable_synthons.png)

The amicability of a synthon

## Other

* Enamine downloads on cluster see https://github.com/matteoferla/Fragment-hit-follow-up-chemistry
* Enamine subsampled with [enamine_subsample.py](library_subsetting/enamine_random_subsample.py)