#####################################
A Complete Worked Example: *Lionepha*
#####################################

.. role:: filepath
.. role:: program

(Software Prerequisites)
========================

For this exercise, in addition to |delineate| you will need the following programs installed on your execution host:

    -   |StarBeast2|_
    -   |BPP|_
    -   |DendroPy|_ (which provides |SumTrees|_)

(Data Files)
============

All data files for this exercise are available for download here:

.. rst-class:: framebox center

    :download:`lionepha.zip </downloads/lionepha.zip>`


Introduction
============

We are going to to carry out a full species delimitation analysis, or, actually, a series of species delimitation analyses.
We are going to use the the *Lionepha* system, using a subset of the data published by David R. Maddison and John S. Sproul:

-   David R Maddison, John S Sproul. *Species delimitation, classical taxonomy and genome skimming: a review of the ground beetle genus Lionepha (Coleoptera: Carabidae)*, Zoological Journal of the Linnean Society, , zlz167, https://doi.org/10.1093/zoolinnean/zlz167.

The original data is available in the following `Dryad <https://datadryad.org>`_ repository:

-   https://datadryad.org/stash/dataset/doi:10.5061/dryad.2jm63xsjq

.. The work by Maddison and Sproul (2020) nicely lays out the basic alpha taxonomy of this system, organizing these 143 lineages into 12 nominal species:

.. -   "*erasa*" group:
..     -   *Lionepha australerasa*
..     -   *Lionepha casta*
..     -   *Lionepha disjuncta*
..     -   *Lionepha erasa*
..     -   *Lionepha kavanaughi*
..     -   *Lionepha lindrothi*
..     -   *Lionepha probata*
.. -   "*osculans*" group:
..     -   *Lionepha osculans*
..     -   *Lionepha pseudoerasa*
..     -   *Lionepha sequoiae*
..     -   *Lionepha tuulukwa*

The work by Maddison and Sproul (2020) nicely lays out the basic alpha taxonomy of this system, organizing these 143 lineages into 12 nominal species.
This alpha taxonomy is based on a thorough analyses integrating morphological as well as genomic information, and thus we already have a very clear understanding of the species boundaries in this group, precluding the need for a species delimitation analysis.

However, for the purposes of this exercise we are going to pretend that this is not the case.
Instead, we are going to consider that while the species identities of *some* of the population lineages are well known, in this sample there are number of lineages for which this is not the case.

Background
==========

Study Design
------------

We can imagine our study sampled individuals from an extensive range of localities, supplementing it with samples from previous collections, to develop a comprehensive data set spanning a broad taxonomic as well as geographic range.
Preliminary phylogenetic analysis, which included genomic data from known, stable, species, as well as comparative morphological analyses, has confidentally established the species identities of some of these lineages.
However, a number of lineages in the sample, while coming out to sister to known species, are sufficiently different in some respects (e.g., genetic distance or morphology), to make us question whether or not we are looking at different populations of the same species or perhaps new species altogether.

Alternatively, we might imagine a study to delineate species units in a particularly problematic clade.
Here, in addition to samples from multiple populations from the clade, we would also include populations from a broader range of species in which the focal clade is a nested subtree.
The analytical workflow for this sort of study will be identical to that discussed here, with only the distribution of population lineages of unknown identity being different --- clustered into a single subtree on the larger phylogeny as opposed to scattered into multiple smaller subtrees (or single lineages).

Data
----

Here we will focus on just a subset of the data from Maddison and Sproul (2020), focusing only on lineages in the genus *Lionepha* itself.
This data consists of 8 gene alignments for 143 samples from multiple individuals from multiple populations from multiple species:

-   :filepath:`lionepha/00-alignments/18s.nex`
-   :filepath:`lionepha/00-alignments/28s.nex`
-   :filepath:`lionepha/00-alignments/argk.nex`
-   :filepath:`lionepha/00-alignments/cad.nex`
-   :filepath:`lionepha/00-alignments/coi.nex`
-   :filepath:`lionepha/00-alignments/msp.nex`
-   :filepath:`lionepha/00-alignments/topo.nex`
-   :filepath:`lionepha/00-alignments/wg.nex`

Stage I: Identification of Population Units using the Multipopulation Coalescent Model
======================================================================================

The first stage of our analysis will be delimitation of *population* units, i.e. assignment of population identities to individual lineages, under the "multipopulation coalescent".
This model was first published under the name "censored coalescent", but has since become more popularly known as the "Multispecies Coalescent" or "MSC" model.
It is a very powerful and elegant model that allows us to identify structure in genomic data arising from disruptions in Wright-Fisher panmixia due to restrictions in gene flow.

We will be using |BPP|_ to identify our population units, and |StarBeast2|_ to obtain a guide tree for the |BPP|_ analysis.
Note that these programs use the term "species" to reference the higher-level units they analyze rather than populations (as in "species delimitation" or "species tree inference", rather than "population delimitation" and "population tree inference").
This can be correct in some contexts: e.g., where the data sampling does, either coincidentally or by design, indeed represent a 1-to-1 correspondence between populations and species, or, equivalently the species being studied do indeed *not* have any within-species lineage structure.
However, more generally, this is, at best, misleading.
These programs use the MSC as their underlying model.
From a statistical perspective, it is uncontroversial that the biological units that these programs model through the "censored" coalescent/MSC that model that underlies them are (Wright-Fisher) populations.
Strictly speaking, when used for species delimitation or model species trees, programs such as |BPP|_ and |StarBeast2|_ are used to identify *population* units under the MSC, and then these population units are interpreted (with or without justification) by the investigator to correspond on a 1-to-1 basis to species units.
Thus, when using the MSC for "species delimitation" or "species tree inference", the difference between "population" and "species" is purely lexical rather than statistical.

Here, however, we are going to interpret higher-level units of organization as exactly what they are as modeled by the MSC --- populations, no more and no less.
Unfortunately, this may result in some confusion as both |BPP|_ and |StarBeast2|_ refer to the higher-level units they target as "species" at least some of the time in the various program options and documentation.
(In fact, |BPP|_, in some of the program documentation as well as n some of the various papers presenting or referencing the theory behind it acutally use the term "species" and "populations" interchangeably).
This is simply the cost of doing business.

Candidate Population Units
--------------------------

A |BPP|_ analysis requires us to identify "population" lineages as input *a priori*, some of which it will then collapse together to form "species" lineages.
We will consider these to be "candidate population" and "actual population" lineages respectively.
That is, we will provide |BPP|_ with the finest-grain units that could possibly be distinct populations as input population lineages, then use the power of the MSC model to accurately merge together our candidate populations into distinct populations ("species", in |BPP|_ terminology).
For this analysis, we will err on the side of caution, not hestitating our over-split our candidate populations, as we can rely on the MSC to collapse them if there is insufficient gene flow restriction between them to form population boundaries.
As such, we will consider every distinct geographical sample to be a distinct candidate population.

.. rst-class:: small-text compressed-table

    +-----------------------------------------------+----------------------------------------+
    | Individual                                    | Candidate Population Assignment        |
    +===============================================+========================================+
    | - L_australerasa_CA_Carson_Spur_3839          | L_australerasa_CA_Carson_Spur          |
    | - L_australerasa_CA_Carson_Spur_3840          |                                        |
    | - L_australerasa_CA_Carson_Spur_3841          |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_australerasa_CA_Homewood_Canyon_5214      | L_australerasa_CA_Homewood_Canyon      |
    +-----------------------------------------------+----------------------------------------+
    | - L_australerasa_CA_Martin_Meadow_3838        | L_australerasa_CA_Martin_Meadow        |
    +-----------------------------------------------+----------------------------------------+
    | - L_australerasa_CA_Mill_Creek_5212           | L_australerasa_CA_Mill_Creek           |
    | - L_australerasa_CA_Mill_Creek_5213           |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_australerasa_CA_Nanny_Creek_3864          | L_australerasa_CA_Nanny_Creek          |
    | - L_australerasa_CA_Nanny_Creek_3896          |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_australerasa_CA_Oyster_Lake_3844          | L_australerasa_CA_Oyster_Lake          |
    | - L_australerasa_CA_Oyster_Lake_3845          |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_australerasa_OR_Crater_Lake_4984          | L_australerasa_OR_Crater_Lake          |
    | - L_australerasa_OR_Crater_Lake_4986          |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_casta_AK_Ketchikan_4894                   | L_casta_AK_Ketchikan                   |
    +-----------------------------------------------+----------------------------------------+
    | - L_casta_AK_Prince_of_Wales_Island_4523      | L_casta_AK_Prince_of_Wales_Island      |
    +-----------------------------------------------+----------------------------------------+
    | - L_casta_CA_Mt_Tamalpais_3830                | L_casta_CA_Mt_Tamalpais                |
    +-----------------------------------------------+----------------------------------------+
    | - L_casta_CA_Soda_Creek_4049                  | L_casta_CA_Soda_Creek                  |
    +-----------------------------------------------+----------------------------------------+
    | - L_casta_CA_West_Branch_Mill_Creek_3703      | L_casta_CA_West_Branch_Mill_Creek      |
    +-----------------------------------------------+----------------------------------------+
    | - L_casta_OR_Lost_Prairie_5204                | L_casta_OR_Lost_Prairie                |
    +-----------------------------------------------+----------------------------------------+
    | - L_casta_OR_Marys_Peak_2545                  | L_casta_OR_Marys_Peak                  |
    +-----------------------------------------------+----------------------------------------+
    | - L_casta_OR_School_Creek_3041                | L_casta_OR_School_Creek                |
    +-----------------------------------------------+----------------------------------------+
    | - L_casta_WA_Taneum_Creek_1400                | L_casta_WA_Taneum_Creek                |
    +-----------------------------------------------+----------------------------------------+
    | - L_disjuncta_BC_Summit_Creek_1896            | L_disjuncta_BC_Summit_Creek            |
    | - L_disjuncta_BC_Summit_Creek_3090            |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_disjuncta_CA_Emerson_Creek_4122           | L_disjuncta_CA_Emerson_Creek           |
    +-----------------------------------------------+----------------------------------------+
    | - L_disjuncta_CA_Lily_Lake_3069               | L_disjuncta_CA_Lily_Lake               |
    +-----------------------------------------------+----------------------------------------+
    | - L_disjuncta_CA_Salmon_Creek_4133            | L_disjuncta_CA_Salmon_Creek            |
    +-----------------------------------------------+----------------------------------------+
    | - L_disjuncta_CA_Trinity_Alps_4115            | L_disjuncta_CA_Trinity_Alps            |
    +-----------------------------------------------+----------------------------------------+
    | - L_disjuncta_ID_Salmon_River_4780            | L_disjuncta_ID_Salmon_River            |
    +-----------------------------------------------+----------------------------------------+
    | - L_disjuncta_MT_Mill_Creek_4716              | L_disjuncta_MT_Mill_Creek              |
    +-----------------------------------------------+----------------------------------------+
    | - L_disjuncta_OR_Lostine_River_3848           | L_disjuncta_OR_Lostine_River           |
    +-----------------------------------------------+----------------------------------------+
    | - L_disjuncta_OR_Mt_Hood_4143                 | L_disjuncta_OR_Mt_Hood                 |
    +-----------------------------------------------+----------------------------------------+
    | - L_erasa_AK_Thompson_Pass_4059               | L_erasa_AK_Thompson_Pass               |
    +-----------------------------------------------+----------------------------------------+
    | - L_erasa_BC_Cherryville_4002                 | L_erasa_BC_Cherryville                 |
    +-----------------------------------------------+----------------------------------------+
    | - L_erasa_OR_Lost_Prairie_5197                | L_erasa_OR_Lost_Prairie                |
    | - L_erasa_OR_Lost_Prairie_5199                |                                        |
    | - L_erasa_OR_Lost_Prairie_5200                |                                        |
    | - L_erasa_OR_Lost_Prairie_5201                |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_erasa_OR_Marys_Peak_2575                  | L_erasa_OR_Marys_Peak                  |
    | - L_erasa_OR_Marys_Peak_2586                  |                                        |
    | - L_erasa_OR_Marys_Peak_2615                  |                                        |
    | - L_erasa_OR_Marys_Peak_2616                  |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_erasa_OR_Mount_Hebo_3013                  | L_erasa_OR_Mount_Hebo                  |
    | - L_erasa_OR_Mount_Hebo_3016                  |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_erasa_OR_Mt_Hood_4144                     | L_erasa_OR_Mt_Hood                     |
    +-----------------------------------------------+----------------------------------------+
    | - L_erasa_OR_Prairie_Peak_2580                | L_erasa_OR_Prairie_Peak                |
    +-----------------------------------------------+----------------------------------------+
    | - L_kavanaughi_MT_Bitterroot_River_4646       | L_kavanaughi_MT_Bitterroot_River       |
    +-----------------------------------------------+----------------------------------------+
    | - L_kavanaughi_MT_Lost_Horse_Creek_4648       | L_kavanaughi_MT_Lost_Horse_Creek       |
    +-----------------------------------------------+----------------------------------------+
    | - L_kavanaughi_OR_Little_Philips_Creek_4998   | L_kavanaughi_OR_Little_Philips_Creek   |
    +-----------------------------------------------+----------------------------------------+
    | - L_kavanaughi_OR_Lostine_River_4996          | L_kavanaughi_OR_Lostine_River          |
    +-----------------------------------------------+----------------------------------------+
    | - L_kavanaughi_OR_Lostine_River_Valley_4990   | L_kavanaughi_OR_Lostine_River_Valley   |
    | - L_kavanaughi_OR_Lostine_River_Valley_4992   |                                        |
    | - L_kavanaughi_OR_Lostine_River_Valley_4993   |                                        |
    | - L_kavanaughi_OR_Lostine_River_Valley_5000   |                                        |
    | - L_kavanaughi_OR_Lostine_River_Valley_5002   |                                        |
    | - L_kavanaughi_OR_Lostine_River_Valley_5006   |                                        |
    | - L_kavanaughi_OR_Lostine_River_Valley_5008   |                                        |
    | - L_kavanaughi_OR_Lostine_River_Valley_5010   |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_lindrothi_CA_Deadman_Creek_4140           | L_lindrothi_CA_Deadman_Creek           |
    +-----------------------------------------------+----------------------------------------+
    | - L_lindrothi_CA_East_Fork_Kaweah_River_4120  | L_lindrothi_CA_East_Fork_Kaweah_River  |
    +-----------------------------------------------+----------------------------------------+
    | - L_lindrothi_CA_Emerald_Lake_4116            | L_lindrothi_CA_Emerald_Lake            |
    | - L_lindrothi_CA_Emerald_Lake_4117            |                                        |
    | - L_lindrothi_CA_Emerald_Lake_4118            |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_lindrothi_CA_Kaiser_Pass_4121             | L_lindrothi_CA_Kaiser_Pass             |
    +-----------------------------------------------+----------------------------------------+
    | - L_lindrothi_CA_Long_Valley_Creek_5072       | L_lindrothi_CA_Long_Valley_Creek       |
    +-----------------------------------------------+----------------------------------------+
    | - L_lindrothi_CA_Sonora_Pass_4134             | L_lindrothi_CA_Sonora_Pass             |
    +-----------------------------------------------+----------------------------------------+
    | - L_lindrothi_CA_South_Fork_Bishop_Creek_3568 | L_lindrothi_CA_South_Fork_Bishop_Creek |
    +-----------------------------------------------+----------------------------------------+
    | - L_lindrothi_CA_Tioga_Lake_4132              | L_lindrothi_CA_Tioga_Lake              |
    +-----------------------------------------------+----------------------------------------+
    | - L_osculans_CA_Carson_Spur_3164              | L_osculans_CA_Carson_Spur              |
    +-----------------------------------------------+----------------------------------------+
    | - L_osculans_CA_Cold_Creek_1387               | L_osculans_CA_Cold_Creek               |
    | - L_osculans_CA_Cold_Creek_1390               |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_osculans_CA_Los_Padres_NF_3162            | L_osculans_CA_Los_Padres_NF            |
    +-----------------------------------------------+----------------------------------------+
    | - L_osculans_CA_Nanny_Creek_3721              | L_osculans_CA_Nanny_Creek              |
    +-----------------------------------------------+----------------------------------------+
    | - L_osculans_CA_Pike_County_Gulch_3846        | L_osculans_CA_Pike_County_Gulch        |
    +-----------------------------------------------+----------------------------------------+
    | - L_osculans_CA_Rainbow_1401                  | L_osculans_CA_Rainbow                  |
    +-----------------------------------------------+----------------------------------------+
    | - L_osculans_CA_Stanislaus_NF_3157            | L_osculans_CA_Stanislaus_NF            |
    +-----------------------------------------------+----------------------------------------+
    | - L_osculans_CA_Strawberry_Creek_3163         | L_osculans_CA_Strawberry_Creek         |
    +-----------------------------------------------+----------------------------------------+
    | - L_osculans_CA_Warner_Range_3161             | L_osculans_CA_Warner_Range             |
    +-----------------------------------------------+----------------------------------------+
    | - L_osculans_OR_Eugene_4593                   | L_osculans_OR_Eugene                   |
    +-----------------------------------------------+----------------------------------------+
    | - L_osculans_OR_Goodman_Creek_3158            | L_osculans_OR_Goodman_Creek            |
    +-----------------------------------------------+----------------------------------------+
    | - L_osculans_OR_Little_Philips_Creek_5001     | L_osculans_OR_Little_Philips_Creek     |
    +-----------------------------------------------+----------------------------------------+
    | - L_osculans_OR_Middle_Fork_Berry_Creek_3095  | L_osculans_OR_Middle_Fork_Berry_Creek  |
    +-----------------------------------------------+----------------------------------------+
    | - L_osculans_OR_School_Creek_2638             | L_osculans_OR_School_Creek             |
    +-----------------------------------------------+----------------------------------------+
    | - L_osculans_OR_Walton_Lake_4743              | L_osculans_OR_Walton_Lake              |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_BC_Summit_Creek_3720              | L_probata_BC_Summit_Creek              |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_CA_Algoma_Camp_3855               | L_probata_CA_Algoma_Camp               |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_CA_Ellery_Lake_4138               | L_probata_CA_Ellery_Lake               |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_CA_Middle_Martis_Creek_1161       | L_probata_CA_Middle_Martis_Creek       |
    | - L_probata_CA_Middle_Martis_Creek_1970       |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_CA_Nanny_Creek_3895               | L_probata_CA_Nanny_Creek               |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_CA_Sherman_Pass_3730              | L_probata_CA_Sherman_Pass              |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_CA_South_Fork_Bishop_Creek_3686   | L_probata_CA_South_Fork_Bishop_Creek   |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_CA_Squaw_Valley_Resort_5211       | L_probata_CA_Squaw_Valley_Resort       |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_CA_Strawberry_Creek_3832          | L_probata_CA_Strawberry_Creek          |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_CA_Tamarack_Lake_4137             | L_probata_CA_Tamarack_Lake             |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_CA_Warner_Range_3863              | L_probata_CA_Warner_Range              |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_CA_White_Mountains_3833           | L_probata_CA_White_Mountains           |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_ID_Baker_Creek_3865               | L_probata_ID_Baker_Creek               |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_ID_Galena_Summit_3722             | L_probata_ID_Galena_Summit             |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_ID_Park_Creek_3866                | L_probata_ID_Park_Creek                |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_MT_Mill_Creek_4713                | L_probata_MT_Mill_Creek                |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_MT_Prospect_Creek_4645            | L_probata_MT_Prospect_Creek            |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_NV_Ruby_Mountains_3684            | L_probata_NV_Ruby_Mountains            |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_OR_Little_Philips_Creek_4995      | L_probata_OR_Little_Philips_Creek      |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_OR_Lonesome_Spring_4744           | L_probata_OR_Lonesome_Spring           |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_OR_Lost_Prairie_3723              | L_probata_OR_Lost_Prairie              |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_OR_Lostine_River_Valley_4991      | L_probata_OR_Lostine_River_Valley      |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_OR_Lostine_River_Valley_5004      | L_probata_OR_Lostine_River_Valley      |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_OR_Mt_Ashland_3165                | L_probata_OR_Mt_Ashland                |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_OR_Odell_Creek_3867               | L_probata_OR_Odell_Creek               |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_OR_Steens_Mountains_2724          | L_probata_OR_Steens_Mountains          |
    | - L_probata_OR_Steens_Mountains_3717          |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_UT_Shingle_Creek_4198             | L_probata_UT_Shingle_Creek             |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_UT_Stansbury_Mtns_3601            | L_probata_UT_Stansbury_Mtns            |
    | - L_probata_UT_Stansbury_Mtns_3685            |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_UT_Tushar_Mountains_5037          | L_probata_UT_Tushar_Mountains          |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_WA_Blue_Mountains_3854            | L_probata_WA_Blue_Mountains            |
    +-----------------------------------------------+----------------------------------------+
    | - L_probata_WA_Taneum_Creek_1320              | L_probata_WA_Taneum_Creek              |
    +-----------------------------------------------+----------------------------------------+
    | - L_pseudoerasa_CA_Kaiser_Pass_4139           | L_pseudoerasa_CA_Kaiser_Pass           |
    +-----------------------------------------------+----------------------------------------+
    | - L_pseudoerasa_CA_Lily_Lake_3073             | L_pseudoerasa_CA_Lily_Lake             |
    +-----------------------------------------------+----------------------------------------+
    | - L_pseudoerasa_CA_Sherman_Pass_3599          | L_pseudoerasa_CA_Sherman_Pass          |
    | - L_pseudoerasa_CA_Sherman_Pass_3688          |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_pseudoerasa_CA_Strawberry_Creek_3072      | L_pseudoerasa_CA_Strawberry_Creek      |
    | - L_pseudoerasa_CA_Strawberry_Creek_3083      |                                        |
    | - L_pseudoerasa_CA_Strawberry_Creek_3086      |                                        |
    | - L_pseudoerasa_CA_Strawberry_Creek_3087      |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_pseudoerasa_CA_Trinity_Alps_4114          | L_pseudoerasa_CA_Trinity_Alps          |
    +-----------------------------------------------+----------------------------------------+
    | - L_sequoiae_CA_Bridal_Veil_Falls_3078        | L_sequoiae_CA_Bridal_Veil_Falls        |
    +-----------------------------------------------+----------------------------------------+
    | - L_sequoiae_CA_Nanny_Creek_3702              | L_sequoiae_CA_Nanny_Creek              |
    +-----------------------------------------------+----------------------------------------+
    | - L_sequoiae_CA_Strawberry_Creek_3075         | L_sequoiae_CA_Strawberry_Creek         |
    | - L_sequoiae_CA_Strawberry_Creek_3085         |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_sequoiae_OR_Oakridge_2647                 | L_sequoiae_OR_Oakridge                 |
    +-----------------------------------------------+----------------------------------------+
    | - L_sequoiae_OR_School_Creek_2614             | L_sequoiae_OR_School_Creek             |
    +-----------------------------------------------+----------------------------------------+
    | - L_tuulukwa_CA_Trinity_Alps_4113             | L_tuulukwa_CA_Trinity_Alps             |
    +-----------------------------------------------+----------------------------------------+
    | - L_tuulukwa_OR_Knowles_Creek_3700            | L_tuulukwa_OR_Knowles_Creek            |
    | - L_tuulukwa_OR_Knowles_Creek_3701            |                                        |
    +-----------------------------------------------+----------------------------------------+
    | - L_tuulukwa_OR_Marys_Peak_2581               | L_tuulukwa_OR_Marys_Peak               |
    | - L_tuulukwa_OR_Marys_Peak_2635               |                                        |
    | - L_tuulukwa_OR_Marys_Peak_2636               |                                        |
    | - L_tuulukwa_OR_Marys_Peak_2637               |                                        |
    | - L_tuulukwa_OR_Marys_Peak_2642               |                                        |
    | - L_tuulukwa_OR_Marys_Peak_2643               |                                        |
    | - L_tuulukwa_OR_Marys_Peak_3782               |                                        |
    +-----------------------------------------------+----------------------------------------+

Generating a Guide Tree for Population Delimitation
---------------------------------------------------

We will provide |BPP|_ with a guide tree for its population delimitation analysis.
We will use |StarBeast2|_ to generate this guide tree.

The full |StarBeast2|_ configuration file, generated using ``BEAUTi``, can be found at in :filepath:`lionepha/01-guidetree-estimation/sb00500M.xml`.
The alignments listed above were imported, and the following "traits" file was used to map alignment sequences to canidate population units: :filepath:`lionepha/01-guidetree-estimation/traits.txt`.

We used a single strict clock model across all genes, and a HKY+G model of substitution.

We ran this four replicates of this analysis for 500 million generations each, sampling every 500,000 generations for a total of 1000 samples from each replicates.
The first 250 samples were discarded from the 100 samples of each replicates.
Convergence was diagnosed through inspection of traces for each parameter as well as the likelihood and posterior.
ESS values for each parameter were established to be more than 250, and distributions of parameter values were compared to a "null" run (i.e., a run without data where just the prior was sampled) to confirm that the analysis learned from the data sufficiently to shift the posterior away from the prior.

The post-burn in samples from the posterior were summarized using |SumTrees|_, with the following command::

    $ sumtrees.py -b 250 \
                -s mcct \
                -e clear \
                -l clear \
                --force-rooted \
                --suppress-annotations \
                -r -o summary.guide.tre \
                run1/species.trees \
                run2/species.trees \
                run3/species.trees \
                run4/species.trees

This selects the Maximum Clade Credibility Tree (MCCT) tree for the summary topolgy, stripping all branch lengths and metadata annotations, to result in the following:

.. rst-class:: framebox center

:filepath:`examples/lionepha/01-guidetree-estimation/guidetree.nex`


Delimitation of Population Units
--------------------------------

Now that we have a guide tree that treats each distinct geographical lineage as a candidate distinct Wright-Fisher population, we will run |BPP|_ in "A10" mode to delimit the true population units under the "multipopulation coalescent" (i.e., the MSC).

For Small Datasets: the Single-Analysis Approach
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The files provided in the :filepath:`lionepha/02a-population-delimitation-pooled` directory set up a fairly straightforward |BPP|_ "A10" analyses using the data that we have collected and the guide tree we have estimated:

.. rst-class:: compressed-table
.. table:: Pooled |BPP|_ Population Delimitation Analysis
    :widths: 50 50

    +------------------------+----------------------------+
    | File                   | Description                |
    +========================+============================+
    | bpprun.input.ctl       | Control file               |
    +------------------------+----------------------------+
    | bpprun.input.chars.txt | Character data             |
    +------------------------+----------------------------+
    | bpprun.input.imap.txt  | Sequence to population map |
    +------------------------+----------------------------+
    | bpp00.job              | Execution job              |
    +------------------------+----------------------------+

The problem is that, when executing the analysis by either running the job file::

    $ bash bpp00.job

or submitting it to an execution host::

    $ qsub bpp00.job

or simply running the command directly::

    $ bpp --cfile bpprun.input.ctl

we might find that the data set is too large to analyze::

    bpp v4.1.4_linux_x86_64, 31GB RAM, 12 cores
    https://github.com/bpp
    .
    .
    .
    Total species delimitations: 23958541050464777
    Unable to allocate enough memory.


For Larger Datasets: the Subtree Decomposition Approach
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The solution is to break the data set up into subtrees and carry out a separate population delimitation analysis on each subtree.
While designing a decomposition scheme, it is important to note that we should *not* separate out candidate population lineages that might potentially belong to the same (actual) population into separate subtrees, thereby preventing |BPP|_ from being able to collapse them if it does not detect any gene flow restriction between them.
In most studies, this should not be too difficult to identify.
Even in cases where we might not know where the *species* boundaries are, with reference to the guide tree phylogeny we should be able to identify subtree clusters that do not disrupt population boundaries.
In this case of this *Lionepha* study, the nominal species clades provide very nice granularity --- small enough to be analyzed by |BPP|_, yet with no danger of breaking up an actual population.

.. figure:: images/lionepha-guidetree.png
    :alt: Guide tree for Lionepha population delimitation, showing subtrees used when decomposing into separate analyses.
    :width: 100%
    :class: figure-image

..  rst-class:: figure-caption

        **Figure**:  We cannot delimit populations for the entire data set simulatneously in |BPP|_ due to the number of lineages. Instead, we decompose the analysis into a set of smaller analyses based on subtrees.

The set up for this set of analyses can be found at:

.. rst-class:: framebox center

    :filepath:`lionepha/02b-population-delimitation-subtrees`

with each subdirectory containing a stand-alone analysis.

Collating Results of the Subtree Approach
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each of the subtree analysis now has the populations delimited under the MSC model.
Having identified these population units across various subtrees, we now need to collate and pool them.
|delineate| helpfully provides a script for you to do this fairly robustly: |bppsum|_.
This script takes as its input two sets of files:

-   the "imap" files you provided to |BPP|_ as input, which maps sequences to candidate population lineages
-   the output log of |BPP|_ analyses (*not* the MCMC log), which provides a tree at the end with posterior probability of internal nodes indicated by labels

These files can be specified in any order, but must collectively span the entire analysis: the set of all candidate population names defined across all "imap" files must be equal to the set of all candidate population names found across all trees in all |BPP|_ output log files.

In this example, we have all the independent subtree analyses packed away in subdirectories, "``01``", "``02``", etc.
Assuming we are in the |BPP|_ analysis subdirectory, :filepath:`lionepha/02b-population-delimitation-subtrees`, we could just type in all the paths::

    delineate-bppsum \
        --imap    01/bpprun.input.imap.txt \
                  02/bpprun.input.imap.txt \
                  .
                  .
                  (etc.)
        --results 01/results.out.txt \
                  02/results.out.txt \
                  .
                  .
                  (etc.)

but because of judicious naming of the files, we can use some basic shell commands to generate the list of input files::

    delineate-bppsum \
        --imap $(find . -name "*imap*") \
        --results $(find . -name "*out.txt")

Executing the above command results in::

    [delineate-bppsum] 11 BPP 'imap' files specified
    [delineate-bppsum] - Reading mapping file   1 of 11: ./00/bpprun.input.imap.txt
    [delineate-bppsum]   - (13 lineages, 7 candidate populations)
    [delineate-bppsum] - Reading mapping file   2 of 11: ./01/bpprun.input.imap.txt
    [delineate-bppsum]   - (9 lineages, 9 candidate populations)
    [delineate-bppsum] - Reading mapping file   3 of 11: ./02/bpprun.input.imap.txt
    [delineate-bppsum]   - (10 lineages, 9 candidate populations)
    [delineate-bppsum] - Reading mapping file   4 of 11: ./03/bpprun.input.imap.txt
    [delineate-bppsum]   - (14 lineages, 7 candidate populations)
    [delineate-bppsum] - Reading mapping file   5 of 11: ./04/bpprun.input.imap.txt
    [delineate-bppsum]   - (12 lineages, 5 candidate populations)
    [delineate-bppsum] - Reading mapping file   6 of 11: ./05/bpprun.input.imap.txt
    [delineate-bppsum]   - (10 lineages, 8 candidate populations)
    [delineate-bppsum] - Reading mapping file   7 of 11: ./06/bpprun.input.imap.txt
    [delineate-bppsum]   - (16 lineages, 15 candidate populations)
    [delineate-bppsum] - Reading mapping file   8 of 11: ./07/bpprun.input.imap.txt
    [delineate-bppsum]   - (34 lineages, 30 candidate populations)
    [delineate-bppsum] - Reading mapping file   9 of 11: ./08/bpprun.input.imap.txt
    [delineate-bppsum]   - (9 lineages, 5 candidate populations)
    [delineate-bppsum] - Reading mapping file  10 of 11: ./09/bpprun.input.imap.txt
    [delineate-bppsum]   - (6 lineages, 5 candidate populations)
    [delineate-bppsum] - Reading mapping file  11 of 11: ./10/bpprun.input.imap.txt
    [delineate-bppsum]   - (10 lineages, 3 candidate populations)
    [delineate-bppsum] 11 BPP output files specified
    [delineate-bppsum] - Reading output file   1 of 11: ./00/results.out.txt
    [delineate-bppsum]   - (7 candidate populations)
    [delineate-bppsum] - Reading output file   2 of 11: ./01/results.out.txt
    [delineate-bppsum]   - (9 candidate populations)
    [delineate-bppsum] - Reading output file   3 of 11: ./02/results.out.txt
    [delineate-bppsum]   - (9 candidate populations)
    [delineate-bppsum] - Reading output file   4 of 11: ./03/results.out.txt
    [delineate-bppsum]   - (7 candidate populations)
    [delineate-bppsum] - Reading output file   5 of 11: ./04/results.out.txt
    [delineate-bppsum]   - (5 candidate populations)
    [delineate-bppsum] - Reading output file   6 of 11: ./05/results.out.txt
    [delineate-bppsum]   - (8 candidate populations)
    [delineate-bppsum] - Reading output file   7 of 11: ./06/results.out.txt
    [delineate-bppsum]   - (15 candidate populations)
    [delineate-bppsum] - Reading output file   8 of 11: ./07/results.out.txt
    [delineate-bppsum]   - (30 candidate populations)
    [delineate-bppsum] - Reading output file   9 of 11: ./08/results.out.txt
    [delineate-bppsum]   - (5 candidate populations)
    [delineate-bppsum] - Reading output file  10 of 11: ./09/results.out.txt
    [delineate-bppsum]   - (5 candidate populations)
    [delineate-bppsum] - Reading output file  11 of 11: ./10/results.out.txt
    [delineate-bppsum]   - (3 candidate populations)
    [delineate-bppsum] Posterior probability threshold of 0.50: 86 populations
    [delineate-bppsum] Posterior probability threshold of 0.75: 72 populations
    [delineate-bppsum] Posterior probability threshold of 0.90: 65 populations
    [delineate-bppsum] Posterior probability threshold of 0.95: 59 populations
    [delineate-bppsum] Posterior probability threshold of 1.00: 46 populations

and produces the following files:

- :filepath:`coalescent-pops.sb2-traits.p050.txt`
- :filepath:`coalescent-pops.sb2-traits.p075.txt`
- :filepath:`coalescent-pops.sb2-traits.p090.txt`
- :filepath:`coalescent-pops.sb2-traits.p095.txt`
- :filepath:`coalescent-pops.sb2-traits.p100.txt`
- :filepath:`coalescent-pops.summary.csv`

The population boundaries and identities of the various individuals are reported at different posterior probabilities (0.50, 0.75, 0.90, 0.95, and 1.00).
A comprehensive overview of all the identities under the different posterior probabilities is provided in the file: :filepath:`coalescent-pops.summary.csv`.
The other files (ending with filenames "...sb2-traits.p0xx.txt") are |StarBeast2|_ "``traits``" files at each of those posterior probability thresholds.
These latter make setting up a |StarBeast2|_ analysis to estimate an ultrametric phylogeny relating the population units delimited at each of the posterior probability thresholds very straightforward.
With the latter, do *not* be confused by the column header "species"!
Again, this is simply due to the misleading terminology adopted for the higher-level units of organization in |StarBeast2|_ and most other programs that use the MSC.
What we are working with here are *population units*, some of which *may* also correspond to species while others may not.

Stage II. Generating the (Multipopulation Coalescent, Ultrametric) Phylogeny of Populations
===========================================================================================

We use |StarBeast2|_ to estimate an ultrametric phylogeny of population lineages.
We use the original set of alignments (found in :filepath:`examples/lionepha/00-alignments`) for as the input data for this, and a "``traits``" file that maps each of the sequence labels in the alignment to population identities assigned in the previous step.

(INCOMPLETE --- WIP)

Stage III. Assignment of Known vs. Unknown Species Identities
=============================================================

We now inspect out phylogeny and determine *a priori* species assignments for as many population lineages as we can.

(INCOMPLETE --- WIP)

Stage IV. Delimitation of Species Units
=======================================

We are now ready to run |delineate|, using the phylogeny obtained in Stage III and the constraints table in Stage IV.

(INCOMPLETE --- WIP)
