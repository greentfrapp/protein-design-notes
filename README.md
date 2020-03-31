# protein-design-notes

Notes while learning about protein design.

**What is the difference between *ab initio* and *de novo*?**

From [here](http://biotech.fyicenter.com/faq/Robetta/What_is_the_difference_between_Ab_Initio_and_De_.html).

> *Ab initio* structure prediction classically refers to structure prediction using nothing more than first-principles (i.e. physics). *De novo* is a more general term that refers to the greater category of methods that do not use templates from homologous PDB structures.

**Homology**

From [Wikipedia](https://en.wikipedia.org/wiki/Homology_(biology)).

> ... similarity due to shared ancestry between a pair of structures or genes in different taxa ...

In the context of proteins, homologous proteins refer to non-identical proteins that perform the same functions in different species, implying a similar ancestry.

From [Biology LibreTexts](https://bio.libretexts.org/Bookshelves/Microbiology/Book%3A_Microbiology_(Boundless)/7%3A_Microbial_Genetics/7.13%3A_Bioinformatics/7.13C%3A_Homologs%2C_Orthologs%2C_and_Paralogs).

> A homologous gene (or homolog) is a gene inherited in two species by a common ancestor. While homologous genes can be similar in sequence, similar sequences are not necessarily homologous.

> Orthologous are homologous genes where a gene diverges after a speciation event, but the gene and its main function are conserved.

> If a gene is duplicated in a species, the resulting duplicated genes are paralogs of each other, even though over time they might become different in sequence composition and function.

**Conservation**

From [Biology LibreTexts](https://bio.libretexts.org/Bookshelves/Microbiology/Book%3A_Microbiology_(Boundless)/7%3A_Microbial_Genetics/7.13%3A_Bioinformatics/7.13C%3A_Homologs%2C_Orthologs%2C_and_Paralogs) (bottom).

> In biology, conserved sequences are similar or identical sequences that  occur within nucleic acid sequences (such as RNA and DNA sequences),  protein sequences, protein structures.

On the other hand for proteins, conservation in amino acid sequences refer to: 

> sequences in which the amino acid at a specific position has been  substituted with a different one with functionally equivalent  physicochemical properties

**Domain**

From [Wikipedia](https://en.wikipedia.org/wiki/Protein_domain).

> A **protein domain** is a conserved part of a given protein sequence and [tertiary structure](https://en.wikipedia.org/wiki/Biomolecular_structure#Tertiary_structure) that can [evolve](https://en.wikipedia.org/wiki/Biological_evolution), function, and exist independently of the rest of the protein chain. Each domain forms a compact three-dimensional structure and often can be independently stable and [folded](https://en.wikipedia.org/wiki/Protein_folding). Many proteins consist of several structural domains. One domain may appear in a variety of different proteins. [Molecular evolution](https://en.wikipedia.org/wiki/Molecular_evolution) uses domains as building blocks and these may be recombined in different arrangements to create [proteins](https://en.wikipedia.org/wiki/Protein) with different functions.

![](https://upload.wikimedia.org/wikipedia/commons/thumb/7/7a/Pyruvate_kinase_protein_domains.png/250px-Pyruvate_kinase_protein_domains.png)

> Pyruvate kinase, a protein with three domains ([PDB](https://en.wikipedia.org/wiki/Protein_Data_Bank): [1PKN](https://www.rcsb.org/structure/1PKN)).

**Motifs**

***Sequence Motifs***

From [Wikipedia](https://en.wikipedia.org/wiki/Sequence_motif).

> In [genetics](https://en.wikipedia.org/wiki/Genetics), a **sequence motif** is a [nucleotide](https://en.wikipedia.org/wiki/Nucleotide) or [amino-acid](https://en.wikipedia.org/wiki/Amino_acid) [sequence](https://en.wikipedia.org/wiki/Sequence) pattern that is widespread and has, or is conjectured to have, a [biological](https://en.wikipedia.org/wiki/Biology) significance. For proteins, a sequence motif is distinguished from a [structural motif](https://en.wikipedia.org/wiki/Structural_motif), a motif formed by the three-dimensional arrangement of amino acids which may or may not be adjacent.

***Structural Motifs***

From [Wikipedia](https://en.wikipedia.org/wiki/Structural_motif).

> In a [chain-like](https://en.wikipedia.org/wiki/Polymer) biological [molecule](https://en.wikipedia.org/wiki/Molecule), such as a [protein](https://en.wikipedia.org/wiki/Protein) or  [nucleic acid](https://en.wikipedia.org/wiki/Nucleic_acid), a **structural motif** is a [supersecondary structure](https://en.wikipedia.org/wiki/Supersecondary_structure), which also appears in a variety of other molecules. Motifs do not allow us to predict the biological functions: they are found in proteins and  enzymes with dissimilar functions.

Some examples of structural motifs in proteins:

> - [Beta hairpin](https://en.wikipedia.org/wiki/Beta_hairpin):
>
>   Extremely common. Two [antiparallel](https://en.wikipedia.org/wiki/Antiparallel_(biochemistry)) beta strands connected by a tight turn of a few amino acids between them.
>
> - [Greek key](https://en.wikipedia.org/wiki/Beta_sheet#Greek_key_motif):
>   Four beta strands, three connected by hairpins, the fourth folded over the top.
> - [Omega loop](https://en.wikipedia.org/wiki/Omega_loop):
> A loop in which the residues that make up the beginning and end of the loop are very close together.
> - [Helix-loop-helix](https://en.wikipedia.org/wiki/Basic-helix-loop-helix):
> Consists of [alpha helices](https://en.wikipedia.org/wiki/Alpha_helix) bound by a looping stretch of amino acids. This motif is seen in transcription factors.
> - [Zinc finger](https://en.wikipedia.org/wiki/Zinc_finger):
> Two beta strands with an alpha helix end folded over to bind a zinc [ion](https://en.wikipedia.org/wiki/Ion). Important in DNA binding proteins.
> - [Helix-turn-helix](https://en.wikipedia.org/wiki/Helix-turn-helix):
> Two α helices joined by a short strand of amino acids and found in many proteins that regulate gene expression.
> - [Nest](https://en.wikipedia.org/wiki/Nest_(protein_structural_motif)):
> Extremely common. Three consecutive amino acid residues form an anion-binding concavity.
> - [Niche](https://en.wikipedia.org/wiki/Niche_(protein_structural_motif)):
> Extremely common. Three or four consecutive amino acid residues form a cation-binding feature.

***Short Linear Motifs***

> In molecular biology **Short Linear Motifs** (also known as **SLiMs**, **Linear Motifs** or **minimotifs**) are short stretches of [protein](https://en.wikipedia.org/wiki/Protein) sequence that mediate [protein–protein interaction](https://en.wikipedia.org/wiki/Protein–protein_interaction).

**Position-specific Scoring Matrix (PSSM)**

*Or Position Weight Matrix (PWM) or Position-specific Weight Matrix (PSWM)*

Given a set of same-length motifs or sequences, the PSSM is a representation of the set and can be used to calculate the probability of a new sequence belonging to that set.

In the PSSM, each row corresponds to an amino acid (or nucleotide for DNA sequences) and each column corresponds to a position in the sequence. Each element in the PSSM represents the log likelihood of the amino acid in that position, after accounting for a background model.

See the [Wikipedia article](https://en.wikipedia.org/wiki/Position_weight_matrix) for a quick example.

PSSMs can also be used to calculate Information Content (IC) which is a measure of how different the motif is from a uniform distribution.

**Multiple Sequence Alignment**

From [Wikipedia](https://en.wikipedia.org/wiki/Multiple_sequence_alignment).

> A **multiple sequence alignment** (**MSA**) is a [sequence alignment](https://en.wikipedia.org/wiki/Sequence_alignment) of three or more [biological sequences](https://en.wikipedia.org/wiki/Biological_sequence), generally [protein](https://en.wikipedia.org/wiki/Protein), [DNA](https://en.wikipedia.org/wiki/DNA), or [RNA](https://en.wikipedia.org/wiki/RNA). In many cases, the input set of query sequences are assumed to have an [evolutionary](https://en.wikipedia.org/wiki/Evolutionary) relationship by which they share a linkage and are descended from a common ancestor. From the resulting MSA, sequence [homology](https://en.wikipedia.org/wiki/Homology_(biology)) can be inferred and [phylogenetic analysis](https://en.wikipedia.org/wiki/Molecular_phylogeny) can be conducted to assess the sequences' shared evolutionary origins.

![](https://upload.wikimedia.org/wikipedia/commons/thumb/7/79/RPLP0_90_ClustalW_aln.gif/800px-RPLP0_90_ClustalW_aln.gif)

> Note the two completely conserved residues [arginine](https://en.wikipedia.org/wiki/arginine) (R) and [lysine](https://en.wikipedia.org/wiki/lysine) (K) marked with an asterisk at the top of the alignment.
>
> (Also note the gaps, represented by dashes, that are inserted for alignment.)

Given a set of sequences (not necessarily of equal length), gaps are added in the sequences to produce a new set of equal-length sequences that are maximally aligned, with minimal insertions and substitutions between sequences. 

This is related to the longest common subsequence problem which is NP-complete.

**PyMol**

Software for visualization of proteins.

See [this](https://www.youtube.com/watch?v=wiKyOF-pGw4) for a tutorial series on using PyMol.

**BLAST**

aka Basic Local Alignment Search Tool

From [Wikipedia](https://en.wikipedia.org/wiki/BLAST_(biotechnology)).

> In [bioinformatics](https://en.wikipedia.org/wiki/Bioinformatics), **BLAST** (**basic local alignment search tool**) is an [algorithm](https://en.wikipedia.org/wiki/Algorithm) and program for comparing [primary](https://en.wikipedia.org/wiki/Primary_structure) biological sequence information, such as the [amino-acid](https://en.wikipedia.org/wiki/Amino_acid) sequences of [proteins](https://en.wikipedia.org/wiki/Protein) or the [nucleotides](https://en.wikipedia.org/wiki/Nucleotide) of [DNA](https://en.wikipedia.org/wiki/DNA_sequence) and/or [RNA](https://en.wikipedia.org/wiki/RNA) sequences. A BLAST search enables a researcher to compare a subject  protein or nucleotide sequence (called a query) with a library or [database](https://en.wikipedia.org/wiki/Database) of sequences, and identify library sequences that resemble the query sequence above a certain threshold.
>
> Different types of BLASTs are available according to the query sequences and the target databases. For example, following the discovery of a  previously unknown gene in the [mouse](https://en.wikipedia.org/wiki/Mus_musculus), a scientist will typically perform a BLAST search of the [human genome](https://en.wikipedia.org/wiki/Human_genome) to see if humans carry a similar gene; BLAST will identify sequences in the human genome that resemble the mouse gene based on similarity of  sequence. 

**Redundancy**

Largely similar protein sequences may be "redundant", where a single sequence can be used to concisely represent a cluster of similar sequences. Non-redundant datasets typically use BLAST to identify similar sequences before selecting a single representative sequence from each cluster.

See [Redundancy in the PDB](http://www.rcsb.org/pdb/statistics/clusterStatistics.do).

**Molecular Dynamics (MD)**

From [Wikipedia](https://en.wikipedia.org/wiki/Molecular_dynamics).

> **Molecular dynamics** (**MD**) is a [computer simulation](https://en.wikipedia.org/wiki/Computer_simulation) method for analyzing the [physical movements](https://en.wikipedia.org/wiki/Motion_(physics)) of [atoms](https://en.wikipedia.org/wiki/Atoms) and [molecules](https://en.wikipedia.org/wiki/Molecules). The atoms and molecules are allowed to interact for a fixed period of time, giving a view of the [dynamic](https://en.wikipedia.org/wiki/Dynamics_(mechanics)) "evolution" of the system. In the most common version, the [trajectories](https://en.wikipedia.org/wiki/Trajectory) of atoms and molecules are determined by [numerically solving](https://en.wikipedia.org/wiki/Numerical_integration) [Newton's equations of motion](https://en.wikipedia.org/wiki/Newton's_laws_of_motion) for a system of interacting particles, where [forces](https://en.wikipedia.org/wiki/Force_(physics)) between the particles and their [potential energies](https://en.wikipedia.org/wiki/Potential_energy) are often calculated using [interatomic potentials](https://en.wikipedia.org/wiki/Interatomic_potential) or [molecular mechanics](https://en.wikipedia.org/wiki/Molecular_mechanics) [force fields](https://en.wikipedia.org/wiki/Force_field_(chemistry)).