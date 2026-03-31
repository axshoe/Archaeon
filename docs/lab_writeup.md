# Archaeon: Mining Extremophile Enzymes from the World's Worst Neighborhoods

**The Xiu Lab | Project 4**

I have zero background in molecular biology. Not even cultural osmosis. When the average person hears "enzyme," they picture a diagram from ninth-grade biology with little pac-man shapes eating a substrate. When I started this project, that was roughly my level. Hyperthermophile was not a word I knew. I thought MGnify was probably a British startup. I had to google what an ORF was three times before it stuck.

So naturally, I decided to build a bioinformatics pipeline to discover novel enzymes from organisms that live in volcanoes.

This is that project.

*Cheers, Angie X.*

---

## Phase 1: The Problem I Didn't Know Was a Problem

It started with a throwaway question I found myself asking after a particularly long night working on DermEquity: "What's the most expensive bottleneck in industrial biotech that hasn't been solved by software yet?"

The answer I kept running into was enzyme discovery.

The industrial enzyme market is enormous. Enzymes are used to make laundry detergent work at cold temperatures, to convert corn starch into fuel, to manufacture drugs, to treat wastewater, to break down cellulose for paper production. The global market is in the tens of billions of dollars annually, and the fundamental rate-limiting step is still: find a protein that does what you want under the conditions you need.

The conditions industrial processes impose are brutal. High temperatures (60-120°C), extreme pH (2 or 12), high salt, organic solvents. Normal enzymes from normal organisms evolved to function at 37°C in a buffered aqueous environment inside a living cell. They fall apart immediately under industrial conditions.

Except for enzymes from extremophiles.

Extremophiles are organisms that thrive in exactly those hostile environments. Sulfolobus acidocaldarius lives in volcanic hot springs at pH 2 and 80°C. Pyrococcus furiosus grows near hydrothermal vents at 100°C. Halobacterium salinarum packs itself into salt crystals with water activity close to zero. Their enzymes evolved under selective pressure to function under conditions that would denature every protein in a normal organism.

The problem is finding them. Culturing these organisms requires specialized lab equipment that most academic labs don't have. Metagenomic databases contain sequences from thousands of environmental samples, including samples from every extreme environment anyone has ever bothered to sequence, but the databases are enormous and unfiltered and navigating them systematically requires infrastructure that has historically been available only to large research groups or paid bioinformatics platforms.

The computational discovery step, the step before any wet lab work happens, is the one that needed to be democratized. That's what Archaeon does.

---

## Phase 2: Learning Molecular Biology from Scratch in Two Weeks

I want to be precise about the knowledge gap I was bridging, because I think it's instructive.

When I started, I knew: proteins are made of amino acids, genes encode proteins, DNA -> RNA -> protein. That was it. I did not know what a reading frame was. I did not know what BLAST stood for or that it was from 1990. I did not know that "metagenomic" meant sequencing DNA directly from an environmental sample without culturing any specific organism first. I did not know what pLDDT meant or why AlphaFold was a big deal.

So I read. For roughly two weeks, in parallel with actual coding, I was consuming a completely unreasonable amount of bioinformatics textbook content, review papers, and EBI training materials. I read the original Kyte-Doolittle hydrophobicity paper (1982) because I wanted to understand where the numbers in the GRAVY index formula actually came from, not just use them. I read the Szilagyi and Zavodszky 2000 paper on charged residue analysis of thermophilic proteins because I wanted to understand the biological mechanism behind the pattern I was implementing, not just implement it.

Three things genuinely surprised me in this process:

First, how much of protein biochemistry reduces to surprisingly simple physical intuition. Thermostable proteins have tightly packed hydrophobic cores (hence aliphatic index), lots of salt bridges (hence charged residue fraction), and rigid loops (hence proline content). The reason these features correlate with stability is basically just: entropy is the enemy. A protein stays folded when folding is energetically favorable. Tight packing and salt bridges lower the energy of the folded state. High proline content reduces the entropy of the unfolded state (proline is conformationally restricted). These are not mysterious black-box correlations; they are consequences of thermodynamics that I could actually reason about.

Second, how much metagenomic diversity is genuinely uncultured. The fraction of microbial diversity that has been characterized in isolation in a laboratory is somewhere between 1% and 10% depending on the environment. The rest is known only from DNA sequences extracted from environmental samples. MGnify alone contains protein sequences from millions of organisms that have never been grown in a lab, some of which have never had any species-level taxonomy assigned to them. Building a pipeline that can query this space directly felt genuinely significant, not just as an exercise, but as an actual capability.

Third, ESMFold (Meta AI's protein structure predictor) is absurdly fast for a structure predictor and has a free public API. AlphaFold2 changed everything in structural biology, but it requires a multiple sequence alignment (MSA) computation that's computationally expensive and not always feasible for novel metagenomic sequences with few known homologs. ESMFold predicts structure directly from sequence using an evolutionary language model. For metagenomic sequences from undersampled biomes, the lack of MSA requirement is a feature, not a limitation. And it runs on Meta's servers for free. I love the internet.

---

## Phase 3: Designing the Industrial Applicability Score

The core intellectual challenge of this project was: given a list of protein sequences, how do you rank them by industrial usefulness without running any experiments?

The naive answer is: you can't. Experimental validation is the only ground truth. And that's correct. Archaeon does not replace experiments. What it does is prioritize which experiments are worth doing.

The IAS (Industrial Applicability Score) combines three independent signal sources into a single 0-100 ranking metric:

**Thermostability features (40% weight).** These are sequence-based predictions of thermostability: aliphatic index (core hydrophobic packing), instability index (dipeptide-based stability prediction), charged residue fraction (salt bridge density), proline content (loop rigidity), and aromaticity (aromatic stacking networks). Each feature is normalized against biological reference ranges derived from the extremophile enzyme literature, then combined with weights calibrated from the literature's consensus on feature importance. The aliphatic index gets the highest weight (30% of the 40%) because it is empirically the strongest single-feature predictor of thermostability in protein databases.

**Sequence quality and length (20% weight).** This is a penalty curve: sequences in the 150-800 amino acid range score full marks; shorter sequences are more likely to be incomplete ORFs (assembly fragments rather than real enzymes); longer sequences are harder to express in recombinant systems. E. coli, the standard industrial workhorse for recombinant enzyme production, notoriously struggles with proteins over 800 aa.

**BLAST identity to reference enzymes (40% weight).** This is the most important component. I searched the literature for the best-characterized thermostable reference enzyme for each of eight major industrial enzyme families (lipase, protease, amylase, cellulase, xylanase, laccase, peroxidase, isomerase), retrieved their NCBI accessions, and use BLAST percent identity to the best-matching reference as the scoring signal. An enzyme with 80% identity to a known Thermotoga maritima cellulase is almost certainly a cellulase variant. Function annotation is the primary bottleneck for industrial enzyme development; BLAST identity is the most efficient computational proxy for it.

An important disclaimer that I built into the report itself: the IAS weights are literature-informed priors, not experimental calibrations. I don't have a training set of labeled (thermostable/not-thermostable) proteins to fit against. If I did, I would use logistic regression or a gradient-boosted tree to learn the weights from data rather than setting them by hand. What I have instead is a theoretically grounded ranking system designed for relative comparison within a candidate set, not absolute prediction.

---

## Phase 4: Building the Pipeline

The architecture is straightforward by design. Four modules feed into a single orchestrating CLI:

**ncbi.py:** Biopython's Entrez wrappers are well-documented and the API is stable. The tricky part was constructing search queries that reliably returned extremophile sequences rather than just any protein from any organism. NCBI's [organism] and environment metadata fields are inconsistently populated. The query strings in EXTREME_BIOME_QUERIES were tuned iteratively against live API responses until they consistently returned sequences with reasonable diversity. The BLAST integration uses NCBIWWW.qblast, which submits to NCBI's remote servers, is free, and is painfully slow (30-60 seconds per query). I flag it as optional in the CLI for exactly this reason.

**mgnify.py:** The MGnify API follows JSON:API specification, which means the response structure is predictable but verbose. The key insight is that MGnify's protein database (MGnify Proteins) is distinct from its sample database; protein search requires hitting a different endpoint with biome lineage strings from their controlled vocabulary. Getting those lineage strings right took more debugging than anything else in this module. The strings are exact and case-sensitive and documented in their API docs but not in any intuitive way.

**features.py:** This is the most academically interesting module. I implemented the GRAVY, aliphatic, instability, aromaticity, proline, and charged residue formulas from scratch by reading the original papers rather than using BioPython's built-in ProteinAnalysis class. I did this partly because I wanted to understand the formulas rather than just call them, and partly because ProteinAnalysis's instability index implementation uses a different dipeptide weight table than the original Guruprasad 1990 paper and I wanted the original. The thermostability_summary function adds a human-readable interpretation layer on top of the raw features, which feeds directly into the HTML report's narrative section for each candidate.

**scorer.py:** The IAS computation itself is not complex; it's a weighted sum of normalized component scores. The interesting design decision was the normalization. For the aliphatic index, I use a linear normalization anchored at biological reference ranges (AI < 60 = mesophilic, AI > 120 = highly thermostable) rather than min-max normalization across the candidate set. Min-max normalization would make the scores depend on which candidates happen to be in the current batch, which would make IAS scores incomparable across runs. Absolute normalization makes scores comparable across any set of candidates.

**structure.py:** The ESMFold API accepts a raw amino acid sequence as the POST body and returns a PDB file as plain text. The pLDDT scores are embedded in the B-factor column of ATOM records (one value per residue). Parsing these out is a matter of reading the fixed-width PDB format. The visualization in the HTML report uses 3Dmol.js's gradient coloring on the B-factor channel to reproduce the standard pLDDT color scheme from AlphaFold2.

**report.py:** The HTML report is self-contained: Chart.js and 3Dmol.js are loaded from CDN, all data is embedded as JavaScript variables, and PDB strings for structure visualization are embedded inline. This means you can email the file, put it on GitHub Pages, or open it locally and it works without a server. The 3D structure viewer renders in the browser using WebGL; clicking a structure card shows an interactive viewer with pLDDT coloring and standard rotate/zoom/pan controls from 3Dmol.js.

---

## Phase 5: Results and Honest Assessment

Running the pipeline on NCBI with 20 sequences per biome query yields roughly 150-200 unique candidates after deduplication. In a typical run, IAS scores range from approximately 25 to 85, with most candidates clustering in the 40-60 range (which is roughly expected given the uncertain BLAST component, which defaults to 50 when BLAST is not run).

With BLAST enabled, the score distribution shifts; candidates with high identity to known enzyme families score significantly higher, and candidates with no recognizable homologs score lower. This is the intended behavior. The top-scoring candidates in full-pipeline runs are typically proteins from Thermus thermophilus, Pyrococcus furiosus, or thermophilic Bacillus species with high identity to known amylases or proteases. This makes sense: those are the best-characterized thermostable enzyme families and their sequences are dense in NCBI.

What's more interesting from a discovery standpoint are the candidates from the MGnify metagenomic sources, which tend to have lower BLAST identity to reference sequences (because they come from poorly characterized organisms) but strong thermostability feature signals. These are the candidates that represent genuine computational novelty: sequences that don't match well-characterized references but show the structural hallmarks of thermostability. They are also the riskiest candidates to pursue experimentally.

The honest limitation is this: I can't tell you whether the top-IAS candidates are actually thermostable without running wet-lab experiments. The IAS is a prioritization metric, not a ground truth. The aliphatic index threshold of >80 for thermostability is a population-level observation from a 1980 paper; individual proteins violate it constantly. A real industrial enzyme discovery pipeline would follow this computational screen with differential scanning fluorimetry (DSF) for rapid experimental thermostability measurement before any larger-scale expression work.

What Archaeon does well: it compresses a multi-step computational pipeline (database query, feature extraction, scoring, BLAST annotation, structure prediction, visualization) into a single reproducible command that any researcher with a laptop and an internet connection can run. The time from `python archaeon.py` to a ranked HTML report with 3D structures is about 5-20 minutes depending on which optional steps you enable. That is a meaningful reduction in barrier to entry for metagenomic enzyme discovery.

---

## Phase 6: What I Actually Learned

The molecular biology was interesting. The pipeline architecture was satisfying to build. But the thing I keep thinking about is the epistemological question underneath the whole project: when is a computational score useful even when it's not calibrated?

The IAS is not calibrated. I said that already. But it is grounded: each component traces back to a mechanism I understand and literature I read. The aliphatic index predicts hydrophobic core packing, which is a real thing that contributes to thermostability via a real physical mechanism. The instability index measures a real pattern in amino acid co-occurrence statistics from proteins with known stability properties. BLAST identity to a characterized homolog is the best proxy we have for function prediction short of running the experiment.

A score built from mechanistic, theoretically-grounded features is not the same as a random number generator even when it isn't fit to experimental data. It is a structured summary of prior knowledge applied to a new observation. That's what most scientific scoring systems are, including many that get used in actual industrial decisions.

What Archaeon can't do is tell you the probability that any given candidate has a thermostability above some threshold. What it can do is tell you, with some confidence grounded in mechanism, that a candidate with AI > 90, II < 35, and CRF > 0.25 is a more promising experimental target than a candidate with AI < 60, II > 55, and no BLAST hits. That narrowing, applied to a space of hundreds of candidates, is exactly the role computational screening is supposed to play before experimental work.

I started this project knowing nothing about metagenomics. I ended it understanding why the organisms living in the worst places on Earth might contain the most useful proteins on Earth. That seems like a worthwhile trade.

The universe keeps its most interesting things in the hardest-to-reach places. This is probably not a coincidence.

---

**Technologies:** Python 3.11, BioPython, Requests, ESMFold (Meta AI), NCBI Entrez, MGnify API, Chart.js, 3Dmol.js

**Data Sources:** NCBI Protein Database (nr), MGnify Metagenomic Protein Database (EBI)

**GitHub:** [github.com/axshoe/archaeon](https://github.com/axshoe/archaeon)
