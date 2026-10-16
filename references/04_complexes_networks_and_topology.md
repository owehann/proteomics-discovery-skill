# Reference Playbook 04: Protein Complexes, Biophysical Signatures, & Network Topology

This playbook specifies the **Stage 4 (`scripts/04_systems_complexes_topology.py`)** rules for pathway enrichment, CORUM macromolecular complex stoichiometry, biophysical/post-translational motif enrichment, and degree-corrected network topology.

---

## 1. Custom Detected-Proteome Background Rule (Mandatory)

Standard ORA (Over-Representation Analysis) against the entire ~20,400 human genome background severely inflates false-positive enrichment for high-abundance housekeeping modules (ribosomes, proteasome, spliceosome, metabolic enzymes, vesicles) that are naturally over-represented in mass spectrometry detections.

**Rule**: Every Fisher's Exact / Hypergeometric ORA test **MUST** set its universe $U$ strictly to the **empirically detected proteins in the current MS experiment** ($U = \text{all quantified proteins passing QC in the current matrix}$, typically $N = 1,500\text{–}9,000$).
For continuous fold-change rankings, run **preranked GSEA (`fgsea` / `gseapy.prerank`)** on the moderated $t$-statistic (or $\text{sign}(\Delta\log_2\text{FC}) \times (-\log_{10} p)$) with $10,000$ permutations.

---

## 2. CORUM Macromolecular Complex Stoichiometry & Co-Shift Audit

Proteins rarely act in isolation; they assemble into stoichiometric macromolecular machines (e.g., 26S Proteasome, CopII coat, Exocyst, Integrin $\alpha/\beta$ heterodimers, ESCRT-0/I/II/III, Respiratory Complex I–V).
For every complex $C_k$ with $\ge 3$ detected subunits in $U$:
1. **Complex Coordinate Shift**:
   $$\bar{\Delta}_{C_k} = \frac{1}{|C_k|}\sum_{g \in C_k} \Delta\log_2\text{FC}_g$$
   tested against $10,000$ random size-matched protein sets drawn from $U$.
2. **Subunit Coherence / Stoichiometry Variance**:
   $$\text{IQR}_{C_k} = Q_{0.75}(\{\Delta_g\}_{g \in C_k}) - Q_{0.25}(\{\Delta_g\}_{g \in C_k})$$
   - **Coordinated Complex Assembly/Loss**: $|\bar{\Delta}_{C_k}| \ge 0.40$, permutation $q < 0.05$, and directional agreement $\ge 75\%$ of subunits.
   - **Orphan Subunit / Stoichiometric Uncoupling**: One subunit shifts strongly ($|\Delta_{g^*}| \ge 1.0$) while the remaining core complex satisfies **`TOST_Equiv`** ($|\bar{\Delta}_{C_k \setminus \{g^*\}}| < 0.25$).

---

## 3. Modality-Scoped Upstream Regulator Inference (`TF` vs `Kinase` vs `E3`)

Never conflate transcriptional regulation with post-translational stability or vesicular sorting:
- **In Global WCL / Tissue Abundance (`M1`)**: Run `CollecTRI` / `DoRothEA` (`decoupleR` ULM/VIPER) ONLY on intracellular/nuclear/metabolic fractions and explicitly label outputs as **`Inferred Transcriptional Program Footprint (Requires RNA-seq/ChIP Confirmation)`**.
- **In EV / Secretome / Organelle Fractions (`M1-EV`)**: Prioritize **Trafficking, Sorting, & Biophysical Motif Enrichment** (Rab GTPases, SNAREs, Tetraspanins, ESCRT, Signal Peptide, Transmembrane helices, GPI-anchor, S-Palmitoylation / SwissPalm, N-glycosylation) rather than TF footprints!
- **In Phosphoproteomics (`M2-STY`)**: Run **The Kinase Library**, `PTM-SEA`, and `KSEA`/`RoKAI` directly on site-level net stoichiometry ($\Delta\text{PTM} - \Delta\text{WCL}$).
- **In Ubiquitinomics (`M2-Ub`)**: Run `UbiBrowser` / E3–substrate enrichment and classify degradative vs non-degradative chains.

---

## 4. Degree-Corrected PPI Topology, Random Walk with Restart (`RWR`), & Key Driver Analysis (`KDA`)

High-degree hubs in `STRING` (`TP53`, `MYC`, `AKT1`, `EGFR`, `HSPA8`, `UBC`) appear in almost any naive subnetwork simply because of literature bias and high node degree.
To prevent hub bias:
1. **Confidence & Co-Abundance Gating**: Restrict `STRING` edges to high physical/functional confidence ($\text{score} \ge 700$) AND require positive empirical co-abundance across samples in the current dataset ($r_{ij} > 0.35$ when $N_{\text{samples}} \ge 8$).
2. **Random Walk with Restart (`RWR`)**:
   $$\mathbf{p}^{(t+1)} = (1 - \alpha) \mathbf{W} \mathbf{p}^{(t)} + \alpha \mathbf{p}^{(0)}, \qquad \alpha = 0.30$$
   where $\mathbf{W}$ is the column-normalized adjacency matrix and $\mathbf{p}^{(0)}$ is seeded with experimental perturbation targets or `E1/E2` core hits.
3. **Degree-Binned Permutation Null**: Each node's RWRstationary score $p_i^{(\infty)}$ or KDA enrichment is compared against $1,000$ random seed permutations preserving the exact degree distribution ($k$-bins) of the seed nodes. Only nodes with degree-adjusted empirical $p_{\text{deg-null}} < 0.01$ are called **True Topological Bridge / Key Driver Nodes**.
