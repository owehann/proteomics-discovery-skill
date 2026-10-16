<div align="center">

**[🏠 Main README](../README.md)** · **[📘 Skill Specification](../SKILL.md)** · **[01 QC & Null Audit](./01_forensic_qc_and_null_calibration.md)** · **[02 Abundance & Archetypes](./02_abundance_and_archetype_models.md)** · **[03 PTM, Interactome & TPP](./03_ptm_interactome_tpp_timecourse.md)** · **[04 Complexes & Networks](./04_complexes_networks_and_topology.md)** · **[05 Cohorts & RRA](./05_literature_cohorts_and_rra.md)** · **[06 BioRender & Figures](./06_biorender_svg_and_npg_figures.md)** · **[07 Editorial PPTX & Gate](./07_editorial_pptx_and_fact_gate.md)**

</div>

---

# Reference Playbook 04: Protein Complexes, Biophysical Topology, & Network KDA

> [!IMPORTANT]
> **Stage Scope (`Stage 4` — Engines `[04]` through `[08]`)**
> This playbook specifies the rules implemented in [`scripts/04_systems_complexes_topology.py`](../scripts/04_systems_complexes_topology.py) for detected-background pathway enrichment, Jaccard de-duplication, CORUM macromolecular complex stoichiometry, UniProt biophysical membrane topology, and degree-corrected STRING Key Driver Analysis (`KDA`).

---

## 1. Custom Detected-Proteome Background Rule (Mandatory)

Standard Over-Representation Analysis (`ORA`) against the entire ~20,400 human genome background severely inflates false-positive enrichment for high-abundance housekeeping modules (ribosomes, proteasomes, spliceosomes, vesicles, and metabolic enzymes) that are naturally over-represented in mass spectrometry detections.

> [!WARNING]
> **Detected-Universe Mandate**
> Every Fisher's Exact / Hypergeometric ORA test **MUST** set its universe $U$ strictly to the **empirically detected proteins in the current MS experiment** ($U = \text{all quantified proteins passing QC in the current matrix}$, typically 1,500–9,000 proteins).
> Pairwise overlapping gene sets (`Jaccard > 0.35`) are clustered into de-duplicated functional modules (`Engine [04]`), and continuous fold-change rankings are evaluated via preranked GSEA (`10,000` permutations).

---

## 2. CORUM Macromolecular Complex Stoichiometry & Co-Shift Audit (`Engine [05]`)

Proteins rarely act in isolation; they assemble into stoichiometric macromolecular machines (e.g., `26S Proteasome`, `COPII Coat`, `Exocyst`, `Integrin–Tetraspanin Web`, `ESCRT-0/I/II/III`, `MCM Replisome`, `ISGF3`, `Ferritin–NCOA4`).

For every complex $C_k$ with $\ge 3$ detected subunits in $U$:

1. **Complex Coordinate Shift**:

$$
\bar{\Delta}_{C_k} = \frac{1}{|C_k|} \sum_{g \in C_k} \Delta\log_2\mathrm{FC}_g
$$

   tested against $10,000$ random size-matched protein sets drawn from the detected universe $U$.

2. **Subunit Coherence / Stoichiometry Dispersion (`IQR`)**:

$$
\mathrm{IQR}_{C_k} = Q_{0.75}\!\left(\{\Delta_g\}_{g \in C_k}\right) - Q_{0.25}\!\left(\{\Delta_g\}_{g \in C_k}\right)
$$

   - **Coordinated Complex Assembly / Depletion**: Complex mean shift $\ge 0.40$ ($\log_2$), permutation $q < 0.05$, and directional agreement across $\ge 75\%$ of detected subunits.
   - **Orphan Subunit / Stoichiometric Uncoupling**: One subunit shifts strongly ($|\Delta_{g^*}| \ge 1.0$) while the remaining core complex satisfies **`TOST_Equiv`** (rest-of-complex shift $< 0.25$).

---

## 3. Biophysical Membrane Topology & Subcellular Redistribution (`Engine [07]`)

Directly on any protein-level matrix (`audit_membrane_and_secretory_topology()`), `ProteoSkill` annotates each protein against **UniProt Structural Features** and **Human Protein Atlas (`HPA`)** subcellular locations:
- **Classical Secretory / 1-TM / ECM Route**: `Signal Peptide (ER/Golgi)`, `Single-Pass Transmembrane (1-TM Integrin/Receptor)`, `Extracellular Matrix (Matrisome)`.
- **Non-Classical Multi-TM / Raft / Vesicular Route**: `Multi-Pass Transmembrane (>=2 TM Helices)`, `4-TM MARVEL / Tetraspanin Web`, `GPI-Anchored (KW-0336)`, `Lipidated / Palmitoylated / Prenylated Raft (KW-0449)`.

For each divergence group or compartment, a $2 \times 2$ contingency table ($a, b, c, d$) is evaluated by two-sided Fisher's exact test with Haldane–Anscombe corrected $\log_2(\text{Odds Ratio})$ and $95\%$ Wald confidence intervals:

$$
\mathrm{OR}_{\mathrm{HA}} = \frac{(a + 0.5)(d + 0.5)}{(b + 0.5)(c + 0.5)}, \qquad \mathrm{SE}(\ln \mathrm{OR}_{\mathrm{HA}}) = \sqrt{\frac{1}{a+0.5} + \frac{1}{b+0.5} + \frac{1}{c+0.5} + \frac{1}{d+0.5}}
$$

---

## 4. Modality-Scoped Upstream Regulator Inference (`Engine [08]`)

Never conflate transcriptional regulation with post-translational stability or vesicular sorting:

| Modality | Allowed Upstream Regulator Engine | Mandatory Labeling Rule |
| :--- | :--- | :--- |
| **Global WCL / Tissue (`M1`)** | `CollecTRI` / `DoRothEA` (`decoupleR` Univariate Linear Model `ULM`) | Label as **`Inferred Transcriptional Program Footprint (Requires RNA/ChIP Validation)`** |
| **EV / Secretome / Biofluid (`M1-EV`)** | Biophysical Membrane Topology + Trafficking Machinery (`Rab`, `SNARE`, `ESCRT`, `4-TM`) | Do **not** attribute cell-free EV cargo shifts directly to nuclear transcription factors |
| **Phosphoproteomics (`M2-STY`)** | **The Kinase Library** PSSMs + `PTM-SEA` + `KSEA` / `RoKAI` | Run strictly on host-protein decoupled net stoichiometry (`ΔPTM − ΔWCL`) |
| **Ubiquitinomics (`M2-Ub`)** | `UbiBrowser 2.0` / `OmniPath` E3 Ligase–Substrate networks | Separate degradative chains (`ΔWCL < −c`) from non-degradative signaling chains |

---

## 5. Degree-Corrected PPI Topology, Random Walk with Restart (`RWR`), & Key Driver Analysis (`Engine [06]`)

High-degree literature hubs in `STRING` (`TP53`, `MYC`, `AKT1`, `EGFR`, `HSPA8`, `UBC`) appear in almost any naive subnetwork simply because of high node degree. To eliminate hub bias:

1. **Confidence & Empirical Co-Abundance Gating**: Restrict `STRING` edges to high physical/functional confidence ($\mathrm{score} \ge 700$, or $\ge 450$ for cross-sector bridge discovery) and weight edges by empirical sample co-abundance when $N_{\mathrm{samples}} \ge 8$.
2. **Random Walk with Restart (`RWR`)**: Given column-normalized adjacency matrix $\mathbf{W}$, restart probability $\alpha = 0.30$, and seed vector $\mathbf{p}^{(0)}$ supported by perturbed targets and `E1/E2` core hits:

$$
\mathbf{p}^{(t+1)} = (1 - \alpha)\mathbf{W}\mathbf{p}^{(t)} + \alpha\mathbf{p}^{(0)}
$$

3. **Degree-Binned Permutation Null ($1,000\times$)**: Each node's stationary RWR score $p_i^{(\infty)}$ and neighborhood enrichment are compared against $1,000$ random seed permutations preserving the exact degree distribution ($k$-bins) of the seed nodes:

$$
Z_{\mathrm{deg}, i} = \frac{p_i^{(\infty)} - \mu_{\mathrm{perm}}(k_i)}{\sigma_{\mathrm{perm}}(k_i)}
$$

   Only nodes with degree-adjusted $Z_{\mathrm{deg}, i} \ge 1.96$ (and cross-sector "Lock-and-Key" edges linking perturbed baits to downstream cargo modules) are designated **True Topological Key Drivers (`KDA`)**.
