<div align="center">
  <img src="./assets/logo.png" alt="ProteoSkill 3D Protein Logo" width="136" height="136" />

  # ProteoSkill

  **Turn Raw Proteomics `.xlsx` Tables into Calibrated Statistics, BioRender-Style Mechanism Figures, and Group-Meeting Slides — Just by Talking to Your AI.**

  [![CI Status](https://github.com/owehann/proteomics-discovery-skill/actions/workflows/ci.yml/badge.svg)](https://github.com/owehann/proteomics-discovery-skill/actions/workflows/ci.yml)
  [![Agent Skills Standard](https://img.shields.io/badge/Agent_Skills-agentskills.io-0F766E?style=flat-square)](https://agentskills.io)
  [![License: MIT](https://img.shields.io/badge/License-MIT-0284C7?style=flat-square)](./LICENSE)
  [![Engine: R + Python](https://img.shields.io/badge/Engine-Bioconductor_%2B_Python-4F46E5?style=flat-square)](#how-to-use-zero-coding-required)
</div>

---

## How It Works

<div align="center">
  <img src="./assets/workflow_banner.png" alt="ProteoSkill 4-Step Workflow" width="900" />
</div>

**ProteoSkill** is a self-contained AI Agent Skill built for **wet-lab biologists, clinicians, and bioinformaticians**. Instead of writing R or Python scripts by hand, you simply drop your mass-spectrometry or affinity proteomics spreadsheet (`.xlsx`, `.tsv`, `.csv`) into your folder and talk to your AI assistant (**Claude Code**, **Antigravity**, **Cursor**, or **Codex**).

---

## How to Use (Zero Coding Required)

You do **not** need to know how to use the command line or write code. Just copy and paste the following prompts to your AI assistant:

### Step 1 — Ask Your AI to Install ProteoSkill
Copy and send this message to your AI:
> **"Please install ProteoSkill from `https://github.com/owehann/proteomics-discovery-skill` into my global skills directory and set up its Python/R dependencies for me."**

### Step 2 — Drop Your `.xlsx` File & Check the Plan (Gate 1)
Put your exported proteomics table (e.g., `results.xlsx`) in your project folder and tell your AI:
> **"I placed my proteomics data in `results.xlsx`. Please use ProteoSkill to inspect the sample groups, check data quality, and show me the Gate 1 analysis plan before running."**

Your AI will automatically read the file headers, infer your control and treatment groups, check sample sizes, and present a clear **Gate 1 Confirmation Plan** for you to review.

### Step 3 — Confirm & Let Your AI Do Everything Else
Once you confirm the sample grouping (or tell the AI how to adjust it), reply:
> **"The grouping looks right. Please run the full ProteoSkill workflow, compare our top hits with published literature, draw a BioRender-style mechanism figure, and generate a group-meeting PowerPoint."**

Your AI will automatically execute the statistical modeling, query PubMed/EuropePMC and public clinical cohorts, generate publication-grade figures, and verify every number before handing you the final report and `.pptx` slide deck.

> **When will the AI pause to ask you?**
> - **Gate 1 (Always)**: To confirm sample groupings and experimental design before running.
> - **Gate 2 (Only if needed)**: If an abnormal outlier sample tube is detected or your target protein is missing in the mass-spec table.
> - **Gate 3 (Only if needed)**: If an external public dataset requires manual download.

---

## Why ProteoSkill? (Solving 4 Common Traps in Proteomics Analysis)

1. **Unreplicated or Small-Sample Overclaiming (`n = 1, 2, or more`)**:
   Standard AI prompts treat single-tube pilot screens or small-sample shared-control designs (`shRNA1/2` sharing one control) as if they were large cohorts. ProteoSkill calculates the exact **Empirical Null False Discovery Rate** and assigns an honest evidence grade (`E1`–`E4`) to every comparison.
2. **Confusing "Not Significant" with "Unchanged"**:
   Having `p > 0.05` does not prove a pathway is unaffected or that two antibodies are equivalent. ProteoSkill runs formal **Two One-Sided Tests (`TOST`) equivalence testing** to prove true biological conservation.
3. **Hidden Spreadsheet & Sample Artifacts**:
   ProteoSkill automatically repairs Excel date-corrupted gene names (`1-Mar` $\rightarrow$ `MARCHF1`, `1-Sep` $\rightarrow$ `SEPTIN1`), strips exogenous antibody/serum carry-over, checks vesicle purity, and quarantines failed sample tubes before normalization.
4. **AI Hallucinations in Reports & Citations**:
   Before delivering any report or slide deck, a built-in **Fact & Citation Gate** cross-checks every gene fold-change against the raw CSV table, flags single-peptide (`1-pep`) detections, and verifies every `PMID`/`DOI` live against PubMed/EuropePMC.

---

## Visual Showcase: BioRender-Style Mechanism Schematics

<div align="center">
  <img src="./assets/biorender_showcase.png" alt="BioRender-Style Cellular Mechanism Showcase" width="880" />
  <p><em>Example cellular signaling and extracellular vesicle (EV) mechanism schematic. High-confidence statistical findings (<code>E1/E2</code>) are rendered as solid arrows, whereas exploratory hypotheses (<code>E3</code>) are automatically distinguished as dashed arrows.</em></p>
</div>

---

## Supported Proteomics Modalities

| Modality | Typical Experiments | Key Capabilities |
| :--- | :--- | :--- |
| **1. Global Abundance** | Gene knockdown/overexpression, multi-shRNA, antibody comparison, clinical tumor vs. adjacent cohorts, fluid vs. tissue secretion | Missing-value recovery, peptide-stratified empirical Bayes (`DEqMS`/`limma`), `TOST` equivalence, `KD × OE` rescue grid, secretion gradient (`SGI`) |
| **2. Post-Translational Modifications** | Phosphorylation, ubiquitination, acetylation, glycosylation | Site-to-protein abundance decoupling and upstream kinase/ligase activity inference |
| **3. Interactomics & Proximity** | AP-MS, Co-IP, BioID, TurboID | True interactor scoring (`SAINT`/`MiST`), background contaminant removal, and protein complex stoichiometry |
| **4. Thermal Stability & Drug Binding** | TPP, CETSA, PISA, LiP-MS | Melting curve fitting ($\Delta T_m$) and functional target engagement testing (`NPARC`) |
| **5. Dose & Time Dynamics** | Multi-dose $\text{EC}_{50}/\text{IC}_{50}$ curves and longitudinal time-courses | Dose-response curve fitting (`CurveCurator`) and temporal trajectory clustering |
| **6. Targeted & Affinity** | PRM/MRM, Single-Cell Proteomics, Olink, SomaScan | Heavy-isotope ratio calibration, limit-of-quantification auditing, and bridge normalization |

---

<details>
<summary><strong>For Bioinformaticians & Developers: Manual CLI Installation</strong></summary>

<br/>

If you prefer to install and test the repository manually from the terminal:

```bash
# 1. Clone into your global Agent Skills directory
mkdir -p ~/.gemini/antigravity/skills ~/.agents/skills
git clone https://github.com/owehann/proteomics-discovery-skill.git ~/.gemini/antigravity/skills/proteomics-discovery-skill
ln -sfn ~/.gemini/antigravity/skills/proteomics-discovery-skill ~/.agents/skills/proteomics-discovery-skill

# 2. Install Python dependencies (or use environment.yml for full R/Bioconductor + Python)
pip install -r ~/.gemini/antigravity/skills/proteomics-discovery-skill/requirements.txt

# 3. Run the 5-suite synthetic validation benchmark
python3 ~/.gemini/antigravity/skills/proteomics-discovery-skill/tests/run_validation_suite.py
```

</details>

---

## Citation & License

Released under the [MIT License](./LICENSE). Third-party biological databases queried at runtime are governed by their respective academic licenses as documented in [`NOTICE.md`](./NOTICE.md). See [`CITATION.cff`](./CITATION.cff) for citation metadata.
