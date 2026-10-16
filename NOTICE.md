# Third-Party Notices & Database Licensing

The **Proteomics Discovery & Mechanism Skill (`proteomics-discovery-skill`)** is distributed under the **MIT License**.

To comply with third-party academic and commercial database licensing terms, this repository **does not redistribute** proprietary or non-commercial bulk database dumps (such as MSigDB Hallmarks, DoRothEA, OmniPath, DrugBank, or PhosphoSitePlus). Instead, `proteomics-discovery-skill` queries public REST APIs or user-installed local packages at runtime.

Users of this skill are responsible for complying with the licenses of external databases and software packages invoked during analysis:

1. **limpa / limma / edgeR** (Bioconductor): GPL (>= 2). Cite: Li, M., Cobbold, S.A., & Smyth, G.K. (2025). Quantification and differential analysis of mass spectrometry proteomics data with probabilistic recovery of information from missing values. *bioRxiv* / *Bioconductor*; Ritchie et al. (2015) *Nucleic Acids Research*.
2. **DEqMS** (Bioconductor): LGPL. Cite: Zhu et al. (2020) *Molecular & Cellular Proteomics*.
3. **msqrob2** (Bioconductor): Artistic-2.0. Cite: Goeminne et al. (2020) *Analytical Chemistry*.
4. **MSstats / MSstatsTMT / MSstatsPTM** (Bioconductor): Artistic-2.0. Cite: Kohler et al. (2023) *Journal of Proteome Research*.
5. **The Kinase Library**: Academic/open tool. Cite: Johnson et al. (2023) *Nature* 613:759–766; Yaron-Barir et al. (2024) *Nature* 629:1174–1181.
6. **CurveCurator**: Apache-2.0. Cite: Bayer et al. (2023) *Nature Communications* 14:7902.
7. **CORUM (Comprehensive Resource of Mammalian Protein Complexes)**: CC BY-NC-SA 4.0 (fetched on-demand from `mips.helmholtz-muenchen.de`).
8. **STRING-db API**: CC BY 4.0.
9. **Reactome / WikiPathways / Gene Ontology**: CC BY 4.0.
10. **MISEV2023 EV Markers**: Adapted from Welsh et al. (2024) *Journal of Extracellular Vesicles* 13:e12404.
11. **cRAP (Common Repository of Adventitious Proteins)**: Public domain / GPM list.
