# Windows Fast-Track Setup & Execution Channel (`win64` / `x64` / `ARM64` / `win32`)

> [!IMPORTANT]
> **Dedicated Isolated Channel for Windows (`Windows 10 / 11`, 64-bit `win64/x86_64`, `ARM64`, and `win32`)**
> When the host operating system is **Windows** (`platform.system() == "Windows"` or PowerShell / CMD), AI assistants **MUST** follow this fast-track channel instead of the macOS/Linux commands in `README.md` or `environment.yml`. This guarantees a **< 30-second installation** without modifying or polluting the repository's core files.

---

## 1. Why Windows Uses This Dedicated Channel

1. **Skip `conda` / `bioconda` & R Compilation Traps**:
   - Do **NOT** run `conda env create -f environment.yml` or attempt to compile R / Bioconductor packages (`MSstats`, `DEqMS`, `fgsea`) via `Rtools` on Windows unless `Rscript.exe` and those packages are *already* installed on `PATH`.
   - `scripts/02_statistical_engine_py.py` implements 100% of the statistical models (`DEqMS` peptide-stratified empirical Bayes, `DPC` missingness, `TOST` equivalence, `NPARC`, and `CurveCurator`) natively in Python with zero loss of functionality.
2. **Skip `cairosvg` (`libcairo-2.dll`) & `pyarrow` C++ Source Builds**:
   - Windows does not ship with the Unix `libcairo-2.dll` / `GTK+` C runtime required by `cairosvg`, and `pyarrow` is only needed for `.parquet` files.
   - Instead of `pip install -r requirements.txt`, install the 9 core pre-compiled binary wheels (`--only-binary=:all:`) + `resvg-py` (a self-contained Rust wheel for `.svg` -> `.png` rendering on Windows with zero external DLL dependencies).

---

## 2. One-Step Fast Installation on Windows (PowerShell)

Run the following PowerShell commands to install `ProteoSkill` and its precompiled binary dependencies in **20–30 seconds**:

```powershell
# 1. Create global Agent Skills directories on Windows
$SkillRoot = "$env:USERPROFILE\.gemini\antigravity\skills\proteomics-discovery-skill"
New-Item -ItemType Directory -Force -Path "$env:USERPROFILE\.gemini\antigravity\skills" | Out-Null
New-Item -ItemType Directory -Force -Path "$env:USERPROFILE\.agents\skills" | Out-Null
New-Item -ItemType Directory -Force -Path "$env:USERPROFILE\.claude\skills" | Out-Null

# 2. Clone repository (or update if already present)
if (-not (Test-Path "$SkillRoot\.git")) {
    git clone https://github.com/owehann/proteomics-discovery-skill.git "$SkillRoot"
}

# 3. Create Windows directory junctions (replaces Unix `ln -sfn` without requiring Admin rights)
if (-not (Test-Path "$env:USERPROFILE\.agents\skills\proteomics-discovery-skill")) {
    cmd /c mklink /J "$env:USERPROFILE\.agents\skills\proteomics-discovery-skill" "$SkillRoot"
}
if (-not (Test-Path "$env:USERPROFILE\.claude\skills\proteomics-discovery-skill")) {
    cmd /c mklink /J "$env:USERPROFILE\.claude\skills\proteomics-discovery-skill" "$SkillRoot"
}

# 4. Fast binary-only pip install (~20 seconds; skips cairosvg DLL & pyarrow C++ compilation)
python -m pip install --upgrade --only-binary=:all: numpy scipy pandas matplotlib openpyxl pyyaml python-pptx pillow networkx resvg-py
```

*(Note: If the user is on a slow network in mainland China, append `-i https://pypi.tuna.tsinghua.edu.cn/simple` to the `pip install` command).*

---

## 3. Windows SVG-to-PNG Rendering Rule (`Stage 6`)

When generating the 7 BioRender schematics (`scripts/07_biorender_canvas_engine.py`) on Windows:
- Do **NOT** attempt to install `GTK+`, `cairo`, or `rsvg-convert` via `conda` or `winget`.
- If `.png` files are not automatically rasterized by `07_biorender_canvas_engine.py` (because `rsvg-convert` and `libcairo-2.dll` are Unix libraries), run this 1-line Python helper using `resvg_py` (installed in Step 4 above) to convert all `.svg` schematics in `figures/` to 300-DPI `.png` files in < 1 second:

```powershell
python -c "import glob, pathlib, resvg_py; [pathlib.Path(p).with_suffix('.png').write_bytes(resvg_py.svg_to_bytes(svg_path=p, width=2560, dpi=300)) for p in glob.glob('figures/*.svg')]"
```

*(Fallback if `resvg-py` is not installed: use pre-installed Windows Microsoft Edge headless via `msedge --headless --disable-gpu --screenshot="figures/s01_biorender_hero.png" --window-size=1400,820 "figures/s01_biorender_hero.svg"`).*

---

## 4. Windows Execution Rules (`Stages 0–7`)

1. **Always set UTF-8 console encoding** in Windows PowerShell before running scripts so bilingual Chinese/English notes and mathematical symbols (`Δ`, `σ`, `≥`) never raise `UnicodeEncodeError` (`cp936` / `cp1252`):
   ```powershell
   $env:PYTHONUTF8 = "1"
   ```
2. **Stage 2 Statistical Engine**: Invoke `python scripts/02_statistical_engine_py.py` directly on Windows (do not block waiting for R/Bioconductor).
3. **Verify Installation**: Run `python tests/run_validation_suite.py` with `$env:PYTHONUTF8 = "1"` to verify all 8 test suites pass.
