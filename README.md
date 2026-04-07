# molecular-simulation-assistant-skill

A Claude Code skill for molecular dynamics simulation—covering AMBER, VMD, PyMOL, Gaussian09, Chimera, and more.

## What it does

This skill activates when you mention anything MD-related. It gives you end-to-end guidance from theory to practice.

## Features

| Category | Coverage |
|----------|----------|
| **AMBER toolchain** | pdb4amber / tleap / antechamber / sander/pmemd / cpptraj / MMPBSA.py / parmed |
| **Visualization** | VMD (Tcl/trajectory) · PyMOL (Python API/ray tracing) · Chimera (MD Movie) |
| **Quantum chemistry** | Gaussian09 input files / basis sets / RESP charges / AMBER coupling |
| **Advanced systems** | Metalloproteins (Zn/12-6-4 LJ) · ligand-protein · membrane proteins · nucleic acids · glycans |
| **Enhanced sampling** | Metadynamics · umbrella sampling · FEP · MSM |
| **AI + MD** | DeepMD · ANI · machine learning force fields |

## File structure

```
molecular-simulation/
├── SKILL.md                          # main skill file (v2.0)
├── integrated_knowledge.md            # compressed knowledge base
├── research_agent1_output.md         # AMBER tutorial research
├── research_agent2_output.md         # software manual research (VMD/PyMOL/Gaussian/Chimera)
├── 3CA2_MD_Simulation_Protocol.md   # 3CA2 carbonic anhydrase MD workflow
├── 3CA2_PreMD_Protocol.md            # tleap preprocessing
└── md-skill_prompt.txt               # original task description
```

## Installation

### npm (recommended)

```bash
npm install -g @shizukuyu0207/molecular-simulation-skill
```

### git clone

```bash
git clone https://github.com/Shizuku0207/molecular-simulation-assistant-skill.git molecular-simulation
```

Or symlink:

```bash
ln -s /path/to/molecular-simulation-assistant-skill ~/.claude/skills/molecular-simulation
```

## Activation

Mention any of these in Claude Code to trigger the skill:

- `分子动力学` / `MD` / `AMBER` / `VMD` / `PyMOL` / `Gaussian` / `Chimera`
- `配体参数化` / `RMSD` / `MM-PBSA` / `金属蛋白` / `GAFF` / `轨迹分析`

## Examples

### AMBER seven-step workflow
```
pdb4amber → tleap → minimize → NVT heat → NPT equilibrate → production MD → cpptraj analysis
```

### VMD loading AMBER trajectories
```tcl
mol new system.prmtop type amberparm
animate read netcdf md.nc waitfor all
```

### Gaussian09 RESP → AMBER
```bash
antechamber -i ligand.mol2 -fi mol2 -o ligand_gaff.mol2 -fo mol2 -c bcc
parmchk2 -i ligand_gaff.mol2 -f mol2 -o ligand.frcmod
```

## Case study

3CA2 carbonic anhydrase MD simulation (ff19SB + OPC3 + Zn²⁺ nonbonded model + AMS ligand + 150mM NaCl):

- `3CA2_PreMD_Protocol.md` — tleap preprocessing
- `3CA2_MD_Simulation_Protocol.md` — full seven-step MD

## Links

- AMBER: https://ambermd.org
- AMBER tutorials: https://ambermd.org/tutorials/
- VMD: https://www.ks.uiuc.edu/Research/vmd/
- PyMOL: https://pymol.org/
- Gaussian: https://gaussian.com/
- Chimera: https://www.cgl.ucsf.edu/chimera/

## Changelog

- **v2.0.0** (2026-04-07) — VMD/PyMOL/Gaussian09/Chimera guides + AI+MD前沿
- **v1.0.0** (2026-04-06) — AMBER toolchain base

## License

MIT
