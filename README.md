# molecular-simulation-assistant-skill

> Claude Code 分子动力学全栈技能包 —— 覆盖 AMBER / VMD / PyMOL / Gaussian09 / Chimera 等主流 MD 软件。

## 简介

本 skill 为 Claude Code (Claude Desktop) 提供完整的**分子动力学模拟**知识支持。当用户提及 MD 相关关键词时自动激活，提供从理论到实战的端到端指导。

## 核心功能

| 类别 | 内容 |
|------|------|
| **AMBER 全工具链** | pdb4amber / tleap / antechamber / sander/pmemd / cpptraj / MMPBSA.py / parmed |
| **可视化软件** | VMD (Tcl脚本/轨迹分析) · PyMOL (Python API/渲染) · Chimera (MD Movie) |
| **量子化学** | Gaussian09 输入文件 / 基组选择 / RESP 电荷 / 与 AMBER 联用 |
| **高阶体系** | 金属蛋白 (Zn/12-6-4 LJ) · 配体-蛋白 · 膜蛋白 · 核酸 · 糖类 |
| **增强采样** | Metadynamics · 伞形采样 · FEP · MSM |
| **AI + MD** | DeepMD · ANI · 机器学习力场 |

## 文件结构

```
molecular-simulation/
├── SKILL.md                          # 主技能文件 (v2.0, 14KB)
├── integrated_knowledge.md            # 整合压缩知识库
├── research_agent1_output.md         # AMBER 官方教程调研
├── research_agent2_output.md         # 软件手册调研 (VMD/PyMOL/Gaussian/Chimera)
├── 3CA2_MD_Simulation_Protocol.md   # 3CA2 碳酸酐酶 MD 实战流程
├── 3CA2_PreMD_Protocol.md            # tleap 前处理详细流程
└── md-skill_prompt.txt               # 原始需求描述
```

## 快速开始

### 安装

将本仓库克隆到 Claude Code skills 目录：

```bash
git clone https://github.com/Shizuku0207/molecular-simulation-assistant-skill.git molecular-simulation
```

或软链接到 skills 目录：

```bash
ln -s /path/to/molecular-simulation-assistant-skill ~/.claude/skills/molecular-simulation
```

### 触发方式

在 Claude Code 中直接提及以下任一关键词即可激活：

- `分子动力学` / `MD模拟` / `AMBER`
- `VMD` / `PyMOL` / `Gaussian` / `Chimera`
- `配体参数化` / `RMSD` / `MM-PBSA` / `金属蛋白`
- `GAFF力场` / `ff19SB` / `轨迹分析`

## Skill 内容示例

### AMBER 七步法
```
pdb4amber → tleap → 最小化 → NVT加热 → NPT平衡 → 成品MD → cpptraj分析
```

### VMD 读取 AMBER 轨迹
```tcl
mol new system.prmtop type amberparm
animate read netcdf md.nc waitfor all
```

### Gaussian09 RESP 电荷 → AMBER
```bash
antechamber -i ligand.mol2 -fi mol2 -o ligand_gaff.mol2 -fo mol2 -c bcc
parmchk2 -i ligand_gaff.mol2 -f mol2 -o ligand.frcmod
```

## 案例

### 3CA2 碳酸酐酶 MD 模拟

完整复现论文级别 MD 流程（ff19SB + OPC3 + Zn²⁺非键合模型 + AMS配体 + 150mM NaCl），详见：
- `3CA2_PreMD_Protocol.md` — tleap 前处理
- `3CA2_MD_Simulation_Protocol.md` — 完整七步 MD

## 参考资料

- AMBER 官网: https://ambermd.org
- AMBER 教程: https://ambermd.org/tutorials/
- VMD: https://www.ks.uiuc.edu/Research/vmd/
- PyMOL: https://pymol.org/
- Gaussian: https://gaussian.com/
- Chimera: https://www.cgl.ucsf.edu/chimera/

## 版本

- v2.0.0 (2026-04-07) — 新增 VMD/PyMOL/Gaussian09/Chimera 完整指南 + AI+MD 前沿
- v1.0.0 (2026-04-06) — AMBER 全工具链基础版

## 许可证

MIT License
