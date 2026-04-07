# Molecular Simulation 全栈技能 v2.0

> **版本**: 2.0.0 | **更新**: 2026-04-07
> **覆盖**: AMBER + VMD + PyMOL + Gaussian09 + Chimera + 高阶体系 + AI+MD

---

## 触发词（全面覆盖）

### 分子动力学理论类
- `分子动力学` / `MD模拟` / `MD` / `Molecular Dynamics`
- `第一性原理` / `势能面` / `力场` / `Verlet积分` / `系综`
- `显式溶剂` / `隐式溶剂` / `GBSA` / `PB`

### AMBER 工具类
- `AMBER` / `AmberTools` / `sander` / `pmemd` / `tleap`
- `antechamber` / `pdb4amber` / `cpptraj` / `parmed`
- `MMPBSA.py` / `pymembrane` / `Li-Merz` / `GAFF` / `ff19SB`

### 可视化软件类
- `VMD` / `PyMOL` / `Chimera` / `UCSF Chimera`
- `分子可视化` / `轨迹可视化` / `制作MD动画` / `发表级图片`
- `Ray tracing` / `渲染` / `RMSD图` / `结构叠图`

### 量子化学计算类
- `Gaussian` / `Gaussian09` / `Gaussian16` / `量子化学`
- `几何优化` / `频率计算` / `NMR计算` / `势能面`
- `RESP电荷` / `QM计算` / `基组` / `B3LYP` / `DFT`

### 力场与参数化类
- `GAFF` / `ff19SB` / `ff14SB` / `GLYCAM` / `ff99bsc0`
- `配体参数化` / `小分子力场` / `金属离子建模`
- `原子类型` / `电荷计算` / `AM1-BCC` / `RESP`

### 轨迹分析与后处理
- `轨迹分析` / `RMSD` / `RMSF` / `PCA` / `聚类分析`
- `氢键分析` / `距离监测` / `主成分分析`
- `MM-PBSA` / `MM-GBSA` / `结合自由能` / `GIST`

### 高阶复杂体系
- `金属蛋白` / `Zn离子建模` / `MCPB.py` / `12-6-4 LJ`
- `配体-蛋白复合物` / `药物设计` / `约束释放`
- `膜蛋白` / `POPC` / `POPE` / `膜构建`
- `糖类` / `GLYCAM` / `核酸` / `RNA MD` / `DNA MD`
- `QM/MM` / `增强采样` / `Metadynamics` / `伞形采样`

### AI + 分子动力学前沿
- `DeepMD` / `DeePMD-kit` / `机器学习力场`
- `ANI` / `分子生成` / `强化学习分子设计`

---

## 第一性原理

### MD 本质
```
F_i = -∇V(r)           ← 牛爵爷第二定律
m_i · a_i = F_i
→ 每步: 求力 → 更新速度 → 更新位置 (Verlet积分)
```

### 势能函数 = 力场
```
V = Σ k_b(r-r_eq)² + Σ k_θ(θ-θ_eq)²
  + Σ V_n/2[1+cos(nφ-γ)]
  + Σ 4ε[(σ/r)^12 - (σ/r)^6]
  + Σ q_iq_j/r
```
项依次为：键伸缩、角弯曲、二面角 torsions、范德华、静电

### 积分算法
| 算法 | 特点 |
|------|------|
| **Verlet** | 简洁、稳定性好 |
| **LEAPFROG** | AMBER 默认，速度快 |
| **RESPA** | 多尺度，可加速 |

### 系综
| 系综 | 参数 | 应用 |
|------|------|------|
| NVT | V,E固定 | 加热、平衡 |
| NPT | N,P,T固定 | 生产常用 |
| NVE | N,V,E固定 | 能量守恒检验 |

---

## AMBER 工具链速查

### pdb4amber — PDB 清洗
```bash
pdb4amber input.pdb -o output.pdb
# 功能：H加、统一命名、处理二硫键、移除杂质
```

### tleap — 系统构建
```bash
tleap -f input.in
# 流程：加载力场 → 加载分子 → 组装 → 加溶剂 → 加离子 → 保存
```

**tleap 典型脚本**:
```bash
source leaprc.protein.ff19SB      # 蛋白力场
source leaprc.water.opc3          # OPC3 水模型
loadamberparams frcmod.ionslm_126_opc3  # 离子参数

PRO = loadpdb protein_his.pdb
addIons PRO Na+ 0
solvateOct PRO OPCBOX 12.0
addIonsRand PRO Na+ 19 Cl- 19

saveamberparm PRO system.prmtop system.rst7
quit
```

### antechamber — 小分子参数化
```bash
# 标准流程: PDB → MOL2 → GAFF + BCC
obabel -ipdb ligand.pdb -omol2 -O ligand.mol2
antechamber -i ligand.mol2 -fi mol2 \
  -o ligand_gaff.mol2 -fo mol2 \
  -c bcc -s 2 -at amber -j 4 -nc 0
parmchk2 -i ligand_gaff.mol2 -f mol2 -o ligand.frcmod
```

### sander / pmemd — 能量最小化与 MD
```bash
# 最小化
pmemd.cuda -O -i min.in -o min.out -p prmtop -c rst7 -r min.rst7

# NPT 平衡/生产
pmemd.cuda -O -i md.in -o md.out -p prmtop -c rst7 -r md.rst7 -x md.nc
```

**标准 mdin 参数**:
```
imin = 0, nstlim = 500000, dt = 0.002
ntf = 2, ntc = 2, ntb = 2, ntp = 1
temp0 = 298.0, ntt = 3, gamma_ln = 1.0
cut = 8.0, iwrap = 1, ioutfm = 1
```

### cpptraj — 轨迹分析
```bash
cpptraj -p system.prmtop -i analysis.in
```

**常用分析命令**:
```
trajin md.nc
autoimage
rms reference :1-261@CA out rmsd_ca.dat
atomicfluct :1-261@CA out rmsf_ca.dat
distance :264@ZN :94@NE2 out zn_his94.dat
hbond avgout hbond.dat
principal :1-261 out principal.dat
```

### MMPBSA.py — 结合自由能
```bash
MMPBSA.py -i mmpbsa.in -o FINAL_RESULTS_MMPBSA -cp complex.prmtop -lp ligand.prmtop -rp receptor.prmtop -y md.nc
```

### parmed — 拓扑操作
```bash
parmed -p system.prmtop << EOF
HMR 1.5
outparm system_hmr.prmtop
quit
EOF
```

---

## VMD 完整指南

> **官网**: https://www.ks.uiuc.edu/Research/vmd/
> **功能**: 分子可视化、轨迹分析、发表级渲染

### 基础操作
```bash
# 启动 VMD
vmd
```

### 读取 AMBER 轨迹
```bash
# 方法1: 手动加载
# File → New Molecule → 加载 .prmtop 和 .nc/.dcd

# 方法2: 命令行
mol new system.prmtop type amberparm
animate read netcdf md.nc waitfor all
```

### 核心命令（Tcl 脚本）
```tcl
# 蛋白显示为 NewCartoon
mol addrep 0
set sel [atomselect 0 "protein and name CA"]
$sel set beta 0
display resize 800 600

# 轨迹播放
animate speed 0.5
animate forward

# RMSD 计算
package require mdff
set ref [atomselect 0 "protein and name CA" frame 0]
set sel [atomselect 0 "protein and name CA"]
set rmsd [measure rmsd $sel $ref]
puts "RMSD: $rmsd"
```

### 轨迹分析
```tcl
# 氢键分析
package require hbonds
set sel [atomselect 0 "protein and resid 10 to 20"]
hbonds -sel $sel -dist 3.0 -angle 120

# 距离监测
set id1 [atomselect 0 "resname ZN and resid 264"]
set id2 [atomselect 0 "resname AMS and resid 265 and name N1"]
measure distance $id1 $id2
```

### 渲染发表级图片
```tcl
# 开启光线追踪
display rendermode GLSL

# 保存图片
render TachyonInternal screenshot.png
# 或
render snapshot screenshot.tga
```

### 与 AMBER 联用
```bash
# VMD 可直接读取 AMBER 输出
# 1. .prmtop 作为结构文件
# 2. .nc (NetCDF) 或 .dcd 作为轨迹
# 3. .rst7 作为坐标文件

# 导出为 PDB 用于其他软件
set sel [atomselect 0 "all"]
$sel writepdb "output.pdb"
```

---

## PyMOL 完整指南

> **官网**: https://pymol.org/
> **功能**: 分子可视化、药物设计、发表级图片动画

### 基础操作
```bash
# 启动
pymol

# 命令行模式
pymol -c script.py
```

### 读取分子
```python
# PyMOL API
from pymol import cmd

cmd.load("protein.pdb")
cmd.load("ligand.mol2", "ligand")
cmd.load("md_traj.nc", "system")

# 叠图
cmd.align("protein", "reference", cycles=0)
```

### 常用命令
```python
# 显示模式
cmd.show("cartoon", "protein")
cmd.show("sticks", "ligand")
cmd.show("spheres", "resname ZN")

# 着色
cmd.color("red", "elem C")
cmd.spectrum("b", "blue_white_red", "protein and name CA")

# 标签
cmd.label("resn AMS", '"%s-%s" % (resn, resi)')

# 制作动画
cmd.mset("1x100")
cmd.mview("store", 1)
cmd.mview("store", 100)
cmd.movie(fps=30)
```

### Ray tracing 渲染
```python
cmd.ray(2400, 1800)
cmd.png("output.png", dpi=300, ray=1)
```

### 与 AMBER 联用
```python
# PyMOL 可读取 PDB/MOL2
# 用于可视化 tleap 输出结构
cmd.load("system.prmtop")  # 只读结构
cmd.load("3CA2_system.pdb")

# RMSD 图
cmd.align("chain A and name CA", "reference and name CA")
```

---

## Gaussian09 完整指南

> **官网**: https://gaussian.com/
> **功能**: 量子化学计算、分子轨道、势能面、RESP 电荷

### 输入文件格式
```bash
# 基本格式 (.gjf)
%mem=4GB
%nprocshared=8
%chk=calc.chk
#p B3LYP/6-31G(d) opt freq

标题行
空行
电荷 多重度
坐标

空行
```

### 常用计算类型

#### 1. 几何优化 (Opt)
```bash
#p B3LYP/6-31G(d) opt

Title
0 1
C  0.0  0.0  0.0
...

```

#### 2. 频率计算 (Freq)
```bash
#p B3LYP/6-31G(d) freq

Title
0 1
...
```

#### 3. NMR 计算
```bash
#p B3LYP/6-311G(d,p) nmr

Title
0 1
...
```

### 基组选择
| 基组 | 用途 |
|------|------|
| 6-31G(d) | 日常优化 |
| 6-311G(d,p) | 高精度 |
| aug-cc-pVDZ | 色散/氢键 |
| LANL2DZ | 过渡金属 |

### 与 AMBER 联用 — RESP 电荷计算
```bash
# 1. 用 Gaussian 计算 ESP
#p b3lyp/6-31G(d) pop=esp

# 2. 用 RED 得到 RESP 电荷
# 需要antechamber的mort：
antechamber -i mol2_file.mol2 -fi mol2 -o mol2_resp.mol2 -fo mol2 -c resp

# 3. 或者用Gaussian + Multiwfn
```

### 常见错误
| 错误 | 原因 | 解决 |
|------|------|------|
| `SCF didn't converge` | 电子结构未收敛 | 增加 DIIS 或改算法 |
| `Opt cycle failed` | 几何优化失败 | 改初始结构或基组 |
| `Basis not found` | 基组名错误 | 检查 Gaussian 文档 |

---

## Chimera 完整指南

> **官网**: https://www.cgl.ucsf.edu/chimera/
> **功能**: 分子可视化、结构建模、能量最小化、MD轨迹处理

### 基础操作
```bash
# 启动
chimera

# 命令行
chimera --nogui script.py
```

### 读取 AMBER 文件
```bash
# Chimera 可直接读取 PDB/MOL2
# 轨迹通过 MD Movie 插件读取

# 常用命令
open system.pdb
open ligand.mol2
```

### MD Movie（轨迹播放）
```bash
# 菜单: Tools → MD/Ensemble → MD Movie
# 支持 AMBER NetCDF 轨迹 (.nc)
```

### 结构编辑
```bash
# 常用命令
# 旋转/移动
turn x 45
move x 5

# 删原子
delete #0:ZN

# 原子着色
color byelement #0
```

### 与 AMBER 联用
```bash
# 1. Chimera 读取 tleap 输出结构
open 3CA2_system.pdb

# 2. 监测 Zn-配体距离
# Measure → Distance

# 3. 突变扫描
# Tools → Sequence → Multalign Viewer
```

---

## 高阶体系 SOP

### 金属蛋白 MD（Zn/Fe/Mg）

#### 非键合模型（推荐快速平衡）
```bash
# Li-Merz 12-6-4 电荷
loadamberparams frcmod.ionslm_1264_opc3

# Zn²⁺ 非键合参数自动从 frcmod 读取
# 关键: restraintmask 覆盖 Zn 配位残基
```

#### Bonded 模型（MCPB.py）
```bash
# MCPB.py 处理金属中心
MCPB.py -i input.pdb -s Zn -o bonded

# 生成新的 frcmod 包含 Zn-配体键角参数
```

#### Zn 距离监测
```bash
# cpptraj
distance :264@ZN :94@NE2 out zn_his94.dat
distance :264@ZN :96@NE2 out zn_his96.dat
# 预期: ~2.0 Å
```

### 配体-蛋白复合物 MD

#### 标准参数化流程
```bash
# 配体 PDB → MOL2 (obabel)
obabel -ipdb ligand.pdb -omol2 -O ligand.mol2

# GAFF + AM1-BCC
antechamber -i ligand.mol2 -fi mol2 \
  -o ligand_gaff.mol2 -fo mol2 -c bcc -s 2 -nc 0

parmchk2 -i ligand_gaff.mol2 -f mol2 -o ligand.frcmod
```

#### 约束释放策略
```
加热: 100 kcal/mol 约束 (1 ns)
NPT1: 100→10 kcal/mol (1 ns)
NPT2: 10→1 kcal/mol (1 ns)
NPT3: 1→0.1 kcal/mol (1 ns)
NPT4: 0.1→0 (1 ns)
生产: 无约束 (50-100 ns)
```

### 膜蛋白 MD

#### 膜构建
```bash
# 方法1: CHARMM-GUI Membrane Builder
# https://www.charmm-gui.org/

# 方法2: tleap 手动构建
source leaprc.lipid14
# 需要 D Membrane Builder 或 PACKMOL
```

#### 膜蛋白取向
```bash
# 检查膜蛋白是否有合理的插入方向
# 使用 OPM 数据库验证
# https://opm.phar.umich.edu/
```

### 增强采样

#### Metadynamics
```bash
# PLUMED 插件
plumed driver --plumed plumed.dat --mf_xtc prod.nc
```

#### 伞形采样
```bash
# 每个窗口运行
pmemd.cuda -O -i umbrella_1.in -o umbrella_1.out \
  -p system.prmtop -c window_0.rst7 -r window_1.rst7

# WHAM 分析
wham -100 100 50 0.01 300 metabias.dat
```

---

## AI + MD 前沿

### 机器学习力场

#### DeepMD
```bash
# DeePMD-kit 训练
dp train input.json

# freeze 模型
dp freeze -o graph.pb

# lammps 运行
dpgen lammps confs
```

#### ANI
```bash
# ANI-1ccx
# 可直接作为 AMBER 的 QM/MM 边界处理
```

### 分子生成与强化学习
```python
# 药物设计
# 常用框架: RDKit + 强化学习

# 示例: 生成类药分子
from rdkit import Chem
from rdkit.Chem import Descriptors

# 结合 MD 采样做引导
# 关键: FEP (Free Energy Perturbation)
```

### MSM (Markov State Models)
```bash
# pyemma 使用
import pyemma

# 聚类
features = pyemma.coordinates.featurizer(prmtop)
features.add_all()
data = pyemma.coordinates.load(traj, features=features)

# K-means 聚类
cluster = pyemma.coordinates.cluster_kmeans(data, k=100)

# 构建 MSM
msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=100)
```

---

## 轨迹分析速查表

| 分析目标 | 工具/命令 | 关键参数 |
|---------|----------|---------|
| RMSD | cpptraj: `rms reference` | :CA, backbone |
| RMSF | cpptraj: `atomicfluct` | :CA |
| 距离监测 | cpptraj: `distance` | 原子选择符 |
| 氢键 | cpptraj: `hbond` | dist=3.0, angle=120 |
| 主成分 | cpptraj: `principal` |  |
| 能量 | `process_mdout.perl` |  |
| 结合自由能 | MMPBSA.py |  |
| 密度 | cpptraj + pytraj |  |

---

## 参考资料

### 官方资源
- AMBER官网: https://ambermd.org
- AMBER教程: https://ambermd.org/tutorials/
- VMD官网: https://www.ks.uiuc.edu/Research/vmd/
- PyMOL官网: https://pymol.org/
- Gaussian官网: https://gaussian.com/
- Chimera官网: https://www.cgl.ucsf.edu/chimera/

### 本地文档
- `3CA2_MD_Simulation_Protocol.md` — 3CA2 碳酸酐酶完整MD模拟流程
- `3CA2_PreMD_Protocol.md` — tleap前处理详细流程
- `integrated_knowledge.md` — 整合知识库（Agent-3产出）

### 核心参考文献
- Case et al., J. Chem. Inf. Model. (2024) — ff19SB力场
- Li et al., JCTC (2020) — Li-Merz 12-6-4离子参数
- Maier et al., J. Chem. Theory Comput. (2015) — OPC水模型
