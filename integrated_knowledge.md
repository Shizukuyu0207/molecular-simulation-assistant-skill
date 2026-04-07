# 分子动力学全栈知识库（整合压缩版）

> 版本: 1.0.0 | 整合日期: 2026-04-07
> 基于: AMBER官方文档 + 3CA2实战协议 + 主流软件手册

---

## 一、MD 第一性原理

### 核心方程
```
F = ma  →  Verlet积分: r(t+dt) = 2r(t) - r(t-dt) + a(t)dt²
```

### 势能函数（五项）
| 项 | 物理含义 | AMBER力场 |
|----|---------|-----------|
| 键伸缩 | k_b(r-r_eq)² | 简谐 |
| 角弯曲 | k_θ(θ-θ_eq)² | 简谐 |
| 二面角 | V_n/2[1+cos(nφ-γ)] | Fourier |
| 范德华 | 4ε[(σ/r)^12-(σ/r)^6] | LJ 12-6 |
| 静电 | q_iq_j/r | Coulomb |

### 关键参数
- dt = 2 fs（SHAKE约束H）
- ntc=2, ntf=2（约束H键）
- ntt=3（Langevin热浴）
- ntp=1（NPT等压）

---

## 二、AMBER 全工具链速查卡

### pdb4amber
```
pdb4amber raw.pdb -o clean.pdb
作用: 清洗/加H/二硫键/命名标准化
```

### tleap
```
tleap -f build.in
流程: source力场 → loadmol2 → loadpdb → bond → addIons → solvate → saveamberparm
```

### antechamber
```
antechamber -i lig.mol2 -fi mol2 -o lig_gaff.mol2 -fo mol2 -c bcc -s 2 -nc 0
parmchk2 -i lig_gaff.mol2 -f mol2 -o lig.frcmod
作用: 小分子GAFF参数化+电荷
```

### sander/pmemd
```
pmemd.cuda -O -i md.in -o md.out -p prmtop -c rst7 -r rst7 -x nc
imin=0(nvt/npt) | imin=1(最小化)
```

### cpptraj
```
trajin md.nc
autoimage
rms reference :1-261@CA out rmsd.dat
atomicfluct :1-261@CA out rmsf.dat
distance :264@ZN :94@NE2 out dist.dat
hbond avgout hb.dat
principal out pc.dat
```

### MMPBSA.py
```
MMPBSA.py -i mmpbsa.in -cp complex.prmtop -lp lig.prmtop -y md.nc
```

### parmed
```
HMR 1.5  # 氢质量重分配
outparm system_hmr.prmtop
```

---

## 三、各软件速查表

### VMD
```
mol new prmtop type amberparm
animate read netcdf md.nc
# Tcl: atomselect, measure rmsd, hbonds
render TachyonInternal img.png
```

### PyMOL
```
cmd.load("protein.pdb")
cmd.load("traj.nc")
cmd.align("target", "reference")
cmd.ray(2400,1800); cmd.png("out.png")
```

### Gaussian09
```
%mem=4GB %nproc=8
#p B3LYP/6-31G(d) opt freq
# 标题
0 1
坐标...
```

### Chimera
```
open system.pdb
# MD Movie: Tools → MD/Ensemble → MD Movie
# Measure → Distance
```

---

## 四、高阶体系处理 SOP

### 金属蛋白（Zn²⁺）
```
非键合: loadamberparams frcmod.ionslm_1264_opc3
约束: restraintmask 覆盖Zn配位残基
监测: distance :264@ZN :94@NE2  预期~2.0Å
```

### 配体-蛋白
```
obabel → antechamber GAFF → tleap加载
平衡: 100→10→1→0.1→0 kcal/mol约束
分析: RMSD(CA) + 配体RMSD + 氢键
```

### 膜蛋白
```
力场: ff19SB + CHARMM36m/POPC
膜构建: CHARMM-GUI 或 D Membrane Builder
取向: OPM数据库验证
```

### 核酸
```
力场: ff99bsc0 / ff14bsc0
离子: Mg²⁺用12-6-4 Li-Merz
分析: 凹槽宽度/盐桥/二级结构
```

---

## 五、前沿技术要点

### 机器学习力场
- DeepMD: 训练→freeze→lammps运行
- ANI-1ccx: 高精度NN势能面

### 增强采样
- Metadynamics: PLUMED驱动
- 伞形采样 + WHAM分析
- FEP: 结合自由能精确计算

### MSM
- pyemma: 聚类→构建MSM→动力学分析

---

## 六、实战案例：3CA2碳酸酐酶MD流程

### 前处理（tleap）
```
Step1: awk提取蛋白/配体/Zn
Step2: obabel AMS→MOL2→GAFF+Hg自定义参数
Step3: pdb4amber→His质子化(HIE/HID)
Step4: tleap组装(ff19SB+OPC3+12-6-4离子)
```

### MD七步法
```
min1(100约束) → heat(1ns) → eq1-4(逐级释放) → prod(50ns+)
```

### 关键监测指标
```
CA RMSD < 3Å | Zn-配位 ~2.0Å | 密度 ~1.0g/cm³ | 温度298K±5K
```
