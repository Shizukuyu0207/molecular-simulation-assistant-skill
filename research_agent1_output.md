# AMBER 教程与 MD 基础知识调研报告

> 调研日期: 2026/04/07
> 数据来源: 本地知识库（Web访问受限，基于已有知识编译）

---

## 一、AMBER 官网教程体系

### 1.1 基础教程 (Basic Tutorials)

#### Tutorial 1: 蛋白体系初级教程
**主题**: 牛心细胞色素C (Crambin) - 最小蛋白体系
**内容**:
- pdb4amber系统清洗
- tleap系统构建(ff19SB力场)
- 能量最小化 → NVT加热 → NPT平衡 → 成品MD
- cpptraj基础分析(RMSD、能量监测)

#### Tutorial 2: 核酸体系 (DNA/RNA)
**主题**: DNA双螺旋结构
**内容**:
- Nucleic Acid Builder构建DNA
- ff19SB力场参数
- Na+/Mg2+离子环境处理

#### Tutorial 3: 糖类 (Glycan) 教程
**主题**: 糖复合物处理
**内容**:
- GLYCAM力场加载
- 糖基化位点识别

---

### 1.2 高级教程 (Advanced Tutorials)

#### Tutorial 10: QM/MM联用
**主题**: 金属蛋白/酶催化体系
**内容**:
- sander QM/MM接口
- ONIOM型分层计算
- 活性位点量子力学处理

#### Tutorial 11: 配体参数化
**主题**: 小分子GAFF参数化
**流程**:
```bash
reduce ligand.pdb > ligand_h.pdb
antechamber -i ligand_h.pdb -fi pdb -o ligand.mol2 -fo mol2 -c bcc
parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod
```

#### Tutorial 12: 金属蛋白
**主题**: Zn2+、Fe2+/Fe3+、Cu2+等金属离子
**方案**:
- 非键合模型: 12-6/12-6-4 LJ + 点电荷
- 键合模型: MCPB.py + QM计算
- Li-Merz离子参数

#### Tutorial 13: 膜蛋白
**主题**: 脂质双层体系
**内容**:
- OPM数据库获取膜定位
- CHARMM-GUI膜构建
- 各向异性压控 (ntp=2)

---

## 二、MD 基础知识

### 2.1 理论基石

#### 核心方程
```
F_i = m_i · a_i = -∇V(r_1, r_2, ..., r_N)
```

#### 势能函数 (力场数学形式)
V = Σk_b(r-r_eq)^2 + Σk_θ(θ-θ_eq)^2 + Σ(V_n/2)[1+cos(nφ-γ)] + Σ4ε[(σ/r)^12-(σ/r)^6] + Σq_iq_j/r

**四项分解**:
| 项 | 类型 | 物理含义 |
|----|------|---------|
| Bond stretch | 键合 | 共价键弹性 |
| Angle bend | 键合 | 三原子角 |
| Dihedral | 键合 | 四原子绕键旋转 |
| LJ + Coulomb | 非键合 | 范德华 + 静电 |

---

### 2.2 力场函数

#### AMBER力场进化
ff99 → ff99SB → ff10 → ff12SB → ff14SB → ff19SB(2019)

#### GAFF (General Amber Force Field)
**定位**: 有机小分子（非标准蛋白/核酸）
**参数化流程**: 原始PDB → reduce(加H) → antechamber → GAFF原子类型 + AM1-BCC/RESP电荷

#### 水模型对比
| 模型 | 位点 | 推荐 |
|------|------|------|
| OPC3 | 3 | 首选 |
| OPC | 4 | 可选 |
| TIP3P | 3 | 老系统 |

#### 推荐组合
| 组件 | 推荐 |
|------|------|
| 蛋白力场 | ff19SB |
| 水模型 | OPC3 |
| 离子参数 | Li-Merz 12-6 |

---

### 2.3 积分算法

#### Verlet算法
r(t+Δt) = 2r(t) - r(t-Δt) + a(t)Δt²

#### Leapfrog算法
r(t+Δt) = r(t) + v(t+Δt/2)Δt; v(t+Δt/2) = v(t-Δt/2) + a(t)Δt

#### RESPA
- 内层: 短程力 (dt=2fs)
- 外层: 长程力 (dt=6-8fs)

---

### 2.4 系综

| 系综 | 条件 | MD实现 |
|------|------|--------|
| NVE | 孤立系统 | 无温压控 |
| NVT | 恒T | +Thermostat |
| NPT | 恒T恒P | +Thermostat+Barostat |

#### 恒温器 (ntt)
| ntt | 方法 | 评价 |
|-----|------|------|
| 0 | 无(NVE) | 能量守恒 |
| 1 | Berendsen | 不产生正确分布 |
| 3 | Langevin | 推荐 |
| 5 | Bussi-Gentile | 正确分布+低噪声 |

#### 恒压器 (ntp/barostat)
| ntp | 方法 |
|-----|------|
| 0 | 无(NVT) |
| 1 | 各向同性 |
| 2 | 各向异性(膜) |

---

### 2.5 溶剂模型

#### 显式溶剂
- OPC3/OPC/TIP3P水模型
- 周期性边界条件
- 计算量大但精度高

#### 隐式溶剂
- Generalized Born (GB) 模型
- 适合快速筛选

---

### 2.6 边界条件与截断

PME (Particle Mesh Ewald): O(NlogN)，实空间(短程)+倒空间(长程)分离
```
ntb=2, cut=8-10 Å, ewald=1
```

---

## 三、AMBER 工具链用法

### 3.1 pdb4amber
```bash
pdb4amber input.pdb -o output.pdb
```
**功能**: PDB清洗、H加、二硫键处理、原子命名标准化

### 3.2 tleap
```bash
tleap -f input.in > output.log
```
**核心命令**:
```bash
source leaprc.protein.ff19SB
source leaprc.water.opc3
loadamberparams frcmod.ionslm_126_opc3
solvateOct PRO OPCBOX 12.0
addIonsRand PRO Na+ 19 Cl- 19
saveamberparm PRO system.prmtop system.rst7
```

### 3.3 antechamber
```bash
antechamber -i ligand.mol2 -fi mol2 -o ligand_gaff.mol2 -fo mol2 -c bcc -s 2
parmchk2 -i ligand_gaff.mol2 -f mol2 -o ligand.frcmod
```

### 3.4 sander/pmemd

**核心mdin参数**:
```fortran
&cntrl
 imin=0, ntf=2, ntc=2, dt=0.002,
 ntb=2, ntt=3, gamma_ln=2.0,
 ntp=1, barostat=2, temp0=300.0, cut=8.0,
 ntwx=5000, ntpr=5000
/
```

**标准MD流程 (7步法)**:
```
pdb4amber → tleap → 最小化 → 加热(NVT) → 平衡(NPT) → 生产 → 分析
```

### 3.5 cpptraj

**RMSD分析**:
```bash
trajin prod.nc
reference initial.pdb
autoimage
rms reference :1-260@CA out rmsd_ca.dat
```

**RMSF分析**:
```bash
atomicfluct :1-260@CA out rmsf_ca.dat
```

**PCA分析**:
```bash
matrix covname myMatrix :1-260@CA out covmat.dat
diagmatrix cvmat myMatrix out evecs.dat vecs 10
projection modes 1-3 :1-260@CA out pc12.dat
```

### 3.6 MMPBSA.py
```bash
MMPBSA.py -O -i mmpbsa.in -o results.dat -sp complex.prmtop -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop
```

### 3.7 parmed
```bash
parmed prmtop
HMassRepartition  # 氢质量重映射
```

---

## 四、金属蛋白建模

### 4.1 Zn2+ 建模策略

| 模型 | 精度 | 复杂度 |
|------|------|--------|
| Nonbonded (12-6-4 LJ) | 中 | 低(首选) |
| Bonded (MCPB.py) | 高 | 高 |

### 4.2 Li-Merz 12-6-4 参数 (OPC3)
| 离子 | 质量 | Rmin | epsilon |
|------|------|------|---------|
| Zn2+ | 65.4 | 1.441 | 0.02343735 |

---

## 五、质子化状态

| 形式 | Amber名称 | 适用场景 |
|------|----------|---------|
| HIE | ε-氮质子化 | 普通中性His |
| HID | δ-氮质子化 | 与金属配位时 |
| HIP | 两者质子化 | 催化His/带正电 |

---

## 六、常见报错与解决方案

| 问题 | 原因 | 解决方案 |
|------|------|---------|
| RMSD持续增长 | 约束释放过快 | 增加eq步骤 |
| Zn-配体距离>3A | 配位键断裂 | 添加bond |
| 能量爆炸 | dt过大/初始冲突 | dt=0.001 |
| 密度过低 | 盒子太小 | 检查OPCBOX |

---

## 七、关键监测指标

| 指标 | 目标 | 工具 |
|------|------|------|
| CA RMSD | <3A | cpptraj |
| Zn-配体距离 | ~2.0A | cpptraj |
| 系统密度 | ~1.0g/cm3 | process_mdout |

---

## 八、参考文献

- AMBER官网: https://ambermd.org
- 官方教程: https://ambermd.org/tutorials/
- AMBER 24 Manual: https://ambermd.org/doc24/
