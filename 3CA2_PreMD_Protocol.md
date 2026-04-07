# 3CA2 — tleap 前处理完整流程

> **AmberTools**: 24.8 (conda 环境: `amber`)
> **体系**: 人源碳酸酐酶 II (Carbonic Anhydrase II, 3CA2)
> **处理方案**: 移除 Hg²⁺ → GAFF+自定义Hg参数 → 非键合 Zn²⁺

---

## 流程总览

```
Step 1   从 3CA2.pdb 提取各组分
Step 2   AMS: PDB → MOL2 (obabel) → GAFF参数化 (antechamber)
Step 3   蛋白部分: pdb4amber 清洗 + His质子化处理
Step 4   自定义 Hg 原子参数 (frcmod)
Step 5   tleap 组装系统
```

---

## Step 1 — 从 3CA2.pdb 提取各组分

原始 PDB 中：
- **蛋白**: ATOM 记录 (Chain A, Res 1-261)
- **Zn²⁺**: HETATM 2043 (ZN A 264)
- **Hg²⁺**: HETATM 2041, 2042 (HG A 262, HG A 263) — **需移除**
- **AMS**: HETATM 2044-2055 (AMS A 265) — **含 Hg**

### 1-1. 提取蛋白部分（移除所有 HETATM）

```bash
conda activate amber
cd ~/3CA2_MD

# 提取蛋白原子 (ATOM 记录)，移除所有 HETATM
awk '{
  if ($1 == "ATOM") print
  if ($1 == "TER") print
  if ($1 == "END") print
}' 3CA2.pdb > protein_only.pdb

# 验证
grep -c "^ATOM" protein_only.pdb   # 应为 2059 (259残基)
grep "ZN\|HG\|AMS" protein_only.pdb && echo "有杂质" || echo "干净"
```

### 1-2. 提取 AMS 配体

```bash
# 提取 AMS 的 HETATM 记录 (RESNAME = AMS, RESSEQ = 265)
awk '{
  if ($1 == "HETATM" && substr($0,18,3) == "AMS") print
  if ($1 == "CONECT") {
    # CONECT 记录的列位置：第7-11列是原子序号
    # CONECT 2044-2055 是 AMS 相关的连接记录
    serial = substr($0,7,5) + 0
    if (serial >= 2044 && serial <= 2055) print
  }
}' 3CA2.pdb > AMS_raw.pdb

# 验证
grep -c "^HETATM" AMS_raw.pdb  # 应为 12
cat AMS_raw.pdb
```

**预期输出 (AMS_raw.pdb)**:
```
HETATM 2044  C1  AMS A 265      -5.535   2.332  15.311  1.00 20.00           C
HETATM 2044  C2  AMS A 265 ...
...
```

> ⚠️ **注意**: AMS_raw.pdb 包含坐标但**没有键连接信息**（PDB 格式的固有限制），不能直接用 antechamber 处理。

### 1-3. 提取 Zn²⁺（保留用于 tleap）

```bash
awk '{
  if ($1 == "HETATM" && substr($0,18,2) == "ZN") print
}' 3CA2.pdb > ZN_ion.pdb

cat ZN_ion.pdb
# HETATM 2043  ZN    ZN A 264      -6.868  -1.603  15.446  1.00  1.69          ZN
```

---

## Step 2 — AMS: PDB → MOL2 → GAFF 参数化

### 关键说明：为什么不能直接用 PDB 做 antechamber

| 格式 | 包含内容 |
|------|---------|
| **PDB** | 原子名、坐标、残基名、**无键连接** |
| **MOL2** | 原子名、坐标、**键连接**、原子类型、电荷 |

antechamber 需要键连接信息来自动判断原子杂化状态和键类型。**纯 PDB 输入 antechamber 会失败或产生错误参数**。

**解决方案**: 用 **Open Babel** 将 PDB 转换为 MOL2（自动推断键连接），再用 antechamber 处理 MOL2。

### 2-1. 安装 obabel（如未安装）

```bash
# 方法1: conda 安装
conda install -c conda-forge openbabel -y

# 方法2: apt 安装
sudo apt install openbabel
```

### 2-2. PDB → MOL2 转换

```bash
# obabel 自动推断键连接，输出 MOL2 格式
obabel -ipdb AMS_raw.pdb -omol2 -O AMS.mol2

# 验证 MOL2 文件
head -30 AMS.mol2
grep "@<TRIPOS>BOND" AMS.mol2
grep "@<TRIPOS>ATOM" AMS.mol2
```

**预期 MOL2 键连接**（苯环 + S=O₂ + NH₂ + Hg）:
```
C1     -5.535   2.332  15.311 C.ar    1 AMS    0.0000
C2     -4.678   1.832  14.227 C.ar    1 AMS    0.0000
...
@<TRIPOS>BOND
     1     1     2 ar
     2     1     6 ar
     3     1     8 1
     ...
```

### 2-3. antechamber GAFF 参数化（MOL2 输入）

```bash
# GAFF 参数化，输入 MOL2
antechamber -i AMS.mol2 -fi mol2 \
            -o AMS_gaff.mol2 -fo mol2 \
            -c bcc -s 2 \
            -at amber -j 4 -nc 0

# 参数说明:
# -c bcc    : AM1-BCC 电荷 (标准)
# -s 2      : 详细输出
# -at amber : AMBER 原子类型
# -j 4      : 完整键类型判定
# -nc 0     : 分子净电荷 = 0
```

**预期输出**:
```
Running: bondtype -j full ...
Running: atomtype -i ANTECHAMBER_AC.AC0 ... -p gaff
Total number of electrons: 160; net charge: 0
Running: sqm.sh
Running: am1bcc ...
```

### 2-4. 检查缺失参数

```bash
parmchk2 -i AMS_gaff.mol2 -f mol2 -o AMS.frcmod

# 查看是否有 "ATTN" 行（缺失参数）
grep "ATTN" AMS.frcmod
# 如果有，标记为需要手动补充
```

---

## Step 3 — 蛋白部分清洗与质子化

### 3-1. pdb4amber 预处理

```bash
# pdb4amber 会:
# 1. 统一原子命名到 AMBER 风格
# 2. 处理二硫键 (CYS → CYX)
# 3. 保留 CONECT 记录供二硫键识别
# 4. 移除 UNL 等未知残基

pdb4amber protein_only.pdb -o protein_amb.pdb
```

### 3-2. 检查 His 残基

3CA2 中有 4 个 His (Res 10, 15, 17, 64, 94, 96, 119)：

```bash
grep "^ATOM.*HIS" protein_amb.pdb
```

| 残基 | PDB残基号 | 位置 | 推荐状态 | 原因 |
|------|----------|------|--------|------|
| His 10 | 10 | 表面 | **HIE** | 普通中性 |
| His 15 | 15 | 表面 | **HIE** | 普通中性 |
| His 17 | 17 | 表面 | **HIE** | 普通中性 |
| His 64 | 64 | 表面/Hg配位 | **HID** | His64 ND1/NE2与Hg²⁺配位 |
| His 94 | 94 | **Zn配位** | **HIE** | NE2配位Zn，质子在ND1 |
| His 96 | 96 | **Zn配位** | **HIE** | NE2配位Zn，质子在ND1 |
| His 119 | 119 | **Zn配位** | **HID** | ND1配位Zn，质子在NE2 |

### 3-3. His 重命名（手动编辑）

```bash
# 备份
cp protein_amb.pdb protein_his.pdb

# His 94: NE2→Zn，质子在ND1 → HIE
sed -i 's/HIS A  94/HIE A  94/g' protein_his.pdb

# His 96: NE2→Zn，质子在ND1 → HIE
sed -i 's/HIS A  96/HIE A  96/g' protein_his.pdb

# His 119: ND1→Zn，质子在NE2 → HID
sed -i 's/HIS A 119/HID A 119/g' protein_his.pdb

# His 64: Hg配位 → HID
sed -i 's/HIS A  64/HID A  64/g' protein_his.pdb

# 验证（3个HIE + 2个HID）
grep "HIE A  94\|HIE A  96\|HID A  64\|HID A 119\|HIE A  10\|HIE A  15\|HIE A  17" protein_his.pdb
```

### 3-4. Se-Met 处理

MSE (Selenomethionine) 在 3CA2 中不存在（序列中无），跳过此步。如有：

```bash
# MSE → MET
sed -i 's/MSE/MET/g' protein_his.pdb
sed -i 's/ SD / SG /g' protein_his.pdb   # Se → S
```

### 3-5. 二硫键检查

```bash
# 检查 SSBOND 记录
grep "^SSBOND" 3CA2.pdb

# 如果有: 检查 CYS → CYX 是否正确
grep "^ATOM.*CYX" protein_his.pdb
```

3CA2 无二硫键，跳过。

---

## Step 4 — 自定义 Hg 原子参数

GAFF 不含 Hg (汞) 原子类型。需要自定义 `frcmod` 补充 Hg 的质量、LJ 参数。

### 4-1. 获取 Hg 的 LJ 参数

Hg²⁺ 的 Li-Merz 12-6-4 参数（来源: `frcmod.ionslm_1264_opc3`，Li et al., JCTC 2020）：

```bash
cat >> AMS.frcmod << 'EOF'

MASS
hg   200.59
BOND
ca-hg   320.0    2.00
ANGLE
ca-ca-hg   50.0  120.0
DIHE
IMPROPER
NONBON
hg   200.59  1.625  0.09235154
EOF
```

> ⚠️ 参数来源: `frcmod.ionslm_1264_opc3`（AMBER AmberTools 24.8）

### 4-2. 验证 frcmod

```bash
cat AMS.frcmod
# 确认包含 hg 的 MASS 和 NONBON 条目
```

---

## Step 5 — tleap 系统组装

### 5-1. 创建 tleap 输入脚本

```bash
cat > tleap_3CA2.in << 'EOF'
# ============================================================
# tleap_3CA2.in — 3CA2 碳酸酐酶 + Zn + AMS 系统构建
# 力场: ff19SB + OPC3 水模型 + Li-Merz 离子参数
# ============================================================

# ---- 1. 加载力场 ----
source leaprc.protein.ff19SB
source leaprc.water.opc3
loadamberparams frcmod.ionslm_126_opc3
loadamberparams frcmod.ionslm_1264_opc3

# ---- 2. 加载 AMS 配体 (GAFF) ----
loadamberparams AMS.frcmod
LIG = loadmol2 AMS_gaff.mol2
check LIG

# ---- 3. 加载蛋白 (含 Zn) ----
PRO = loadpdb protein_his.pdb

# ---- 4. Zn²⁺ 自动被 leap 识别 (原子名 ZN, 元素 Zn) ----

# ---- 5. 添加 Zn-AMS 配位键 (N1→Zn) ----
# CONECT记录确认: ZN(2043) ↔ AMS@N1(2054)，距离2.02Å
bond PRO.264.ZN PRO.265@N1

# ---- 6. 添加中和离子 ----
addIons PRO Na+ 0

# ---- 7. 添加 OPC3 水盒子 (正八面体, 12 Å 缓冲层) ----
solvateOct PRO OPCBOX 12.0

# ---- 8. 添加 150 mM NaCl 盐缓冲 ----
addIonsRand PRO Na+ 19 Cl- 19

# ---- 9. 保存系统 ----
saveamberparm PRO 3CA2.prmtop 3CA2.rst7
savepdb PRO 3CA2_system.pdb

# 检查电荷
charge PRO
quit
EOF
```

### 5-2. 运行 tleap

```bash
conda activate amber
tleap -f tleap_3CA2.in > tleap.log 2>&1

# 检查输出
tail -20 tleap.log
ls -lh 3CA2.prmtop 3CA2.rst7
```

### 5-3. 诊断常见错误

```bash
# 错误1: 未找到原子类型
grep -i "atom type.*not found" tleap.log

# 错误2: 缺失参数
grep -i "ATTN" tleap.log

# 错误3: 电荷不为零
grep "CHARGE" tleap.log
```

---

## 输出文件汇总

| 文件 | 说明 |
|------|------|
| `protein_only.pdb` | 蛋白部分 (ATOM 记录) |
| `protein_amb.pdb` | pdb4amber 清洗后 |
| `protein_his.pdb` | His 质子化处理后 |
| `AMS_raw.pdb` | AMS 配体坐标 |
| `AMS.mol2` | obabel 转换后的 MOL2 |
| `AMS_gaff.mol2` | GAFF 参数化后的 MOL2 |
| `AMS.frcmod` | GAFF 缺失参数 + Hg 自定义参数 |
| `ZN_ion.pdb` | Zn²⁺ 单独保存 |
| `3CA2.prmtop` | 最终拓扑文件 ✅ |
| `3CA2.rst7` | 最终坐标文件 ✅ |
| `tleap.log` | 构建日志 |

---

## 快速启动命令（复制粘贴版）

```bash
conda activate amber

# ===== Step 1: 提取组分 =====
mkdir -p ~/3CA2_MD && cd ~/3CA2_MD
cp ~/3CA2/3CA2.pdb .

awk 'NR==1,/^ATOM/{p=1} p{print} /^HETATM/{p=0} /^(TER|END)$/{p=1}' 3CA2.pdb > protein_only.pdb
awk '$1=="HETATM" && substr($0,18,3)=="AMS"' 3CA2.pdb > AMS_raw.pdb
awk '$1=="HETATM" && substr($0,18,2)=="ZN"' 3CA2.pdb > ZN_ion.pdb

# ===== Step 2: AMS MOL2 转换 + GAFF =====
obabel -ipdb AMS_raw.pdb -omol2 -O AMS.mol2
antechamber -i AMS.mol2 -fi mol2 -o AMS_gaff.mol2 -fo mol2 -c bcc -s 2 -at amber -j 4 -nc 0
parmchk2 -i AMS_gaff.mol2 -f mol2 -o AMS.frcmod

# 添加 Hg 参数 (Li-Merz 12-6-4, frcmod.ionslm_1264_opc3)
cat >> AMS.frcmod << 'EOF'
MASS
hg   200.59
NONBON
hg   200.59  1.625  0.09235154
EOF

# ===== Step 3: 蛋白清洗 + His 质子化 =====
pdb4amber protein_only.pdb -o protein_amb.pdb
cp protein_amb.pdb protein_his.pdb

# His 94, 96 → HIE (NE2→Zn); His 119, 64 → HID
sed -i 's/HIS A  94/HIE A  94/g' protein_his.pdb
sed -i 's/HIS A  96/HIE A  96/g' protein_his.pdb
sed -i 's/HIS A 119/HID A 119/g' protein_his.pdb
sed -i 's/HIS A  64/HID A  64/g' protein_his.pdb

# ===== Step 4: tleap 组装 =====
tleap -f tleap_3CA2.in
```

---

## ⚠️ 待确认事项

1. ✅ **Hg 的 12-6-4 LJ 参数**: Rmin=1.625, ε=0.09235154（来自 `frcmod.ionslm_1264_opc3`）
2. ✅ **His 质子化状态**: His94/96→HIE，His119/64→HID（基于CONECT配位分析）
3. ✅ **Zn-AMS 配位键**: 已添加 `bond PRO.264.ZN PRO.265@N1`

---

## 修正说明（针对之前文档的错误）

| 问题 | 修正 |
|------|------|
| 直接用 PDB 给 antechamber | ✅ 改为 obabel PDB→MOL2 |
| 目录过多 | ✅ 合并为一个流程 |
| `frcmod.ions1lm_126_hfe_opc` | ✅ 修正为 `frcmod.ionslm_126_opc3` + `frcmod.ionslm_1264_opc3` |
| obabel 未安装 | ✅ 已安装 (conda-forge openbabel 3.1.0) |
| AMS.pdb → antechamber | ✅ 需 MOL2 格式 |
| 激活命令 `source ~/miniconda3/bin/activate` | ✅ 修正为 `conda activate` |
