# 3CA2 — MD 模拟完整流程

> **AmberTools**: 24.8 (conda 环境: `amber`)
> **体系**: 人源碳酸酐酶 II (Carbonic Anhydrase II, 3CA2)
> **体系组成**: ff19SB + OPC3 + Zn²⁺ (nonbonded) + AMS (GAFF+Hg) + 150mM NaCl
> **GPU**: CUDA 加速 (pmemd.cuda)

---

## 流程总览

```
Step 0    tleap 系统构建 (3CA2_PreMD_Protocol.md)
Step 1    溶剂最小化 (100 kcal/mol约束)
Step 2    加热 0→298K (1 ns, 100 kcal/mol约束)
Step 3    NPT 298K (1 ns, 100→10 kcal/mol约束)
Step 4    NPT 298K (1 ns, 10→1 kcal/mol约束)
Step 5    NPT 298K (1 ns, 1→0.1 kcal/mol约束)
Step 6    NPT 298K (1 ns, 0.1→无约束)
Step 7    成品MD (NVT/NPT, 298K, 50-100 ns+)
```

---

## 通用参数约定

| 参数 | 值 | 说明 |
|------|-----|------|
| 时间步长 | 2 fs | dt = 0.002 ps |
| 约束常数 | 100/10/1/0.1 kcal/mol·Å² | 逐步释放 |
| 温度 | 298 K | |
| 压力 | 1 atm (NPT阶段) | |
| 电荷中和 | Na⁺ | addIons |
| 盐浓度 | 150 mM | addIonsRand |
| 水模型 | OPC3 | 正八面体盒子 |

---

## mdin 文件通用模板

```bash
# 最小化mdin模板
minimize = 1,           # 开启最小化
maxcyc = 10000,         # 最大循环数
ncyc = 5000,            # 前ncyc步使用最陡下降
ntpr = 100,             # 每ntpr步输出能量
ntr = 1,                # 开启位置约束
restraintmask = ':1-261@CA',  # CA原子约束
restraint_wt = 100,     # 约束力常數 (kcal/mol)

# MDmdin模板
ntpr = 10000,           # 每10000步(~20ps)输出能量
ntwx = 5000,            # 每5000步(~10ps)输出轨迹
ntwr = 50000,           # 每50000步重启文件
ntwprt = 0,             # 输出所有原子

imin = 0,               # MD模拟
dt = 0.002,             # 2 fs步长
nstlim = 500000,        # 500000步 = 1 ns (@2fs+SHAKE)
ntf = 2,                # 不计算受约束的键
ntc = 2,                # 启用SHAKE
temp0 = 298.0,          # 目标温度
ntt = 3,                # Langevin thermostat
gamma_ln = 1.0,         # 碰撞频率
ig = -1,                # 随机种子
ntr = 1,                # 开启约束
restraintmask = ':1-261@CA',
restraint_wt = 100,
iwrap = 1,              # 成像轨迹
ioutfm = 1,             # NetCDF格式
```

---

## Step 1 — 溶剂最小化

**目标**: 消除溶剂碰撞产生的极端冲突

```bash
conda activate amber

cat > min1.in << 'EOF'
Minimization — 溶剂 (100 kcal/mol CA约束)
 &cntrl
  imin = 1, maxcyc = 10000, ncyc = 5000,
  ntpr = 100, ntr = 1,
  restraintmask = ':1-261@CA', restraint_wt = 100
 /
EOF

pmemd.cuda -O -i min1.in -o min1.out -p 3CA2.prmtop -c 3CA2.rst7 -r min1.rst7 -ref 3CA2.rst7

# 检查能量下降
tail -30 min1.out
grep "minimization" min1.out
```

**预期**: EAMBER 从数万降到 < 1000 kcal/mol

---

## Step 2 — 加热 (0 → 298 K)

**目标**: 逐步升温，避免突然热冲击导致结构崩溃

```bash
cat > heat.in << 'EOF'
Heating — 0→298K, 1ns, 100 kcal/mol CA约束
 &cntrl
  imin = 0, nstlim = 500000,
  dt = 0.002, ntf = 2, ntc = 2,
  ntpr = 10000, ntwx = 5000, ntwr = 500000,
  temp0 = 298.0, ntt = 3, gamma_ln = 1.0, ig = -1,
  ntr = 1, restraintmask = ':1-261@CA', restraint_wt = 100,
  iwrap = 1, ioutfm = 1,
  nmropt = 1,          # 启用温度渐变
 /
 &wt type = 'TEMP0', istep1 = 0, istep2 = 500000,
   value1 = 0.0, value2 = 298.0 /
 &wt type = 'END' /
EOF

pmemd.cuda -O -i heat.in -o heat.out -p 3CA2.prmtop -c min1.rst7 -r heat.rst7 -x heat.nc

# 检查温度曲线
grep "^TEMP" heat.out | awk '{print $2}' | sort -u
```

**预期**: 温度在500ps内升至298K并保持稳定

---

## Step 3 — NPT 平衡 (100→10 kcal/mol约束)

**目标**: 在恒压下收紧溶剂，缓慢释放CA约束

```bash
cat > eq1.in << 'EOF'
NPT平衡 — 298K, 100→10 kcal/mol, 1ns
 &cntrl
  imin = 0, nstlim = 500000,
  dt = 0.002, ntf = 2, ntc = 2,
  ntpr = 10000, ntwx = 5000, ntwr = 500000,
  temp0 = 298.0, ntt = 3, gamma_ln = 1.0, ig = -1,
  ntp = 1, taup = 2.0,           # 等压，松弛时间2ps
  ntr = 1,
  restraintmask = ':1-261@CA',
  restraint_wt = 100,             # 从100开始
  iwrap = 1, ioutfm = 1,
 /
EOF

pmemd.cuda -O -i eq1.in -o eq1.out -p 3CA2.prmtop -c heat.rst7 -r eq1.rst7 -x eq1.nc

# 下一步用10 kcal/mol
cp eq1.rst7 eq2.rst7
```

---

## Step 4 — NPT 平衡 (10→1 kcal/mol约束)

```bash
cat > eq2.in << 'EOF'
NPT平衡 — 298K, 10→1 kcal/mol, 1ns
 &cntrl
  imin = 0, nstlim = 500000,
  dt = 0.002, ntf = 2, ntc = 2,
  ntpr = 10000, ntwx = 5000, ntwr = 500000,
  temp0 = 298.0, ntt = 3, gamma_ln = 1.0, ig = -1,
  ntp = 1, taup = 2.0,
  ntr = 1,
  restraintmask = ':1-261@CA',
  restraint_wt = 10,              # 降至10
  iwrap = 1, ioutfm = 1,
 /
EOF

pmemd.cuda -O -i eq2.in -o eq2.out -p 3CA2.prmtop -c eq2.rst7 -r eq2.rst7 -x eq2.nc
```

---

## Step 5 — NPT 平衡 (1→0.1 kcal/mol约束)

```bash
cat > eq3.in << 'EOF'
NPT平衡 — 298K, 1→0.1 kcal/mol, 1ns
 &cntrl
  imin = 0, nstlim = 500000,
  dt = 0.002, ntf = 2, ntc = 2,
  ntpr = 10000, ntwx = 5000, ntwr = 500000,
  temp0 = 298.0, ntt = 3, gamma_ln = 1.0, ig = -1,
  ntp = 1, taup = 2.0,
  ntr = 1,
  restraintmask = ':1-261@CA',
  restraint_wt = 1,               # 降至1
  iwrap = 1, ioutfm = 1,
 /
EOF

pmemd.cuda -O -i eq3.in -o eq3.out -p 3CA2.prmtop -c eq2.rst7 -r eq3.rst7 -x eq3.nc
```

---

## Step 6 — NPT 无约束平衡

```bash
cat > eq4.in << 'EOF'
NPT平衡 — 298K, 无约束, 1ns
 &cntrl
  imin = 0, nstlim = 500000,
  dt = 0.002, ntf = 2, ntc = 2,
  ntpr = 10000, ntwx = 5000, ntwr = 500000,
  temp0 = 298.0, ntt = 3, gamma_ln = 1.0, ig = -1,
  ntp = 1, taup = 2.0,
  ntr = 0,                # 完全放开约束
  iwrap = 1, ioutfm = 1,
 /
EOF

pmemd.cuda -O -i eq4.in -o eq4.out -p 3CA2.prmtop -c eq3.rst7 -r eq4.rst7 -x eq4.nc

# 验证系统稳定性
grep "EAMBER" eq4.out | tail -5
```

---

## Step 7 — 成品MD

**目标**: 无约束生产运行

```bash
cat > prod.in << 'EOF'
成品MD — 298K, NPT, 50ns (可延长至100ns+)
 &cntrl
  imin = 0, nstlim = 25000000,      # 50ns @ 2fs
  dt = 0.002, ntf = 2, ntc = 2,
  ntpr = 25000,                      # 每50ps输出能量
  ntwx = 25000,                      # 每50ps输出轨迹
  ntwr = 2500000,                    # 每5ns保存重启
  temp0 = 298.0, ntt = 3, gamma_ln = 1.0, ig = -1,
  ntp = 1, taup = 2.0,
  ntr = 0,
  iwrap = 1, ioutfm = 1,
 /
EOF

pmemd.cuda -O -i prod.in -o prod.out -p 3CA2.prmtop -c eq4.rst7 -r prod.rst7 -x prod.nc

# 如需延长，直接续跑
cp prod.rst7 prod_restart.rst7
# 修改nstlim继续运行
```

---

## 轨迹分析

### RMSD (稳定性监测)

```bash
conda activate amber

cat > rmsd_cpptraj.in << 'EOF'
trajin prod.nc
reference eq4.rst7
autoimage
rms reference :1-261@CA out rmsd_ca.dat
rms reference :1-261@C,N,O,CA out rmsd_backbone.dat
EOF

cpptraj -p 3CA2.prmtop -i rmsd_cpptraj.in -o rmsd.out
```

**判读标准**:
- CA RMSD < 3 Å（稳态）
- > 3 Å 需检查是否未平衡

### RMSF (柔性区域)

```bash
cat > rmsf_cpptraj.in << 'EOF'
trajin prod.nc
autoimage
rms first :1-261@CA
atomicfluct :1-261@CA out rmsf_ca.dat
EOF

cpptraj -p 3CA2.prmtop -i rmsf_cpptraj.in -o rmsf.out
```

**判读标准**:
- > 2 Å：高度柔性 (loops/termimi)
- < 1 Å：刚性区域

### Zn-AMS 配位距离监测

```bash
cat > distance_cpptraj.in << 'EOF'
trajin prod.nc
autoimage
distance :264@ZN :265@N1 out zn_ams_n1_dist.dat
distance :264@ZN :94@NE2 out zn_his94_dist.dat
distance :264@ZN :96@NE2 out zn_his96_dist.dat
distance :264@ZN :119@ND1 out zn_his119_dist.dat
EOF

cpptraj -p 3CA2.prmtop -i distance_cpptraj.in -o distance.out
```

**预期**: Zn-N1 ≈ 2.0 Å，Zn-His ≈ 2.0 Å

### 氢键分析

```bash
cat > hbond_cpptraj.in << 'EOF'
trajin prod.nc
autoimage
hbond avgout hbond_avg.dat writecontacts hb_cont.dat
EOF

cpptraj -p 3CA2.prmtop -i hbond_cpptraj.in -o hbond.out
```

### 能量与热力学监测

```bash
# 使用 process_mdout.perl
$AMBERHOME/bin/process_mdout.perl heat.out eq1.out eq2.out eq3.out eq4.out prod.out

# 生成文件
ls summary.*
```

---

## 输出文件汇总

| 阶段 | 输出文件 | 说明 |
|------|---------|------|
| Step 1 | `min1.rst7` | 最小化坐标 |
| Step 2 | `heat.nc`, `heat.rst7` | 加热轨迹+坐标 |
| Step 3 | `eq1.nc`, `eq1.rst7` | NPT 100约束 |
| Step 4 | `eq2.nc`, `eq2.rst7` | NPT 10约束 |
| Step 5 | `eq3.nc`, `eq3.rst7` | NPT 1约束 |
| Step 6 | `eq4.nc`, `eq4.rst7` | NPT 无约束 |
| Step 7 | `prod.nc`, `prod.rst7` | 成品MD ✅ |

---

## 快速启动脚本

```bash
conda activate amber

# ===== 运行全部平衡阶段 =====
for step in min1 heat eq1 eq2 eq3 eq4; do
  inp="${step}.in"
  out="${step}.out"
  rst="${step}.rst7"
  prv="${step%_*}.rst7"  # 前一步的rst

  if [ "$step" == "min1" ]; then
    prv="3CA2.rst7"
  fi

  echo "Running $step..."
  pmemd.cuda -O -i $inp -o $out -p 3CA2.prmtop -c $prv -r $rst ${step}.nc 2>&1 | tail -5
done

# ===== 运行成品MD =====
pmemd.cuda -O -i prod.in -o prod.out -p 3CA2.prmtop -c eq4.rst7 -r prod.rst7 -x prod.nc

# ===== 分析 =====
cpptraj -p 3CA2.prmtop -i rmsd_cpptraj.in
cpptraj -p 3CA2.prmtop -i rmsf_cpptraj.in
cpptraj -p 3CA2.prmtop -i distance_cpptraj.in
```

---

## 关键监测指标

| 指标 | 目标值 | 工具 | 警告阈值 |
|------|--------|------|---------|
| CA RMSD | < 3 Å | cpptraj | > 4 Å |
| Zn-AMS(N1) | ~2.0 Å | cpptraj | > 3 Å |
| Zn-His | ~2.0 Å | cpptraj | > 3 Å |
| 系统密度 | ~1.0 g/cm³ | process_mdout | < 0.95 |
| 温度 | 298 K ± 5 K | mdout | ± 20 K |
| 压力 | ~1 atm | mdout | ± 100 |
| 总能量 | 稳定 | mdout | 持续漂移 |

---

## 已知问题与解决方案

| 问题 | 原因 | 解决 |
|------|------|------|
| RMSD持续增长 | 约束释放过快 | 增加eq步骤，延长约束时间 |
| Zn-AMS距离>3Å | 配位键断裂 | 检查是否需要bond，或增加约束 |
| 能量爆炸 | dt过大/初始冲突 | dt=0.001，延长最小化 |
| 密度过低 | 盒子太小/水模型问题 | 检查OPCBOX大小 |

---

## 修改记录

| 日期 | 修改内容 |
|------|---------|
| 2026-04-06 | 初次编写，基于Tutorial 13体系参数 |
