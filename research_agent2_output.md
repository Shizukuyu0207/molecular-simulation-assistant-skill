# 分子动力学软件手册调研报告

> 调研日期: 2026/04/07
> 数据来源: 本地知识库（Web访问受限，基于已有知识编译）
> 覆盖: VMD + PyMOL + Gaussian09 + Chimera

---

## 一、VMD (Visual Molecular Dynamics)

> **官网**: https://www.ks.uiuc.edu/Research/vmd/
> **最新版本**: 1.9.4+ | **平台**: Linux/Mac/Windows
> **语言**: Tcl脚本 + Python插件

### 1.1 安装
```bash
# conda
conda install -c conda-forge vmd

# 或下载二进制
wget https://www.ks.uiuc.edu/Research/vmd/alpha/vmd-1.9.4.bin.LINUXAMD64.tar.gz
tar xzf vmd-1.9.4.bin.LINUXAMD64.tar.gz
cd vmd-1.9.4
./configure
make install
```

### 1.2 AMBER 轨迹读取
```tcl
# 方法1: GUI 加载
# File → New Molecule → Browse
# File type: AMBER Parmtop (选择.prmtop)
# 然后 Add: NetCDF Grid (选择.nc)

# 方法2: Tcl命令行
mol new system.prmtop type amberparm
animate read netcdf md.nc waitfor all

# 方法3: Python
package require mdff
mol new system.prmtop type amberparm
animate read netcdf md.nc
```

### 1.3 核心 Tcl 命令

#### 显示模式
```tcl
# 蛋白Cartoon
mol addrep 0
set sel [atomselect 0 protein]
$sel set beta 0
display renames

# 常用显示模式
# NewCartoon: 蛋白/核酸(默认)
# Lines: 仅显示化学键
# CPK: 球棍模型
# Licorice: 棍模型(ball&stick)
mol modstyle 0 0 NewCartoon
mol modstyle 0 0 Licorice
```

#### RMSD 计算
```tcl
package require mdff

# 蛋白CA原子RMSD
set ref [atomselect 0 "protein and name CA" frame 0]
set sel [atomselect 0 "protein and name CA"]

set nframes [molinfo 0 get numframes]
for {set i 1} {$i < $nframes} {incr i} {
    $sel frame $i
    set rmsd [measure rmsd $sel $ref]
    puts "Frame $i: RMSD = $rmsd"
}
```

#### 氢键分析
```tcl
package require hbonds

set sel [atomselect 0 "protein"]
hbonds -sel $sel -dist 3.5 -angle 120 -iterations 1000

# 距离监测
set id1 [atomselect 0 "resname ZN and resid 264"]
set id2 [atomselect 0 "resname AMS and resid 265 and name N1"]
measure distance $id1 $id2
```

#### 轨迹操作
```tcl
# 播放速度
animate speed 0.5

# 前进/后退
animate forward
animate backward

# 跳到指定帧
animate goto 100
```

### 1.4 渲染发表级图片
```tcl
# 背景色
color Display Background white

# 光线追踪模式
display rendermode GLSL

# Tachyon渲染器(高质量)
render TachyonInternal screenshot.png

# 或 snapshot
render snapshot screenshot.tga

# 设置光源
light 0 on
light 1 on
light 2 off
```

### 1.5 与 AMBER 联用实战

#### 读取 AMBER 输出
```tcl
# 读取拓扑+轨迹
mol new 3CA2.prmtop type amberparm
animate read netcdf prod.nc waitfor all

# 自动成像
package require autoimagemd
autoimage prod.nc
```

#### 分析 Zn-配位距离
```tcl
set nframes [molinfo 0 get numframes]
set outfile [open "zn_dist.dat" w]

for {set i 0} {$i < $nframes} {incr i} {
    set sel1 [atomselect 0 "resid 264 and name ZN" frame $i]
    set sel2 [atomselect 0 "resid 265 and name N1" frame $i]
    set dist [measure distance $sel1 $sel2]
    puts $outfile "$i $dist"
}

close $outfile
```

---

## 二、PyMOL

> **官网**: https://pymol.org/
> **版本**: 2.5+ (Python 3.8+) | **授权**: 开源 (Python API)
> **功能**: 分子可视化、药物设计、动画制作

### 2.1 安装
```bash
# conda (推荐)
conda install -c conda-forge pymol

# pip
pip install pymol

# 启动
pymol
pymol -c script.py  # 无GUI模式
```

### 2.2 Python API 基础

#### 读取分子
```python
from pymol import cmd

# 读取PDB
cmd.load("protein.pdb", "protein")
cmd.load("ligand.mol2", "ligand")

# 读取轨迹
cmd.load("md_traj.nc", "traj")

# 叠图对齐
cmd.align("protein", "reference", cycles=0)
```

#### 显示控制
```python
# 显示模式
cmd.show("cartoon", "protein")      # -cartoon
cmd.show("sticks", "ligand")         # 棍
cmd.show("spheres", "resname ZN")   # 球
cmd.show("surface", "protein")       # 表面

# 隐藏
cmd.hide("lines", "protein")

# 着色
cmd.color("red", "elem C")              # 按元素
cmd.color("cyan", "chain A")            # 按链
cmd.spectrum("b", "blue_white_red", "protein and name CA")  # B因子着色
```

#### 标签
```python
cmd.label("resn AMS", '"%s-%s" % (resn, resi)')
cmd.set("label_size", 14)
cmd.set("label_font_id", 7)
```

### 2.3 动画制作
```python
# 设置场景序列
cmd.mset("1x100")  # 100帧

# 旋转动画
cmd.rotate("x", 360, "protein")

# 缩放
cmd.zoom("protein", 2.0)

# 保存动画
cmd.movie(fps=30)
cmd.png("frame.png", dpi=300, ray=1)

# 批量渲染
for i in range(100):
    cmd.frame(i)
    cmd.png(f"frame_{i:03d}.png", ray=1)
```

### 2.4 Ray Tracing 高质量渲染
```python
# 设置尺寸
cmd.ray(3200, 2400)

# 渲染
cmd.png("output.png", dpi=300, ray=1)

# 调整光照
cmd.set("light", [0.5, -1.0, 1.0])
cmd.set("ambient", 0.5)
cmd.set("specular", 1.0)
```

### 2.5 与 AMBER 联用
```python
# PyMOL 可直接读取 PDB/MOL2/NC
cmd.load("3CA2_system.pdb")
cmd.load("prod.nc")

# 分析RMSD
cmd.align("chain A and name CA", "reference and name CA")

# 制作RMSD曲线
cmd.rama("all")  # Ramachandran图
```

---

## 三、Gaussian09

> **官网**: https://gaussian.com/
> **功能**: 量子化学计算、几何优化、频率分析、势能面
> **授权**: 商业软件 (学术许可证可用)

### 3.1 输入文件格式 (.gjf/.com)
```bash
%mem=4GB                    # 内存
%nprocshared=8              # CPU核数
%chk=calc.chk               # 检查点文件

#p B3LYP/6-31G(d) opt freq  # 方法/基组 计算类型

标题行
空行
电荷  多重度
原子坐标...

空行
其他指令
```

### 3.2 常用计算类型

#### 几何优化 (Opt)
```bash
%mem=8GB
%nproc=16
%chk=mycalc.chk
#p B3LYP/6-31G(d) opt

Title
0 1
C  0.0  0.0  0.0
H  0.5  0.5  0.5
...

```

#### 频率计算 (Freq)
```bash
#p B3LYP/6-31G(d) freq

Title
0 1
...
```
**输出**: 能量、Hessian矩阵、热力学修正、振动频率
**判读**: 无负频率 = 真正极小点

#### NMR 计算
```bash
#p B3LYP/6-311G(d,p) nmr=giao

Title
0 1
...
```

#### 单点能
```bash
#p CCSD(T)/aug-cc-pVTZ

Title
0 1
...
```

### 3.3 基组选择

| 基组 | 描述 | 用途 |
|------|------|------|
| 6-31G(d) | DZP | 日常优化/频率 |
| 6-311G(d,p) | TZP | 高精度 |
| aug-cc-pVDZ | aug-cc-pVnZ | 氢键/色散 |
| LANL2DZ | ECP | 过渡金属 |
| def2-TZVP | 极化 | 高精度+效率 |

### 3.4 与 AMBER 联用 — RESP 电荷

```bash
# 步骤1: Gaussian计算ESP
#p b3lyp/6-31G(d) pop=esp

Title
0 1
... (配体坐标) ...

# 步骤2: 用 Multiwfn 拟合RESP
# 或用antechamber的resp程序

# 步骤3: GAFF参数化
antechamber -i ligand.mol2 -fi mol2 \
  -o ligand_resp.mol2 -fo mol2 -c resp -s 2
```

### 3.5 常见错误与解决

| 错误 | 原因 | 解决 |
|------|------|------|
| SCF didn't converge | 电子未收敛 | 加 DIIS 或改算法 |
| Opt cycle failed | 优化失败 | 改初始结构 |
| Basis not found | 基组名错误 | 检查文档 |
| File size limit | 文件过大 | 分块计算 |

---

## 四、Chimera (UCSF Chimera)

> **官网**: https://www.cgl.ucsf.edu/chimera/
> **最新版本**: 1.17+ | **平台**: Mac/Windows/Linux
> **功能**: 分子可视化、结构编辑、能量最小化、MD轨迹

### 4.1 安装
```bash
# 下载
wget https://www.cgl.ucsf.edu/chimera/files/1.17/chimera.bin.tar.gz
tar xzf chimera.bin.tar.gz
cd chimera/bin
./chimera
```

### 4.2 基础命令

#### 文件操作
```bash
# 启动GUI
chimera

# 无GUI模式
chimera --nogui script.py

# 读取文件
open system.pdb
open ligand.mol2
open md.nc  # 需要MD Movie插件
```

#### 常用操作
```bash
# 旋转/移动
turn x 45
move x 5
translate [0,0,5] model #0

# 删原子
delete #0:ZN

# 原子着色
color byelement #0
color red #0:ALA

# 显示模式
show cartoon #0
show sticks #0:ligand
```

### 4.3 MD Movie (轨迹播放)
```
菜单: Tools → MD/Ensemble → MD Movie
支持格式: NetCDF (.nc), DCD, AMBER轨迹
功能:
- 播放/暂停/停止
- 跳帧
- RMSD/RMSF绑定图
- 距离监测图
```

### 4.4 结构编辑
```bash
# 旋转侧链
rotsl #0:ALA

# 突变
swapaa #0:ALA:CA  # 突变为任意氨基酸

# 加H
addh #0
```

### 4.5 与 AMBER 联用
```bash
# 读取AMBER系统
open 3CA2_system.pdb

# 监测Zn-配位
# Tools → MD/Ensemble → Distance Monitor
# 选择Zn原子和配位原子

# 制作动画
# MD Movie → File → Record → 保存为MP4
```

---

## 五、软件对比速查表

| 软件 | 擅长 | 渲染质量 | AMBER联用 | 学习曲线 |
|------|------|---------|-----------|---------|
| VMD | 轨迹分析/大体系 | ★★★★ | ★★★★★ 原生支持 | 中 |
| PyMOL | 发表图片/药物设计 | ★★★★★ | ★★★ 需转换 | 低 |
| Gaussian09 | 量子化学/RESP | N/A | ★★★★ ESP电荷 | 高 |
| Chimera | 结构编辑/MD动画 | ★★★★ | ★★★★ 轨迹播放 | 低 |

---

## 六、参考链接

- VMD: https://www.ks.uiuc.edu/Research/vmd/
- PyMOL: https://pymol.org/
- Gaussian: https://gaussian.com/
- Chimera: https://www.cgl.ucsf.edu/chimera/
- AMBER: https://ambermd.org/
