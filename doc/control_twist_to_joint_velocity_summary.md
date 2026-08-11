# 控制部分总结：由末端 Twist 到关节速度（含约束）

本文档总结该库在控制阶段如何将末端速度指令（twist）映射为关节速度，并明确列出优化问题及约束项。

---

## 1. 总体流程

在 `SerialLinkMPCC` 中，优化器得到局部坐标系下的最优控制：

\[
\mathbf{u}_{opt}=
\begin{bmatrix}
\mathbf{v}_T\\
\boldsymbol{\omega}_T\\
v_s
\end{bmatrix}
\in\mathbb{R}^7
\]

其中前 6 维是末端线速度/角速度，最后 1 维是进度速度。随后将前 6 维旋回基坐标系：

\[
\mathbf{v}_B = R_{TW}^\top \mathbf{v}_T,
\qquad
\boldsymbol{\omega}_B = R_{TW}^\top \boldsymbol{\omega}_T
\]

构成末端 twist：

\[
\mathbf{t}_B=
\begin{bmatrix}
\mathbf{v}_B\\
\boldsymbol{\omega}_B
\end{bmatrix}
\in\mathbb{R}^6
\]

然后调用
`SerialLinkKinematic::resolve_endpoint_twist()`，并进一步进入
`resolve_endpoint_motion()` 求解关节速度 \(\dot{\mathbf{q}}\)。

---

## 2. Twist 到关节速度的优化模型

记末端雅可比为 \(J\in\mathbb{R}^{6\times n}\)，目标末端速度为 \(\mathbf{t}\in\mathbb{R}^6\)，则核心关系是：

\[
J\dot{\mathbf{q}}\approx \mathbf{t}
\]

### 2.1 非冗余/欠驱动（\(n\le 6\)）

实现中求解二次规划：

\[
\min_{\dot{\mathbf{q}}}
\frac{1}{2}\dot{\mathbf{q}}^\top H\dot{\mathbf{q}} + \mathbf{f}^\top\dot{\mathbf{q}}
\]

其中

\[
H = J^\top J,
\qquad
\mathbf{f} = -J^\top\mathbf{t}
\]

等价于最小化跟踪误差 \(\|J\dot{\mathbf{q}}-\mathbf{t}\|^2\)（在约束下）。

### 2.2 冗余（\(n>6\)）

实现中使用约束最小二乘：

\[
\min_{\dot{\mathbf{q}}}
(\dot{\mathbf{q}}_d-\dot{\mathbf{q}})^\top W(\dot{\mathbf{q}}_d-\dot{\mathbf{q}})
\]

\[
\text{s.t.}
\quad J\dot{\mathbf{q}} = \mathbf{t},
\quad B\dot{\mathbf{q}}\le \mathbf{z}
\]

其中 \(\dot{\mathbf{q}}_d\) 是冗余任务（默认由可操作度梯度构造），\(W\) 使用关节惯性矩阵。

### 2.3 奇异邻域（Damped Least Squares）

当接近奇异位形时，采用阻尼项：

\[
H = J^\top J + \lambda I,
\qquad
\lambda = \Big(1-\frac{\mu}{\mu_{min}}\Big)^2\cdot 0.01
\]

其中 \(\mu\) 是当前可操作度，\(\mu_{min}\) 为阈值。然后仍解带约束 QP。

---

## 3. 约束（Constraints）构造

优化统一采用不等式形式：

\[
B\dot{\mathbf{q}}\le \mathbf{z}
\]

包含两类约束：

### 3.1 关节速度上下界约束

对每个关节 \(i\)：

\[
\dot{q}_i^{min} \le \dot{q}_i \le \dot{q}_i^{max}
\]

上下界不是固定常数，而是基于当前位置实时计算，综合以下因素取最紧边界：

1. 距离关节位置极限的“一步可达”限制（与控制频率相关）
2. 关节自身速度上限
3. 最大关节加速度诱导的速度限制

实现公式（与代码一致）可写为：

\[
\dot{q}_i^{min}=
\max\left(
-\Delta_i^- f_c,
-v_i^{lim},
-2\sqrt{a_{max}\Delta_i^-}
\right)
\]

\[
\dot{q}_i^{max}=
\min\left(
\Delta_i^+ f_c,
 v_i^{lim},
 2\sqrt{a_{max}\Delta_i^+}
\right)
\]

其中：

- \(\Delta_i^- = q_i-q_i^{low}\)
- \(\Delta_i^+ = q_i^{up}-q_i\)
- \(f_c\) 为控制频率
- \(v_i^{lim}\) 为关节速度极限
- \(a_{max}\) 为设定的关节最大加速度

### 3.2 可操作度（Manipulability）屏障约束

实现中还加入一条基于可操作度梯度 \(\nabla\mu\) 的线性不等式（控制屏障思想）：

\[
-\nabla\mu(\mathbf{q})^\top\dot{\mathbf{q}}
\le
(\mu-\mu_{min})\,100\sqrt{f_c}
\]

该项用于在优化中抑制向低可操作度区域运动，降低奇异风险。

---

## 4. 可直接放论文的方法描述（简版）

本文在控制层采用“末端 twist 到关节速度”的约束优化逆解。首先由上层控制器（如 MPCC）产生末端线/角速度指令，并转换至机器人基坐标系，形成 \(\mathbf{t}\in\mathbb{R}^6\)。随后基于雅可比关系 \(J\dot{\mathbf{q}}\approx\mathbf{t}\) 求解关节速度：对非冗余系统使用带不等式约束的二次规划，对冗余系统使用带等式（任务一致性）与不等式（安全约束）的约束最小二乘。约束项包括关节位置/速度/加速度联合诱导的瞬时速度上下界，以及基于可操作度梯度的屏障约束；当接近奇异位形时，引入阻尼最小二乘项提高数值鲁棒性。

---

## 5. 对应源码位置

- MPCC 输出 twist 并映射回基坐标系：`Control/src/SerialLinkMPCC.cpp`
- twist 入口与求解：`Control/src/SerialLinkKinematic.cpp`
- 关节限幅与可操作度约束：`Control/src/SerialLinkKinematic.cpp`
- 控制基类接口与雅可比访问：`Control/include/Control/SerialLinkBase.h`

---

## 6. 技术细节补充（建议补充在方法说明后）

### 6.1 目标函数项的完整说明

建议在文中把控制优化目标明确拆分为：

1. 末端 twist 跟踪误差项（主任务）
2. 输入幅值正则项（抑制过大关节速度）
3. 变化率正则项（若使用 \(\Delta u\) 或相邻步平滑）
4. 冗余任务项（可操作度提升或其他次任务）
5. 奇异区阻尼项（DLS）

这样比仅写“QP求解逆运动学”更完整。

### 6.2 约束激活与可行性

建议补充：

- 约束激活比例（实验统计中有说服力）；
- 当约束导致不可行时的处理策略（例如放松项、降级控制或保持上一时刻解）；
- 起始点（warm-start）如何选取（当前实现使用当前关节速度并裁剪到可行域附近）。

### 6.3 奇异切换机制

建议写清楚奇异判据（\(\mu<\mu_{min}\)）与阻尼系数公式，避免读者认为是经验切换。可补一句：阻尼随可操作度下降而增大，保证数值稳定与连续性。

### 6.4 约束来源的物理意义

建议在正文中点出三类硬限制对应物理来源：

- 位置极限避免撞限位；
- 速度极限遵循执行器额定能力；
- 加速度诱导限制用于限制瞬时速度变化，改善控制平滑性与可实现性。

### 6.5 计算代价与实时性

建议补充每周期求解时间统计（均值/最大值/95分位），并与控制频率对照，说明实时性是否满足。

---

## 7. 可补充的一段简述

在实现层面，本方法不仅给出 twist-to-joint 的约束优化形式，还通过瞬时速度上下界、可操作度屏障和奇异区阻尼机制保证可行性与数值鲁棒性。为确保工程可复现，建议同时报告权重参数、约束触发统计、求解耗时分布以及不可行工况下的回退策略。
