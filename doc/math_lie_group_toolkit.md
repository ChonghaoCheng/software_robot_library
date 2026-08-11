# Math 李群工具集（SE(3) / SO(3)）

本文件记录 `Math` 模块中用于刚体运动学/控制的李群函数：约定、数学定义、数值稳定性处理，以及已验证的精度。所有函数声明在 [`MathFunctions.h`](../Math/include/Math/MathFunctions.h)、实现在 [`MathFunctions.cpp`](../Math/src/MathFunctions.cpp)，仅依赖 Eigen 与本模块的 `SkewSymmetric`（**不依赖 Model 层**）。

---

## 1. 全局约定

- **twist 排序**：统一为 \(\boldsymbol{\xi}=[\mathbf{v};\ \boldsymbol{\omega}]\)（线速度在前、角速度在后）。
- **旋转向量**：\(\mathbf{r}=\theta\hat{\mathbf{a}}\)，\(\theta\in[0,\pi]\)，取最短路径（\(q_w\ge0\)）。
- **指数/对数互逆**：`so3_exponential`↔`so3_logarithm`、`se3_exponential`↔`se3_logarithm` 在数值上互为逆。
- **平移块约定**：`se3_logarithm` 用 \(V^{-1}=J_l^{-1}\) 还原线速度（\(\mathbf{v}=J_l^{-1}(\boldsymbol{\omega})\,\mathbf{p}\)），故 `se3_exponential` 用 \(V=J_l\)（\(\mathbf{p}=J_l(\boldsymbol{\omega})\,\mathbf{v}\)）。两者据此严格互逆。
- \([\cdot]_\times\) 表示 3 维向量的反对称矩阵（由 `SkewSymmetric` 提供）。

---

## 2. 函数与数学定义

### 2.1 指数映射

**`so3_exponential(r) -> R`** —— Rodrigues 公式
\[
R=\exp([\mathbf{r}]_\times)=I+\frac{\sin\theta}{\theta}[\mathbf{r}]_\times+\frac{1-\cos\theta}{\theta^2}[\mathbf{r}]_\times^2 .
\]

**`se3_exponential([v;w]) -> T`**
\[
T=\begin{bmatrix}R & J_l(\boldsymbol{\omega})\,\mathbf{v}\\ \mathbf{0}^\top & 1\end{bmatrix},\qquad R=\exp([\boldsymbol{\omega}]_\times).
\]

### 2.2 SO(3) 左/右雅可比（旋转向量速率 ↔ 角速度）

这组函数用于把旋转向量导数 \(\dot{\mathbf{r}}\) 与角速度互相转换——正是修正「\(\dot{\mathbf{r}}\) 被当成 \(\boldsymbol{\omega}\)」类问题的关键。
\[
\boldsymbol{\omega}_{\text{spatial}}=J_l(\mathbf{r})\,\dot{\mathbf{r}},\qquad
\boldsymbol{\omega}_{\text{body}}=J_r(\mathbf{r})\,\dot{\mathbf{r}},
\]
\[
J_l=I+\frac{1-\cos\theta}{\theta^2}[\mathbf{r}]_\times+\frac{\theta-\sin\theta}{\theta^3}[\mathbf{r}]_\times^2,\qquad
J_r(\mathbf{r})=J_l(-\mathbf{r})=J_l(\mathbf{r})^\top .
\]
逆映射（已知角速度求 \(\dot{\mathbf{r}}\)）：
\[
J_l^{-1}=I-\tfrac12[\mathbf{r}]_\times+c(\theta)\,[\mathbf{r}]_\times^2,\qquad
c(\theta)=\frac{1}{\theta^2}-\frac{\cot(\theta/2)}{2\theta}\ \xrightarrow{\theta\to0}\ \frac{1}{12}.
\]
- `so3_left_jacobian(r)` / `so3_right_jacobian(r)`
- `so3_left_jacobian_inverse(r)` / `so3_right_jacobian_inverse(r)`

### 2.3 SE(3) 伴随

**`adjoint(T) -> 6x6`**（\([\mathbf{v};\boldsymbol{\omega}]\) 排序）：把 body twist 映射为 spatial twist。
\[
\mathrm{Ad}_T=\begin{bmatrix}R & [\mathbf{p}]_\times R\\ \mathbf{0} & R\end{bmatrix},\qquad
\boldsymbol{\xi}_{\text{spatial}}=\mathrm{Ad}_T\,\boldsymbol{\xi}_{\text{body}} .
\]
等价于 \(\mathbf{v}_s=R\mathbf{v}_b+\mathbf{p}\times(R\boldsymbol{\omega}_b),\ \boldsymbol{\omega}_s=R\boldsymbol{\omega}_b\)。注意它输出的是 screw-theory spatial twist；若目标 Jacobian 的线速度语义是末端点速度 \(\dot p\)，不能直接使用该伴随映射。

### 2.4 四元数辅助

- **`quaternion_to_rotation_vector(q) -> r`**：\(\theta=2\,\mathrm{atan2}(\|\mathbf{q}_v\|,q_w)\)，\(\mathbf{r}=\theta\,\mathbf{q}_v/\|\mathbf{q}_v\|\)；对非归一/非有限输入与近零向量做了保护，强制 \(q_w\ge0\)。
- **`quaternion_orientation_error(current, desired) -> r`**：世界系姿态误差 \(\mathbf{e}=\log(q_{des}\,q_{cur}^{-1})\)。

> 这些函数统一了控制器原先各自维护的旋转与位姿数学副本。

---

## 3. 数值稳定性

每个有 \(0/0\) 形式的系数都按区间处理：

1. **\(\theta\to0\)**：`so3_exponential`、`so3_left_jacobian(_inverse)` 的系数均改用 Taylor 展开（如 \(\sin\theta/\theta\approx1-\theta^2/6\)、\(c(\theta)\approx 1/12+\theta^2/720\)），阈值 `1e-6`。
2. **\(\theta\to\pi\)**：`J_l^{-1}` 的系数原写法 \((1+\cos\theta)/(2\theta\sin\theta)\) 在 \(\theta=\pi\) 处分子分母**同时趋零**、发生灾难性抵消。改用半角恒等式
   \[
   \frac{1+\cos\theta}{\sin\theta}=\cot\!\frac{\theta}{2},
   \]
   其中 \(\sin(\theta/2)\) 在 \(\theta=\pi\) 处为 1，远离零点，因此数值稳定。
   `se3_logarithm` 现在**直接复用** `so3_left_jacobian_inverse`，同步获得该改进（影响 RMPCC 与 `CartesianSpline::tangent_at_progress`）。

---

## 4. 已验证精度

随机 2×10⁵ 次、覆盖一般角 / 近 π / 纳弧度三种区间（验证脚本见提交说明）。代数恒等式直接反映函数本身精度：

| 校验（代数恒等式，无有限差分） | 最大误差 |
|---|---|
| `J_r = J_lᵀ` | 2.8e-16 |
| `J_l · r = r`（轴不变） | 9.5e-16 |
| `adjoint(T)` vs 共轭 `log(T·exp(ξ)·T⁻¹)` | 1.9e-15 |
| `J_l⁻¹ · J_l = I` | 7.2e-13 |
| `J_l = R · J_r` | 1.4e-12 |

| 往返 `se3` log(exp)（按角度区间） | 最大误差 |
|---|---|
| 一般角（<2.5 rad） | 4.9e-13 |
| 近 π | 1.7e-15 |
| 纳弧度（~1e-9 rad） | 2.5e-10（四元数 log 在亚纳弧度下的理论下限，无实际意义） |

**`cot(θ/2)` 修复前后对比**：`J_l⁻¹·J_l=I` 从 7.4e-9 → 7.2e-13；近 π 往返从 7.8e-9 → 1.7e-15。

---

## 5. 控制器复用情况

已完成的下沉/去重：

- ✅ **RMPCC**：明确区分 spatial twist 与机器人 Jacobian 的端点速度语义，出口只旋转 \([v;\omega]\)，不加入 \(\mathbf p\times\boldsymbol\omega\)。
- ✅ **MPC / MovingFrameMPC**：已改用 Math 的 `quaternion_to_rotation_vector`，删除 `SerialLinkMPC` 的静态副本。
- ✅ **MPCC**：已改用 Math 的 `quaternion_to_rotation_vector`，删除静态四元数副本。

仍待处理：

- ⬜ **MPC 姿态 bug**：仍是「绝对旋转向量相减」。修复方向是用 `so3_left_jacobian_inverse` / `so3_right_jacobian_inverse` 在 \(\dot{\mathbf{r}}\) 与 \(\boldsymbol{\omega}\) 间正确换算（工具已就位，算法尚未改）。

详见 [`mpc_controllers_review.md`](mpc_controllers_review.md)。

---

## 6. 附：MPC 凝聚装配（CondensedMPC）

与李群工具一同下沉的还有线性 MPC 的「凝聚装配」助手，见 [`CondensedMPC.h`](../Math/include/Math/CondensedMPC.h) / [`.cpp`](../Math/src/CondensedMPC.cpp)，同样纯 Eigen、不依赖 Model/Control。

- **`condense_prediction(A, B, N) -> {stateTransition, inputResponse}`**：为 \(x_{k+1}=Ax_k+Bu_k\) 生成堆叠预测 \(X=A_x x_0+B_u U\)（\(X=[x_1;\dots;x_N]\)），其中 \((A_x)_k=A^{k+1}\)、\((B_u)_{k,j}=A^{k-j}B\ (j\le k)\)。\(A\) 的幂只预计算一次。
- **`block_diagonal(block, N)`**：把方块沿对角复制 \(N\) 次（用于堆叠 \(Q,R\)）。
- **`box_constraint(lower, upper) -> {constraintMatrix, constraintVector}`**：把 \(\text{lower}\le z\le\text{upper}\) 写成 `QPSolver` 期望的 \(Bz\le c\)（\(B=[I;-I],\ c=[\text{upper};-\text{lower}]\)）。

**已替换**：`SerialLinkMPC::solveMPC` 与 `SerialLinkMPCC::solve_mpcc` 中各自手写的 \(A_x/B_u\) 循环、\(Q_b/R_b\) 块对角、盒约束装配（约 3 处 × 各 ~40 行）现已统一调用这三个函数。

**正确性**：随机算例下 `stateTransition*x0 + inputResponse*U` 与逐步前向仿真一致（误差 ~1e-16，含 \(A\neq I\)）；`box_constraint` 的 \(c=[\text{upper};-\text{lower}]\) 与原实现逐元素一致。RMPCC 的预测雅可比是误差动力学专有结构，不走该通用通道，保持不变。
