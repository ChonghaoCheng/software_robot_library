# 轨迹生成与控制器调用流程

本文档说明本库**笛卡尔轨迹是如何生成的**，以及三个预测控制器（`SerialLinkMPC` / `SerialLinkMPCC` / `SerialLinkRMPCC`）**分别如何调用轨迹**，最后给出从「一组路点」到「关节速度指令」的完整调用链。控制器内部算法细节见配套文档 [`mpc_controllers_review.md`](mpc_controllers_review.md)。

---

## 1. 总体流程

```text
路点 poses[] + 时间 times[]
        │  CartesianSpline 构造：四元数 → 旋转向量
        ▼
  6 维向量 [p; r] 序列
        │  SplineTrajectory：每维独立建样条
        ▼
  Spline（分段）→ Polynomial（每段三次多项式）
        │  solve_cubic_spline_derivatives：拟合中间点速度（C² 连续）
        ▼
  时间索引轨迹  T(t), ṫ(t)
        │
        ├── 时间接口：query_state(t)            ← SerialLinkMPC 使用
        └── 进度接口：pose_at_progress(s),       ← SerialLinkMPCC / RMPCC 使用
                     tangent_at_progress(s)
        ▼
  控制器求解最优末端 twist  t_B ∈ ℝ⁶
        │  SerialLinkVelocityBase::resolve_endpoint_twist()
        ▼
  关节速度 q̇（QP，含关节限位/奇异/冗余）
```

---

## 2. 轨迹是如何生成的

轨迹生成自底向上分四层：`Polynomial` → `Spline` → `SplineTrajectory` → `CartesianSpline`。

### 2.1 Polynomial：单段多项式

[`Math/Polynomial.cpp`](../Math/src/Polynomial.cpp) 表示一段
\[
f(x)=\sum_{i=0}^{n} a_i\,x^{i}
\]
给定段起点 \(x_0\)、终点 \(x_1\) 处的值与导数（位置、速度、二阶导），用线性方程组求系数。构造一个 \((n{+}1)\times(n{+}1)\) 的范德蒙德型矩阵 \(X\)（行依次为 0/1/2 阶导在两端的取值），解
\[
X\,\mathbf{a}=\mathbf{b}\quad\Rightarrow\quad \mathbf{a}=X^{-1}\mathbf{b}.
\]
阶数必须为奇数（边界条件成对）：三次（order=3）用两端的「值 + 一阶导」，五次（order=5）再加二阶导。`evaluate_point(x)` 同时返回 \(f,\,f',\,f''\)。

### 2.2 Spline：分段拼接

[`Math/Spline.cpp`](../Math/src/Spline.cpp)：\(n\) 个支撑点之间生成 \(n{-}1\) 段 `Polynomial`，相邻段共享端点的值/导数从而保证连续。`evaluate_point(x)` 按区间查找落在第几段并求值；超出支撑范围则**钳位**到端点（轨迹在支撑区间外无定义）。

### 2.3 中间点速度拟合：C² 连续

三次样条只能保证每段内部光滑，**段与段之间的二阶导连续**需要先求出各中间点的速度。[`solve_cubic_spline_derivatives`](../Math/src/MathFunctions.cpp#L202) 解三对角系统
\[
A\,\dot{\mathbf{y}} = B\,\mathbf{y}\quad\Rightarrow\quad \dot{\mathbf{y}}=A^{-1}B\,\mathbf{y},
\]
其中中间行（令 \(h_1=x_i-x_{i-1},\ h_2=x_{i+1}-x_i\)）为
\[
\frac{1}{h_1}\dot y_{i-1}+2\Big(\frac{1}{h_1}+\frac{1}{h_2}\Big)\dot y_i+\frac{1}{h_2}\dot y_{i+1}
=-\frac{3}{h_1^2}y_{i-1}+3\Big(\frac{1}{h_1^2}-\frac{1}{h_2^2}\Big)y_i+\frac{3}{h_2^2}y_{i+1},
\]
首末行固定边界速度（起始速度由用户给，末端速度默认 0）。求得的 \(\dot y_i\) 作为每段 `Polynomial` 的端点一阶导，从而整条样条 C² 连续。

### 2.4 SplineTrajectory：多维 + 时间

[`Trajectory/SplineTrajectory.cpp`](../Trajectory/src/SplineTrajectory.cpp)：对 \(d\) 维向量轨迹，**每一维独立**建一条 `Spline`（共享同一组时间结点）。当路点数 > 2 时调用 §2.3 拟合中间速度；只有 2 个点时直接用起末速度。`query_state(t)` 在区间外返回起/末状态，区间内对每维求值，组装成
\[
\text{State}=\{\mathbf{y}(t),\ \dot{\mathbf{y}}(t),\ \ddot{\mathbf{y}}(t)\}.
\]

### 2.5 CartesianSpline：SE(3) 封装 + 进度接口

[`Trajectory/CartesianSpline.cpp`](../Trajectory/src/CartesianSpline.cpp) 把 SE(3) 位姿序列变成可插值的 6 维实向量轨迹。

**构造（四元数 → 旋转向量）。** 每个路点位姿拆成
\[
\mathbf{y}_i=\big[\underbrace{\mathbf{p}_i}_{3},\ \underbrace{\theta_i\,\hat{\mathbf{a}}_i}_{3}\big],\qquad
\theta_i=2\arccos(w_i),\ \hat{\mathbf{a}}_i=\frac{\mathbf{q}_{v,i}}{\|\mathbf{q}_{v,i}\|},
\]
交给 `SplineTrajectory` 在时间上插值。

**时间接口 `query_state(t)`（→ MPC 使用）。** 把插值得到的 6 维 \(\mathbf{y}(t)\) 还原成位姿：位置直接取前 3 维，姿态由旋转向量 \(\mathbf{r}=\mathbf{y}_{3:6}\) 重建四元数 \(\big(\cos\tfrac{\theta}{2},\ \sin\tfrac{\theta}{2}\hat{\mathbf{a}}\big)\)。返回 `CartesianState{pose, twist, acceleration}`。
> 注意：`twist` 角分量是 \(\dot{\mathbf{r}}\)（旋转向量导数）**而非角速度** \(\boldsymbol{\omega}\)；位置分量 \(\dot{\mathbf{p}}\) 才是真正的线速度。

**进度接口（→ MPCC / RMPCC 使用）。** 统一几何接口把归一化进度 \(s\in[0,1]\) 当作路径参数：
\[
t(s)=t_{start}+s\,(t_{end}-t_{start}),\qquad T(s)=\text{query\_state}(t(s)).\text{pose},
\]
并提供 **body 系 SE(3) 切向**（有限差分，接近 \(s=1\) 时退化为后向差分）：
\[
\boldsymbol{\tau}(s)=\frac{1}{\delta}\log_{SE(3)}\!\big(T(s)^{-1}T(s+\delta)\big)\in\mathbb{R}^6.
\]
`tangent_at_progress` 因为先重建位姿再取 SE(3) log，所以对旋转向量表示法稳健（不受 \(\dot{\mathbf{r}}\neq\boldsymbol{\omega}\) 影响）。

---

## 3. 控制器分别如何调用轨迹

三者的差异本质是「当前进度从哪里来 / 用什么接口采样参考」。

### 3.1 SerialLinkMPC —— 时间索引（`query_state`）

- 通过 `set_trajectory(CartesianSpline)` **持有**完整轨迹。
- 每拍 `track_endpoint_trajectory_at_time(t)` 在 \(t+k\Delta t\) 处对 `query_state()` 采样 \(N\) 个未来参考，构成状态/控制参考栈。
- 参考姿态用旋转向量、控制参考用 `CartesianState.twist`。
- 也保留兼容入口 `track_endpoint_trajectory(pose,vel,acc)`：无完整轨迹时用单点 + 速度线性外推。

### 3.2 SerialLinkMPCC —— 外部进度 + 局部欧氏误差

- 通过 `set_trajectory(CartesianSpline)` 持有完整轨迹。
- 每拍 `step(dt, estimatedProgress)` 接收外部最近点估计 \(\hat s\)，它是当前进度的唯一真值。
- QP 用 warm start 产生名义未来进度，并在每个 stage 查询 \(T(s_k)\)、\(\tau(s_k)\)。
- 局部误差预测为 \(e_{k+1}\approx e_k+\Delta t(u_k-\tau(s_k)\dot s_k)\)，因此进度与路径几何直接耦合。

### 3.3 SerialLinkRMPCC —— 外部进度 + SE(3) 误差

- 同样持有完整轨迹，并由 `step(dt, estimatedProgress)` 接收外部最近点估计 \(\hat s\)。
- `reference_transform(s)` → `pose_at_progress(s)` 取参考位姿（叠加刚性扰动 \(D\)）。
- `body_twist_reference_at_progress(s)` 与 QP 内每个 stage → `tangent_at_progress(s)` 取 body 切向 \(\boldsymbol{\tau}(s)\)，用于前馈与误差预测雅可比。
- \(\dot s_k\) 只负责预测域内的未来进度；控制器不再把积分结果当作下一拍真实进度。

| | 持有轨迹 | 当前进度来源 | 采样接口 | 主入口 |
|---|---|---|---|---|
| MPC   | ✔ | ✔（以时间为进度） | `query_state(t)` | `track_endpoint_trajectory_at_time(t)` |
| MPCC  | ✔ | 外部最近点估计 \(\hat s\) | `pose_at_progress(s)` / `tangent_at_progress(s)` | `step(dt, ŝ)` |
| RMPCC | ✔ | 外部最近点估计 \(\hat s\) | `pose_at_progress(s)` / `tangent_at_progress(s)` | `step(dt, ŝ)` |

---

## 4. 从轨迹到关节速度的完整调用链

以 RMPCC 为例（MPC/MPCC 结构一致，仅参考采样方式不同）：

```text
set_trajectory(spline)
   └─ reset(); 自动推导 ṡ_ref
step(dt, estimatedProgress)
   ├─ 用 estimatedProgress 同步当前进度 s
   ├─ SerialLinkBase::update()                     // 刷新 J、末端位姿
   ├─ solve_rmpcc(T_cur, dt)
   │     ├─ reference_transform(s) → pose_at_progress(s)   // 参考位姿
   │     ├─ tangent_at_progress(s) × 每 stage             // body 切向 τ(s)
   │     └─ QPSolver::solve(H,f,Bineq,zineq,z0)           // 解 7N 维 QP
   ├─ 旋转：body twist → 基座系端点速度 [ṗ; ω]
   ├─ （可选）+ 位姿反馈 K·e_pose
   └─ SerialLinkVelocityBase::resolve_endpoint_twist(t_B)
         └─ SerialLinkVelocityBase::resolve_endpoint_motion(t_B)
               ├─ compute_control_limits(i)  ∀ 关节        // 关节限速/限位
               └─ QPSolver::solve(...)  → q̇               // 末端 twist → 关节速度
```

末端 twist → 关节速度这一层（求解 \(J\dot{\mathbf{q}}\approx\mathbf{t}\) 的 QP 与约束）见 [`control_twist_to_joint_velocity_summary.md`](control_twist_to_joint_velocity_summary.md)。

> 关节限位由共同父类 `SerialLinkVelocityBase` 的速度解析层统一施加。
