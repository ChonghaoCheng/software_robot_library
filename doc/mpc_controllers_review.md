# 笛卡尔空间 MPC 系列控制器：算法、函数与缺陷总结

本文档系统总结 `Control` 模块中三个基于预测优化的笛卡尔空间控制器：

| 控制器 | 文件 | 一句话定位 |
|--------|------|-----------|
| `SerialLinkMPC`   | [`SerialLinkMPC.cpp`](../Control/src/SerialLinkMPC.cpp)   | 任务空间速度级 MPC，状态为绝对位姿，控制为末端 twist |
| `SerialLinkMPCC`  | [`SerialLinkMPCC.cpp`](../Control/src/SerialLinkMPCC.cpp) | 局部路径对齐的轮廓控制（contouring），含进度变量 |
| `SerialLinkRMPCC` | [`SerialLinkRMPCC.cpp`](../Control/src/SerialLinkRMPCC.cpp) | 直接在 \(SE(3)\) 上的黎曼 MPCC，控制器自身拥有参考路径与进度 |

三者都遵循同一条出口链路：外层在任务空间求得最优末端 twist \(\mathbf{t}_B\in\mathbb{R}^6\)，再委托内层
`SerialLinkKinematic::resolve_endpoint_twist()` 把 twist 映射为关节速度 \(\dot{\mathbf{q}}\)（关节限位、奇异、冗余都在内层处理）。

---

## 0. 公共依赖：被修改后的 `CartesianSpline`

三个控制器都依赖 [`CartesianSpline`](../Trajectory/include/Trajectory/CartesianSpline.h)。该类近期新增了一套**统一的「路径进度」几何接口**，供拥有进度的控制器使用：

- `progress_to_time(s)`：把归一化进度 \(s\in[0,1]\) 线性映射到轨迹时间
  \[
  t(s) = t_{start} + s\,(t_{end}-t_{start})
  \]
- `pose_at_progress(s)`：进度 \(s\) 处的参考位姿 \(T(s)\)。
- `tangent_at_progress(s, \delta)`：**body 系下**的 \(SE(3)\) 路径切向（6 维 body twist），用有限差分定义
  \[
  \boldsymbol{\tau}(s)=\frac{1}{\delta}\,\log_{SE(3)}\!\big(T(s)^{-1}\,T(s+\delta)\big)
  \]
  在 \(s\) 接近 1 时自动退化为后向差分。

> 关键点：`CartesianSpline` 内部用 6 维向量 \([\,\mathbf{p};\ \mathbf{r}\,]\)（位置 + **旋转向量** rotation-vector）插值。因此 `query_state().twist` 的角分量是 \(\dot{\mathbf{r}}\)（旋转向量导数），**不是**真正的角速度 \(\boldsymbol{\omega}\)。该细节在下文 MPC 缺陷中再次出现。而 `tangent_at_progress` 因为先重建位姿再取 \(SE(3)\) log，所以不受该表示法影响，是数值稳健的。

---

## 1. SerialLinkMPC

### 1.1 算法

**状态 / 控制 / 模型。** 状态为基座系下的绝对位姿，控制为末端 twist：
\[
\mathbf{x}=\begin{bmatrix}\mathbf{p}\\ \mathbf{r}\end{bmatrix}\in\mathbb{R}^6,
\qquad
\mathbf{u}=\begin{bmatrix}\mathbf{v}\\ \boldsymbol{\omega}\end{bmatrix}\in\mathbb{R}^6,
\]
其中 \(\mathbf{r}\) 为姿态的旋转向量。离散模型为单积分器：
\[
\mathbf{x}_{k+1}=A\,\mathbf{x}_k+B\,\mathbf{u}_k,
\qquad A=I_6,\quad B=\Delta t\,I_6 .
\]

**凝聚（condensed）预测。** 在时域 \(N\) 上堆叠 \(X=[\mathbf{x}_1;\dots;\mathbf{x}_N]\)、\(U=[\mathbf{u}_0;\dots;\mathbf{u}_{N-1}]\)：
\[
X = A_x\,\mathbf{x}_0 + B_u\,U,
\qquad
(A_x)_k = A^{k+1},\quad
(B_u)_{k,j}=A^{\,k-j}B\ \ (j\le k) .
\]
由于 \(A=I\)，实际上 \(A_x=[\,I;\dots;I\,]\)、\(B_u\) 为分块下三角且每块为 \(\Delta t\,I\)。

**代价。** 取分块对角权重 \(Q_b=\mathrm{blkdiag}(Q)\)、\(R_b=\mathrm{blkdiag}(R)\)，
\[
Q=\mathrm{diag}(w_p,w_p,w_p,\ w_o,w_o,w_o),\quad
R=\mathrm{diag}(w_v,w_v,w_v,\ w_\omega,w_\omega,w_\omega),
\]
最小化
\[
J=\tfrac12\,(X-X_{ref})^\top Q_b (X-X_{ref}) + \tfrac12\,(U-U_{ref})^\top R_b (U-U_{ref}).
\]
代入预测式得标准 QP \( \min_U \tfrac12 U^\top H U + f^\top U\)：
\[
H = B_u^\top Q_b B_u + R_b + \epsilon I,\qquad
f = B_u^\top Q_b (A_x\mathbf{x}_0 - X_{ref}) - R_b\,U_{ref}.
\]

**约束。** 逐轴盒约束 \(\mathbf{u}_{\min}\le \mathbf{u}_k\le \mathbf{u}_{\max}\)，写成 \(B_{ineq}U\le z_{ineq}\) 交给 `QPSolver`。求得 \(U^\star\) 后取首步 \(\mathbf{u}_0=U^\star_{0:6}\)，再做一次安全限幅。

### 1.2 函数功能与调用（逐函数）

**`track_endpoint_trajectory_at_time(const double &time)`** —— 时间索引轨迹的主入口（推荐用法）。
- *前置*：必须先 `set_trajectory()`，否则抛 `runtime_error`。
- *输入*：当前控制时刻 `time`（会被 `clamp` 到 `[start_time, end_time]` 得到 \(t_0\)）。
- *步骤*：① 调 `update()` 刷新雅可比与末端位姿；② 取当前末端位姿构造 \(\mathbf{x}_0=[\mathbf{p};\ \mathbf{r}]\)；③ 对 \(k=0\dots N-1\)，在 `stateTime=t_0+(k{+}1)\Delta t` 处 `query_state()` 得状态参考 \(\mathbf{x}_{ref,k}\)、在 `controlTime=t_0+k\Delta t` 处得控制参考 \(\mathbf{u}_{ref,k}=\)`controlReference.twist`；④ 调 `solveMPC()`；⑤ 对首步 \(\mathbf{u}_0\) 逐轴限幅到 \(\pm v_{\max},\pm\omega_{\max}\)。
- *输出*：关节速度 \(\dot{\mathbf{q}}\)（经继承的 `resolve_endpoint_twist(u0)`）。
- *调用链*：`update()` → `query_state()`×2N → `solveMPC()` → `resolve_endpoint_twist()`。

**`track_endpoint_trajectory(desiredPose, desiredVelocity, desiredAcceleration)`** —— 兼容（单点）入口，覆写自基类。
- *输入*：单个参考位姿 + 期望 twist（`desiredAcceleration` 被忽略，速度级控制不用）。
- *步骤*：因为只有一个采样点，用 `desiredVelocity` 把该点**线性外推**成整段本地时域：\(\mathbf{x}_{ref,k}=\mathbf{x}_{ref,0}+(k{+}1)\Delta t\,\mathbf{v}_{des}\)，\(\mathbf{u}_{ref,k}=\mathbf{v}_{des}\)；其余同上。
- *用途*：当上层没有完整 `CartesianSpline`、只能逐拍喂参考时的降级模式。

**`solveMPC(x0, xRefStack, uRefStack)`** —— 核心求解器（堆叠版）。
- *输入*：初始状态 \(\mathbf{x}_0\)、长度 \(N\) 的状态参考栈（对应 \(\mathbf{x}_1\dots\mathbf{x}_N\)）与控制参考栈（对应 \(\mathbf{u}_0\dots\mathbf{u}_{N-1}\)）；栈不足 \(N\) 时用末元素填充。
- *步骤*：构造 \(A_x,B_u\) → 组 \(Q_b,R_b\) → 组 \(H,f\) → 组盒约束 \(B_{ineq},z_{ineq}\) → 用 `_warmStart` 调 `_qpSolver.solve()`；求解或返回尺寸异常时**回退**为 \(U_{ref}\)；更新 `_warmStart`。
- *输出*：首步 \(\mathbf{u}_0=U^\star_{0:6}\)。

**`solveMPC(x0, xRef, uRef)`** —— 单点重载：把同一 \((\mathbf{x}_{ref},\mathbf{u}_{ref})\) 铺满 \(N\) 步后转调堆叠版。

**`set_trajectory(traj)` / `clear_trajectory()` / `has_trajectory()`** —— 存/清除时间索引轨迹；任一改动都 `resize(0)` 清空 `_warmStart`（避免跨轨迹的脏热启动）。

**`resolve_endpoint_motion()` / `resolve_endpoint_twist()` / `track_joint_trajectory()`** —— 继承自共同的 `SerialLinkVelocityBase`。

**`quaternion_to_rotation_vector(q)`**（静态）—— 把四元数转成最短旋转向量 \(\theta\hat{\mathbf{a}}\)：归一化、强制 \(w\ge0\)（取最短路径）、对 NaN/近零向量返回 0、\(\theta=2\,\mathrm{atan2}(\|\mathbf{q}_v\|,w)\)。仅用于构造状态的姿态分量。

### 1.3 缺陷

1. **【严重】姿态误差用「绝对旋转向量直接相减」，几何不正确。** 代价中的姿态误差是 \(\mathbf{r}_{ref}-\mathbf{r}\)（两个绝对旋转向量之差，[SerialLinkMPC.cpp:136](../Control/src/SerialLinkMPC.cpp#L136) / [153](../Control/src/SerialLinkMPC.cpp#L153)），这**不是**测地姿态误差。在 \(\pm\pi\) 附近会跳变并产生巨大伪误差（例如 \(\mathbf{r}=[0,0,3]\) 与 \(\mathbf{r}_{ref}=[0,0,-3]\) 实际只差 \(\approx0.28\,\mathrm{rad}\)，相减却得 \(6.0\)）。应改为在误差坐标系中建模，姿态误差取
   \[
   \mathbf{e}_R=\log_{SO(3)}\!\big(R_{ref}R^\top\big),
   \]
   即 MPCC / RMPCC 已采用的做法。

2. **【中】角速度语义不一致。** 模型把 \(\mathbf{u}\) 的角分量当作 \(\dot{\mathbf{r}}\)（与轨迹 `twist` 的角分量一致，自洽），但最终 \(\mathbf{u}_0\) 又被当作真实角速度 \(\boldsymbol{\omega}\) 送入 `resolve_endpoint_twist`。小误差下 \(\dot{\mathbf{r}}\approx\boldsymbol{\omega}\) 成立，大姿态偏差时不一致。位置项 \(\dot{\mathbf{p}}=\mathbf{v}\) 精确，无此问题。

3. **【小】预测矩阵冗余。** [solveMPC](../Control/src/SerialLinkMPC.cpp#L324-L345) 完整展开 \(A^{k}\)、\(\Phi\)，但 \(A=I\) 时全部折叠为 \(\Delta t\,I\)，是无谓计算。

---

## 2. SerialLinkMPCC

### 2.1 算法

**思路。** 控制器通过 `set_trajectory()` 持有完整 `CartesianSpline`。每个控制周期由 `step(dt, estimatedProgress)` 接收外部最近点估计 \(\hat s\)，以 \(T(\hat s)\) 建立局部坐标系，并沿真实轨迹曲率预测未来误差。

**局部坐标变换。** 记参考姿态 \(R_{ref}\)、位置 \(\mathbf{p}_{ref}\)，则世界→局部：
\[
R_{TW}=R_{ref}^\top,\qquad \mathbf{p}_{TW}=-R_{TW}\,\mathbf{p}_{ref}.
\]
当前末端在局部系的位置与姿态误差：
\[
\mathbf{p}_T = R_{TW}\mathbf{p}_{cur}+\mathbf{p}_{TW},\qquad
\mathbf{r}_T = \log_{SO(3)}\big(R_{TW}R_{cur}\big).
\]

**误差 / 控制。**
\[
\mathbf{e}=\begin{bmatrix}\mathbf{p}_T\\ \mathbf{r}_T\end{bmatrix}\in\mathbb{R}^6,
\qquad
\mathbf{u}=\begin{bmatrix}\mathbf{v}_T\\ \boldsymbol{\omega}_T\\ v_s\end{bmatrix}\in\mathbb{R}^7,
\]
其中 \(v_s=\dot s\) 是虚拟进度控制。利用 warm start 得到名义 \(s_k\)，在每个 stage 查询轨迹切向 \(\tau(s_k)\)，线性化模型为
\[
\mathbf e_{k+1}\approx\mathbf e_k+\Delta t\big(\mathbf u_{k,1:6}-\tau(s_k)\dot s_k\big).
\]

**轮廓 / 滞后加权。** 用参考线速度方向 \(\hat{\mathbf{t}}\) 构造正交基 \([\,\mathbf{n}_1,\mathbf{n}_2,\hat{\mathbf{t}}\,]=R_{path}\)，位置权重各向异性：
\[
Q_{pos}=R_{path}\,\mathrm{diag}(w_c,\,w_c,\,w_l)\,R_{path}^\top,
\]
即垂直路径方向（contour）重罚 \(w_c\)、沿路径方向（lag）轻罚 \(w_l\)。姿态用各向同性 \(w_o\)。

**代价项。** 含状态项、输入正则 \(R_k\)、增量正则
\(\;\|\,E_d U\,\|^2_{R_{\Delta u}}\)（\(E_d\) 为一阶差分算子，偏置围绕上次控制 \(\mathbf{u}_{last}\)）、路径速度一致性
\(\;\|\mathbf u_{1:6}-\tau(s)\dot s\|^2_{R_v}\)，以及线性进度奖励 \(f_{s}\mathrel{-}=q_{rew}\Delta t\)。
若堆叠预测写成 \(E=A_uU+b\)，则误差项为 \(A_u^\top Q A_u\) 与 \(A_u^\top Qb\)，再叠加输入、增量、速度参考和进度奖励。盒约束对线/角速度对称，对进度速度为 \(v_{s,\min}\le v_s\le v_{s,\max}\)，并约束预测域累计进度不超过 \(1-\hat s\)。求解后取首步前 6 维旋回基座系：
\[
\mathbf{v}_B=R_{TW}^\top\mathbf{v}_T,\qquad \boldsymbol{\omega}_B=R_{TW}^\top\boldsymbol{\omega}_T,
\]
下一拍当前进度重新由外部 \(\hat s\) 校准，QP 内的 \(\dot s_k\) 只用于未来预测。

### 2.2 函数功能与调用（逐函数）

**`set_trajectory(traj)`** —— 复制统一轨迹，清空 warm start，并由 \(1/(t_{end}-t_{start})\) 初始化名义进度速率。

**`step(dt, estimatedProgress)`** —— 正式入口。钳位外部 \(\hat s\)，刷新机器人状态，在 \(T(\hat s)\) 的局部系构造当前误差，求解后把首步局部 twist 旋回基座系并交给 `resolve_endpoint_twist()`。

**`solve_mpcc(error0, referenceRotation)`** —— 用 warm start 生成名义进度序列；逐 stage 查询真实轨迹姿态和切向；组装 \(e_{k+1}\) 预测、contour/lag 权重、输入/平滑/前馈代价以及速度和终点约束；求解并滚动 warm start。

**`track_endpoint_trajectory(...)`** —— 仅为满足基类接口，始终抛 `logic_error`，避免退化回单点速度跟踪。

### 2.3 当前边界

MPCC 仍采用局部欧氏一阶误差模型；曲率较大或姿态误差较大时，RMPCC 的 \(SE(3)\) 模型更合适。这是两种算法的预期差别，不再是轨迹接口差异。

---

## 3. SerialLinkRMPCC

### 3.1 算法

**定位。** 控制器持有参考路径，当前进度由外部最近点估计 \(\hat s\) 提供，直接在 \(SE(3)\) 上做黎曼 MPCC。

**误差。** 在带刚性扰动 \(D\) 的参考帧下，定义 body 系误差
\[
T_{ref}(s)=D\,T_{traj}(s),\qquad
\mathbf{e}_0=\log_{SE(3)}\!\big(T_{ref}(s)^{-1}\,T_{cur}\big)\in\mathbb{R}^6.
\]

**预测（误差动力学一阶线性化）。** 决策变量 \(\mathbf{z}=[\mathbf{u}_1,\dots,\mathbf{u}_N,\ \dot s_1,\dots,\dot s_N]\in\mathbb{R}^{7N}\)。参考自身以 body 切向 \(\boldsymbol{\tau}(s)\) 运动，故相对误差速率 \(\approx \mathbf{u}-\boldsymbol{\tau}(s)\dot s\)：
\[
\mathbf{e}_j \approx \mathbf{e}_0 + \Delta t\sum_{i\le j}\big(\mathbf{u}_i-\boldsymbol{\tau}(s_i)\,\dot s_i\big).
\]
写成 \(\mathbf{e}=A\mathbf{z}+\mathbf{b}\)（\(\mathbf{b}\) 各段为 \(\mathbf{e}_0\)），其中 \(\boldsymbol{\tau}\) 在 warm-start 预测的进度点 \(s_i\) 处求值。

**黎曼 contour / lag 分解。** 给定度量 \(M\)，沿切向的（斜）投影与正交投影：
\[
P_{lag}=\frac{\boldsymbol{\tau}\,\boldsymbol{\tau}^\top M}{\boldsymbol{\tau}^\top M\boldsymbol{\tau}},\qquad
P_{c}=I-P_{lag},\qquad
\ell^\top=\frac{\boldsymbol{\tau}^\top M}{\sqrt{\boldsymbol{\tau}^\top M\boldsymbol{\tau}}}.
\]
每段权重（末段乘终端系数 \(\gamma\)）：
\[
Q_{stage}=P_c^\top W_c P_c + w_{lag}\,\ell\,\ell^\top .
\]

**代价。** 误差跟踪 \( \mathbf{e}^\top Q\,\mathbf{e}\) 加上：控制正则 \(W_u\)、路径速度跟踪 \(\|\mathbf{u}-\boldsymbol{\tau}\dot s_{ref}\|^2_{W_{pv}}\)、控制增量平滑 \(W_{\Delta u}\)、进度速率跟踪/平滑 \(w_{\dot s},w_{\Delta\dot s}\) 以及线性进度奖励 \(-r\,\dot s\)。统一为
\[
H=2A^\top QA+(\text{各正则块})+\lambda I,\qquad
f=2A^\top Q\,\mathbf{b}+(\text{各参考项}).
\]

**约束。** 线/角速度盒约束、进度速率 \(\dot s\in[\dot s_{\min},\dot s_{\max}]\)，外加两条首步约束：
\[
\Delta t\,\dot s_1\le (1-s)+\text{slack}\quad(\text{不越过终点}),\qquad
\Delta t\,\dot s_1\le s_{limit}+\text{slack}-s\quad(\text{时钟调度上限}).
\]

**出口。** 取首步 body twist \(\mathbf{u}_1\)。由于机器人 Jacobian 的线速度语义是末端点速度 \(\dot p\)，这里只旋转线/角分量到基座系：
\[
\boldsymbol{\omega}_W=R\,\boldsymbol{\omega}_b,\qquad
\dot{\mathbf{p}}_W=R\,\mathbf{v}_b,
\]
可选叠加外环位姿反馈 \(K\,\mathbf{e}_{pose}\)。QP 内的进度速率只预测未来；下一拍当前进度由新的外部 \(\hat s\) 校准。

### 3.2 函数功能与调用（逐函数）

**`step(dt, estimatedProgress)`** —— 单步主入口。

- **前置**：必须先 `set_trajectory()`，否则抛 `runtime_error`。
- **输入**：距上一步的时间 `dt` 与外部最近点估计 `estimatedProgress`。
- **步骤**：① 用外部估计同步当前 \(s\)；② `update()`；③ 取当前末端位姿 \(T_{cur}\)；④ 调 `solve_rmpcc()` 得首步 body twist 与进度速率；⑤ 将 body twist 的线/角分量旋转到基座系端点速度；⑥ 可选叠加外环位姿反馈并更新 warm-start 记忆。
- **输出**：关节速度 \(\dot{\mathbf{q}}\)。
- **调用链**：`update()` → `solve_rmpcc()`（内含 `reference_transform`/`tangent_at_progress`/`_qpSolver.solve`）→ `reference_transform()` → `pose_feedback_error()` → `resolve_endpoint_twist()`。

**`solve_rmpcc(currentTransform, dt)`** —— 核心求解器（\(7N\) 维 QP），无返回值，结果落在 `_warmStart` 与 `_diagnostics`。

- **步骤**：① 算初始 body 误差 \(\mathbf{e}_0=\log_{SE(3)}(T_{ref}^{-1}T_{cur})\)；② 用 warm-start 的 \(\dot s\) 猜测把进度沿时域**外推**为 `predictedProgress`（线性化点）；③ 逐 stage 在预测进度处取切向 \(\boldsymbol{\tau}(s)\)，组黎曼 contour/lag 权重 \(Q_{stage}\)（末段乘终端系数），并填预测雅可比 \(A\)（\(+\Delta t I\) 对 \(\mathbf{u}\)，\(-\Delta t\boldsymbol{\tau}\) 对 \(\dot s\)）；④ \(H=2A^\top QA,\ f=2A^\top Q\mathbf{b}\)，再叠加控制/路径速度/控制增量/进度速率与平滑/进度奖励各项；⑤ 组盒约束 + 两条首步进度约束（不越终点、不越调度上限）；⑥ 初始化/裁剪 warm-start → `_qpSolver.solve()`，失败回退到裁剪后的种子；⑦ 滚动 warm-start（整体前移一格），填充全部诊断量。

**`set_trajectory(traj)`** —— 存路径并 `reset()`；若 `autoProgressRate` 则由轨迹时长推导 \(\dot s_{ref}=\mathrm{clamp}(1/(T_{end}-T_{start}),\min,\max)\)，并据此确定 \(\dot s_{\max}\)。

**`set_disturbance(D)`** —— 设刚性工作空间扰动，\(T_{ref}(s)=D\,T_{traj}(s)\)；刚性 \(D\) 只平移 \(\mathbf{e}_0\)，不改 body 切向。

**`set_schedule_limit(s)`** —— 设时钟调度上限（首步进度不得越过 \(s+\)slack），`clamp` 到 \([0,1]\)；传 1.0 关闭。

**`reset()`** —— 进度、warm-start、上次 twist/速率、调度上限、诊断全部复位（换目标时调用）。

**`reference_transform(s)` / `reference_pose(s)`** —— 进度 \(s\) 处（含扰动）的参考变换/位姿；内部调 `_trajectory.pose_at_progress(s)`。

**`body_twist_reference_at_progress(s)`** —— 前馈参考 body twist \(\boldsymbol{\tau}(s)\,\dot s_{ref}\)；调 `_trajectory.tangent_at_progress(s)`。

**`clipped_warm_start(seed, lower, upper, dt, remaining, scheduleRemaining)`** —— 把种子逐元素夹到盒约束，并对**首步进度速率**额外施加「不越终点 / 不越调度上限」两个上限。

**`path_progress()` / `diagnostics()`** —— 查询当前进度 / 最近一步诊断（含 contour/lag 误差、QP 状态、是否回退等）。

**`compute_control_limits(j)`** —— 见 §4（同式，实际不被调用）。

> 注：`pose_to_matrix` / `matrix_to_pose` / `pose_feedback_error` 是该 `.cpp` 内的匿名命名空间辅助函数，分别做位姿↔4×4 矩阵互转、世界系 PD 位姿误差 \([\,\mathbf{p}_{des}-\mathbf{p};\ \log_{SO(3)}(R_{des}R^\top)\,]\)。

### 3.3 当前边界

预测仍是一阶线性化，名义 \(s_i\) 来自 warm start；轨迹曲率或误差很大时可能需要迭代线性化。这是当前实时 QP 实现的明确近似。

---

## 4. 三者共同点

1. **出口统一。** 所有任务空间指令都经继承的 `SerialLinkVelocityBase::resolve_endpoint_twist()` 出口，twist→\(\dot{\mathbf{q}}\) 的 QP 与关节限位由共同速度控制层处理（参见 [`control_twist_to_joint_velocity_summary.md`](control_twist_to_joint_velocity_summary.md)）。

2. **QP 形式统一。** 三者都化为 \(\min_z \tfrac12 z^\top H z+f^\top z\) s.t. \(B_{ineq}z\le z_{ineq}\)，统一交给 `QPSolver<double>`，并带 warm-start 与求解失败回退。

3. **轨迹语义统一。** MPCC 与 RMPCC 都持有 `CartesianSpline`，都以外部 \(\hat s\) 为当前进度，并通过 `pose_at_progress()` / `tangent_at_progress()` 查询未来几何。

---

## 5. 修改优先级建议

> 更新（2026-06）：支撑这些修复的 SE(3)/SO(3) 工具已下沉到 Math，见 [`math_lie_group_toolkit.md`](math_lie_group_toolkit.md)（含 `so3_left/right_jacobian(_inverse)`、`adjoint`、统一的 `quaternion_*`，并经数值验证）。

1. 修 **MPC 的姿态误差**：用 Math 的 SO(3) Jacobian 正确处理旋转向量速率与角速度。
2. 为外部最近点估计增加独立测试，验证噪声、路径自交和进度单调策略；MPCC/RMPCC 都只消费其输出的 \(\hat s\)。
3. 若大曲率下单次线性化精度不足，再考虑对 MPCC/RMPCC 增加迭代线性化；没有实测问题前不增加复杂度。
