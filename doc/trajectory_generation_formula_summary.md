# 轨迹生成方法（公式化总结）

本文档给出该库轨迹生成模块的公式化描述，重点包括：

1. 统一轨迹定义（`TrajectoryBase`）
2. 多维样条轨迹（`SplineTrajectory`）
3. 笛卡尔位姿样条（`CartesianSpline`）
4. 梯形速度轨迹（`TrapezoidalVelocity`）
5. 基于进度/弧长变量 \(\mathbf{r}(s)\) 的统一写法（用于 MPCC 语境）

---

## 1. 统一轨迹定义（TrajectoryBase）

库中所有轨迹统一为连续时间映射：

\[
\mathbf{x}(t)=
\begin{bmatrix}
\mathbf{p}(t)\\
\dot{\mathbf{p}}(t)\\
\ddot{\mathbf{p}}(t)
\end{bmatrix},\quad t\in[t_0,t_f]
\]

其中状态结构为：

\[
\texttt{State}=\{\texttt{position},\ \texttt{velocity},\ \texttt{acceleration}\}
\]

并满足统一的时间与维度约束：

\[
t_0 < t_f,
\]

\[
\dim(\mathbf{p})=\dim(\dot{\mathbf{p}})=\dim(\ddot{\mathbf{p}}).
\]

---

## 2. 多维样条轨迹（SplineTrajectory）

对 \(m\) 维轨迹按维度独立构造分段样条：

\[
p_j(t)=\mathcal{S}_j(t),\quad j=1,\dots,m
\]

整体写为：

\[
\mathbf{p}(t)=\begin{bmatrix}\mathcal{S}_1(t)&\cdots&\mathcal{S}_m(t)\end{bmatrix}^\top,
\]

\[
\dot{\mathbf{p}}(t)=\begin{bmatrix}\dot{\mathcal{S}}_1(t)&\cdots&\dot{\mathcal{S}}_m(t)\end{bmatrix}^\top,
\quad
\ddot{\mathbf{p}}(t)=\begin{bmatrix}\ddot{\mathcal{S}}_1(t)&\cdots&\ddot{\mathcal{S}}_m(t)\end{bmatrix}^\top.
\]

若输入是“路点位置 + 起始速度”（常见三次样条场景），中间路点速度通过线性方程组估计：

\[
\mathbf{A}\,\mathbf{d}=\mathbf{B}\,\mathbf{y},
\]

其中 \(\mathbf{y}\) 为路点位置，\(\mathbf{d}\) 为路点一阶导（速度）估计，再据此构造各段三次多项式，保证到二阶导连续。

---

## 3. 笛卡尔位姿样条（CartesianSpline）

`CartesianSpline` 先把位姿映射到 6 维向量，再调用 `SplineTrajectory` 插值：

\[
\mathbf{z}(t)=
\begin{bmatrix}
\mathbf{r}(t)\\
\boldsymbol{\phi}(t)
\end{bmatrix}\in\mathbb{R}^6,
\]

其中 \(\mathbf{r}\) 是平移，\(\boldsymbol{\phi}=\theta\mathbf{u}\) 是旋转向量（axis-angle）。

四元数 \(q=[w,\mathbf{v}]\) 到旋转向量：

\[
\theta=2\arccos(\mathrm{clamp}(w,-1,1)),\qquad
\mathbf{u}=\frac{\mathbf{v}}{\|\mathbf{v}\|},\qquad
\boldsymbol{\phi}=\theta\mathbf{u}.
\]

样条查询后再反变换回四元数：

\[
\theta=\|\boldsymbol{\phi}\|,\quad
\mathbf{u}=\frac{\boldsymbol{\phi}}{\|\boldsymbol{\phi}\|},\quad
q=\left[\cos\frac{\theta}{2},\ \mathbf{u}\sin\frac{\theta}{2}\right].
\]

因此，SE(3) 轨迹被转化为 \(\mathbb{R}^6\) 样条问题处理。

---

## 4. 梯形速度轨迹（TrapezoidalVelocity）

对两点 \(\mathbf{p}_0\to\mathbf{p}_f\)，定义归一化进度 \(s(t)\in[0,1]\)：

\[
\mathbf{p}(t)=(1-s(t))\mathbf{p}_0+s(t)\mathbf{p}_f,
\]

\[
\dot{\mathbf{p}}(t)=\dot{s}(t)(\mathbf{p}_f-\mathbf{p}_0),
\quad
\ddot{\mathbf{p}}(t)=\ddot{s}(t)(\mathbf{p}_f-\mathbf{p}_0).
\]

其中 \(s,\dot{s},\ddot{s}\) 采用三段式梯形速度律：

1. 加速段 \((0\le\tau<t_r)\)
\[
s=\frac12 a_n\tau^2,\quad \dot{s}=a_n\tau,\quad \ddot{s}=a_n.
\]

2. 匀速段 \((t_r\le\tau<t_r+t_c)\)
\[
s=s_r+v_n(\tau-t_r),\quad \dot{s}=v_n,\quad \ddot{s}=0.
\]

3. 减速段 \((t_r+t_c\le\tau\le t_f)\)
\[
t'=\tau-t_r-t_c,
\]
\[
s=s_r+s_c+v_n t'-\frac12 a_n t'^2,\quad
\dot{s}=v_n-a_n t',\quad
\ddot{s}=-a_n.
\]

多路点轨迹通过多个两点梯形段串接生成。

---

## 5. 基于弧长/进度变量 \(\mathbf{r}(s)\) 的统一描述（MPCC 语境）

将几何路径写为：

\[
\mathbf{r}(s)\in\mathbb{R}^m,\quad s\in[0,L],
\]

其中 \(s\) 是进度（可取弧长），\(L\) 为总路径长度。

时间演化通过：

\[
\dot{s}=v_s,\qquad \ddot{s}=a_s.
\]

于是：

\[
\mathbf{r}(t)=\mathbf{r}(s(t)),
\]

\[
\dot{\mathbf{r}}(t)=\mathbf{r}'(s)\dot{s},
\qquad
\ddot{\mathbf{r}}(t)=\mathbf{r}''(s)\dot{s}^2+\mathbf{r}'(s)\ddot{s}.
\]

在 MPCC 中常将进度作为扩展状态：

\[
\mathbf{x}=
\begin{bmatrix}
\mathbf{x}_{pose}\\
s
\end{bmatrix},
\quad
\mathbf{u}=
\begin{bmatrix}
\mathbf{u}_{motion}\\
v_s
\end{bmatrix},
\]

离散化：

\[
s_{k+1}=s_k+\Delta t\,v_{s,k},
\]

并约束：

\[
0\le s_k\le L,
\qquad
v_s^{\min}\le v_{s,k}\le v_s^{\max}.
\]

---

## 6. 可直接放论文的一段话

本文轨迹生成采用“统一状态接口 + 几何与时标解耦”的框架：在底层以 `TrajectoryBase` 统一表示 \((\mathbf{p},\dot{\mathbf{p}},\ddot{\mathbf{p}})\)；在几何层分别采用多维分段样条（`SplineTrajectory`）、SE(3) 位姿样条（`CartesianSpline`）与梯形速度参数化（`TrapezoidalVelocity`）构造参考路径；在时标层引入进度/弧长变量 \(s\)，将参考写为 \(\mathbf{r}(s(t))\)，从而可在 MPCC 中直接优化路径推进速度并保持轨迹几何不变。

---

## 7. 技术细节补充（建议写入正文或附录）

### 7.1 坐标系与变量语义

建议在文中显式声明：

- \(W\)：世界/基坐标系（或你定义的主参考系）
- \(T\)：路径切向局部坐标系（path-aligned frame）
- \(B\)：机械臂控制基坐标系（若与 \(W\) 不同需明确）

并给出速度变量所在坐标系，例如 \(\mathbf{v}_T,\boldsymbol{\omega}_T\) 与 \(\mathbf{v}_B,\boldsymbol{\omega}_B\) 的关系。

### 7.2 姿态参数化数值细节

在 `CartesianSpline` 中，四元数到旋转向量转换存在小角度阈值处理：当 \(|\theta|\) 足够小时近似零旋转，避免除零和方向抖动。建议在文中注明该阈值策略用于数值稳定。

### 7.3 样条边界条件说明

建议明确：

- 样条支持“给定路点状态（位置/速度/加速度）”的直接构造；
- 仅给位置时，速度通过线性系统估计；
- 常用边界为给定初速度、末速度设为零（机械臂停靠更安全）。

### 7.4 进度变量的可行域与饱和

虽然写作 \(s\in[0,L]\)，实现时建议说明是否采用：

- 硬约束截断（clamp），
- 软惩罚，
- 或两者结合。

并注明 \(v_s\) 的上下界与标称值（如 \(v_s^{nom}\)）关系。

### 7.5 参数复现性

建议在文末或附录提供“轨迹参数表”：

- 插值阶次、关键阈值；
- \(v_{max},a_{max}\)（梯形速度）；
- 采样周期 \(\Delta t\) 与控制频率。

没有参数表会显著影响复现性。
