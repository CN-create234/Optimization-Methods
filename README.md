# 最优化方法考研复习讲义·详细版

## 第一章 基本概念

### 1.1 最优化问题简介

**最优化问题的定义**：从所有可能方案中选择使某一指标达到最优的方案，即在满足一定约束的条件下，寻找使目标函数达到最小值（或最大值）的决策变量取值。

**最优化问题的数学模型**：

$$
\begin{aligned}
& \min_{\boldsymbol{x} \in \mathbb{R}^n} \quad f(\boldsymbol{x}) \\
& \text{s.t.} \quad c_i(\boldsymbol{x}) = 0, \quad i = 1, 2, \dots, p \\
& \qquad c_i(\boldsymbol{x}) \leq 0, \quad i = p+1, \dots, p+q
\end{aligned}
$$

其中：
- **目标函数**：\( f(\boldsymbol{x}) \)，需要最小化（或最大化）的函数
- **决策变量**：\( \boldsymbol{x} = (x_1, x_2, \dots, x_n)^T \)，待确定的未知量
- **约束条件**：等式约束和不等式约束，共同定义**可行域** \( D \)

**分类标准**：

| 分类维度 | 类型 | 说明 |
|----------|------|------|
| 有无约束 | **无约束优化** | 决策变量无任何限制，在整个 \( \mathbb{R}^n \) 上搜索 |
| | **约束优化** | 存在等式或不等式约束 |
| 函数线性性 | **线性规划（LP）** | 目标函数和约束条件均为线性函数 |
| | **非线性规划（NLP）** | 至少有一个函数是非线性的 |
| 凸性 | **凸规划** | 目标函数为凸函数，可行域为凸集 |
| | **非凸规划** | 不满足凸性条件 |
| 变量类型 | **连续优化** | 决策变量取连续值 |
| | **整数规划** | 部分或全部变量取整数值 |
| | **混合整数规划** | 连续变量与整数变量并存 |

**全局最优与局部最优**：
- **全局最小值点**：对可行域内任意 \( \boldsymbol{x} \)，都有 \( f(\boldsymbol{x}) \geq f(\boldsymbol{x}^*) \)
- **局部极小值点**：存在邻域 \( U(\boldsymbol{x}^*, \delta) \)，使邻域内任意可行点 \( \boldsymbol{x} \) 满足 \( f(\boldsymbol{x}) \geq f(\boldsymbol{x}^*) \)
- **重要性质**：全局最小值点一定是局部极小值点，但反之不成立

**最优化的数学基础需求**：
- 线性代数：向量空间、矩阵运算、特征分解、范数
- 多元微积分：梯度、海色矩阵、泰勒展开、链式法则

### 1.2 凸集和凸函数

凸性是最优化理论中最核心的概念，凸规划具有良好性质——局部最优即为全局最优。

**凸集的定义**：设 \( C \subseteq \mathbb{R}^n \)，若对任意 \( \boldsymbol{x}, \boldsymbol{y} \in C \) 和任意 \( \theta \in [0,1] \)，都有：

$$
\theta \boldsymbol{x} + (1-\theta) \boldsymbol{y} \in C
$$

则称 \( C \) 为凸集。

**几何理解**：凸集中任意两点的连线仍在该集合内。

**常见凸集**：
- 超平面：\( \{ \boldsymbol{x} \mid \boldsymbol{a}^T \boldsymbol{x} = b \} \)
- 半空间：\( \{ \boldsymbol{x} \mid \boldsymbol{a}^T \boldsymbol{x} \leq b \} \)
- 多面体：有限个半空间的交集
- 球体：\( \{ \boldsymbol{x} \mid \|\boldsymbol{x} - \boldsymbol{x}_0\| \leq r \} \)
- 椭球体

**凸函数的定义**：设 \( f: C \to \mathbb{R} \)，\( C \) 为凸集。若对任意 \( \boldsymbol{x}, \boldsymbol{y} \in C \)，\( \theta \in [0,1] \)，有：

$$
f(\theta \boldsymbol{x} + (1-\theta) \boldsymbol{y}) \leq \theta f(\boldsymbol{x}) + (1-\theta) f(\boldsymbol{y})
$$

则称 \( f \) 为凸函数。若上述不等式严格成立（<），称为**严格凸函数**。

**几何理解**：函数图像上任意两点的连线位于函数图像之上。

**凸函数的判定条件**（可微情形）：

一阶条件（\( f \) 可微）：
$$
f(\boldsymbol{y}) \geq f(\boldsymbol{x}) + \nabla f(\boldsymbol{x})^T (\boldsymbol{y} - \boldsymbol{x}), \quad \forall \boldsymbol{x}, \boldsymbol{y} \in C
$$

二阶条件（\( f \) 二阶可微）：
- 对一维函数：\( f''(x) \geq 0 \)
- 对多元函数：海色矩阵 \( \nabla^2 f(\boldsymbol{x}) \) 半正定

**海色矩阵（Hessian Matrix）** ：
$$
\nabla^2 f(\boldsymbol{x}) = \begin{bmatrix}
\frac{\partial^2 f}{\partial x_1^2} & \frac{\partial^2 f}{\partial x_1 \partial x_2} & \cdots & \frac{\partial^2 f}{\partial x_1 \partial x_n} \\
\frac{\partial^2 f}{\partial x_2 \partial x_1} & \frac{\partial^2 f}{\partial x_2^2} & \cdots & \frac{\partial^2 f}{\partial x_2 \partial x_n} \\
\vdots & \vdots & \ddots & \vdots \\
\frac{\partial^2 f}{\partial x_n \partial x_1} & \frac{\partial^2 f}{\partial x_n \partial x_2} & \cdots & \frac{\partial^2 f}{\partial x_n^2}
\end{bmatrix}
$$

若所有二阶偏导数连续，海色矩阵对称。

**强凸函数**：存在 \( m > 0 \)，使 \( \nabla^2 f(\boldsymbol{x}) - mI \) 半正定，即海色矩阵的最小特征值不小于 \( m \)。强凸函数有更好的收敛性质。

**凸规划的性质**：
1. 凸规划的任意局部极小值点必为全局极小值点
2. 凸规划的最优值点集为凸集
3. 若目标函数严格凸，则最优解唯一

### 1.3 最优性条件

最优性条件是判断一个点是否为最优解的理论依据。

#### 无约束优化问题的最优性条件

**一阶必要条件**：若 \( \boldsymbol{x}^* \) 是 \( f \) 的局部极小值点，且 \( f \) 在 \( \boldsymbol{x}^* \) 处可微，则：

$$
\nabla f(\boldsymbol{x}^*) = \boldsymbol{0}
$$

满足该条件的点称为**驻点**（平稳点）。

**二阶必要条件**：若 \( \boldsymbol{x}^* \) 是 \( f \) 的局部极小值点，且 \( f \) 在 \( \boldsymbol{x}^* \) 处二阶可微，则海色矩阵 \( \nabla^2 f(\boldsymbol{x}^*) \) **半正定**。

**二阶充分条件**：若 \( \nabla f(\boldsymbol{x}^*) = \boldsymbol{0} \) 且 \( \nabla^2 f(\boldsymbol{x}^*) \) **正定**，则 \( \boldsymbol{x}^* \) 是严格局部极小值点。

**无约束最优化方法的分类**：

| 方法类型 | 特点 | 典型算法 |
|----------|------|----------|
| **直接法** | 仅利用目标函数值 | 模式搜索、Nelder-Mead单纯形法 |
| **梯度法** | 利用一阶导数（梯度）信息 | 最速下降法、共轭梯度法 |
| **二阶法** | 利用二阶导数（海色矩阵）信息 | 牛顿法 |
| **拟牛顿法** | 用迭代信息近似海色矩阵 | DFP、BFGS算法 |

#### 约束优化问题的最优性条件

**等式约束情况**：考虑 \( \min f(\boldsymbol{x}) \) s.t. \( c_i(\boldsymbol{x}) = 0, i = 1,\dots,p \)

**拉格朗日乘子法**：构造拉格朗日函数：
$$
L(\boldsymbol{x}, \boldsymbol{\lambda}) = f(\boldsymbol{x}) + \sum_{i=1}^p \lambda_i c_i(\boldsymbol{x})
$$

一阶必要条件（在 \( \boldsymbol{x}^* \) 处）：
$$
\nabla_{\boldsymbol{x}} L(\boldsymbol{x}^*, \boldsymbol{\lambda}^*) = \nabla f(\boldsymbol{x}^*) + \sum_{i=1}^p \lambda_i^* \nabla c_i(\boldsymbol{x}^*) = \boldsymbol{0}
$$
$$
c_i(\boldsymbol{x}^*) = 0, \quad i = 1,\dots,p
$$

#### KKT条件（Karush-Kuhn-Tucker）

对于同时包含等式和不等式约束的优化问题：
$$
\min f(\boldsymbol{x}) \quad \text{s.t.} \quad c_i(\boldsymbol{x}) = 0, \quad c_j(\boldsymbol{x}) \leq 0
$$

**KKT条件的完整形式**（假设 \( \boldsymbol{x}^* \) 处满足约束规范，如线性独立约束规范LICQ）：

1. **平稳性**：\( \nabla f(\boldsymbol{x}^*) + \sum_{i} \lambda_i \nabla c_i(\boldsymbol{x}^*) + \sum_{j} \mu_j \nabla c_j(\boldsymbol{x}^*) = \boldsymbol{0} \)

2. **原始可行性**：\( c_i(\boldsymbol{x}^*) = 0, \quad c_j(\boldsymbol{x}^*) \leq 0 \)

3. **对偶可行性**：\( \mu_j \geq 0 \)

4. **互补松弛性**：\( \mu_j c_j(\boldsymbol{x}^*) = 0 \)

互补松弛性的含义：若不等式约束松驰（\( c_j(\boldsymbol{x}^*) < 0 \)），则对应的 \( \mu_j = 0 \)；若 \( \mu_j > 0 \)，则对应的约束必须活跃（\( c_j(\boldsymbol{x}^*) = 0 \)）。

**对偶理论**：原问题的对偶问题在凸优化中具有强对偶性——原问题最优值等于对偶问题最优值。


## 第二章 线性规划

线性规划（Linear Programming）是最基本、应用最广泛的一类最优化问题，目标函数和约束条件均为线性函数。

### 2.1 线性规划问题的数学模型

**标准形式**：

$$
\begin{aligned}
\min \quad & \boldsymbol{c}^T \boldsymbol{x} \\
\text{s.t.} \quad & A\boldsymbol{x} = \boldsymbol{b} \\
& \boldsymbol{x} \geq \boldsymbol{0}
\end{aligned}
$$

其中 \( A \in \mathbb{R}^{m \times n} \)，\( \boldsymbol{b} \in \mathbb{R}^m \)，\( \boldsymbol{c} \in \mathbb{R}^n \)。

**一般形式化标准形的方法**：
- 最大化目标 → 乘以 -1 变为最小化
- 不等式约束 \( \leq \) → 添加松弛变量变为等式
- 不等式约束 \( \geq \) → 减去剩余变量变为等式
- 自由变量 → 用两个非负变量之差表示

**图解法**（二维情形）：在平面上画出约束边界，找出可行域的顶点，目标函数等值线平移至顶点即得最优解。

**线性规划的基本性质**：
- 可行域是凸多面体（有限个半空间交集）
- 若最优解存在，则最优解可在可行域的某个顶点（基本可行解）取到
- 最优解不唯一时，最优解集为凸多面体

### 2.2 单纯形法

单纯形法是求解线性规划问题的核心算法。

**基本思想**：从可行域的一个顶点出发，沿可行域的边移动到使目标函数值下降的相邻顶点，直至到达最优顶点。

**算法步骤**：
1. 将问题化为标准形，找出初始基本可行解（单位矩阵对应）
2. 计算检验数 \( \sigma_j = c_j - \boldsymbol{c}_B^T B^{-1} \boldsymbol{a}_j \)
3. 若所有 \( \sigma_j \geq 0 \)，当前解最优，终止
4. 选取最负的检验数对应的非基变量 \( x_k \) 作为入基变量
5. 按最小比值规则确定离基变量（\( \theta = \min\{ b_i / a_{ik} \mid a_{ik} > 0 \} \)）
6. 用入基变量替换离基变量，更新基矩阵，转步骤2

**退化情形**：当出现退化基本可行解时，可能出现**循环**，可用Bland法则避免（选择最小下标入基/离基）。

### 2.3 对偶理论

每个线性规划问题都存在一个对应的对偶问题。

**原问题与对偶问题的关系**：

| 原问题（最小化） | 对偶问题（最大化） |
|------------------|-------------------|
| \( \min \boldsymbol{c}^T \boldsymbol{x} \) | \( \max \boldsymbol{b}^T \boldsymbol{y} \) |
| s.t. \( A\boldsymbol{x} \geq \boldsymbol{b} \) | s.t. \( A^T \boldsymbol{y} \leq \boldsymbol{c} \) |
| \( \boldsymbol{x} \geq \boldsymbol{0} \) | \( \boldsymbol{y} \geq \boldsymbol{0} \) |

**对偶定理**：
- **弱对偶定理**：对任意可行解，原问题目标值 \( \geq \) 对偶问题目标值
- **强对偶定理**：若原问题有最优解，则对偶问题也有最优解，且最优值相等
- **互补松弛性**：对最优解 \( \boldsymbol{x}^*, \boldsymbol{y}^* \)，有 \( (\boldsymbol{A}\boldsymbol{x}^* - \boldsymbol{b})^T \boldsymbol{y}^* = 0 \) 和 \( (\boldsymbol{A}^T\boldsymbol{y}^* - \boldsymbol{c})^T \boldsymbol{x}^* = 0 \)

**对偶单纯形法**：从对偶可行解出发，逐步达到原始可行解，适用于需要多次求解的情况。

**影子价格**：对偶变量的最优值表示资源的边际价值——单位资源变化引起目标函数最优值的变化量。

### 2.4 整数规划

当决策变量要求取整数值时称为整数规划，这是一类NP难问题。

**0-1规划**：变量只能取0或1，常用于方案选择问题。

**求解方法**：
- **分枝定界法**：将问题分解为子问题求解
- **割平面法**：逐步添加约束排除非整数最优解
- **匈牙利法**：专门求解分配问题


## 第三章 一维搜索

一维搜索（线性搜索）是求解多维最优化问题的基础——在确定搜索方向后，需要沿该方向寻找最优步长。

### 3.1 基本概念

**一维搜索问题**：
$$
\min_{\alpha > 0} \phi(\alpha) = f(\boldsymbol{x}_k + \alpha \boldsymbol{d}_k)
$$

其中 \( \boldsymbol{d}_k \) 为搜索方向，\( \alpha \) 为步长。

**精确线性搜索**：沿方向精确求解最优步长 \( \alpha^* \)

**不精确线性搜索**：只要求步长使目标函数充分下降，以减少计算量。

### 3.2 精确线性搜索方法

**（1）0.618法（黄金分割法）**

适用于单峰函数（在区间上先单调递减后单调递增）。

- 区间缩比：每次保留 \( 0.618 \) 的区间长度
- 黄金比例：\( \sqrt{5} - 1 \approx 0.618 \)
- 每次只需计算一个新的函数值

**（2）Fibonacci法**

- 使用Fibonacci数列确定试探点位置
- 对于给定迭代次数，区间缩比最优
- 缩比：\( 1 / F_n \)（\( F_n \) 为Fibonacci数）

**（3）二分法**

利用导数的符号确定区间：
- 计算中点 \( \alpha_m = (a+b)/2 \)
- 若 \( \phi'(\alpha_m) > 0 \)，搜索区间缩为 \( (a, \alpha_m) \)
- 若 \( \phi'(\alpha_m) < 0 \)，搜索区间缩为 \( (\alpha_m, b) \)

### 3.3 不精确线性搜索

精确搜索计算量大，实际中常采用不精确搜索。

**（1）Armijo准则**：
$$
f(\boldsymbol{x}_k + \alpha \boldsymbol{d}_k) \leq f(\boldsymbol{x}_k) + c_1 \alpha \nabla f(\boldsymbol{x}_k)^T \boldsymbol{d}_k
$$

其中 \( c_1 \in (0,1) \)，通常取 \( 10^{-4} \)。

**（2）Goldstein准则**：
$$
f(\boldsymbol{x}_k) + (1-c)\alpha \nabla f(\boldsymbol{x}_k)^T \boldsymbol{d}_k \leq f(\boldsymbol{x}_k + \alpha \boldsymbol{d}_k) \leq f(\boldsymbol{x}_k) + c \alpha \nabla f(\boldsymbol{x}_k)^T \boldsymbol{d}_k
$$

其中 \( 0 < c < 1/2 \)。

**（3）Wolfe准则**：
Armijo条件 + 曲率条件 \( \nabla f(\boldsymbol{x}_k + \alpha \boldsymbol{d}_k)^T \boldsymbol{d}_k \geq c_2 \nabla f(\boldsymbol{x}_k)^T \boldsymbol{d}_k \)，\( 0 < c_1 < c_2 < 1 \)。

### 3.4 信赖域方法

与线性搜索不同，信赖域方法先确定步长范围（信赖域半径），再在该范围内选择方向和步长。

**基本框架**：
1. 在当前点 \( \boldsymbol{x}_k \) 构造模型 \( m_k(\boldsymbol{s}) \)（通常是二次模型）
2. 在信赖域 \( \|\boldsymbol{s}\| \leq \Delta_k \) 内极小化模型
3. 计算实际下降量与预测下降量的比值 \( \rho_k \)
4. 若 \( \rho_k \) 接近1，扩大信赖域；若接近0或为负，缩小信赖域


## 第四章 无约束最优化方法

### 4.1 最速下降法（梯度下降法）

**思想**：沿负梯度方向 \( \boldsymbol{d}_k = -\nabla f(\boldsymbol{x}_k) \) 搜索，这是局部下降最快的方向。

**迭代公式**：
$$
\boldsymbol{x}_{k+1} = \boldsymbol{x}_k - \alpha_k \nabla f(\boldsymbol{x}_k)
$$

**优缺点**：
- 优点：算法简单，每步只需计算梯度
- 缺点：收敛速度慢（线性收敛），尤其在接近最优时容易出现“锯齿现象”

**收敛性分析**：
- 对于强凸函数，收敛速度为线性
- 条件数越大，收敛越慢

### 4.2 牛顿法

**思想**：利用目标函数的二阶信息（海色矩阵）构造二次模型，一步求得二次模型的最小点。

**牛顿步**：
$$
\boldsymbol{d}_k = -[\nabla^2 f(\boldsymbol{x}_k)]^{-1} \nabla f(\boldsymbol{x}_k)
$$

迭代公式：\( \boldsymbol{x}_{k+1} = \boldsymbol{x}_k + \boldsymbol{d}_k \)

**优缺点**：
- 优点：收敛速度快（二次收敛）
- 缺点：需计算和求逆海色矩阵（计算量大），海色矩阵不正定（非凸函数）时可能发散

**牛顿法的变体**：若海色矩阵不正定，可用修正牛顿法（加正则化 \( \nabla^2 f + \mu I \) 使其正定）。

### 4.3 共轭梯度法

**思想**：构造一组共轭方向，使在有限步内达到二次函数的最优解。

**共轭方向**：对对称正定矩阵 \( A \)，若非零向量 \( \boldsymbol{d}_i, \boldsymbol{d}_j \) 满足 \( \boldsymbol{d}_i^T A \boldsymbol{d}_j = 0 \)，则称它们关于 \( A \) 共轭。

**共轭梯度法（CG）** ：

**初始方向**：\( \boldsymbol{d}_0 = -\nabla f(\boldsymbol{x}_0) \)

**步长**（精确搜索时）：\( \alpha_k = -\frac{\boldsymbol{g}_k^T \boldsymbol{d}_k}{\boldsymbol{d}_k^T A \boldsymbol{d}_k} \)

**方向更新**：
$$
\boldsymbol{d}_{k+1} = -\boldsymbol{g}_{k+1} + \beta_k \boldsymbol{d}_k
$$

其中 \( \beta_k \) 的计算公式：
- **Fletcher-Reeves**：\( \beta_k = \frac{\|\boldsymbol{g}_{k+1}\|^2}{\|\boldsymbol{g}_k\|^2} \)
- **Polak-Ribière**：\( \beta_k = \frac{\boldsymbol{g}_{k+1}^T (\boldsymbol{g}_{k+1} - \boldsymbol{g}_k)}{\|\boldsymbol{g}_k\|^2} \)

**重要性质**：
- 对n维二次函数，最多n步收敛
- 对非二次函数，可结合重启动策略

### 4.4 拟牛顿法

拟牛顿法的核心思想：用迭代信息近似海色矩阵（或其逆），既获得超线性收敛速度，又避免了计算二阶导数。

**拟牛顿条件**（割线条件）：
$$
\boldsymbol{H}_{k+1} \boldsymbol{s}_k = \boldsymbol{y}_k
$$

其中 \( \boldsymbol{s}_k = \boldsymbol{x}_{k+1} - \boldsymbol{x}_k \)，\( \boldsymbol{y}_k = \nabla f(\boldsymbol{x}_{k+1}) - \nabla f(\boldsymbol{x}_k) \)。

**DFP公式**（Davidon-Fletcher-Powell）：

$$
\boldsymbol{H}_{k+1} = \boldsymbol{H}_k - \frac{\boldsymbol{H}_k \boldsymbol{y}_k \boldsymbol{y}_k^T \boldsymbol{H}_k}{\boldsymbol{y}_k^T \boldsymbol{H}_k \boldsymbol{y}_k} + \frac{\boldsymbol{s}_k \boldsymbol{s}_k^T}{\boldsymbol{y}_k^T \boldsymbol{s}_k}
$$

**BFGS公式**（Broyden-Fletcher-Goldfarb-Shanno）——最广泛使用的拟牛顿法：

$$
\boldsymbol{H}_{k+1} = \left( I - \frac{\boldsymbol{s}_k \boldsymbol{y}_k^T}{\boldsymbol{y}_k^T \boldsymbol{s}_k} \right) \boldsymbol{H}_k \left( I - \frac{\boldsymbol{y}_k \boldsymbol{s}_k^T}{\boldsymbol{y}_k^T \boldsymbol{s}_k} \right) + \frac{\boldsymbol{s}_k \boldsymbol{s}_k^T}{\boldsymbol{y}_k^T \boldsymbol{s}_k}
$$

**DFP与BFGS的比较**：
- DFP是第一个拟牛顿公式，性能良好
- BFGS通常表现更优，数值稳定性更好

拟牛顿法具有**超线性收敛**速度，是无约束优化中最实用的方法之一。


## 第五章 约束最优化方法

### 5.1 约束优化问题的分类

| 问题类型 | 特点 | 典型方法 |
|----------|------|----------|
| 线性规划（LP） | 目标和约束均为线性 | 单纯形法、内点法 |
| 二次规划（QP） | 目标为二次函数，约束线性 | 有效集法、内点法 |
| 一般非线性规划 | 至少有一个非线性函数 | 罚函数法、SQP |
| 凸规划 | 目标凸，可行域凸 | KKT条件有效 |

### 5.2 罚函数法

**基本思想**：将约束优化问题转化为一系列无约束优化问题，对违反约束的行为施加惩罚。

**（1）二次罚函数法（外点法）** ：

将约束以惩罚项形式加入目标函数：
$$
P(\boldsymbol{x}; \rho) = f(\boldsymbol{x}) + \frac{\rho}{2} \sum_{i=1}^p c_i(\boldsymbol{x})^2 + \frac{\rho}{2} \sum_{j=p+1}^{p+q} [\max(0, c_j(\boldsymbol{x}))]^2
$$

其中 \( \rho > 0 \) 为惩罚参数。当 \( \rho \to \infty \) 时，罚问题的解趋近于原问题的解。

**（2）内点障碍函数法** ：

**基本思想**：在可行域内部搜索，用障碍项阻止迭代点接近可行域边界：

$$
B(\boldsymbol{x}; \mu) = f(\boldsymbol{x}) - \mu \sum_{j} \ln(-c_j(\boldsymbol{x}))
$$

其中 \( \mu > 0 \)，当 \( \mu \to 0^+ \) 时，障碍问题的解趋近于原问题解。适用于不等式约束问题。

### 5.3 序列二次规划（SQP）

SQP是求解约束非线性规划最有效的算法之一。

**基本思想**：在每次迭代中，用二次规划子问题近似原问题：
- 目标函数：二次近似（使用拉格朗日函数的海色矩阵）
- 约束条件：线性近似（一阶泰勒展开）

**SQP子问题**（在 \( \boldsymbol{x}_k \) 处）：
$$
\min_{\boldsymbol{d}} \quad \frac{1}{2} \boldsymbol{d}^T \boldsymbol{B}_k \boldsymbol{d} + \nabla f(\boldsymbol{x}_k)^T \boldsymbol{d}
$$
$$
\text{s.t.} \quad c_i(\boldsymbol{x}_k) + \nabla c_i(\boldsymbol{x}_k)^T \boldsymbol{d} = 0, \quad i \in \mathcal{E}
$$
$$
\qquad c_j(\boldsymbol{x}_k) + \nabla c_j(\boldsymbol{x}_k)^T \boldsymbol{d} \leq 0, \quad j \in \mathcal{I}
$$

其中 \( \boldsymbol{B}_k \) 是拉格朗日函数海色矩阵的近似（用拟牛顿法更新）。

**SQP的优势**：
- 收敛速度快（超线性）
- 能有效处理大规模非线性约束问题
- 可直接利用KKT条件构造子问题


## 第六章 特殊问题

### 6.1 线性与非线性最小二乘

**问题描述**：
$$
\min_{\boldsymbol{x}} \quad \frac{1}{2} \sum_{i=1}^m r_i(\boldsymbol{x})^2 = \frac{1}{2} \| \boldsymbol{r}(\boldsymbol{x}) \|^2
$$

其中 \( \boldsymbol{r}(\boldsymbol{x}) = (r_1(\boldsymbol{x}), \dots, r_m(\boldsymbol{x}))^T \) 为残差向量。

**线性最小二乘**：\( r_i(\boldsymbol{x}) = \boldsymbol{a}_i^T \boldsymbol{x} - b_i \)，可用正规方程或QR分解求解。

**非线性最小二乘**：残差非线性，常用特殊方法。

**Gauss-Newton法**：

利用残差结构避免计算二阶导数。每次迭代，用线性最小二乘近似非线性残差：
$$
\boldsymbol{x}_{k+1} = \boldsymbol{x}_k + \boldsymbol{d}_k
$$
其中 \( \boldsymbol{d}_k \) 解线性最小二乘问题 \( \min_{\boldsymbol{d}} \| J(\boldsymbol{x}_k) \boldsymbol{d} + \boldsymbol{r}(\boldsymbol{x}_k) \| \)。

**Levenberg-Marquardt方法**：

Gauss-Newton的改进版，当雅可比矩阵 \( J \) 秩不足或接近奇异时加入正则化：
$$
(J_k^T J_k + \mu I) \boldsymbol{d}_k = -J_k^T \boldsymbol{r}_k
$$

\( \mu \) 根据信赖域策略自适应调整，是处理病态问题的有效方法。

### 6.2 二次规划（QP）

**标准形式**：
$$
\min_{\boldsymbol{x}} \quad \frac{1}{2} \boldsymbol{x}^T Q \boldsymbol{x} + \boldsymbol{c}^T \boldsymbol{x}
$$
$$
\text{s.t.} \quad A\boldsymbol{x} = \boldsymbol{b}, \quad G\boldsymbol{x} \leq \boldsymbol{h}
$$

其中 \( Q \) 为对称矩阵。

**凸二次规划**：\( Q \) 半正定，易于求解。

**非凸二次规划**：\( Q \) 不定，可能有多局部最优（NP难）。

**等式约束QP**（\( Q \) 正定）：
$$
\min \frac{1}{2} \boldsymbol{x}^T Q \boldsymbol{x} + \boldsymbol{c}^T \boldsymbol{x}, \quad \text{s.t.} \quad A\boldsymbol{x} = \boldsymbol{b}
$$

KKT条件：
$$
\begin{bmatrix} Q & A^T \\ A & 0 \end{bmatrix} \begin{bmatrix} \boldsymbol{x} \\ \boldsymbol{\lambda} \end{bmatrix} = \begin{bmatrix} -\boldsymbol{c} \\ \boldsymbol{b} \end{bmatrix}
$$

求解该线性方程组即得最优解。

### 6.3 凸优化在现代AI中的应用

**凸优化的地位**：许多现代AI问题可转化为凸优化问题，包括线性回归、逻辑回归、支持向量机（SVM）、主成分分析（PCA）等。

**支持向量机（SVM）** ：

SVM本质上是凸二次规划问题：
$$
\min_{\boldsymbol{w}, b} \quad \frac{1}{2} \|\boldsymbol{w}\|^2 + C \sum_{i=1}^n \xi_i
$$
$$
\text{s.t.} \quad y_i(\boldsymbol{w}^T \boldsymbol{x}_i + b) \geq 1 - \xi_i, \quad \xi_i \geq 0
$$

其对偶形式是典型的凸二次规划。

**主成分分析（PCA）** ：

PCA的目标是最大化方差，等价于求解协方差矩阵的特征值分解问题。

**随机优化**：机器学习中数据规模巨大，常用随机梯度下降（SGD）及其变体（Adam、RMSprop）处理大规模优化问题。


## 附录：参考书目与备考建议

### 推荐参考书目

1. **《最优化方法》（第三版）** ，孙文瑜、徐成贤、朱德通、孙海琳，高等教育出版社，2025——“十二五”国家级规划教材，内容涵盖基本概念、线性规划、无约束/约束优化、最小二乘、二次规划等

2. **《最优化：建模、算法与理论》** ，刘浩洋等，高教出版社，2025——侧重AI应用，强调向量化思维，含主成分分析、SVM、随机优化等现代内容

3. **《数值最优化方法》** ，高立，北京大学出版社——国内常用教材

4. **《最优化理论与方法》** ，北京交通大学考研指定参考书——适用于备考

5. **《Practical Optimization》** ，陆吾生——华东师大短期课程配套教材

### 课程内容结构（参考多校大纲）

| 章节 | 内容 | 重要度 |
|------|------|--------|
| 第一章 | 基本概念：优化问题、凸集、凸函数、最优性条件 | ★★★ |
| 第二章 | 线性规划：数学模型、单纯形法、对偶理论 | ★★★ |
| 第三章 | 一维搜索：0.618法、Fibonacci法、不精确搜索准则 | ★★ |
| 第四章 | 无约束优化：最速下降法、牛顿法、共轭梯度法、拟牛顿法 | ★★★ |
| 第五章 | 约束优化：KKT条件、罚函数法、SQP | ★★★ |
| 第六章 | 特殊问题：最小二乘、二次规划、凸优化 | ★★ |

### 复习重点提醒

- **重中之重**：凸性判断（凸集+凸函数）、KKT条件推导、BFGS算法框架、SQP子问题构造
- **线性规划必考**：单纯形法表格运算、对偶问题构造、互补松弛性应用
- **无约束优化高频**：最速下降法与牛顿法的收敛性对比、共轭梯度法原理
- **工具技能**：MATLAB/Python的优化工具箱使用（fmincon、linprog、quadprog）

---

> 以上为最优化方法完整知识框架。如需配套例题解析或华东师大历年考题，可进一步告知。
