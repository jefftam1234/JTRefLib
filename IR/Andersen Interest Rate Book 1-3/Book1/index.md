---
tags: [quant/knowledge, quant/books, interest-rates, pricing, modeling]
---

[Master Index](../index.md)

# Book 1: Chapters 1–9 Index

This document provides a highly detailed, mathematically self-sufficient, and rigorously structured index for Chapters 1–9 of Andersen & Piterbarg's *Interest Rate Modeling* (Book 1).

---

## Chapter-by-Chapter Summaries

### Chapter 1: Introduction to Arbitrage Pricing Theory

Chapter 1 establishes the mathematical foundations of continuous-time arbitrage pricing theory under the real-world probability measure $\mathrm{P}$ and its risk-neutral equivalent martingale measure $\mathrm{Q}$. The chapter begins by specifying the stochastic evolution of a vector of primary dividend-free asset prices $X(t)$ driven by a multi-dimensional standard Brownian motion $W(t)$:

This opening orientation identifies how the chapter's assumptions, state variables, and computational objective fit into the broader arbitrage-pricing framework before the first displayed formula.

The Radon-Nikodym likelihood ratio process is the Doléans-Dade exponential that implements the chapter's change of probability measure.

$$
\begin{equation*}
d X(t)=\mu(t, \omega) d t+\sigma(t, \omega) d W(t) \tag{1.2}
\end{equation*}
$$
[[Andersen_Ch01#^eq1-2|Eq. (1.2)]]

where the drift vector $\mu(t, \omega)$ represents instantaneous expected returns and the volatility matrix $\sigma(t, \omega)$ models systematic risk exposure.

To define arbitrage-free pricing without relying on specific utility assumptions, the framework introduces the state price deflator (SDF) or stochastic discount factor $D(t)$:

$$
d D(t) = -r(t) D(t) d t - \theta(t)^{\top} D(t) d W(t)
$$
[[Andersen_Ch01#^eq-custom-1-sdf|Eq. (Custom 1.SDF)]]

which acts as a continuous Arrow-Debreu state-contingent pricing tool, discounting payoffs for the time value of money via the risk-free short rate $r(t)$ and adjusting for risk via the market price of risk vector $\theta(t)$. Under stochastic calculus, the product rule on the deflated asset price process $X_i(t) D(t)$ reveals that the drift is given by:



For this deflated process to be a martingale, its drift must vanish, yielding the fundamental no-arbitrage return relation:



which states that the excess return of each asset above the risk-free rate must equal the covariance of the asset's return with the market price of risk $\theta(t)$ — the continuous-time CAPM pricing relation. The market price of risk $\theta(t)$ is the unique process that simultaneously prices all traded assets under no-arbitrage.

The absence of arbitrage is mathematically equivalent to the existence of a positive deflator $D(t)$ such that the deflated price process $V(t) D(t)$ of any tradeable security is a martingale. Under a chosen numeraire, this martingale representation allows derivative pricing to be formulated as the conditional expectation of deflated terminal payoffs:

$$
\begin{equation*}
V(t)=D(t) \mathrm{E}_{t}^{\mathrm{Q}^{D}}\left(\frac{V(T)}{D(T)}\right) \tag{1.15}
\end{equation*}
$$
[[Andersen_Ch01#^eq1-15|Eq. (1.15)]]

This pricing relationship forms the core foundation of valuation across all interest rate models in the book.

For the Radon-Nikodym density $\varsigma^\theta(t)$ to define a valid equivalent martingale measure rather than merely a local martingale, the market price of risk must satisfy the Novikov integrability condition:

$$
\begin{equation*}
\mathrm{E}^{\mathrm{P}}\left[\exp \left(\frac{1}{2} \int_{0}^{t} \theta(s)^{\top} \theta(s) d s\right)\right]<\infty \tag{1.17}
\end{equation*}
$$
[[Andersen_Ch01#^eq1-17|Eq. (1.17)]]

If this condition fails, $\varsigma^\theta(t)$ is only a strict local martingale and $\mathrm{E}^P[\varsigma^\theta(T)] < 1$, so the candidate measure $\mathrm{Q}$ would not integrate to one — a failure mode that arises in certain models with explosion risk. Novikov's condition is verified by bounding the squared market price of risk, which is always finite under linear growth conditions on $\sigma$.

The transition between probability measures is handled by Girsanov's theorem, which shifts the drift of the Brownian motion by the market price of risk $\theta(t)$ using the Radon-Nikodym density process:

$$
\begin{align*}
\varsigma^{\theta}(t) & =\exp \left(-\int_{0}^{t} \theta(s)^{\top} d W(s)-\frac{1}{2} \int_{0}^{t} \theta(s)^{\top} \theta(s) d s\right) \\
& \triangleq \mathcal{E}\left(-\int_{0}^{t} \theta(s)^{\top} d W(s)\right) \tag{1.16}
\end{align*}
$$
[[Andersen_Ch01#^eq1-16|Eq. (1.16)]]

to define a new Brownian motion $W^{\theta}(t)$ under the equivalent martingale measure $\mathrm{Q}$:

$$
W^{\theta}(t)=W(t)+\int_{0}^{t} \theta(s) d s
$$
[[Andersen_Ch01#^eq-custom-1-girsanov-bm|Eq. (Custom 1.Girsanov BM)]]

where $\theta(s)$ is the market price of risk vector.

The uniqueness of this equivalent martingale measure (EMM) is tied to market completeness. In a complete market, any derivative payoff can be perfectly replicated by a dynamic trading portfolio, meaning the market price of risk vector must satisfy the system of linear equations:

$$
\begin{equation*}
\sigma^{D}(t) \theta(t)=\mu^{D}(t), \quad t \in[0, T], \tag{1.18}
\end{equation*}
$$
[[Andersen_Ch01#^eq1-18|Eq. (1.18)]]

If the asset volatility matrix $\sigma^D(t)$ is invertible, the market price of risk is uniquely determined, and the market is complete.

To connect no-arbitrage pricing to PDE methods, the framework constructs a self-financing replicating portfolio whose incremental value equals the gain from asset positions:

$$
\begin{equation*}
d V(t)=\phi(t)^{\top} d X(t)=\sum_{i=1}^{p} \phi_{i}(t) \mu_{i}(t) d t+\sum_{i=1}^{p} \phi_{i}(t) \sigma_{i}(t) d W(t) \tag{1.25}
\end{equation*}
$$
[[Andersen_Ch01#^eq1-25|Eq. (1.25)]]

For perfect replication, the holdings must equal the delta of the derivative (the hedging relation):

$$
\begin{equation*}
\phi_{i}(t)=\frac{\partial V(t, X(t))}{\partial x_{i}}, \quad i=1, \ldots, p \tag{1.26}
\end{equation*}
$$
[[Andersen_Ch01#^eq1-26|Eq. (1.26)]]

Substituting these delta positions into Itô's lemma and requiring the replicating portfolio to earn exactly the risk-free rate $r(t)$ eliminates all stochastic Brownian terms and yields the pricing PDE. This is the fundamental mechanism by which Itô's calculus converts an SDE problem into a boundary value problem.

For pricing derivatives that are smooth functions of the underlying state processes, the chapter applies Itô's lemma to expand the stochastic differentials of the pricing function $V(t, X(t))$:

$$
\begin{align*}
& d V(t)=V_{t}(t) d t+\sum_{i=1}^{p} V_{x_{i}}(t) \mu_{i}(t) d t \\
& \quad+\frac{1}{2} \sum_{i=1}^{p} \sum_{j=1}^{p} V_{x_{i} x_{j}}(t) \Sigma_{i, j}(t) d t+\sum_{i=1}^{p} V_{x_{i}}(t) \sigma_{i}(t) d W(t) \tag{1.24}
\end{align*}
$$
[[Andersen_Ch01#^eq1-24|Eq. (1.24)]]

where the first-order terms represent directional drift and volatility, and the second-order terms capture the covariance $\Sigma_{i, j}(t) = \sum_{k=1}^d \sigma_{i,k}(t) \sigma_{j,k}(t)$ of the underlying asset paths. Matching the drift of the replicating portfolio to the risk-free rate yields the general Feynman-Kac pricing partial differential equation (PDE):

$$
\begin{equation*}
\frac{\partial u(t, x)}{\partial t}+\mathcal{A} u(t, x)+h(t, x)=r(t, x) u(t, x) \tag{1.34}
\end{equation*}
$$
[[Andersen_Ch01#^eq1-34|Eq. (1.34)]]

This parabolic PDE governs the value of the derivative backward in time. By the Feynman-Kac theorem, the unique classical solution to this PDE admits the stochastic representation:

$$
\begin{equation*}
u(t, x)=\mathrm{E}^{\mathrm{P}}\left(\psi(t, T) g(X(T))+\int_{t}^{T} \psi(t, s) h(s, X(s)) d s \mid X(t)=x\right) \tag{1.35}
\end{equation*}
$$
[[Andersen_Ch01#^eq1-35|Eq. (1.35)]]

where $\psi(t,s) = \exp(-\int_t^s r(u,X(u))\,du)$ is the path-dependent discount factor. This PDE-expectation duality is the architectural foundation of the entire book: it licenses the use of both PDE solvers (Chapters 2, 9) and Monte Carlo simulation (Chapter 3) as interchangeable pricing engines for the same derivative.

### Logical Bridge to Chapter 2

Chapter 2 introduces two-dimensional ADI stencils and predictor-corrector time-stepping schemes to resolve multi-dimensional boundary value problems when analytic solutions do not exist.

---
A supplementary checkpoint from 1.9.1 basics makes another step of the source argument directly citable without changing the existing derivation:

$$
\begin{equation*}
\frac{\partial c}{\partial t}+\frac{1}{2} \sigma^{2} S^{2} \frac{\partial^{2} c}{\partial S^{2}}=0 \tag{1.45}
\end{equation*}
$$
[[Andersen_Ch01#^eq1-45|Eq. (1.45)]]

### Chapter 2: Arbitrage Theory Numerical Solutions

Chapter 2 transitions from the continuous-time theoretical PDEs of Chapter 1 to practical numerical finite difference methods (FDM) for solving parabolic boundary value problems. The general Backward Kolmogorov equation is written as:

This opening orientation identifies how the chapter's assumptions, state variables, and computational objective fit into the broader arbitrage-pricing framework before the first displayed formula.

$$
\begin{equation*}
\frac{\partial V}{\partial t}+\mathcal{L} V=0 \tag{2.1}
\end{equation*}
$$
[[Andersen_Ch02#^eq2-1|Eq. (2.1)]]

where $\mathcal{L}$ is the infinitesimal generator of the underlying state variable, defined for a general process as:



This operator captures the expected rate of change of the option value due to drift and volatility.

Under a one-dimensional state variable $x(t)$ governed by the drift $\mu$ and volatility $\sigma$:

$$
\begin{equation*}
d x(t)=\mu(t, x(t)) d t+\sigma(t, x(t)) d W(t) \tag{2.2}
\end{equation*}
$$
[[Andersen_Ch02#^eq2-2|Eq. (2.2)]]

the spatial differential operator is replaced by finite difference stencils. However, if the grid Peclet condition is violated (meaning the local drift dominates the local volatility over the grid spacing), central differences can introduce spurious oscillations, requiring a switch to first-order upwind differences.

For a Black-Scholes-Merton setting with constant drift $r$ and volatility $\sigma$, the exact terminal value in log-coordinates ($x = \ln S$) follows from Itô's lemma applied to $\ln S$:

$$
\begin{equation*}
x(T)=x(0)+\left(r-\frac{1}{2} \sigma^{2}\right) T+\sigma(W(T)-W(0)) \tag{2.4}
\end{equation*}
$$
[[Andersen_Ch02#^eq2-4|Eq. (2.4)]]

Since $x(T)$ is Gaussian, prices are lognormal and the distribution is fully characterized analytically — this is the motivation for the log-coordinate transform. Converting the PDE into log-coordinates ($x = \ln S$) yields:

$$
\begin{equation*}
\frac{\partial V}{\partial t}+\left(r-\frac{1}{2} \sigma^{2}\right) \frac{\partial V}{\partial x}+\frac{1}{2} \sigma^{2} \frac{\partial^{2} V}{\partial x^{2}}-r V=0 \tag{2.3}
\end{equation*}
$$
[[Andersen_Ch02#^eq2-3|Eq. (2.3)]]

which eliminates the state-dependence of the coefficients. Central difference operators are applied to approximate the first- and second-order spatial derivatives:

$$
\begin{gather*}
\delta_{x} V\left(t, x_{j}\right) \triangleq \frac{V\left(t, x_{j+1}\right)-V\left(t, x_{j-1}\right)}{2 \Delta_{x}}  \tag{2.5}\\
\delta_{x x} V\left(t, x_{j}\right) \triangleq \frac{V\left(t, x_{j+1}\right)+V\left(t, x_{j-1}\right)-2 V\left(t, x_{j}\right)}{\Delta_{x}^{2}} \tag{2.6}
\end{gather*}
$$
[[Andersen_Ch02#^eq2-5|Eq. (2.5)-(2.6)]]

achieving a second-order truncation error of $O(\Delta_x^2)$.

Replacing the spatial derivatives with their grid approximations transforms the PDE into a semi-discretized system of ordinary differential equations (ODEs):

$$
\begin{equation*}
\widehat{\mathcal{L}} \mathbf{V}(t)=\mathbf{A}(t) \mathbf{V}(t)+\mathbf{\Omega}(t) \tag{2.10}
\end{equation*}
$$
[[Andersen_Ch02#^eq2-10|Eq. (2.10)]]

where $\mathbf{V}(t)$ is the vector of grid values and $\mathbf{A}(t)$ is the tridiagonal coefficient matrix. The explicit entries at each interior grid node $j$ are derived by substituting the central difference operators into the spatial operator $\widehat{\mathcal{L}}$:

$$
\begin{align*}
c_{j}(t) & \triangleq-\sigma\left(t, x_{j}\right)^{2} \Delta_{x}^{-2}-r\left(t, x_{j}\right)  \tag{2.7}\\
u_{j}(t) & \triangleq \frac{1}{2} \mu\left(t, x_{j}\right) \Delta_{x}^{-1}+\frac{1}{2} \sigma\left(t, x_{j}\right)^{2} \Delta_{x}^{-2}  \tag{2.8}\\
l_{j}(t) & \triangleq-\frac{1}{2} \mu\left(t, x_{j}\right) \Delta_{x}^{-1}+\frac{1}{2} \sigma\left(t, x_{j}\right)^{2} \Delta_{x}^{-2} \tag{2.9}
\end{align*}
$$
[[Andersen_Ch02#^eq2-7|Eq. (2.7)-(2.9)]]

Here $c_j < 0$ is the on-diagonal entry (combining diffusion and discounting), while $u_j, l_j > 0$ are the super- and sub-diagonal entries. The sign structure $c_j < 0$, $u_j,l_j > 0$, $|c_j| \geq u_j + l_j$ is the M-matrix (monotone matrix) property that guarantees the discrete scheme preserves the maximum principle — option prices cannot exceed their boundary values.

To solve the system in time, a generalized $\theta$-scheme is introduced, which linearizes the backward step using a weighted average of the spatial operator at a parameterized time level:

$$
\begin{equation*}
t_{i}^{i+1}(\theta)=(1-\theta) t_{i+1}+\theta t_{i} \tag{2.14}
\end{equation*}
$$
[[Andersen_Ch02#^eq2-14|Eq. (2.14)]]

Substituting these weighted spatial operators into the backward time step leads to the theta-scheme finite difference system and makes the implicit-explicit balance explicit:

$$
\begin{equation*}
\left(\mathbf{I}-\theta \Delta_{t} \mathbf{A}\left(t_{i}^{i+1}(\theta)\right)\right) \mathbf{V}\left(t_{i}\right)=\left(\mathbf{I}+(I-\theta) \Delta_{t} \mathbf{A}\left(t_{i}^{i+1}(\theta)\right)\right) \mathbf{V}\left(t_{i+1}\right)+e_{i}^{i+1} \tag{2.16}
\end{equation*}
$$
[[Andersen_Ch02#^eq2-16|Eq. (2.16)]]

The local truncation error $e_i^{i+1}$ of the $\theta$-scheme encodes both spatial and temporal discretization errors:

$$
\begin{equation*}
e_{i}^{i+1}=\Delta_{t} O\left(\Delta_{x}^{2}\right)+1_{\left\{\theta \neq \frac{1}{2}\right\}} O\left(\Delta_{t}^{2}\right)+O\left(\Delta_{t}^{3}\right) \tag{2.17}
\end{equation*}
$$
[[Andersen_Ch02#^eq2-17|Eq. (2.17)]]

For $\theta = \frac{1}{2}$ (Crank-Nicolson), the first-order indicator $\mathbf{1}_{\{\theta \neq 1/2\}}$ vanishes, upgrading temporal accuracy from $O(\Delta_t)$ to $O(\Delta_t^2)$. In matrix form, the backward propagation from $t_{i+1}$ to $t_i$ is $\widehat{\mathbf{V}}(t_i) = \mathbf{B}_i^{i+1}\widehat{\mathbf{V}}(t_{i+1})$, where the amplification matrix is:

$$
\begin{equation*}
\widehat{\mathbf{V}}\left(t_{i}\right)=\mathbf{B}_{i}^{i+1} \widehat{\mathbf{V}}\left(t_{i+1}\right), \tag{2.20}
\end{equation*}
$$
[[Andersen_Ch02#^eq2-20|Eq. (2.20)]]

with $\mathbf{B}_i^{i+1} = (\mathbf{I}-\theta\Delta_t\mathbf{A})^{-1}(\mathbf{I}+(1-\theta)\Delta_t\mathbf{A})$. Setting $\theta = 0$ yields the fully explicit scheme (conditionally stable), $\theta = 1.0$ yields the fully implicit scheme (unconditionally stable, tridiagonal inversion required), and $\theta = 0.5$ is Crank-Nicolson achieving $O(\Delta_t^2 + \Delta_x^2)$.

Setting $\theta < 0.5$ restricts the maximum stable time step size by the grid ratio $\alpha = \Delta_t / (\Delta_x)^2$. The Von Neumann stability bound requires:

$$
\begin{equation*}
1 \geq \theta \geq \frac{1}{2}-\frac{1}{\alpha}\left(\frac{\sigma_{k, j}^{2}}{\sigma_{k, j}^{4}+\mu_{k, j}^{2} \Delta_{x}^{2}+\left|\mu_{k, j}^{2} \Delta_{x}^{2}-\sigma_{k, j}^{4}\right|}\right) \tag{2.24}
\end{equation*}
$$
[[Andersen_Ch02#^eq2-24|Eq. (2.24)]]

For Crank-Nicolson ($\theta = 0.5$) and fully implicit ($\theta = 1.0$), the scheme is unconditionally stable in the $L^2$ norm, though Crank-Nicolson can produce spurious oscillations for large time steps when the initial payoff is non-smooth, necessitating Rannacher startup steps.

### Logical Bridge to Chapter 3

Finite difference grids are efficient for low-dimensional options (1D or 2D). For multi-asset contracts or complex path-dependent structures, the curse of dimensionality makes FDM impractical. Chapter 3 introduces Monte Carlo methods, which scale linearly with the number of dimensions.

---
A supplementary checkpoint from 2.8.1 mesh refinements for multiple events makes another step of the source argument directly citable without changing the existing derivation:

$$
\begin{equation*}
\Delta_{t}^{k}=T_{k} / n, \quad \Delta_{x}^{k}=2 \alpha \sigma \sqrt{T_{k}} /(m+1) . \tag{2.52}
\end{equation*}
$$
[[Andersen_Ch02#^eq2-52|Eq. (2.52)]]

### Chapter 3: Monte Carlo Methods

Chapter 3 introduces Monte Carlo simulation as the primary numerical method for pricing high-dimensional and path-dependent interest rate derivatives. The foundation rests on the martingale pricing representation: option value is the expected value of deflated cash flows under the equivalent martingale measure $\mathrm{Q}^N$ associated with a chosen numeraire asset $N(t)$:

This opening orientation identifies how the chapter's assumptions, state variables, and computational objective fit into the broader arbitrage-pricing framework before the first displayed formula.

$$
\begin{equation*}
V(t)=N(t) \mathrm{E}_{t}^{\mathrm{Q}^{N}}(g(T) / N(T)) \tag{3.1}
\end{equation*}
$$
[[Andersen_Ch03#^eq3-1|Eq. (3.1)]]

To evaluate this expectation, $n$ independent paths are simulated and payoffs are averaged. The Monte Carlo price estimator is:

$$
\begin{equation*}
\bar{V}(t)=N(t) \frac{1}{n} \sum_{i=1}^{n} g_{i} / N_{i} . \tag{3.3}
\end{equation*}
$$
[[Andersen_Ch03#^eq3-3|Eq. (3.3)]]

where $g_i$ is the discounted payoff on path $i$ and $N_i$ is the numeraire on that path. By the CLT, the error scales as $O(1/\sqrt{n})$ regardless of dimension, with the two-sided confidence interval:

$$
\begin{equation*}
\left[\bar{V}(t)-u_{\gamma / 2} \cdot s_{n} / \sqrt{n}, \bar{V}(t)+u_{\gamma / 2} \cdot s_{n} / \sqrt{n}\right] \tag{3.4}
\end{equation*}
$$
[[Andersen_Ch03#^eq3-4|Eq. (3.4)]]

where $s_n$ is the sample standard deviation of the discounted payoffs. This dimension-free $O(1/\sqrt{n})$ convergence is the decisive advantage of Monte Carlo over finite difference methods for problems with more than three state variables.

The underlying vector state variables $X(t)$ are modeled using systems of stochastic differential equations (SDEs):

$$
\begin{equation*}
d X(t)=\mu(t, X(t)) d t+\sigma(t, X(t)) d W(t) \tag{3.26}
\end{equation*}
$$
[[Andersen_Ch03#^eq3-26|Eq. (3.26)]]

where $\mu$ is the drift vector and $\sigma$ is the volatility matrix. The simplest scheme is the Euler-Maruyama discretization, which evaluates the drift $\mu$ and volatility $\sigma$ at the beginning of each step:

$$
\begin{align*}
\widehat{X}_{i+1} & =\widehat{X}_{i}+\mu\left(i \Delta, \widehat{X}_{i}\right) \Delta+\sigma\left(i \Delta, \widehat{X}_{i}\right)(W(i \Delta+\Delta)-W(i \Delta))  \tag{3.34}\\
i & =0,1, \ldots, m-1
\end{align*}
$$
[[Andersen_Ch03#^eq3-34|Eq. (3.34)]]

where the Brownian increments are simulated by generating independent standard Gaussian random variables $Z_i \sim \mathcal{N}(0, \mathbf{I})$ scaled by $\sqrt{\Delta}$. For assets following geometric Brownian motion with constant drift and volatility, an exact (zero-discretization-error) simulation is available by integrating the SDE analytically:

$$
\begin{equation*}
S\left(t_{i}\right)=S\left(t_{i-1}\right) \exp \left(\left(r-\frac{1}{2} \sigma^{2}\right) \Delta_{i}\right) \exp \left(\sigma \sqrt{\Delta_{i}} Z_{i}\right), \quad i=1, \ldots, m \tag{3.23}
\end{equation*}
$$
[[Andersen_Ch03#^eq3-23|Eq. (3.23)]]

This step exploits the lognormality of GBM, bypassing discretization bias entirely. When drift or volatility is state-dependent, this exact step is unavailable and one must revert to Euler-Maruyama or its variants.

Euler-Maruyama has a weak convergence order of $1.0$ (errors in expectations scale as $O(\Delta)$) and a strong convergence order of $0.5$ (errors in individual paths scale as $O(\sqrt{\Delta})$). In the presence of highly mean-reverting (stiff) drift coefficients, explicit schemes can become unstable. For example, testing the explicit scheme on a simple test SDE:

$$
\begin{equation*}
\widehat{X}_{i+1}=\widehat{X}_{i}(1+\lambda \Delta)+\sqrt{\Delta} Z_{i} \tag{3.35}
\end{equation*}
$$
[[Andersen_Ch03#^eq3-35|Eq. (3.35)]]

reveals that the solution can oscillate and blow up unless the time step $\Delta$ is extremely small.

To handle mean-reverting interest rate models, the textbook considers the linear-drift SDE:

$$
\begin{equation*}
d X(t)=\varkappa(\theta(t)-X(t)) d t+\sigma(t, X(t)) d W(t) \tag{3.36}
\end{equation*}
$$
[[Andersen_Ch03#^eq3-36|Eq. (3.36)]]

where $\varkappa$ is the speed of mean reversion. The implicit Euler-Maruyama scheme resolves this instability by evaluating the drift term implicitly at the end of the step:

$$
\begin{equation*}
\widehat{X}_{i+1}=\widehat{X}_{i}+\mu\left(i \Delta+\Delta, \widehat{X}_{i+1}\right) \Delta+\sigma\left(i \Delta, \widehat{X}_{i}\right)(W(i \Delta+\Delta)-W(i \Delta)) \tag{3.37}
\end{equation*}
$$
[[Andersen_Ch03#^eq3-37|Eq. (3.37)]]

To ensure the scheme remains consistent with Itô calculus, the volatility term must remain explicit.

To balance the truncation errors of both the drift and volatility terms, the generalized drift-corrected theta-eta schemes are formulated:

$$
\begin{align*}
\widehat{X}_{i+1}=\widehat{X}_{i}+ & {\left[(1-\theta) \bar{\mu}_{\eta}\left(i \Delta, \widehat{X}_{i}\right)+\theta \bar{\mu}_{\eta}\left(i \Delta+\Delta, \widehat{X}_{i+1}\right)\right] \Delta } \\
& +\left[(1-\eta) \sigma\left(i \Delta, \widehat{X}_{i}\right)+\eta \sigma\left(i \Delta+\Delta, \widehat{X}_{i+1}\right)\right] Z_{i} \sqrt{\Delta} \tag{3.39}
\end{align*}
$$
[[Andersen_Ch03#^eq3-39|Eq. (3.39)]]

where a drift correction must be applied if the volatility has state-dependent derivatives, preventing any discretization bias.

For multi-factor models, the Brownian increments across correlated factors must be generated consistently. Given a target correlation matrix $\Sigma$, the Cholesky factorization:

$$
\begin{equation*}
C C^{\top}=\Sigma \tag{3.14}
\end{equation*}
$$
[[Andersen_Ch03#^eq3-14|Eq. (3.14)]]

produces a lower-triangular factor $C$. Multiplying uncorrelated standard normals $X \sim \mathcal{N}(0,\mathbf{I})$ by $C$ yields correlated increments $Z = CX \sim \mathcal{N}(0, \Sigma)$. When the full correlation matrix is not required, a principal component approximation retains only the $r$ largest eigenvalue directions, reducing simulation cost while preserving most of the variance structure.

For more complex models, the predictor-corrector scheme is introduced. This scheme runs an explicit predictor step to estimate the state variable at the end of the interval, then uses this estimate to evaluate the drift implicitly in the corrector step:

$$
\begin{align*}
& \bar{X}_{i+1}=\widehat{X}_{i}+\mu\left(i \Delta, \widehat{X}_{i}\right) \Delta+\sigma\left(i \Delta, \widehat{X}_{i}\right)(W(i \Delta+\Delta)-W(i \Delta))  \tag{3.41}\\
& \widehat{X}_{i+1}=\widehat{X}_{i}+\left[(1-\theta) \bar{\mu}_{\eta}\left(i \Delta, \widehat{X}_{i}\right)+\theta \bar{\mu}_{\eta}\left(i \Delta+\Delta, \bar{X}_{\imath+1}\right)\right] \Delta \\
& \quad+\left[(1-\eta) \sigma\left(i \Delta, \widehat{X}_{i}\right)+\eta \sigma\left(i \Delta+\Delta, \bar{X}_{i+1}\right)\right](W(i \Delta+\Delta)-W(i \Delta)) \tag{3.42}
\end{align*}
$$
[[Andersen_Ch03#^eq3-41|Eq. (3.41)-(3.42)]]

which can be written step-by-step as:



This scheme preserves the stability of implicitness without requiring non-linear solver iterations. The chapter also discusses variance reduction. The antithetic variate estimator pairs each Gaussian draw $Z_i$ with its mirror $-Z_i$, forming the antithetic pair estimator:



which halves variance whenever $\mathrm{Cov}(Y(Z), Y(-Z)) < 0$ — satisfied for monotone payoffs. Low-discrepancy (quasi-random) sequences such as Sobol or Halton sequences fill the unit hypercube more uniformly than pseudo-random numbers, improving effective convergence toward $O((\log n)^d/n)$ for $d$-dimensional integrands at the cost of losing the CLT-based standard error estimate.

### Logical Bridge to Chapter 4

With numerical engines established, we need a realistic framework to model interest rates. Chapter 4 introduces the fundamentals of term structure modeling, transitioning from single-asset dynamics to infinite-dimensional forward rate curves.

---
A supplementary checkpoint from 3.3.1.1 black-scholes delta makes another step of the source argument directly citable without changing the existing derivation:

$$
\begin{equation*}
\frac{d V}{d S_{0}}=\lim _{h \rightarrow 0} \frac{V\left(S_{0}+h\right)-V\left(S_{0}\right)}{h} \tag{3.68}
\end{equation*}
$$
[[Andersen_Ch03#^eq3-68|Eq. (3.68)]]

### Chapter 4: Fundamentals of Interest Rate Modeling

Chapter 4 introduces the core mathematical concepts and frameworks of interest rate modeling. It starts by defining the zero-coupon bond price $P(t, T)$ as the price at time $t$ of a contract that pays 1 (normalized notional) at maturity $T \geq t$. The zero-coupon yield and the instantaneous forward rate $f(t, T)$ are defined, with the latter representing the rate contracted at time $t$ for lending at maturity $T$:

This opening orientation identifies how the chapter's assumptions, state variables, and computational objective fit into the broader arbitrage-pricing framework before the first displayed formula.

$$
\begin{equation*}
f(t, T)=-\frac{\partial \ln P(t, T)}{\partial T} \tag{4.4}
\end{equation*}
$$
[[Andersen_Ch04#^eq4-4|Eq. (4.4)]]

which allows bond prices to be reconstructed from the full forward curve via continuous compounding. The partial bond spanning the forward period $[T, T+\tau]$ is equivalently:

$$
\begin{equation*}
P(t, T, T+\tau)=\exp \left(-\int_{T}^{T+\tau} f(t, u) d u\right) \tag{4.3}
\end{equation*}
$$
[[Andersen_Ch04#^eq4-3|Eq. (4.3)]]

These reconstruction formulas establish that specifying the entire forward curve $\{f(t,u)\}_{u \geq t}$ is equivalent to specifying all discount factors — a key insight for term structure modeling. The instantaneous short rate $r(t)$ is the limit of the forward rate curve as the maturity matches the observation time, representing the rate of return on the money market account:

$$
\begin{equation*}
r(t) \triangleq f(t, t) \tag{4.5}
\end{equation*}
$$
[[Andersen_Ch04#^eq4-5|Eq. (4.5)]]

Under the risk-neutral pricing measure $\mathrm{Q}$, the numeraire is the continuously compounded money market account $\beta(t)$, satisfying the locally deterministic SDE:

$$
\begin{equation*}
d \beta(t)=r(t) \beta(t) d t, \quad \beta(0)=1 \tag{4.12}
\end{equation*}
$$
[[Andersen_Ch04#^eq4-12|Eq. (4.12)]]

Solving this equation yields the stochastic discount factor.

Under the risk-neutral pricing measure $\mathrm{Q}$, the zero-coupon bond price is the expected value of the stochastic discount factor:

$$
\begin{equation*}
P(t, T)=\mathrm{E}_{t}^{\mathrm{Q}}\left(e^{-\int_{t}^{T} r(u) d u}\right) \tag{4.15}
\end{equation*}
$$
[[Andersen_Ch04#^eq4-15|Eq. (4.15)]]

It follows from this lemma that specification of the dynamics of $r(t)$ under $\mathrm{Q}$ suffices to determine the prices of discount bonds at all times and maturities. The Heath-Jarrow-Morton (HJM) framework specifies the dynamics of the instantaneous forward rate curve directly, where zero-coupon bond prices are governed by:

$$
\begin{equation*}
d P(t, T) / P(t, T)=r(t) d t-\sigma_{P}(t, T)^{\top} d W(t) \tag{4.31}
\end{equation*}
$$
[[Andersen_Ch04#^eq4-31|Eq. (4.31)]]

where $\sigma_P(t, T)$ is the integrated volatility vector.

Under the $T$-forward measure $\mathrm{Q}^T$, the zero-coupon bond $P(t,T)$ serves as the numeraire. The Radon-Nikodym density relating $\mathrm{Q}^T$ to the risk-neutral measure $\mathrm{Q}$ is:

$$
\begin{equation*}
\mathrm{E}_{t}^{\mathrm{Q}}\left(\frac{d \mathrm{Q}^{T}}{d \mathrm{Q}}\right)=\frac{P(t, T) / P(0, T)}{\beta(t)} \tag{4.21}
\end{equation*}
$$
[[Andersen_Ch04#^eq4-21|Eq. (4.21)]]

Under $\mathrm{Q}^T$, forward Libor rates $L(t,T,T+\tau) = \mathrm{E}_t^{T+\tau}(L(T,T,T+\tau))$ are martingales — the key property exploited by the LIBOR Market Model. More broadly, the change from risk-neutral measure to any forward measure removes the stochastic discounting from the expectation, which is the mechanism that makes forward-measure pricing tractable. Enforcing this martingale property on zero-coupon bonds restricts the drift of forward rates, leading to the HJM drift condition:

$$
\begin{align*}
d f(t, T) & =\sigma_{f}(t, T)^{\top} \sigma_{P}(t, T) d t+\sigma_{f}(t, T)^{\top} d W(t) \\
& =\sigma_{f}(t, T)^{\top} \int_{t}^{T} \sigma_{f}(t, u) d u d t+\sigma_{f}(t, T)^{\top} d W(t) \tag{4.38}
\end{align*}
$$
[[Andersen_Ch04#^eq4-38|Eq. (4.38)]]

This condition shows that the drift is completely determined by the volatility structure of forward rates, preventing independent specification of drift and volatility. If the forward rate volatility is separable, the short rate can be represented using a finite number of Markovian state variables. For the Gaussian HJM model, the short rate is written as:

$$
\begin{align*}
r(t) & =f(0, t)+h(t) \int_{0}^{t} g(u)^{\top} g(u)\left(\int_{u}^{t} h(s) d s\right) d u+h(t) \int_{0}^{t} g(u)^{\top} d W(u)  \tag{4.45}\\
& \triangleq f(0, t)+h(t) \int_{0}^{t} m_{f}(t, u) d u+h(t) \int_{0}^{t} g(u)^{\top} d W(u)
\end{align*}
$$
[[Andersen_Ch04#^eq4-45|Eq. (4.45)]]

This representation allows pricing via finite difference grids.

Finally, the swap rate is introduced as a weighted average of simple forward rates on a specified tenor structure. Under no-arbitrage, the par swap rate $S_{k,m}(t)$ is represented by:

$$
\begin{equation*}
S_{k, m}(t)=\frac{\sum_{n=k}^{k+m-1} \tau_{n} P\left(t, T_{n+1}\right) L_{n}(t)}{A_{k, m}(t)}, \quad t \leq T_{k} \tag{4.11}
\end{equation*}
$$
[[Andersen_Ch04#^eq4-11|Eq. (4.11)]]

Under the annuity measure associated with the annuity factor, the par swap rate is a martingale, which simplifies options pricing. Due to the daily mark-to-market feature of futures contracts, the futures rate incorporates a convexity adjustment over the forward rate, making futures rates systematically higher than forward rates in upward-sloping yield environments.

### Logical Bridge to Chapter 5

The continuous-time HJM framework defines the dynamics of the yield curve. To apply this in the market, we must analyze the discrete cash-flow instruments traded in the fixed income markets, which is the focus of Chapter 5.

---
A supplementary checkpoint from 4.3.2 risk neutral measures makes another step of the source argument directly citable without changing the existing derivation:

$$
\begin{equation*}
V_{f}(t)=\beta_{f}(t) \mathrm{E}_{t}^{f}\left(g(T) \beta_{f}(T)^{-1}\right) \tag{4.27}
\end{equation*}
$$
[[Andersen_Ch04#^eq4-27|Eq. (4.27)]]

### Chapter 5: Fixed Income Instruments

Chapter 5 covers the definitions, payoffs, and pricing formulas of standard fixed income instruments traded in interest rate markets. The primary benchmark rate is the spot Libor rate $L(T, T, T+\tau)$, which represents the simple interest rate earned on a deposit over the accrual period $[T, T+\tau]$:

This opening orientation identifies how the chapter's assumptions, state variables, and computational objective fit into the broader arbitrage-pricing framework before the first displayed formula.

$$
\begin{equation*}
L=L(T, T, T+\tau)=\frac{1}{\tau}\left(\frac{1}{P(T, T+\tau)}-1\right) . \tag{5.1}
\end{equation*}
$$
[[Andersen_Ch05#^eq5-1|Eq. (5.1)]]

Libor is quoted using simple compounding, meaning that a deposit of 1 at time $T$ returns $1 + \tau L$ at maturity. Under no-arbitrage, this zero-coupon bond price relation must hold:

$$
\begin{equation*}
1+\tau L(t, T, T+\tau)=1 / P(t, T, T+\tau) . \tag{4.2}
\end{equation*}
$$
[[Andersen_Ch04#^eq4-2|Eq. (4.2)]]

which relates simple rates directly to discount factors.

A Forward Rate Agreement (FRA) is a contract that allows a counterparty to lock in a borrowing or lending rate for a future period. The price of the FRA prior to its start date is the difference between the discount factors:

$$
\begin{align*}
V_{\mathrm{FRA}}(t) & =P(t, T)-P(t, T+\tau)-\tau k P(t, T+\tau) \\
& =\tau P(t, T+\tau)\left(\frac{P(t, T)-P(t, T+\tau)}{\tau P(t, T+\tau)}-k\right) \tag{5.2}
\end{align*}
$$
[[Andersen_Ch05#^eq5-2|Eq. (5.2)]]

Under no-arbitrage, the par FRA rate that makes the contract value equal to zero is the forward Libor rate.

Interest rate swaps exchange a series of fixed-rate payments for floating-rate Libor payments. The annuity factor $A(t)$ and par swap rate $S(t)$ over the entire schedule are defined as:

$$
\begin{align*}
A(t) \triangleq A_{0, N}(t) & =\sum_{n=0}^{N-1} \tau_{n} P\left(t, T_{n+1}\right)  \tag{5.4}\\
S(t) \triangleq S_{0, N}(t) & =\frac{\sum_{n=0}^{N-1} \tau_{n} P\left(t, T_{n+1}\right) L_{n}(t)}{\sum_{n=0}^{N-1} \tau_{n} P\left(t, T_{n+1}\right)} \tag{5.5}
\end{align*}
$$
[[Andersen_Ch05#^eq5-4|Eq. (5.4)-(5.5)]]

Prior to maturity, the swap value can be expressed as a linear function of the par swap rate and the annuity factor:

$$
\begin{equation*}
V_{\mathrm{swap}}(t)=A(t)(S(t)-k) \tag{5.6}
\end{equation*}
$$
[[Andersen_Ch05#^eq5-6|Eq. (5.6)]]

where the par swap rate is the rate that makes the swap value equal to zero.

A swaption is an option to enter into an interest rate swap. A payer swaption gives the holder the right to pay fixed, and its payoff at maturity $T_0$ is defined as:

$$
\begin{equation*}
V_{\text {swaption }}\left(T_{0}\right)=\left(V_{\text {swap }}\left(T_{0}\right)\right)^{+}=\left(\sum_{n=0}^{N-1} \tau_{n} P\left(T_{0}, T_{n+1}\right)\left(L_{n}\left(T_{0}\right)-k\right)\right)^{+} \tag{5.10}
\end{equation*}
$$
[[Andersen_Ch05#^eq5-10|Eq. (5.10)]]

Because the swap value can be expressed in terms of the par swap rate, the swaption is a call option on the swap rate. Under the risk-neutral measure $\mathrm{Q}$, the value of the swaption is the expected value of the discounted payoff:

$$
\begin{equation*}
V_{\text {swaption }}(t)=\beta(t) \mathrm{E}_{t}\left(\beta\left(T_{0}\right)^{-1} A\left(T_{0}\right)\left(S\left(T_{0}\right)-k\right)^{+}\right) . \tag{5.11}
\end{equation*}
$$
[[Andersen_Ch05#^eq5-11|Eq. (5.11)]]

Under the swap measure $\mathrm{Q}^A$ (where the annuity factor is the numeraire), the forward swap rate is a martingale, allowing European swaptions to be priced using Black-style models:

$$
\begin{equation*}
V_{\text {swaption }}(t)=A(t) \mathrm{E}_{t}^{A}\left(S\left(T_{0}\right)-k\right)^{+}, \tag{5.12}
\end{equation*}
$$
[[Andersen_Ch05#^eq5-12|Eq. (5.12)]]

This formulation decouples the swaption price from the annuity factor dynamics.

For general forward-starting swaps starting at a future date $T_n$, the forward annuity factor $A_{n,m}(t)$ and forward swap rate $S_{n,m}(t)$ are represented by:

$$
\begin{align*}
A_{n, m}(t) & =\sum_{i=n}^{n+m-1} \tau_{i} P\left(t, T_{i+1}\right)  \tag{5.13}\\
S_{n, m}(t) & =\frac{\sum_{i=n}^{n+m-1} \tau_{i} P\left(t, T_{i+1}\right) L_{i}(t)}{\sum_{i=n}^{n+m-1} \tau_{i} P\left(t, T_{i+1}\right)} \tag{5.14}
\end{align*}
$$
[[Andersen_Ch05#^eq5-13|Eq. (5.13)-(5.14)]]

These forward rates are the underlying variables for swaptions. Finally, average rate options and average rate swaps are covered. Their pricing requires modeling the distribution of the average rate over the schedule, leading to the valuation formula:

$$
\begin{equation*}
V_{\text {average }}(t)=\beta(t) \sum_{n=0}^{N-1} \tau_{n} \mathrm{E}_{t}\left(\beta\left(T_{n+1}\right)^{-1}\left(\bar{L}_{n}-k\right)\right) . \tag{5.9}
\end{equation*}
$$
[[Andersen_Ch05#^eq5-9|Eq. (5.9)]]

which is typically solved using moment-matching approximations or Monte Carlo simulation.

### Logical Bridge to Chapter 6

Pricing market instruments requires a smooth, consistent yield curve calibrated to benchmark prices. Chapter 6 covers curve construction techniques and risk management frameworks for these yield curves.

---

### Chapter 6: Yield Curve Construction and Risk Management

Chapter 6 focuses on yield curve construction and risk management in fixed income markets. Constructing a yield curve involves finding a continuous discount curve $P(t)$ that fits the prices of $N$ benchmark instruments. A general fixed income security with coupon payments $c_{i,j}$ is priced as:

This opening orientation identifies how the chapter's assumptions, state variables, and computational objective fit into the broader arbitrage-pricing framework before the first displayed formula.

$$
\begin{equation*}
V_{\imath}=\sum_{j=1}^{M} c_{\imath, j} P\left(t_{j}\right), \quad i=1, \ldots, N, \tag{6.1}
\end{equation*}
$$
[[Andersen_Ch06#^eq6-1|Eq. (6.1)]]

Because the number of cash flow dates $M$ is typically much larger than the number of benchmark instruments $N$, the system is underdetermined, and infinite curves can satisfy the pricing equation.

Curve construction algorithms resolve this by specifying a parametric form or optimizing a smoothness criterion. In spline interpolation, the curve is built in the yield domain:

$$
\begin{equation*}
e^{-y(T) T}=P(T) \quad \Rightarrow \quad y(T)=-T^{-1} \ln P(T), \tag{6.5}
\end{equation*}
$$
[[Andersen_Ch06#^eq6-5|Eq. (6.5)]]

where $y(T)$ is the zero-coupon yield and $P(T)$ is the discount factor. The equivalent representation in forward rate space reconstructs the discount factor from the entire forward curve:

$$
\begin{equation*}
P(T)=e^{-\int_{0}^{T} f(u) d u} \tag{6.6}
\end{equation*}
$$
[[Andersen_Ch06#^eq6-6|Eq. (6.6)]]

Differentiating $\ln P(T) = -y(T)\cdot T$ with respect to $T$ gives the identity connecting the instantaneous forward rate to the zero-coupon yield and its slope:

$$
\begin{equation*}
f(T)=y(T)+\frac{d y(T)}{d T} T \tag{6.7}
\end{equation*}
$$
[[Andersen_Ch06#^eq6-7|Eq. (6.7)]]

This identity is critical for practitioners: a steep (upward-sloping) yield curve implies forward rates substantially above spot yields, while an inverted curve implies forward rates below the current short rate. A cubic spline segment polynomial is used to represent the yield curve between two adjacent knot points:

$$
\begin{equation*}
y(T)=a_{3, i}\left(T-T_{i}\right)^{3}+a_{2, i}\left(T-T_{i}\right)^{2}+a_{1, i}\left(T-T_{i}\right)+a_{0, i}, \quad T \in\left[T_{i}, T_{i+1}\right] \tag{6.11}
\end{equation*}
$$
[[Andersen_Ch06#^eq6-11|Eq. (6.11)]]

To ensure the curve is smooth, we enforce continuity of the yield and its first and second derivatives at the knots. The simplest degenerate case is piecewise-linear interpolation, which connects adjacent knot values with straight lines:

$$
\begin{equation*}
y(T)=y\left(T_{i}\right) \frac{T_{i+1}-T}{T_{i+1}-T_{i}}+y\left(T_{i+1}\right) \frac{T-T_{i}}{T_{i+1}-T_{i}}, \quad T \in\left[T_{i}, T_{i+1}\right] \tag{6.9}
\end{equation*}
$$
[[Andersen_Ch06#^eq6-9|Eq. (6.9)]]

Linear yield interpolation implies piecewise-flat forward rates with kinks at each knot — distorting short-dated option prices. This motivates the cubic and tension spline methods that produce $C^2$-smooth forward curves.

Tension splines generalize cubic splines by introducing a tension parameter $\sigma$ that pulls the curve toward a linear interpolant, mitigating the oscillatory behavior common in natural cubic splines:

$$
\begin{align*}
y(T)=\left(\frac{\sinh \left(\sigma\left(T_{i+1}-T\right)\right)}{\sinh \left(\sigma h_{i}\right)}-\frac{T_{i+1}-T}{h_{i}}\right) \frac{y_{i}^{\prime \prime}}{\sigma^{2}} & \\
& +\left(\frac{\sinh \left(\sigma\left(T-T_{i}\right)\right)}{\sinh \left(\sigma h_{i}\right)}-\frac{T-T_{i}}{h_{i}}\right) \frac{y_{i+1}^{\prime \prime}}{\sigma^{2}} \\
& +y_{i} \frac{T_{i+1}-T}{h_{i}}+y_{i+1} \frac{T-T_{i}}{h_{i}} \tag{6.16}
\end{align*}
$$
[[Andersen_Ch06#^eq6-16|Eq. (6.16)]]

This tension parametrisation allows control over curve oscillations.

Portfolio risk is decomposed into sensitivities with respect to the prices of the benchmark instruments used to build the curve (par-point deltas). The first-order risk exposure is represented as:

$$
\begin{equation*}
\Delta V_{0} \approx \sum_{i=1}^{N} \frac{\partial V_{0}}{\partial V_{i}} \Delta V_{i}+\sum_{i} \frac{\partial V_{0}}{\partial \theta_{i}} \Delta \theta_{i} \tag{6.25}
\end{equation*}
$$
[[Andersen_Ch06#^eq6-25|Eq. (6.25)]]

To compute these, the benchmark price is bumped, the yield curve is rebuilt, and the portfolio is repriced. A typical benchmark instrument is a Forward Rate Agreement (FRA) whose present value is:

$$
\begin{equation*}
V_{\mathrm{FRA}}=\tau(L(0, T, T+\tau)-k) P(T+\tau) \tag{6.3}
\end{equation*}
$$
[[Andersen_Ch06#^eq6-3|Eq. (6.3)]]

Setting $V_{\mathrm{FRA}}=0$ determines the par FRA rate equal to the forward Libor $L(0,T,T+\tau)$. When the FRA price is bumped by $\delta V_i$, the curve is re-bootstrapped and the portfolio is repriced — the resulting $\partial V_0 / \partial V_i$ is the par-point delta.

Forward rate deltas apply localized perturbation shifts directly to the forward rate curve, represented as Gâteaux derivatives:

$$
\begin{equation*}
\partial_{k} V_{0}=\left.\frac{d V_{0}\left(f(t)+\varepsilon \mu_{k}(t)\right)}{d \varepsilon}\right|_{\varepsilon=0}, \quad k=1, \ldots, K . \tag{6.26}
\end{equation*}
$$
[[Andersen_Ch06#^eq6-26|Eq. (6.26)]]

This allows localized risk sensitivities to be calculated directly. Two standard choices for the bump function $\mu_k(t)$ are used in practice:

$$
\begin{equation*}
+\frac{t_{k+1}-t}{t_{k+1}-t_{k}} 1_{\left\{t \in\left\{t_{k}, t_{k+1}\right)\right\}}, \tag{6.27}
\end{equation*}
$$
[[Andersen_Ch06#^eq6-27|Eq. (6.27)]]

Tent bumps produce forward rate deltas with linear tapering into adjacent buckets (smooth but with cross-bucket leakage), while piecewise-flat bumps produce sharp localized shifts with no contamination across buckets. Piecewise-flat deltas are the standard choice for books with distinct tenor buckets. The Jacobian least-squares hedging problem finds a set of hedging weights that minimizes the weighted difference between the sensitivities of the hedge portfolio and the target portfolio:

$$
\begin{equation*}
\left(\partial \mathbf{H}^{\top} \mathbf{p}-\partial \mathbf{V}_{0}\right)^{\top} \mathbf{W}^{2}\left(\partial \mathbf{H}^{\top} \mathbf{p}-\partial \mathbf{V}_{0}\right)+\mathbf{p}^{\top} \mathbf{U}^{2} \mathbf{p} \rightarrow \min . \tag{6.30}
\end{equation*}
$$
[[Andersen_Ch06#^eq6-30|Eq. (6.30)]]

This minimizes tracking error under variance constraints. The explicit penalized weighted least-squares formulation for the optimal hedge weights $\widehat{\mathbf{p}}$ is:

$$
\begin{equation*}
\widehat{\mathbf{p}}=\underset{\mathbf{p}}{\operatorname{argmin}}\left(\sum_{k=1}^{K} W_{k}^{2}\left(\partial_{k} H_{0}(\mathbf{p})-\partial_{k} V_{0}\right)^{2}+\sum_{l=1}^{L} U_{l}^{2} p_{l}^{2}\right) . \tag{6.29}
\end{equation*}
$$
[[Andersen_Ch06#^eq6-29|Eq. (6.29)]]

The first term penalizes sensitivity mismatch between the hedge portfolio and target; the second term ($U_l^2 p_l^2$) is Tikhonov regularization that prevents excessively large positions. The balance between accuracy and position sizing is governed by the weight matrices $\mathbf{W}$ and $\mathbf{U}$.

Furthermore, the concept of cross-currency interest rate parity is covered. Under the no-arbitrage assumption, the foreign discount factor is linked to the domestic discount factor via the spot and forward exchange rates:

$$
\begin{equation*}
P_{\S}(T)=X(0) P_{\Psi}(T) Y(T) \quad \Rightarrow \quad P_{\Psi}(T)=\frac{P_{\S}(T)}{Y(T) X(0)} . \tag{6.38}
\end{equation*}
$$
[[Andersen_Ch06#^eq6-38|Eq. (6.38)]]

This parity relation must hold to prevent cross-currency arbitrage. In modern markets, the curve used for discounting cash flows (often the OIS curve) differs from the index curve used to project floating cash flows (such as Libor), requiring separate construction.

### Logical Bridge to Chapter 7

A calibrated yield curve allows for the pricing of linear instruments. To price options, we need a volatility model. Chapter 7 covers local volatility models, which model the volatility as a deterministic function of the asset price.

---

### Chapter 7: Vanilla Models with Local Volatility

Chapter 7 introduces local volatility models, also known as Dupire's model. The asset price (or forward rate) process is governed by the SDE:

This opening orientation identifies how the chapter's assumptions, state variables, and computational objective fit into the broader arbitrage-pricing framework before the first displayed formula.

$$
\begin{equation*}
d S(t)=\lambda \varphi(S(t)) d W(t) \tag{7.1}
\end{equation*}
$$
[[Andersen_Ch07#^eq7-1|Eq. (7.1)]]

where the volatility coefficient $\varphi(S(t))$ is a deterministic function of the asset price. By choosing an appropriate state-dependent function, local volatility models can generate a variety of implied volatility smile shapes. For any model in this class, the undiscounted European call price is the conditional expectation of the terminal payoff:

$$
\begin{equation*}
c(t, S(t) ; T, K)=\mathrm{E}_{t}\left((S(T)-K)^{+}\right) \tag{7.3}
\end{equation*}
$$
[[Andersen_Ch07#^eq7-3|Eq. (7.3)]]

The Black-Scholes model (constant $\varphi(S)=S$, lognormal dynamics) gives the closed-form baseline:

$$
\begin{gather*}
c(t, S ; T, K)=S \Phi\left(d_{+}\right)-K \Phi\left(d_{-}\right)  \tag{7.6}\\
d_{ \pm}=\frac{\ln (S / K) \pm \frac{1}{2} \sigma_{\mathrm{B}}(t, S ; T, K)^{2}(T-t)}{\sigma_{\mathrm{B}}(t, S ; T, K) \sqrt{T-t}}
\end{gather*} .
$$
[[Andersen_Ch07#^eq7-6|Eq. (7.6)]]

with $d_{\pm} = (\ln(S/K) \pm \frac{1}{2}\lambda^2(T-t))/(\lambda\sqrt{T-t})$. Local volatility models are calibrated to match the market surface of call prices $c(t,S;T,K)$ across all $(K,T)$, which generically requires $\varphi$ to be a non-trivial function of $S$.

The implied probability density of the terminal asset price is determined by the second derivative of the European call option price surface with respect to the strike $K$:

$$
\begin{align*}
\mathrm{P}_{t}(S(T) \in d K) / d K & =\mathrm{E}_{t}(\delta(S(T)-K))  \tag{7.4}\\
& =\mathrm{E}_{t}\left(\frac{\partial^{2} c(T, S(T) ; T, K)}{\partial K^{2}}\right)=\frac{\partial^{2} c(t, S(t) ; T, K)}{\partial K^{2}} \tag{7.5}
\end{align*}
$$
[[Andersen_Ch07#^eq7-5|Eq. (7.5)]]

This relation, derived using Tanaka's formula, forms the basis of Dupire's local volatility construction. Under this formulation, the local volatility $\sigma_{\text{local}}(K, T)$ is resolved from the option price surface via Dupire's PDE:



which expresses local volatility entirely in terms of market-observable option price sensitivities.

The CEV (Constant Elasticity of Variance) model specializes $\varphi(S) = S^p$ with $p < 1$, leading to the SDE:

$$
\begin{equation*}
d S(t)=\lambda S(t)^{p} d W(t) \tag{7.9}
\end{equation*}
$$
[[Andersen_Ch07#^eq7-9|Eq. (7.9)]]

For $p < 1$, the diffusion coefficient shrinks as $S$ falls, generating negative implied volatility skew consistent with equity and rate markets. The CEV model ($\varphi(S) = S^{\beta}$) generates a skew in the implied volatility smile. The call option price under the CEV model is written as:

$$
\begin{equation*}
c_{\mathrm{CEV}}(t, S ; T, K)=S(1-\Upsilon(a, b+2, c))-K \Upsilon(c, b, a) \tag{7.13}
\end{equation*}
$$
[[Andersen_Ch07#^eq7-13|Eq. (7.13)]]

where $\Upsilon(\cdot;\nu,\gamma)$ is the non-central chi-squared CDF with $\nu$ degrees of freedom and non-centrality $\gamma$. The arguments $(a,b,c)$ are derived from the CEV variance transformation:

$$
\begin{align*}
\vartheta & =\frac{1}{2(p-1)}  \tag{7.10}\\
X(t) & =\frac{S(t)^{2(1-p)}}{(1-p)^{2}} \tag{7.11}
\end{align*}
$$
[[Andersen_Ch07#^eq7-10|Eq. (7.10)-(7.11)]]

with $a = n(T-t)K^{2(1-p)}$, $b = 2|\vartheta|+1$, $c = n(T-t)S(0)^{2(1-p)}$ and $n(t) = \lambda^2/(2\vartheta^2 t)$. The exponent $p$ controls the smile skew: $p=1$ recovers Black-Scholes, $p=0$ gives the normal model, and $p<0$ allows arbitrarily steep skew.

The quadratic volatility model represents another analytically tractable local volatility framework. The local volatility coefficient function is represented by:

$$
\begin{equation*}
\varphi(x)=\alpha+\beta x+\gamma x^{2} \tag{7.22}
\end{equation*}
$$
[[Andersen_Ch07#^eq7-22|Eq. (7.22)]]

This forms a quadratic curve in the asset price.

Under the quadratic volatility model, European option payoffs are split into linear combinations of simpler options to obtain closed-form pricing formulas:

$$
\begin{align*}
p(T)= & \frac{1}{u-l}((K-u)(S(T)-l)-(S(T)-u)(K-l))^{+} \\
= & \frac{(K-u)(S(T)-l)}{u-l} 1_{\{(K-u)(S(T)-l)-(S(T)-u)(K-l)>0\}} \\
& -\frac{(S(T)-u)(K-l)}{u-l} 1_{\{(K-u)(S(T)-l)-(S(T)-u)(K-l)>0\}} \\
\triangleq & p_{1}(T)-p_{2}(T) \tag{7.25}
\end{align*}
$$
[[Andersen_Ch07#^eq7-25|Eq. (7.25)]]

This splitting simplifies the integration of the PDF.

The option pricing function can be transformed to use the time-to-maturity parameter $\tau = T-t$, converting the pricing PDE into a time-homogeneous form:

$$
\begin{equation*}
-\frac{\partial g(\tau, x)}{\partial \tau}+\frac{1}{2} \varphi(x)^{2} \frac{\partial^{2} g(\tau, x)}{\partial x^{2}}=0 \tag{7.35}
\end{equation*}
$$
[[Andersen_Ch07#^eq7-35|Eq. (7.35)]]

This PDE governs the value of the option backward in time. The dual counterpart is the forward Kolmogorov (Fokker-Planck) equation governing the transition density $p(\tau,x)$ of the underlying:



Substituting the undiscounted put price as the Green's function produces Dupire's forward PDE, which propagates option values forward in time-to-maturity $\tau$ for a fixed initial spot. The forward PDE solves the option prices for all strikes and maturities in a single finite difference run, improving calibration efficiency. The initial condition for this PDE is the option payoff:

$$
\begin{equation*}
h(0, x)=(S-x)^{+} \tag{7.40}
\end{equation*}
$$
[[Andersen_Ch07#^eq7-40|Eq. (7.40)]]

which acts as the boundary condition at expiry. The time-dependent skew averaging formula defines the effective time-independent skew parameter for a time-dependent model:

$$
\begin{equation*}
\bar{b}=\int_{0}^{T} b(t) w_{T}(t) d t \tag{7.64}
\end{equation*}
$$
[[Andersen_Ch07#^eq7-64|Eq. (7.64)]]

where $w_T(t)$ is a volatility weight function. However, the model predicts that the local volatility smile moves in the opposite direction of the asset price, which contradicts empirical observations and leads to mispriced options dynamics.

### Logical Bridge to Chapter 8

Local volatility models capture the static smile but predict that the smile moves in the opposite direction of the asset price (the "sticky local volatility" problem). Chapter 8 introduces stochastic volatility models to resolve this.

---
A supplementary checkpoint from 7.5.1 expansion around displaced log-normal process makes another step of the source argument directly citable without changing the existing derivation:

$$
\begin{align*}
(S & \left.+\frac{\beta}{\zeta}\right)^{2} \Omega \frac{\partial \Omega}{\partial \tau} \\
& =\frac{1}{2} \varphi(S)^{2}\left[\left(S+\frac{\beta}{\zeta}\right)^{2} \Omega \frac{\partial^{2} \Omega}{\partial S^{2}}+\left(1-h_{-3}\right)-h_{1}\left(1-h_{-1}\right)\right] \tag{7.42}
\end{align*}
$$
[[Andersen_Ch07#^eq7-42|Eq. (7.42)]]

### Chapter 8: Vanilla Models with Stochastic Volatility I

Chapter 8 covers stochastic volatility (SV) models, which introduce a separate stochastic process for the variance parameter. The general stochastic volatility SDE system combines a state-dependent local volatility function with a mean-reverting variance process:

This opening orientation identifies how the chapter's assumptions, state variables, and computational objective fit into the broader arbitrage-pricing framework before the first displayed formula.

$$
\begin{align*}
d S(t) & =\lambda \varphi(S(t)) \sqrt{z(t)} d W(t)  \tag{8.1}\\
d z(t) & =\theta(m(t)-z(t)) d t+\eta \psi(z(t)) d Z(t), \quad z(0)=z_{0} \tag{8.2}
\end{align*}
$$
[[Andersen_Ch08#^eq8-1|Eq. (8.1)-(8.2)]]

where $W(t)$ and $Z(t)$ are correlated Brownian motions with $d\langle W,Z\rangle_t = \rho\,dt$. The correlation $\rho$ controls the asymmetry of the smile (negative $\rho$ produces negative skew), while $\eta$ controls the curvature (convexity). The specific displaced-lognormal CEV model with CIR variance — the primary pricing model developed in this chapter — is:

$$
\begin{align*}
d S(t) & =\lambda(b S(t)+(1-b) L) \sqrt{z(t)} d W(t)  \tag{8.3}\\
d z(t) & =\theta\left(z_{0}-z(t)\right) d t+\eta \sqrt{z(t)} d Z(t), \quad z(0)=z_{0}=1 \tag{8.4}
\end{align*}
$$
[[Andersen_Ch08#^eq8-3|Eq. (8.3)-(8.4)]]

Here $b\in[0,1]$ interpolates between normal ($b=0$) and lognormal ($b=1$) dynamics; $L$ is the displacement level. Under the Feller condition $2\theta z_0 > \eta^2$, the CIR variance process $z(t)$ remains strictly positive almost surely.

In stochastic volatility models (such as the Heston model), the characteristic function of the log-asset price can be solved analytically by solving a system of Riccati ordinary differential equations:

$$
\begin{equation*}
\Psi_{X}(u ; t) \triangleq \mathrm{E}\left(e^{u \ln X(t)}\right)=\mathrm{E}\left(X(t)^{u}\right) \tag{8.10}
\end{equation*}
$$
[[Andersen_Ch08#^eq8-10|Eq. (8.10)]]

This function completely determines the probability distribution of the log-asset price, allowing European option prices to be computed via Fourier inversion. To evaluate $\Psi_X(u;T)$ analytically, a Girsanov measure change $d\widetilde{\mathrm{P}} \propto X(T)^u\,d\mathrm{P}$ twists the drift of the variance process:

$$
\begin{equation*}
d z(t)=\left(\theta\left(z_{0}-z(t)\right)+\rho \eta \lambda b u z(t)\right) d t+\eta \sqrt{z(t)} d \widetilde{Z}(t), \quad z(0)=z_{0} \tag{8.12}
\end{equation*}
$$
[[Andersen_Ch08#^eq8-12|Eq. (8.12)]]

The twisted process remains a CIR process with modified mean-reversion speed $\theta - \rho\eta\lambda bu$. Because CIR characteristic functions are available in closed form (via Riccati ODEs), this observation makes the entire characteristic function analytically tractable.

To evaluate the joint characteristic function, we must also solve the characteristic function of the integrated variance process:

$$
\begin{equation*}
\Psi_{\bar{z}}(v, u ; t) \triangleq \mathrm{E}^{\widetilde{\mathrm{P}}}\left(e^{v \bar{z}(t)}\right), \quad \bar{z}(t) \triangleq \int_{0}^{t} z(s) d s \tag{8.11}
\end{equation*}
$$
[[Andersen_Ch08#^eq8-11|Eq. (8.11)]]

where $\bar{z}(t)$ is the integrated variance over the path. The CIR conditional moments needed for simulation and calibration are:

$$
\begin{align*}
\mathrm{E}(z(T) \mid z(t)) & =z_{0}+\left(z(t)-z_{0}\right) e^{-\theta(T-t)}  \tag{8.7}\\
\operatorname{Var}(z(T) \mid z(t)) & =\frac{z(t) \eta^{2} e^{-\theta(T-t)}}{\theta}\left(1-e^{-\theta(T-t)}\right)+\frac{z_{0} \eta^{2}}{2 \theta}\left(1-e^{-\theta(T-t)}\right)^{2} \tag{8.8}
\end{align*}
$$
[[Andersen_Ch08#^eq8-7|Eq. (8.7)-(8.8)]]

The mean reverts exponentially to $z_0$ at rate $\theta$, while the variance converges to the stationary level $z_0\eta^2/(2\theta)$. These moments underpin the Quadratic-Exponential (QE) scheme for exact simulation of the CIR process without discretization bias.

The Carr-Madan Fourier integration formula expresses the option price as an integral in the frequency domain, avoiding the need to solve the pricing PDE numerically:

$$
\begin{align*}
c_{\mathrm{SV}}(0, S ; T, K)= & \frac{1}{b} c_{\mathrm{B}}\left(0, S^{\prime} ; T, K^{\prime}, \lambda b\right) \\
& -\frac{K^{\prime}}{2 \pi b} \int_{-\infty}^{\infty} \frac{e^{(1 / 2+i \omega) \ln \left(S^{\prime} / K^{\prime}\right)} q(1 / 2+i \omega)}{\omega^{2}+1 / 4} d \omega, \tag{8.17}
\end{align*}
$$
[[Andersen_Ch08#^eq8-17|Eq. (8.17)]]

where $q$ is the normalized characteristic function and $b$ is a scaling parameter.

The integration is performed numerically using quadrature rules. The discretized Fourier inversion integration is written as:

$$
\begin{align*}
& \operatorname{Re}\left(\int_{0}^{\infty} \frac{e^{(1 / 2+i \omega) \ln \left(S^{\prime} / K^{\prime}\right)} q(1 / 2+i \omega)}{\omega^{2}+1 / 4} d \omega\right) \\
& \quad \approx \frac{\omega_{\max }}{N_{\omega}} \sum_{n=0}^{N_{\omega}-1} \frac{e^{\ln \left(S^{\prime} / K^{\prime}\right) / 2}}{\omega_{n}^{2}+1 / 4} \operatorname{Re}\left(e^{i \omega_{\prime \prime} \ln \left(S^{\prime} / K^{\prime}\right)} q\left(1 / 2+i \omega_{n}\right)\right) \tag{8.25}
\end{align*}
$$
[[Andersen_Ch08#^eq8-25|Eq. (8.25)]]

where $N_\omega$ is the number of grid points. To improve numerical accuracy for long-term options, an adaptive oscillatory quadrature algorithm is applied to the tail integration, yielding the approximation for the second component:

$$
\begin{align*}
& I_{2} \approx e^{q_{0}} \sqrt{\frac{S^{\prime}}{K^{\prime}}} \sum_{n=0}^{N_{\omega}^{\prime \prime}-1} \frac{e^{r\left(\omega_{n}^{\prime \prime}\right)-r_{\infty}\left(\omega_{n}^{\prime \prime}\right)}-1}{\left(\omega_{n}^{\prime \prime}\right)^{2}+1 / 4} \int_{\omega_{n}^{\prime \prime}}^{\omega_{n+1}^{\prime \prime}} e^{\omega\left(-q_{\infty}+i \ln \left(S^{\prime} / K^{\prime}\right)\right)} d \omega \\
&= e^{q_{0}} \sqrt{\frac{S^{\prime}}{K^{\prime}}} \sum_{n=0}^{N_{\omega}^{\prime \prime}-1} \frac{e^{r\left(\omega_{n}^{\prime \prime}\right)-r_{\infty}\left(\omega_{n}^{\prime \prime}\right)}-1}{\left(\omega_{n}^{\prime \prime}\right)^{2}+1 / 4} \\
& \times \frac{e^{\omega_{n+1}^{\prime \prime}}\left(-q_{\infty}+i \ln \left(S^{\prime} / K^{\prime}\right)\right)}{-e^{\omega_{n}^{\prime \prime}}\left(-q_{\infty}+i \ln \left(S^{\prime} / K^{\prime}\right)\right)}  \tag{8.33}\\
&-q_{\infty}+i \ln \left(S^{\prime} / K^{\prime}\right)
\end{align*}
$$
[[Andersen_Ch08#^eq8-33|Eq. (8.33)]]

This algorithm concentrates nodes in areas of high oscillation.

The payoff replication formula replicates any twice-continuously differentiable payoff using a portfolio of cash, the underlying asset, and a strip of European put and call options:

$$
\begin{align*}
& \mathrm{E}(f(S(T)))=f\left(K^{*}\right)+f^{\prime}\left(K^{*}\right)\left(S(0)-K^{*}\right) \\
& \quad+\int_{-\infty}^{K^{*}} p(0, S(0) ; T, K) f^{\prime \prime}(K) d K+\int_{K^{*}}^{\infty} c(0, S(0) ; T, K) f^{\prime \prime}(K) d K \tag{8.35}
\end{align*}
$$
[[Andersen_Ch08#^eq8-35|Eq. (8.35)]]

This replication is model-independent. Finally, transformed variance processes are introduced to stabilize the simulation. By applying a coordinate transformation to the variance process, numerical errors near zero are reduced, leading to the SDE:

$$
\begin{equation*}
d z^{*}(t)=\theta\left(z_{0}^{*}-z^{*}(t)\right) d t+\eta \sqrt{z^{*}(t)} d Z(t), \quad z_{0}^{*}=z_{0}-v>0 \tag{8.44}
\end{equation*}
$$
[[Andersen_Ch08#^eq8-44|Eq. (8.44)]]

where $z^*(t)$ represents the transformed variance.

### Logical Bridge to Chapter 9

Fourier pricing methods are powerful but limited to models with analytical characteristic functions. Chapter 9 covers calibration techniques for time-dependent parameters and numerical Monte Carlo simulation methods for stochastic volatility models.

---
A supplementary checkpoint from 8.4.4 refinements of numerical implementation makes another step of the source argument directly citable without changing the existing derivation:

$$
\begin{equation*}
\int_{0}^{\infty} \frac{e^{(1 / 2+i \omega) \ln \left(S^{\prime} / K^{\prime}\right) q(1 / 2+i \omega)}}{\omega^{2}+1 / 4} d \omega=I_{1}+I_{2}+I_{3} \tag{8.29}
\end{equation*}
$$
[[Andersen_Ch08#^eq8-29|Eq. (8.29)]]

### Chapter 9: Vanilla Models with Stochastic Volatility II

Chapter 9 extends the stochastic volatility models to include time-dependent parameters. The SDE system incorporates time-dependency in the volatility scaling $\lambda(t)$ to match the term structure of implied volatilities observed in the market:

This opening orientation identifies how the chapter's assumptions, state variables, and computational objective fit into the broader arbitrage-pricing framework before the first displayed formula.

$$
\begin{align*}
d S(t) & =\lambda(t)(b S(t)+(1-b) L) \sqrt{z(t)} d W(t),  \tag{9.1}\\
d z(t) & =\theta\left(z_{0}-z(t)\right) d t+\eta \sqrt{z(t)} d Z(t), \tag{9.2}
\end{align*}
$$
[[Andersen_Ch09#^eq9-1|Eq. (9.1)-(9.2)]]

where $L$ is the CEV parameter.

To price these models efficiently, parameter averaging techniques are introduced. The effective volatility parameter is recovered via an expansion in powers of the vol-of-vol $\eta$:

$$
\begin{equation*}
v^{*}=\bar{v}+\eta^{2} v_{1}^{*}+\eta^{4} v_{2}^{*}+\ldots \tag{9.20}
\end{equation*}
$$
[[Andersen_Ch09#^eq9-20|Eq. (9.20)]]

where $\bar{v}$ is the leading-order average volatility.

The volatility averaging weight function defines the weighting scheme used to average time-dependent parameters, mapping the time-dependent model to an effective time-independent model:

$$
\begin{align*}
w_{T}(t) & =\frac{\widehat{v}(t)^{2} \lambda(t)^{2}}{\int_{0}^{T} \widehat{v}(t)^{2} \lambda(t)^{2} d t}  \tag{9.36}\\
\widehat{v}(t)^{2} & =z_{0}^{2} \int_{0}^{t} \lambda(s)^{2} d s+z_{0} e^{-\theta t} \int_{0}^{t} \lambda(s)^{2} e^{-\theta s} \int_{0}^{s} \eta(u)^{2} e^{2 \theta u} d u d s
\end{align*}
$$
[[Andersen_Ch09#^eq9-36|Eq. (9.36)]]

where $\widehat{v}(t)$ represents the local volatility weight.

For numerical pricing via finite difference methods, the two-dimensional pricing PDE is solved over a grid. The general time-dependent SV dynamics are written as:

$$
\begin{align*}
d S(t) & =\varphi(t, S(t)) \sqrt{z(t)} d W(t)  \tag{9.54}\\
d z(t) & =\theta\left(z_{0}-z(t)\right) d t+\eta(t) \psi(z(t)) d Z(t) \tag{9.55}
\end{align*}
$$
[[Andersen_Ch09#^eq9-54|Eq. (9.54)-(9.55)]]

where $\psi(z(t))$ represents the variance volatility function. A quadratic grid spacing is introduced for the discretized variance domain, concentrating the grid density near zero to improve numerical resolution in the most statistically relevant region:

$$
\begin{equation*}
\zeta_{n}=\left(\sqrt{z_{\min }}+\frac{n}{N_{z}}\left(\sqrt{z_{\max }}-\sqrt{z_{\min }}\right)\right)^{2}, \quad n=0, \ldots, N_{z} . \tag{9.60}
\end{equation*}
$$
[[Andersen_Ch09#^eq9-60|Eq. (9.60)]]

where $N_z$ is the number of grid points in the variance dimension.

For numerical simulation via Monte Carlo, the standard stochastic volatility model is discretized:

$$
\begin{align*}
d S(t) & =\lambda(b S(t)+(1-b) L) \sqrt{z(t)} d W(t)  \tag{9.65}\\
d z(t) & =\theta\left(z_{0}-z(t)\right) d t+\eta \sqrt{z(t)} d Z(t) \tag{9.66}
\end{align*}
$$
[[Andersen_Ch09#^eq9-65|Eq. (9.65)-(9.66)]]

where $W(t)$ and $Z(t)$ are correlated Brownian motions.

The QE scheme matches the conditional mean and variance of the true non-central chi-square distribution. By representing the simulated variance as a squared Gaussian variable in the high-state regime, the variance update is written as:

$$
\begin{equation*}
\widehat{z}(t+\Delta)=a(b+Z)^{2} \tag{9.75}
\end{equation*}
$$
[[Andersen_Ch09#^eq9-75|Eq. (9.75)]]

where $Z$ is a standard normal variable. Furthermore, the Broadie-Kaya framework is covered for bias-free simulation of stochastic volatility models. The log-asset price SDE in this framework is written as:

$$
\begin{equation*}
d \ln X(t)=-\frac{1}{2} \lambda^{2} b^{2} z(t) d t+\lambda b \sqrt{z(t)} d W(t) \tag{9.84}
\end{equation*}
$$
[[Andersen_Ch09#^eq9-84|Eq. (9.84)]]

which is integrated step-by-step using Fourier inversion.

### Logical Bridge to Chapter 10 (Book 2)

Book 1 focuses on single-currency models. Book 2 (Chapter 10 onwards) extends these concepts to multi-factor short-rate models and Libor market models, which are required to price complex interest rate exotics.

---

A supplementary checkpoint from 9.1 fourier integration with time-dependent parameters makes another step of the source argument directly citable without changing the existing derivation:

$$
\begin{equation*}
\Psi \overline{z \lambda^{2}}(v, u ; t) \triangleq \mathrm{E}^{\widetilde{\mathrm{P}}}\left(e^{v \overline{z \lambda^{2}}(t)}\right), \quad \overline{z \lambda^{2}}(t) \triangleq \int_{0}^{t} z(s) \lambda(s)^{2} d s \tag{9.3}
\end{equation*}
$$
[[Andersen_Ch09#^eq9-3|Eq. (9.3)]]

A supplementary checkpoint from 9.3.1 volatility averaging makes another step of the source argument directly citable without changing the existing derivation:

$$
\begin{equation*}
h(x) \approx a+b e^{c x} \tag{9.25}
\end{equation*}
$$
[[Andersen_Ch09#^eq9-25|Eq. (9.25)]]

A supplementary checkpoint from 9.3.4 calibration by parameter averaging makes another step of the source argument directly citable without changing the existing derivation:

$$
\begin{equation*}
\bar{\lambda}_{n}(\mathcal{X})=F\left(\lambda(\cdot) ; \bar{b}_{n}(\mathcal{X}), \bar{\eta}_{n}(\mathcal{X})\right), \quad n=1, \ldots, N, \tag{9.47}
\end{equation*}
$$
[[Andersen_Ch09#^eq9-47|Eq. (9.47)]]

A supplementary checkpoint from 9.5.3.1 log-normal approximation makes another step of the source argument directly citable without changing the existing derivation:

$$
\begin{align*}
e^{\mu+\frac{1}{2} \sigma^{2}} & =\mathrm{E}(z(t+\Delta) \mid z(t)=\widehat{z}(t))  \tag{9.72}\\
e^{2\left(\mu+\frac{1}{2} \sigma^{2}\right)}\left(e^{\sigma^{2}}-1\right) & =\operatorname{Var}(z(t+\Delta) \mid z(t)=\widehat{z}(t)) \tag{9.73}
\end{align*}
$$
[[Andersen_Ch09#^eq9-72|Eq. (9.72)]]

A supplementary checkpoint from 9.a appendix: proof of proposition 9.3.4 makes another step of the source argument directly citable without changing the existing derivation:

$$
\begin{align*}
v(t)^{2}= & \mathrm{E}\left(\left(S(t)-x_{0}\right)^{2} z(t)\right)  \tag{9.101}\\
= & \int_{0}^{t} \lambda(s)^{2}\left(e^{-\theta(t-s)} \mathrm{E}\left(z(s)^{2}\right)+\left(1-e^{-\theta(t-s)}\right) z_{0}^{2}\right) d s \\
= & \int_{0}^{t} \lambda(s)^{2}\left(e^{-\theta(t-s)} z_{0}^{2}+e^{-\theta(t-s)} z_{0} \eta^{2} \frac{1-e^{-2 \theta s}}{2 \theta}\right.  \tag{9.102}\\
& \left.+z_{0}^{2}\left(1-e^{-\theta(t-s)}\right)\right) d s \\
= & \int_{0}^{t} \lambda(s)^{2}\left(z_{0}^{2}+z_{0} \eta^{2} e^{-\theta(t-s)} \frac{1-e^{-2 \theta s}}{2 \theta}\right) d s \\
= & z_{0}^{2} \int_{0}^{t} \lambda(s)^{2} d s+z_{0} \eta^{2} e^{-\theta t} \int_{0}^{t} \lambda(s)^{2} \frac{e^{\theta s}-e^{-\theta s}}{2 \theta} d s
\end{align*}
$$
[[Andersen_Ch09#^eq9-101|Eq. (9.101)]]

## Mathematical Deep Dives

### Deep Dive 1: Heath-Jarrow-Morton (HJM) Drift Condition Derivation

The HJM framework models the continuous evolution of the entire instantaneous forward rate curve $f(t, T)$ directly. We begin by specifying the SDE for $f(t, T)$ under the risk-neutral equivalent martingale measure $\mathrm{Q}$ as:

where $W(t)$ is a $d$-dimensional standard Brownian motion, $\sigma_f(t, T)$ is a $d$-dimensional volatility vector process, and $\alpha_f(t, T)$ is the drift rate. Under no-arbitrage, the price of a zero-coupon bond $P(t, T)$ must satisfy the relation:

We define the integrated forward rate process $I(t, T) = \int_t^T f(t, u) du$, so that $P(t, T) = e^{-I(t, T)}$. Applying Itô's lemma to find the differential of $I(t, T)$ with respect to $t$:

For convenience, let us define the integrated volatility vector $\sigma_P(t, T) \triangleq \int_t^T \sigma_f(t, u) du$ and the integrated drift $\alpha_P(t, T) \triangleq \int_t^T \alpha_f(t, u) du$. The differential of $I(t, T)$ simplifies to:

Now, we apply Itô's lemma to the function $P(t, T) = e^{-I(t, T)}$. Let $g(y) = e^{-y}$. The derivatives are $g'(y) = -e^{-y}$ and $g''(y) = e^{-y}$. The stochastic differential of $P(t, T)$ is:

Under the risk-neutral pricing measure $\mathrm{Q}$ with the money market account $\beta(t) = \exp(\int_0^t r(u) du)$ as the numeraire, the deflated bond price process $P(t, T) / \beta(t)$ must be a $\mathrm{Q}$-martingale. The differential of the deflated bond price is:

For this process to be a martingale, its drift term must be zero almost surely. This yields the integrated HJM drift condition:
$$
\begin{align*}
d f(t, T) & =\sigma_{f}(t, T)^{\top} \sigma_{P}(t, T) d t+\sigma_{f}(t, T)^{\top} d W(t) \\
& =\sigma_{f}(t, T)^{\top} \int_{t}^{T} \sigma_{f}(t, u) d u d t+\sigma_{f}(t, T)^{\top} d W(t) \tag{4.38}
\end{align*}
$$
[[Andersen_Ch04#^eq4-38|Eq. (4.38)]]
Differentiating this relation with respect to maturity $T$:

This is the HJM drift condition. It shows that in an arbitrage-free market, the drift of the instantaneous forward rate is completely determined by its volatility structure.

---

### Deep Dive 2: The Quadratic-Exponential (QE) Simulation Scheme

The CIR/Heston stochastic variance process satisfies the SDE:
$$
\begin{align*}
d S(t) & =\lambda(b S(t)+(1-b) L) \sqrt{z(t)} d W(t)  \tag{9.65}\\
d z(t) & =\theta\left(z_{0}-z(t)\right) d t+\eta \sqrt{z(t)} d Z(t) \tag{9.66}
\end{align*}
$$
[[Andersen_Ch09#^eq9-65|Eq. (9.65)-(9.66)]]
Due to the square root, standard Euler-Maruyama discretization can yield negative variance values. The QE scheme (Andersen [2008]) resolves this by matching the conditional moments of the true non-central chi-square distribution using two different regimes.
Let $m \triangleq \mathrm{E}(z(t+\Delta) \mid z(t))$ and $s^2 \triangleq \operatorname{Var}(z(t+\Delta) \mid z(t))$. Under the CIR process, these moments are:

We define the switching variable $\psi \triangleq s^2 / m^2$.

#### Regime 1: High-State Scheme ($\psi \leq \psi_c$)
When the variance is high (typically $\psi_c = 1.5$), the distribution of $z(t+\Delta)$ is approximated by a squared Gaussian variable:

We choose the constants $a$ and $b$ to match the mean $m$ and variance $s^2$. The mean and variance of $\widehat{z}(t+\Delta)$ are:

Setting these equal to $m$ and $s^2$ yields the system:

Dividing the second equation by the square of the first:

Solving this quadratic equation for $b^2$ yields:
$$
\begin{equation*}
b^{2}=2 \psi^{-1}-1+\sqrt{2 \psi^{-1}} \sqrt{2 \psi^{-1}-1} \geq 0 \tag{9.80}
\end{equation*}
$$
[[Andersen_Ch09#^eq9-80|Eq. (9.80)]]
Once $b^2$ is computed, $a$ is recovered by:
$$
\begin{equation*}
a=\frac{m}{1+b^{2}} \tag{9.81}
\end{equation*}
$$
[[Andersen_Ch09#^eq9-81|Eq. (9.81)]]

#### Regime 2: Low-State Scheme ($\psi > \psi_c$)
When the variance is low, the distribution of $z(t+\Delta)$ is approximated by a combination of a Dirac delta mass at zero and an exponential tail:
$$
\begin{equation*}
\mathrm{P}(\widehat{z}(t+\Delta) \in[x, x+d x])=\left(p \delta(x)+\beta(1-p) e^{-\beta x}\right) d x, \quad x \geq 0 \tag{9.76}
\end{equation*}
$$
[[Andersen_Ch09#^eq9-76|Eq. (9.76)]]
The cumulative distribution function is:
$$
\begin{equation*}
\Psi(x)=\mathrm{P}(\widehat{z}(t+\Delta) \leq x)=p+(1-p)\left(1-e^{-\beta x}\right), \quad x \geq 0 \tag{9.77}
\end{equation*}
$$
[[Andersen_Ch09#^eq9-77|Eq. (9.77)]]
The mean and variance of this distribution are:

Setting these equal to $m$ and $s^2$:

Dividing the second equation by the square of the first:

Solving for $p$ yields:
$$
\begin{equation*}
p=\frac{\psi-1}{\psi+1} \in[0,1) \tag{9.82}
\end{equation*}
$$
[[Andersen_Ch09#^eq9-82|Eq. (9.82)]]
Once $p$ is computed, $\beta$ is recovered by:
$$
\begin{equation*}
\beta=\frac{1-p}{m}=\frac{2}{m(\psi+1)}>0 . \tag{9.83}
\end{equation*}
$$
[[Andersen_Ch09#^eq9-83|Eq. (9.83)]]
To simulate $\widehat{z}(t+\Delta)$, we draw a uniform random variable $U_z \sim \mathcal{U}[0, 1]$ and invert the CDF:
$$
\begin{equation*}
\widehat{z}(t+\Delta)=\Psi^{-1}\left(U_{z} ; p, \beta\right) \tag{9.79}
\end{equation*}
$$
[[Andersen_Ch09#^eq9-79|Eq. (9.79)]]
This moment-matched scheme ensures strictly non-negative variance paths and exhibits attractive convergence properties.
