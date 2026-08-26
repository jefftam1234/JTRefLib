---
tags: [quant/knowledge, quant/books, interest-rates, pricing, modeling]
---

[Master Index](../index.md)

# Book 2: Chapters 10-15 Index

This index is being executed in chapter batches from the approved roadmap.

- Batch 1 complete: Chapter 10
- Batch 2 complete: Chapter 11
- Batch 3 complete: Chapter 12
- Batch 4 complete: Chapter 13
- Batch 5 complete: Chapter 14
- Batch 6 complete: Chapter 15

---

## Chapter-by-Chapter Summaries

### Chapter 10: One-Factor Short Rate Models I

Chapter 10 starts from the risk-neutral discount-bond valuation identity and then narrows the model class to one-factor short-rate dynamics that are analytically or numerically tractable. The core structural move is to treat the short rate as the unique state variable from which the full discount curve can be reconstructed, so pricing reduces to understanding the state dynamics and the conditional expectation operator:

The chapter's modeling logic is recursive: start from a generic pricing identity, impose structure on the state process, then recover tractable bond and option maps that can be calibrated repeatedly. This is why each new equation family is evaluated not only for economic plausibility but also for algorithmic reuse in calibration and risk engines.

$$
\begin{equation*}
P(t, T)=\mathrm{E}_{t}^{\mathrm{Q}}\left(e^{-\int_{t}^{T} r(u) d u}\right) \tag{10.1}
\end{equation*}
$$
[[Andersen_Ch10#^eq10-1|Eq. (10.1)]]

The first concrete specification is the Vasicek mean-reverting Gaussian short-rate model. Its importance in the chapter is not realism per se, but that it provides a transparent benchmark where state moments are explicit and all later approximations/calibration heuristics can be compared against a closed-form anchor:

$$
\begin{equation*}
d r(t)=\varkappa(\vartheta-r(t)) d t+\sigma_{r} d W(t) \tag{10.10}
\end{equation*}
$$
[[Andersen_Ch10#^eq10-10|Eq. (10.10)]]

The closed-form state solution makes explicit how mean reversion and diffusion accumulate over time.

Operationally, this explicit solution is the bridge from qualitative mean reversion intuition to quantitative variance-term calculations and parameter sensitivity analysis. It also clarifies which parameters can be estimated from historical time-series behavior and which are primarily calibration knobs.

$$
\begin{equation*}
r(t)=\vartheta+(r(0)-\vartheta) e^{-\varkappa t}+\sigma_{r} \int_{0}^{t} e^{-\varkappa(t-s)} d W(s) \tag{10.11}
\end{equation*}
$$
[[Andersen_Ch10#^eq10-11|Eq. (10.11)]]

The model is then generalized to time-dependent parameters to fit the initial term structure while preserving one-factor Markov tractability. The key technical bridge is the shifted state representation and bond reconstitution formula, which turns the entire curve into an explicit function of one state and deterministic kernels:

$$
\begin{equation*}
d r(t)=\varkappa(t)(\vartheta(t)-r(t)) d t+\sigma_{r}(t) d W(t) \tag{10.14}
\end{equation*}
$$
[[Andersen_Ch10#^eq10-14|Eq. (10.14)]]

Using the shifted coordinate $x(t)=r(t)-f(0,t)$ then yields a direct curve reconstitution map.

$$
\begin{align*}
& P(t, T)=\frac{P(0, T)}{P(0, t)} \exp \left(-x(t) G(t, T)-\frac{1}{2} y(t) G(t, T)^{2}\right)  \tag{10.18}\\
& G(t, T)=\int_{t}^{T} e^{-\int_{t}^{u} \varkappa(s) d s} d u
\end{align*}
$$
[[Andersen_Ch10#^eq10-18|Eq. (10.18)]]

With this state/bond representation in hand, the chapter moves to swaption valuation. It first formulates the terminal payoff directly as a nonlinear function of discount bonds, then introduces the Jamshidian critical-rate equation that decomposes the payoff into bond-option components. This is the exact route and remains an important numerical benchmark even when later approximations are used in production:

$$
\begin{equation*}
V_{\text {swaption }}\left(T_{0}\right)=\left(1-P\left(T_{0}, T_{N}\right)-c \sum_{i=0}^{N-1} \tau_{i} P\left(T_{0}, T_{i+1}\right)\right)^{+}, \quad \tau_{i}=T_{i+1}-T_{i} . \tag{10.21}
\end{equation*}
$$
[[Andersen_Ch10#^eq10-21|Eq. (10.21)]]

Jamshidian decomposition introduces a critical state $x^*$ that nulls the underlying swap value at expiry.

$$
\begin{equation*}
P\left(T_{0}, T_{N}, x^{*}\right)+c \sum_{i=0}^{N-1} \tau_{i} P\left(T_{0}, T_{\imath+1}, x^{*}\right)=1 \tag{10.22}
\end{equation*}
$$
[[Andersen_Ch10#^eq10-22|Eq. (10.22)]]

The chapter then switches to annuity-measure valuation to derive a calibration-friendly route based on approximate swap-rate dynamics. This representation is central because it exposes measure choice and numeraire structure directly, and it scales better than repeated Jamshidian decompositions in large calibration loops:

$$
\begin{equation*}
V_{\text {swaption }}(0)=A(0) \mathrm{E}^{A}\left(\left(S\left(T_{0}\right)-c\right)^{+}\right) \tag{10.25}
\end{equation*}
$$
[[Andersen_Ch10#^eq10-25|Eq. (10.25)]]

The practical consequence is a clear separation of concerns: payoff nonlinearity remains in the terminal mapping, while dynamics and measure effects are absorbed into the swap-rate process under a convenient numeraire. This decomposition is what makes large calibration batches numerically manageable.

For numerical pricing of general payoffs, the one-factor PDE is derived in the shifted state coordinate. This PDE is the practical hinge from analytical model definitions to finite-difference implementation and is repeatedly reused later when closed-form transform methods are not available:

$$
\begin{equation*}
\frac{\partial V}{\partial t}+(y(t)-\varkappa(t) x) \frac{\partial V}{\partial x}+\frac{1}{2} \sigma_{r}(t)^{2} \frac{\partial^{2} V}{\partial x^{2}}=(x+f(0, t)) V \tag{10.29}
\end{equation*}
$$
[[Andersen_Ch10#^eq10-29|Eq. (10.29)]]

The affine one-factor extension then broadens volatility structure from constant/deterministic forms to state-dependent forms while retaining transform-based tractability. This is where Chapter 10 moves from the Gaussian benchmark into skew-capable dynamics and sets up the Riccati machinery that will recur across later affine and quasi-Gaussian constructions:

$$
\begin{equation*}
d r(t)=\varkappa(\vartheta-r(t)) d t+\sigma v(r(t)) d W(t) \tag{10.44}
\end{equation*}
$$
[[Andersen_Ch10#^eq10-44|Eq. (10.44)]]

Transform tractability is recovered by solving coupled Riccati equations for transform coefficients.

These coupled Riccati equations are the true computational core of the affine block: once solved robustly, they compress path-dependent pricing objects into low-dimensional deterministic coefficient functions that can be reused across instruments and maturities.

$$
\begin{array}{r}
\frac{d A}{d t}-\varkappa(t) \vartheta(t) B+\frac{1}{2} \sigma(t)^{2} \alpha B^{2}=0 \\
-\frac{d B}{d t}+\varkappa(t) B+\frac{1}{2} \sigma(t)^{2} \beta B^{2}=c_{2} \tag{10.54}
\end{array}
$$
[[Andersen_Ch10#^eq10-54|Eq. (10.54)]]

In curve-fitting affine form, discount-bond pricing is rewritten through the extended transform in shifted coordinates. This gives the operational calibration interface of the chapter: solve Riccati/Volterra components once, then evaluate curve-consistent bond and option quantities repeatedly in calibration and simulation loops:

$$
\begin{equation*}
P(t, T)=g(t, T ; 0,1)=\frac{P(0, T)}{P(0, t)} e^{c(t, T)-x(t) b(t, T)} \tag{10.63}
\end{equation*}
$$
[[Andersen_Ch10#^eq10-63|Eq. (10.63)]]

Finally, swaption pricing is restated in annuity-normalized form in the affine setting, preserving continuity with the Gaussian approximation logic while allowing richer volatility structures. This closes the chapter by linking analytical model extension to practical calibration targets:

$$
\begin{equation*}
V_{\text {swaption }}\left(T_{0}\right)=A\left(T_{0}\right)\left(S\left(T_{0}\right)-c\right)^{+} \tag{10.66}
\end{equation*}
$$
[[Andersen_Ch10#^eq10-66|Eq. (10.66)]]

#### Key Topical Extensions and Implementation Details

* **Notations and First Steps (Section 10.1.1.1)**: The development of the Ho-Lee model begins by establishing basic notations and first steps, setting the short rate process as a drift-adjusted Brownian motion $dr(t) = \theta(t) dt + \sigma_r dW(t)$ without mean reversion, forming the simplest Gaussian interest rate framework.
* **Analysis and Comparison with HJM Approach**: The bottoms-up SDE development of the Ho-Lee model can be contrasted with the top-down Heath-Jarrow-Morton (HJM) approach. In the HJM framework, specifying a constant forward rate volatility $\sigma_f(t, T) = \sigma_r$ immediately dictates the risk-neutral forward rate drift $\mu_f(t, T) = \sigma_r^2(T-t)$. Integrating this yields the discount bond volatilities $\sigma_P(t, T) = \sigma_r(T-t)$ and the bond reconstitution formula directly, without evaluating any expectations.
* **Time-Stationarity and Caplet Hump**: A Gaussian HJM model is time-stationary if $\sigma_f(t, T)$ depends only on $T-t$. In a one-factor GSR model, this forces an exponentially decaying volatility structure $\sigma_f(t, T) = \sigma_r e^{-\varkappa (T-t)}$. Attempting to calibrate to a humped caplet volatility structure at time 0 by introducing a time-dependent mean reversion $\varkappa(t)$ violates time-stationarity: the volatility hump decays and eventually disappears for calendar times $t$ past the hump peak $t_p$. Thus, freezing $\varkappa$ as a constant and using a time-dependent volatility $\sigma_r(t)$ is recommended to preserve the shape of the volatility structure through time.
* **PDE and Spatial Boundary Conditions**: Solving the one-factor PDE requires spatial boundary conditions at the grid boundaries $x_{\min}$ and $x_{\max}$. Using a zero-curvature boundary condition $\partial^2 V / \partial x^2 = 0$ requires large boundaries set by a high confidence multiplier (e.g., $\alpha \approx 5\text{--}6$ standard deviations) to prevent boundary errors from propagating back to the origin, which wastes grid points. Determining spatial boundary conditions from the PDE itself resolves this by evaluating the SDE drift at the boundary using one-sided differences.
* **Upwinding**: When the grid cell size is large, the finite difference scheme can suffer from instability and spatial oscillations unless the grid spacing $\Delta_x$ satisfies $|y(t) - \varkappa(t)x|\Delta_x \le \sigma_r(t)^2$. Since $\sigma_r(t)^2$ is small, this is often violated at the edges of the grid where mean reversion is strong. Applying an upwinding scheme (using first-order one-sided differences in the drift term matching its sign) restores stability.
* **Regularity Issues in Affine Models**: The affine volatility function $v(r(t)) = \sqrt{\alpha + \beta r(t)}$ requires regularity conditions to ensure the square root is well-defined. If $\beta \neq 0$, the drift must satisfy $\varkappa \beta(\vartheta + \alpha / \beta) \geq 0$ to prevent $r(t)$ from leaving its domain $[-\alpha/\beta, \infty)$. To keep the rate strictly positive, the stronger Feller condition $\varkappa \beta(\vartheta + \alpha/\beta) \ge \frac{1}{2}\beta^2 \sigma^2$ is required (which simplifies to the classical CIR condition $2\varkappa \vartheta \ge \sigma^2$ when $\alpha = 0$).
* **Basic Problem in Affine Calibration**: In the Gaussian model ($\beta = 0$), the state loading $G(t, T)$ is deterministic and independent of $\sigma(t)$. However, in the general affine model ($\beta \neq 0$), the coefficient $b(t, T)$ depends on the future volatility path $\{\sigma(s)\}_{t \le s \le T}$. This represents the basic problem of calibration: the price of a short-expiry option depends on long-term volatilities, preventing the use of standard non-iterative bootstrapping algorithms and requiring iterative multi-dimensional optimization.
* **Monte Carlo Simulation**: To value path-dependent exotic options in the GSR model, we perform Monte Carlo simulation. Because the state variable SDE is linear and Gaussian, it can be integrated exactly over any time step without introducing discretization bias.

### Logical Bridge to Chapter 11

Chapter 11 keeps the one-factor short-rate setting but changes core state transformations and admissibility emphasis, introducing log-normal and shadow-rate style constructions where positivity and integrability/tail behavior become central technical constraints.

---

### Chapter 11: One-Factor Short Rate Models II

Chapter 11 stress-tests one-factor short-rate architectures by comparing modeling choices that preserve positivity against choices that preserve calibration robustness and numerical stability. The opening log-normal specification is attractive for non-negativity, but the chapter emphasizes that this does not automatically guarantee acceptable tail behavior or tractable calibration dynamics:

Conceptually, the chapter is about admissibility versus usability: a model may look structurally appealing because it enforces positivity, yet still become impractical if tail moments explode or if calibration surfaces become unstable. The narrative therefore tracks not just model form, but model failure modes.

$$
\begin{equation*}
d \ln r(t)=\varkappa(t)(\vartheta(t)-\ln r(t)) d t+\sigma_{r}(t) d W(t) \tag{11.3}
\end{equation*}
$$
[[Andersen_Ch11#^eq11-3|Eq. (11.3)]]

The chapter then broadens the diffusion family to power-type dynamics, making skew and level-dependence explicit through the exponent parameter $p$:

$$
\begin{equation*}
d r(t)=\varkappa(t)(\vartheta(t)-r(t)) d t+\sigma(t) r(t)^{p} d W(t), \quad p>0 \tag{11.18}
\end{equation*}
$$
[[Andersen_Ch11#^eq11-18|Eq. (11.18)]]

To capture richer volatility structure, the Fong-Vasicek-style two-state setup is represented through bond pricing loadings in $(r,z)$, separating curve fit from extra factor effects:

$$
\begin{equation*}
P(t, t+\delta)=e^{A(\delta)+r(t) B(\delta)+z(t) C(\delta)} \tag{11.21}
\end{equation*}
$$
[[Andersen_Ch11#^eq11-21|Eq. (11.21)]]

This decomposition clarifies what is spanned by the term structure and what is carried by latent volatility state variables. That distinction matters directly for hedge design, because spanned and unspanned components contribute differently to observable instrument sensitivities.

For implementation, the chapter recasts generic one-factor models in shifted coordinates and derives the valuation PDE kernel used in finite-difference workflows:

$$
\begin{equation*}
\frac{\partial V}{\partial t}+\mu_{x}(t, x) \frac{\partial V}{\partial x}+\frac{1}{2} \sigma_{x}(t, x)^{2} \frac{\partial^{2} V}{\partial x^{2}}=(x+f(0, t)) V, \tag{11.25}
\end{equation*}
$$
[[Andersen_Ch11#^eq11-25|Eq. (11.25)]]

The Arrow-Debreu transport equation then becomes the organizing object for forward/backward induction and term-structure fit recursion:

$$
\begin{equation*}
\frac{\partial G}{\partial t}+\mu_{x}(t, x ; \vartheta(t)) \frac{\partial G}{\partial x}+\frac{1}{2} \sigma_{x}(t, x)^{2} \frac{\partial^{2} G}{\partial x^{2}}=(x+f(0, t)) G, \quad(s, y) \text { fixed } \tag{11.29}
\end{equation*}
$$
[[Andersen_Ch11#^eq11-29|Eq. (11.29)]]

Calibration-oriented parameterizations are then expressed in the generic state-dependent-volatility short-rate form:

$$
\begin{equation*}
d r(t)=\mu_{r}(t, r(t) ; \vartheta(t)) d t+\sigma_{r}(t) \psi(r(t)) d W(t) \tag{11.40}
\end{equation*}
$$
[[Andersen_Ch11#^eq11-40|Eq. (11.40)]]

At this point the chapter effectively turns into an implementation manual: drift and level-map choices are treated as calibration architecture decisions, not just modeling decoration. The downstream PDE and Monte Carlo steps depend on this parameterization quality.

The appendix moves to Markov-functional construction, where a one-dimensional state process is mapped to numeraires and terminal payoffs:

$$
\begin{equation*}
d x(t)=\mu_{x}(t, x(t)) d t+\sigma_{x}(t, x(t)) d W(t), \quad x(0)=0 \tag{11.45}
\end{equation*}
$$
[[Andersen_Ch11#^eq11-45|Eq. (11.45)]]

Monte Carlo discounting and implementation closure are then expressed through decomposition of the discount factor into deterministic and stochastic components:

$$
\begin{equation*}
\exp \left(-\int_{0}^{t} r(u) d u\right)=P(0, t) \exp \left(-\int_{0}^{t} x(u) d u\right) \triangleq P(0, t) Y(t) \tag{11.53}
\end{equation*}
$$
[[Andersen_Ch11#^eq11-53|Eq. (11.53)]]

To make the chapter's admissibility and implementation logic explicit, five bridge equations are central. First, the shadow-rate positive-part construction shows how positivity can be imposed without rewriting the entire Gaussian core:

$$
\begin{equation*}
r(t)=\left(r^{*}(t)\right)^{+} \tag{11.19}
\end{equation*}
$$
[[Andersen_Ch11#^eq11-19|Eq. (11.19)]]

Second, unspanned-volatility structure is characterized by bond prices being insensitive to latent volatility states directly:

$$
\begin{equation*}
\frac{\partial P(t, T)}{\partial z_{i}(t)}=0, \quad i=1, \ldots, n, \tag{11.22}
\end{equation*}
$$
[[Andersen_Ch11#^eq11-22|Eq. (11.22)]]

Third, the state-shift map from short-rate to shifted-coordinate coefficients is what turns heterogeneous model families into one reusable PDE template:

$$
\begin{equation*}
\mu_{x}(t, x)=\mu_{r}(t, x+f(0, t))-\frac{\partial f(0, t)}{\partial t}, \quad \sigma_{x}(t, x)=\sigma_{r}(t, x+f(0, t)) \tag{11.26}
\end{equation*}
$$
[[Andersen_Ch11#^eq11-26|Eq. (11.26)]]

Fourth, forward transport of the Arrow-Debreu kernel provides the dual evolution equation needed for forward-from-backward induction workflows:

$$
\begin{align*}
-\frac{\partial G}{\partial s}-\frac{\partial}{\partial y}\left(\mu_{x}(s, y ; \vartheta(s)) G\right)+ & \frac{1}{2} \frac{\partial^{2}}{\partial y^{2}}\left(\sigma_{x}(t, y)^{2} G\right) \\
& =(y+f(0, t)) G, \quad(t, x) \text { fixed } \tag{11.30}
\end{align*}
$$
[[Andersen_Ch11#^eq11-30|Eq. (11.30)]]

Finally, simulation is organized through the coupled state-discount-factor SDE system, which is the implementation object underlying discounted Monte Carlo valuation:

$$
\begin{equation*}
d\binom{x(t)}{Y(t)}=\binom{\mu_{x}(t, x(t))}{-x(t) Y(t)} d t+\binom{\sigma_{x}(t, x(t))}{0} d W(t) \tag{11.54}
\end{equation*}
$$
[[Andersen_Ch11#^eq11-54|Eq. (11.54)]]

Taken together, these five bridges complete a coherent path from admissibility fixes to PDE transport and then to implementable simulation dynamics, so the chapter reads as a single modeling workflow rather than disconnected technical fragments.

[[Andersen_Ch11#^eq11-54|Eq. (11.54)]]

#### Key Topical Extensions and Implementation Details

* **The Black-Derman-Toy (BDT) Model**: The continuous-time limit of the BDT model is given by $r(t) = U(t) e^{\sigma_r(t) W(t)}$. By applying Ito's lemma, the dynamics are represented as $d \ln r(t) = (\vartheta(t) + \frac{\sigma'_r(t)}{\sigma_r(t)} \ln r(t)) dt + \sigma_r(t)dW(t)$. Since the mean reversion speed is locked by volatility as $\varkappa(t) = -\sigma'_r(t)/\sigma_r(t)$, BDT can exhibit undesirable "mean-fleeing" behavior (negative mean reversion) during periods of rising short-rate volatility.
* **Black-Karasinski (BK) Model**: Decoupling mean reversion from volatility, the BK model exogenously specifies mean reversion speed via the SDE $d \ln r(t) = \varkappa(t)(\vartheta(t) - \ln r(t)) dt + \sigma_r(t)dW(t)$. The short rate $r(t) = e^{x(t)}$ remains positive since its logarithm is a Gaussian process $x(t)$, but it lacks a closed-form bond reconstitution formula.
* **Sandmann-Sondermann Transformation**: A major drawback of continuously compounded log-normal short rate models is that the expectation of $e^{\int_t^{t'} r du}$ is infinite, leading to infinite expected returns on the money market account and infinite Libor futures rates. The Sandmann-Sondermann transformation resolves this by writing the log-normal dynamics for a discretely compounded rate $r_d(t)$ defined by $e^{r(t)\delta} = 1 + r_d(t)\delta$ for a finite compounding interval $\delta > 0$, guaranteeing finite expectations.
* **Calibration to Initial Yield Curve**: Finding the time-dependent drift parameter $\vartheta(t)$ to match the initial yield curve is traditionally done using slow, brute-force binomial tree backward induction. The forward Kolmogorov (Fokker-Planck) equation makes this efficient: a single forward sweep computes the Arrow-Debreu green's function densities, allowing $\vartheta(t)$ to be solved recursively.
* **The Dybvig Parameterization**: Parameterizing the short rate as $r(t) = s(t) + \vartheta(t)$ (where $s(t)$ is a tractably modeled state variable and $\vartheta(t)$ is fitted to initial bond prices) works well in Gaussian models. In non-Gaussian models (such as CIR where $s(t)$ is non-negative), this parameterization is highly problematic because the state-dependent volatility of $s(t)$ becomes $\sigma \sqrt{r(t) - \vartheta(t)}$, shifting the boundary and distorting volatility as $\vartheta(t)$ evolves.
* **The Hagan and Woodward Parameterization**: This parameterization shifts the yield-curve fitting parameter from the drift to the deflator $N(t) = \frac{1}{P(0,t)} e^{h(t, x(t)) + a(t)}$. Calibrating $a(t)$ requires only a single forward PDE pass and is independent of the initial yield curve. For the subclass $h(t,x) = h(t)x$, the scaling function $h(t)$ controls the effective mean reversion speed via $\varkappa_r(t) = -h''(t)/h'(t)$.
* **SDE Discretization**: While the Euler-type scheme is simple, it has first-order weak convergence. The second-order Milstein scheme for the coupled state-discount system uses operator expansions $\mathcal{L}_0, \mathcal{L}_1$ to expand the drift and diffusion terms, which significantly reduces discretization bias in Monte Carlo simulation.
* **Markov-Functional Models (Appendix 11.A)**:
    * **State Process and Numeraire Mapping**: MF models specify the state process $dx(t) = \sigma(t)e^{\varkappa t} dW^*(t)$ as a Gaussian process under the terminal measure $\mathrm{Q}^*$. The final discount factor is mapped as $P(t, T^*) = H(t, x(t))$, from which other bonds are recovered via $P(t, T) = H(t, x(t)) \mathrm{E}_t^* [1/H(T, x(T))]$.
    * **Libor and Swap MF Parameterizations**: The mapping function $H(t,x)$ is established inductively. Under the Libor MF parameterization, forward rates are mapped as $L_n(T_n) = c_n L_n(0) \exp(k_n x(T_n) - \frac{1}{2}k_n^2 v_n^2)$, where the constants $c_n$ are calibrated to fit the yield curve. Under the Swap MF parameterization, core swap rates $S_n(T_n) = s_n(x(T_n))$ are mapped to approximate log-normal swap rate distributions, which is especially useful for pricing Bermudan swaptions.
    * **Non-Parametric Calibration**: By equating the theoretical pricing integral to observed option prices, $V(0; K, T) = P(0, T^*) \int p(0,0; z, T) \frac{g(z; K)}{H(T, z)} dz$, the numeraire function $H(T, \cdot)$ can be calibrated non-parametrically. However, this integral equation requires numerical grid inversion, which is prone to high instability and can imply unrealistic volatility smiles through time, making low-dimensional parametric forms preferable.
    * **Comparisons**: Compared to the quasi-Gaussian models of Chapter 13, MF models require numerical integration to reconstitute discount bonds at each step, whereas QG models have fully analytical bond reconstitution formulas.

### Logical Bridge to Chapter 12

Chapter 12 lifts the one-factor constructions into truly multi-factor state systems, where matrix-valued mean reversion and cross-factor correlation become first-class modeling objects.

---

### Chapter 12: Multi-Factor Short Rate Models

Chapter 12 is the transition from one-dimensional rate states to genuinely multi-factor term-structure dynamics. It starts from HJM consistency and shows how finite-dimensional Markov representations recover tractable forward and bond formulas in vector form:

The core theme is dimensional control: add factors to capture realistic curve-shape and correlation effects, but keep state transformations structured enough that pricing remains computationally predictable. Every block in the chapter can be read as a trade-off between expressiveness and solvability.

$$
\begin{align*}
d f(t, T) & =\sigma_{f}(t, T)^{\top} \sigma_{P}(t, T) d t+\sigma_{f}(t, T)^{\top} d W(t) \\
& =\sigma_{f}(t, T)^{\top} \int_{t}^{T} \sigma_{f}(t, u) d u d t+\sigma_{f}(t, T)^{\top} d W(t) \tag{12.1}
\end{align*}
$$
[[Andersen_Ch12#^eq12-1|Eq. (12.1)]]

Diagonalization and state reduction produce a compact vector SDE backbone for two-factor Gaussian models:

$$
\begin{equation*}
d x(t)=(y(t) 1-\varkappa(t) x(t)) d t+\sigma_{x}^{*}(t) d W^{*}(t) \tag{12.21}
\end{equation*}
$$
[[Andersen_Ch12#^eq12-21|Eq. (12.21)]]

Curve reconstruction is then made explicit through deterministic loadings acting on the state vector:

$$
\begin{align*}
f(t, T) & =f(0, T)+M(t, T)^{\top}(x(t)+y(t) G(t, T))  \tag{12.22}\\
P(t, T) & =\frac{P(0, T)}{P(0, t)} \exp \left(-G(t, T)^{\top} x(t)-\frac{1}{2} G(t, T)^{\top} y(t) G(t, T)\right)
\end{align*}
$$
[[Andersen_Ch12#^eq12-22|Eq. (12.22)]]

These reconstruction identities are the bridge between latent-state modeling and observable term-structure outputs. They are also the reason calibration can still be organized instrument-by-instrument even when the underlying dynamics are vector-valued.

In the two-factor block, state dynamics are represented component-wise to expose factor-specific drift and covariance structure:

$$
\begin{align*}
& d x_{1}(t)=\left(\vartheta_{1}(t)-\varkappa_{1}(t) x_{1}(t)\right) d t+\sigma_{11}(t) d W_{1}(t)+\sigma_{12}(t) d W_{2}(t)  \tag{12.38}\\
& d x_{2}(t)=\left(\vartheta_{2}(t)-\varkappa_{2}(t) x_{2}(t)\right) d t+\sigma_{21}(t) d W_{1}(t)+\sigma_{22}(t) d W_{2}(t)
\end{align*}
$$
[[Andersen_Ch12#^eq12-38|Eq. (12.38)-(12.39)]]

For swaptions, annuity-measure pricing remains the central computational interface. This keeps the model aligned with liquid market quoting conventions and simplifies direct calibration mapping:

$$
\begin{equation*}
V_{\text {swaption }}(0)=A(0) \mathrm{E}^{A}\left(\left(S\left(T_{0}\right)-c\right)^{+}\right) \tag{12.45}
\end{equation*}
$$
[[Andersen_Ch12#^eq12-45|Eq. (12.45)]]

This keeps continuity with one-factor chapters: despite higher state dimension, the valuation interface remains expectation under a convenient numeraire. That continuity is important for reusing calibration infrastructure.

The chapter then generalizes to affine multi-factor diffusions, retaining transform tractability:

$$
\begin{equation*}
d x(t)=\mu(x(t)) d t+\sigma(x(t)) d W(t) \tag{12.53}
\end{equation*}
$$
[[Andersen_Ch12#^eq12-53|Eq. (12.53)]]

From this point onward, admissible state-domain and regularity constraints become first-order concerns. Without them, the formal affine machinery exists on paper but fails in robust production calibration.

$$
\begin{equation*}
P(t, T, x)=e^{A(T-t)-B(T-t)^{\top} x} \tag{12.61}
\end{equation*}
$$
[[Andersen_Ch12#^eq12-61|Eq. (12.61)]]

The exponential-affine bond form is the practical payoff of the affine setup: it converts the pricing map into deterministic Riccati objects plus current state, making repeated valuation significantly cheaper.

Quadratic-Gaussian structure is introduced to generate smile while keeping Gaussian state evolution:

$$
\begin{equation*}
r(t)=x(t)^{\top} \widetilde{\gamma}(t) x(t)+\mathbf{1}^{\top} x(t)+a(t) \tag{12.67}
\end{equation*}
$$
[[Andersen_Ch12#^eq12-67|Eq. (12.67)]]

Quadratic short-rate mapping then introduces smile capacity while preserving Gaussian state evolution, which is exactly why the model is attractive in practice: nonlinear payoff behavior with relatively stable simulation and PDE machinery.

Approximation and calibration closure in annuity measure are built on mixture-state terminal densities and quadratic swap-rate expansions:

$$
\begin{align*}
\psi^{A}(z) & =\sum_{i=0}^{N-1} w_{i}^{A} \phi\left(z ; m^{T_{i+1}}\left(0, T_{0}, 0\right), \nu^{T_{\iota+1}}\left(0, T_{0}, 0\right)\right)  \tag{12.86}\\
w_{i}^{A} & =\frac{\tau_{i} P\left(0, T_{i+1}\right)}{A(0)}, \quad i=0, \ldots, N-1
\end{align*}
$$
[[Andersen_Ch12#^eq12-86|Eq. (12.86)]]

This annuity-mixture characterization is what enables fast approximation hierarchies: once terminal-state distribution is explicit enough, higher-order swap-rate expansions can be calibrated with controlled computational cost.

$$
\begin{equation*}
h_{S}=\nabla^{\top} S\left(T_{0}, z^{*}\right), \quad \gamma_{S}=\frac{1}{2} \nabla^{\top} \nabla S\left(T_{0}, z^{*}\right), \tag{12.91}
\end{equation*}
$$
[[Andersen_Ch12#^eq12-91|Eq. (12.91)]]

The gradient and Hessian objects are the local geometry of the swap-rate map. They determine how much curvature must be preserved in approximation to retain pricing fidelity near the chosen expansion point.

$$
\begin{equation*}
\widehat{\alpha}_{S}=-\mathrm{E}^{A}\left(z\left(T_{0}\right)^{\top} R^{\top} \widehat{\gamma}_{S} R z\left(T_{0}\right)+h_{S}^{\top} z\left(T_{0}\right)\right) \tag{12.96}
\end{equation*}
$$
[[Andersen_Ch12#^eq12-96|Eq. (12.96)]]

The correction term enforces repricing consistency, preventing approximation drift in the forward level. This is a key requirement if approximation-based calibration is to remain stable across maturities.

$$
\begin{aligned}
\widehat{m} & =m-\nu(I+2 Q \nu)^{-1}(2 Q m+u)  \\
\widehat{\nu} & =\nu(I+2 Q \nu)^{-1} \tag{12.98}
\end{aligned}
$$
[[Andersen_Ch12#^eq12-98|Eq. (12.98)]]

Five inner-node equations complete the derivation bridges between correlation structure, reconstruction, affine admissibility, measure change, and approximation closure. The correlation bridge is:

$$
\begin{align*}
\operatorname{Var}_{t}(d f(t, T)) & =\sigma_{1}(t)^{2} e^{-2 \int_{t}^{T} \varkappa_{1}(u) d u} b(t, T, T) \\
\rho\left(t, T_{1}, T_{2}\right) & =\operatorname{Corr}_{t}\left(d f\left(t, T_{1}\right), d f\left(t, T_{2}\right)\right)=\frac{b\left(t, T_{1}, T_{2}\right)}{\sqrt{b\left(t, T_{1}, T_{1}\right) b\left(t, T_{2}, T_{2}\right)}} \tag{12.24}
\end{align*}
$$
[[Andersen_Ch12#^eq12-24|Eq. (12.24)]]

This correlation identity clarifies how cross-maturity comovement is constrained by the same loading structure that controls instantaneous variance, so calibration cannot tune one object independently of the other.

The two-factor bond reconstruction used in option decomposition is:

$$
\begin{equation*}
P(T, T+\Delta)=\frac{P(0, T+\Delta)}{P(0, T)} e^{A(T, T+\Delta)-x_{1}(T) G_{1}(T, T+\Delta)-x_{2}(T) G_{2}(T, T+\Delta)} \tag{12.40}
\end{equation*}
$$
[[Andersen_Ch12#^eq12-40|Eq. (12.40)]]

Because this reconstruction is explicit in both factors, it becomes possible to translate market-observable bond prices back into latent states before applying option payoff mappings.

Affine regularity and factor loading structure are encoded by the diagonal diffusion arrangement:

$$
d x(t)=(a-b x(t)) d t+\Sigma\left(\begin{array}{cccc}
\sqrt{v_{1}(x(t))} & 0 & \ddots & 0  \tag{12.55}\\
0 & \sqrt{v_{2}(x(t))} & \ddots & \ddots \\
\ddots & \ddots & \ddots & 0 \\
0 & \ddots & 0 & \sqrt{v_{d}(x(t))}
\end{array}\right) d W(t)
$$
[[Andersen_Ch12#^eq12-55|Eq. (12.55)]]

The diagonalized form is not merely cosmetic: it is the bridge that makes moment and transform calculations modular factor-by-factor in practical implementation pipelines.

Forward-measure dynamics are then transported via Girsanov shift:

$$
\begin{equation*}
d W^{T}(t)=d W(t)+g(t)(2 \gamma(t, T) z(t)+h(t, T)) d t \tag{12.70}
\end{equation*}
$$
[[Andersen_Ch12#^eq12-70|Eq. (12.70)]]

Once this measure change is made explicit, the remaining pricing machinery can stay in forward measure without repeatedly re-deriving drift adjustments for each instrument class.

And quadratic swap-rate approximation closure is expressed through the repricing-adjusted quadratic form:

$$
\begin{gather*}
S\left(T_{0}, z\right) \approx S_{q}\left(T_{0}, z\right) \\
S_{q}\left(T_{0}, z\right)=z^{\top} \gamma_{S} z+h_{S}^{\top} z-\mathrm{E}^{A}\left(z\left(T_{0}\right)^{\top} \gamma_{S} z\left(T_{0}\right)+h_{S}^{\top} z\left(T_{0}\right)\right)+S(0) \tag{12.92}
\end{gather*}
$$
[[Andersen_Ch12#^eq12-92|Eq. (12.92)]]

#### Key Topical Extensions and Implementation Details

* **Development from Separability Condition**: The multi-factor HJM model can be reduced to a finite-dimensional Markovian model by imposing a separability condition on instantaneous forward rate volatilities: $\sigma_f(t, T) = g(t) h(T)$ where $g(t)$ is a stochastic matrix process and $h(T)$ is a deterministic vector. This separates the state variables into mean-reverting components $x(t) \in \mathbb{R}^d$ and a symmetric path-dependent matrix component $y(t) \in \mathbb{R}^{d \times d}$. Further changes of variables scale the state process to isolate mean-reversion effects.
* **Classical Development and Diagonalization**: The classical approach constructs the vector short rate process directly. If the mean reversion matrix $K(t)$ is diagonalizable, we can define a coordinate transformation that diagonalizes $K(t)$ to decouple the SDEs of the state variables, facilitating analytical tractability.
* **Two-Dimensional Jamshidian Decomposition (Section 12.1.6.1)**: While standard Jamshidian decomposition applies only to one-factor models, in a two-factor Gaussian model, swaptions can be priced by defining a two-dimensional Jamshidian decomposition boundary in the state space that separates the exercise region from the non-exercise region.
* **Calibration and Parameterization via Benchmark Rates**: Rather than parameterizing volatility vectors $\lambda_k(t)$ directly, we can calibrate them to a covariance matrix of benchmark rates. Defining $d$ benchmark tenors $\{\delta_i\}$ and a target benchmark covariance $C^f(t)$, we construct the benchmark matrix $H^f(t)$ and solve the matrix equation $H^f(t) g(t)^\top = C^f(t)^\top$ to link the benchmark rate covariance directly to the factor loading matrix $g(t)$.
* **Some Concrete Models**:
    * **Fong-Vasicek Model**: A two-factor affine stochastic volatility model where short rate $r(t)$ satisfies $d r(t) = \varkappa_r(\vartheta_r - r(t))dt + \sqrt{z(t)} dW_1(t)$ and the variance process $z(t)$ follows a CIR-type SDE driven by correlated Brownian motions.
    * **Longstaff-Schwartz Model**: A two-factor CIR model where the short rate is a sum of two independent state variables $r(t) = x_1(t) + x_2(t)$. Because the state variables are independent, discount bond prices are recovered in closed form as a product of two independent one-factor CIR bond pricing formulas.
* **Quadratic Gaussian (QG) Model Basics**: In the QG model, zero-coupon bonds are exponentials of quadratic forms: $-\ln P(t,T) = z(t)^\top \gamma(t,T) z(t) + h(t,T)^\top z(t) + a(t,T) - \ln(P(0,T)/P(0,t))$, where the matrix $\gamma(t,T)$ and the vector $h(t,T)$ satisfy Riccati ordinary differential equations (ODEs), and the scalar $a(t,T)$ is solved from expectations under the $t$-forward measure.
* **QG Numerical Methods**: Because the state variables $z(t)$ in a QG model follow the same Gaussian processes as in linear Gaussian models, numerical pricing is highly efficient. Multi-dimensional PDE schemes and Monte Carlo simulation carry over unchanged from the linear Gaussian case, allowing bias-free Monte Carlo simulation over any step length.
* **Appendix 12.A: Quadratic Forms of Gaussian Vectors**: The moment-generating function for a quadratic form of a Gaussian vector $Z \sim \mathcal{N}(m, \nu)$ is evaluated as $\Psi(u, Q; m, \nu) = \mathrm{E}(\exp(Z^\top Q Z + u^\top Z))$. Taking derivatives of this function yields closed-form formulas for the mean $\mathrm{E}(Z^\top Q Z + u^\top Z) = m^\top Q m + u^\top m + \operatorname{tr}(Q\nu)$ and variance $\operatorname{Var}(Z^\top Q Z + u^\top Z) = (2 m^\top Q + u^\top) \nu (2Qm + u) + 2\operatorname{tr}(Q\nu Q\nu)$.

### Logical Bridge to Chapter 13

Chapter 13 preserves tractable state reconstruction but replaces deterministic/Gaussian volatility structure by local and stochastic volatility layers projected to swap-rate dynamics.

---

### Chapter 13: Quasi-Gaussian Local and Stochastic Volatility

Chapter 13 develops quasi-Gaussian models as a practical compromise between explicit term-structure reconstruction and richer volatility/skew behavior. The baseline two-state backbone keeps the same reconstruction philosophy as earlier Gaussian models:

The chapter's organizing idea is controlled deformation: keep the reconstruction skeleton that made Gaussian models tractable, then deform the volatility structure just enough to match skew and smile without losing numerical control.

$$
\begin{align*}
d x(t) & =(y(t)-\varkappa(t) x(t)) d t+\sigma_{r}(t, \omega) d W(t)  \tag{13.4}\\
d y(t) & =\left(\sigma_{r}(t, \omega)^{2}-2 \varkappa(t) y(t)\right) d t \\
x(0) & =y(0)=0
\end{align*}
$$
[[Andersen_Ch13#^eq13-4|Eq. (13.4)]]

Given this state backbone, the next step is to reconstruct bond prices explicitly from $(x,y)$ so that payoff pricing can be written as a direct state-functional map.

$$
\begin{equation*}
P(t, T, x, y)=\frac{P(0, T)}{P(0, t)} \exp \left(-G(t, T) x-\frac{1}{2} G(t, T)^{2} y\right) \tag{13.5}
\end{equation*}
$$
[[Andersen_Ch13#^eq13-5|Eq. (13.5)]]

This reconstitution formula is the structural invariant across model variants. It allows local-vol and stochastic-vol enrichments to be layered without rewriting the entire pricing interface.

Local-volatility enrichment is introduced by letting diffusion strength depend on current state variables:

$$
\begin{align*}
& d x(t)=(y(t)-\varkappa(t) x(t)) d t+\sigma_{r}(t, x(t), y(t)) d W(t)  \tag{13.10}\\
& d y(t)=\left(\sigma_{r}(t, x(t), y(t))^{2}-2 \varkappa(t) y(t)\right) d t
\end{align*}
$$
[[Andersen_Ch13#^eq13-10|Eq. (13.10)-(13.11)]]

State-dependent volatility here is not just a cosmetic extension; it is the mechanism that transports level effects into option surface asymmetry. The challenge is to keep that mechanism compatible with fast calibration routines.

Pricing focus is shifted to projected swap-rate dynamics under annuity measure, which is the calibration-facing representation:

$$
\begin{equation*}
d S(t)=\varphi(t, S(t)) d W^{A}(t) \tag{13.18}
\end{equation*}
$$
[[Andersen_Ch13#^eq13-18|Eq. (13.18)]]

Projection is the computational compression step: instead of carrying full state complexity directly into every option valuation, the model pushes complexity into an effective local volatility map on swap rates.

Linearized local-vol parameterizations provide fast approximations for strip-level calibration:

$$
\begin{equation*}
\sigma_{r}(t, x, y)=\lambda_{r}(t)\left(\alpha_{r}(t)+b_{r}(t) x\right) \tag{13.29}
\end{equation*}
$$
[[Andersen_Ch13#^eq13-29|Eq. (13.29)]]

This parameterization introduces explicit handles for level and skew control. In practice, those handles are what make calibration explainable rather than purely black-box numerical fitting.

$$
\begin{equation*}
d S(t) \approx \lambda_{S}(t)\left(b_{S}(t) S(t)+\left(1-b_{S}(t)\right) S(0)\right) d W^{A}(t) \tag{13.33}
\end{equation*}
$$
[[Andersen_Ch13#^eq13-33|Eq. (13.33)]]

The linearized swap-rate form is the speed-accuracy compromise used for production calibration cycles. It sacrifices exactness in dynamics to gain stable and interpretable calibration updates.

Stochastic-volatility scaling then augments the framework with a latent variance factor while retaining projected swap-rate interfaces:

$$
\begin{align*}
& d x(t)=(y(t)-\varkappa(t) x(t)) d t+\sqrt{z(t)} \sigma_{r}(t, x(t), y(t)) d W(t)  \tag{13.64}\\
& d y(t)=\left(z(t) \sigma_{r}(t, x(t), y(t))^{2}-2 \varkappa(t) y(t)\right) d t \\
& d z(t)=\theta\left(z_{0}-z(t)\right) d t+\eta(t) \sqrt{z(t)} d Z(t)
\end{align*}
$$
[[Andersen_Ch13#^eq13-64|Eq. (13.64)]]

Adding a stochastic variance factor lifts the model from static skew fitting to dynamic smile behavior. The key question becomes whether this extra state dimension improves fit enough to justify simulation and calibration overhead.

$$
\begin{equation*}
d S(t)=\sqrt{z(t)} \varphi(t, S(t)) d U^{A}(t) \tag{13.79}
\end{equation*}
$$
[[Andersen_Ch13#^eq13-79|Eq. (13.79)]]

Projected stochastic-vol swap-rate dynamics preserve the same pricing interface while carrying richer variance dynamics. This continuity is what allows gradual model upgrades in existing pricing systems.

Numerically, the chapter combines full two-state PDE formulations with reduced one-state approximations and appendix asymptotic machinery:

$$
\begin{align*}
\frac{\partial V}{\partial t}(t, x, u)+(\mathcal{L} V)(t, x, u) & =(f(0, t)+x) V(t, x, u), \quad 0 \leq t<T  \tag{13.55}\\
V(T, x, u) & =V(x, u)
\end{align*}
$$
[[Andersen_Ch13#^eq13-55|Eq. (13.55)]]

The two-state PDE is the fidelity benchmark. Reduced and asymptotic methods should be understood as controlled approximations against this benchmark, not replacements without error accounting.

$$
\begin{align*}
d x(t) & =(\widetilde{y}(t, x(t))-\varkappa(t) x(t)) d t+\widetilde{\sigma}_{r}(t, x(t)) d W(t)  \tag{13.60}\\
P(t, T) & =\frac{P(0, T)}{P(0, t)} \exp \left(-G(t, T) x(t)-\frac{1}{2} G(t, T)^{2} \widetilde{y}(t, x(t))\right) \\
r(t) & =f(0, t)+x(t)
\end{align*}
$$
[[Andersen_Ch13#^eq13-60|Eq. (13.60)]]

The one-state reduction is valuable because it makes calibration loops dramatically faster when market regimes are moderate. The chapter is explicit that this speedup has a regime-of-validity boundary.

$$
\begin{align*}
& d x(t)=(\bar{y}(t)-\varkappa x(t)) d t+\sigma v(x(t)) d W(t)  \tag{13.84}\\
& d \bar{y}(t)=\left(\sigma^{2}-2 \varkappa \bar{y}(t)\right) d t
\end{align*}
$$
[[Andersen_Ch13#^eq13-84|Eq. (13.84)]]

Five additional bridge equations complete the local-vol and stochastic-vol calibration chain. The local-vol baseline anchor is:

$$
\begin{equation*}
\sigma_{r}^{0}(t) \triangleq \sigma_{r}(t, 0,0) \tag{13.12}
\end{equation*}
$$
[[Andersen_Ch13#^eq13-12|Eq. (13.12)]]

This baseline object is the reference scale against which all later state-dependent perturbations are measured, keeping the calibration narrative tied to a concrete deterministic anchor.

Projected swap-rate local variance is linked by the conditional expectation relation:

$$
\begin{equation*}
\varphi(t, s)^{2}=\mathrm{E}^{A}\left(\left.\left(\frac{\partial S}{\partial x}(t, x(t), y(t)) \sigma_{r}(t, x(t), y(t))\right)^{2} \right\rvert\, S(t)=s\right) \tag{13.19}
\end{equation*}
$$
[[Andersen_Ch13#^eq13-19|Eq. (13.19)]]

The projection formula explains where computational difficulty lives: in the conditional expectation map from latent states to observed swap-rate levels under annuity measure.

Mean-reversion calibration bridge through variance ratios is:

$$
\begin{equation*}
\frac{\operatorname{Var}\left(S_{n, m_{1}}\left(T_{n}\right)\right)}{\operatorname{Var}\left(S_{n, m_{2}}\left(T_{n}\right)\right)} \approx \frac{\int_{0}^{T_{n}}\left(\sigma_{r}^{0}(t) \frac{\partial S_{n, m_{1}}}{\partial x}(t, 0,0)\right)^{2} d t}{\int_{0}^{T_{n}}\left(\sigma_{r}^{0}(t) \frac{\partial S_{n, m_{2}}}{\partial x}(t, 0,0)\right)^{2} d t} \tag{13.40}
\end{equation*}
$$
[[Andersen_Ch13#^eq13-40|Eq. (13.40)]]

Using variance ratios reduces absolute level sensitivity and therefore provides a more stable calibration handle for mean-reversion than direct volatility matching alone.

Stochastic-volatility extension is anchored by the CIR-type variance process:

$$
\begin{equation*}
d z(t)=\theta\left(z_{0}-z(t)\right) d t+\eta(t) \sqrt{z(t)} d Z(t), \quad z(0)=z_{0}=1 \tag{13.74}
\end{equation*}
$$
[[Andersen_Ch13#^eq13-74|Eq. (13.74)]]

This additional factor introduces persistent variance clustering, which is precisely the behavior missing in purely local-vol specifications during stressed smile regimes.

Appendix asymptotic closure then uses simplified forward-measure dynamics:

$$
\begin{equation*}
d \widetilde{x}(t)=-\varkappa \widetilde{x}(t) d t+\sigma v(\widetilde{x}(t)) d W^{T}(t) \tag{13.86}
\end{equation*}
$$
[[Andersen_Ch13#^eq13-86|Eq. (13.86)]]

#### Key Topical Extensions and Implementation Details

* **One-Factor QG Model Definition**: The one-factor QG model is defined in shifted coordinates as a system of SDEs for $x(t) \in \mathbb{R}$ and $y(t) \in \mathbb{R}$, where the volatility of instantaneous forward rates satisfies the separability condition $\sigma_f(t, T, \omega) = g(t, \omega) h(T)$.
* **Calibrating Mean Reversion to Inter-Temporal Correlations**: The mean reversion parameter $\varkappa$ determines the decay rate of forward rate correlations over time. In calibration, $\varkappa$ can be fitted to historical or market-implied inter-temporal correlation curves (the correlation between forward rates of the same maturity observed at different calendar times), ensuring the model reproduces realistic curve dynamics.
* **Finite Difference and Monte Carlo Methods**:
    * **Finite Difference Methods**: Solving the pricing PDE requires a grid in three variables $(t, x, y)$. Upwinding schemes are implemented to handle the strong drift term in the boundary zones.
    * **Monte Carlo Simulation**: Because the auxiliary variable $y(t)$ has a deterministic SDE given the path of $x(t)$, it can be integrated exactly over each step, eliminating any discretization bias in $y$ and making Monte Carlo simulation highly robust and efficient.
* **Single-State Approximations**: To avoid solving a two-dimensional spatial PDE in $(x, y)$, the QG model can be reduced to a one-factor state space. By replacing the auxiliary process $y(t)$ with a deterministic function of $x(t)$ or freezing it at its expectation $\mathrm{E}(y(t)) = \bar{y}(t)$, the PDE dimension is reduced to one, accelerating calibration loops.
* **Stochastic Volatility QG Definition and Correlation**: Layering a CIR variance process $z(t)$ yields the SV QG model. Allowing non-zero correlation between the state and variance shocks introduces skew/smile asymmetry. Under the forward measure, this correlation induces a measure-adjusted drift multiplier $\alpha(t)$ in the variance SDE.
* **General Multi-Factor Model (d > 1)**: The multi-factor QG model uses a separable matrix volatility $\sigma_f(t, T, \omega) = g(t, \omega) h(T)$, which yields a state vector $x(t) \in \mathbb{R}^d$ and a symmetric path-dependent matrix $y(t) \in \mathbb{R}^{d \times d}$. Discount bonds are given by $P(t, T) = \frac{P(0,T)}{P(0,t)} \exp(-G(t,T)^\top x(t) - \frac{1}{2} G(t,T)^\top y(t) G(t,T))$.
* **Multi-Factor Mean Reversions, Correlations, and Numerical Schemes**:
    * **Mean Reversions**: We recommend choosing $d$ distinct fixed mean reversion parameters (e.g., spanning the interval $[0, 1.2]$) corresponding to key benchmark tenors.
    * **Correlations**: The correlation matrix between benchmark rates can be calibrated to market-observed prices of CMS spread options to capture the correlation smile.
    * **Numerical Schemes**: For $d > 1$, finite difference methods become impractical because a $d$-factor model requires $d(d+1)/2$ auxiliary $y$ variables (e.g., a 2-factor model requires a 5D PDE). Therefore, Monte Carlo simulation is the preferred numerical scheme.
* **Appendix 13.A: Proof of Theorem 13.1.14**: The proof of the density approximation formula for the short-rate transition density uses asymptotic expansion techniques under simplified forward-measure dynamics.

### Logical Bridge to Chapter 14

Chapter 14 moves from short-rate-centric state reconstructions to tenor-grid Libor state dynamics, where no-arbitrage across measures drives drift and calibration structure.

---

### Chapter 14: Libor Market Model I

Chapter 14 sets up the Libor Market Model directly on tenor-forward rates and then derives consistent measure transforms across forward, spot, and annuity numeraires. The starting point is discount-bond reconstruction from tenor Libors:

Compared with short-rate models, the modeling primitive changes here from a single instantaneous rate to a tenor vector. This shift increases dimensional realism but also makes measure-consistency and drift construction the central technical burden.

$$
\begin{equation*}
P\left(t, T_{n}\right)=P\left(t, T_{q(t)}\right) \prod_{i=q(t)}^{n-1}\left(1+L_{i}(t) \tau_{i}\right)^{-1} \tag{14.3}
\end{equation*}
$$
[[Andersen_Ch14#^eq14-3|Eq. (14.3)]]

Forward-rate SDEs are then specified with measure-dependent drifts. This step is where no-arbitrage consistency is converted into a simulation-ready dynamic system across tenor points:

$$
\begin{equation*}
d L_{n}(t)=\sigma_{n}(t)^{\top}\left(\mu_{n}(t) d t+d W(t)\right) \tag{14.4}
\end{equation*}
$$
[[Andersen_Ch14#^eq14-4|Eq. (14.4)]]

and Girsanov links between adjacent forward measures provide the drift-transport mechanism that carries no-arbitrage constraints consistently across tenor dates during model evolution:

$$
\begin{equation*}
d W^{n}(t)=d W^{n+1}(t)-\frac{\tau_{n} \sigma_{n}(t)}{1+\tau_{n} L_{n}(t)} d t \tag{14.6}
\end{equation*}
$$
[[Andersen_Ch14#^eq14-6|Eq. (14.6)]]

This measure chain is the no-arbitrage backbone of LMM. Most practical implementation issues in calibration and simulation can be traced back to how well this drift transport is approximated numerically.

Volatility structure is parameterized through deterministic loadings and level maps:

$$
\begin{equation*}
\sigma_{n}(t)=\lambda_{n}(t) \varphi\left(L_{n}(t)\right) \tag{14.12}
\end{equation*}
$$
[[Andersen_Ch14#^eq14-12|Eq. (14.12)]]

with stochastic-volatility extension introducing a global variance factor and coupled tenor drifts, thereby transmitting smile dynamics across maturities while preserving a tractable calibration interface:

$$
\begin{align*}
d L_{n}(t) & =\sqrt{z(t)} \varphi\left(L_{n}(t)\right) \lambda_{n}(t)^{\top}\left(\sqrt{z(t)} \mu_{n}(t) d t+d W^{B}(t)\right)  \tag{14.16}\\
\mu_{n}(t) & =\sum_{j=q(t)}^{n} \frac{\tau_{j} \varphi\left(L_{j}(t)\right) \lambda_{j}(t)}{1+\tau_{j} L_{j}(t)}
\end{align*}
$$
[[Andersen_Ch14#^eq14-16|Eq. (14.16)]]

The global variance factor couples smile dynamics across maturities while preserving a manageable volatility parameterization. This is where model flexibility and calibration identifiability begin to pull in opposite directions.

Swaption approximations are derived by projecting high-dimensional forward-rate dynamics to effective swap-rate diffusion:

$$
\begin{equation*}
d S(t)=\sqrt{z(t)} \varphi(S(t)) \sum_{n=j}^{k-1} w_{n}(t) \lambda_{n}(t)^{\top} d W^{A}(t) \tag{14.30}
\end{equation*}
$$
[[Andersen_Ch14#^eq14-30|Eq. (14.30)]]

This projected swap-rate diffusion defines the effective variance channel used in practical swaption approximation. Pricing is then completed under annuity measure through a standard expectation form.

$$
\begin{equation*}
V_{\text {swaption }}(0)=A(0) \mathrm{E}^{A}\left(\left(S\left(T_{j}\right)-c\right)^{+}\right) \tag{14.32}
\end{equation*}
$$
[[Andersen_Ch14#^eq14-32|Eq. (14.32)]]

Annuity-measure swaption valuation remains the practical anchor for calibration. Even when high-dimensional forward dynamics are complex, calibration is still organized around this projected pricing interface.

The numerical section then focuses on robust simulation drift treatment under spot measure:

$$
\begin{equation*}
d L_{n}(t)=\sigma_{n}(t)^{\top}\left(\mu_{n}(t) d t+d W^{B}(t)\right), \quad \mu_{n}(t)=\sum_{j=q(t)}^{n} \frac{\tau_{j} \sigma_{j}(t)}{1+\tau_{j} L_{j}(t)} \tag{14.56}
\end{equation*}
$$
[[Andersen_Ch14#^eq14-56|Eq. (14.56)]]

and martingale-preserving discretization through deflated bond states, which stabilize long-horizon Monte Carlo behavior and reduce drift-discretization bias:

$$
\begin{equation*}
U\left(t, T_{n+1}\right)=\frac{P\left(t, T_{n+1}\right)}{B(t)} \tag{14.74}
\end{equation*}
$$
[[Andersen_Ch14#^eq14-74|Eq. (14.74)]]

Five additional bridge equations make the LMM derivation chain explicit from measure choice to discretization consistency. Under spot measure, drift coupling is:

$$
\begin{equation*}
d L_{n}(t)=\sigma_{n}(t)^{\top}\left(\sum_{j=q(t)}^{n} \frac{\tau_{j} \sigma_{j}(t)}{1+\tau_{j} L_{j}(t)} d t+d W^{B}(t)\right) \tag{14.7}
\end{equation*}
$$
[[Andersen_Ch14#^eq14-7|Eq. (14.7)]]

Expressing spot-measure drift in this explicit summed form makes clear why numerical freezing or predictor-corrector choices materially affect long-horizon Monte Carlo stability.

HJM-to-LMM volatility mapping is:

$$
\begin{equation*}
\sigma_{n}(t)=\tau_{n}^{-1}\left(1+\tau_{n} L_{n}(t)\right) \int_{T_{n}}^{T_{n+1}} \sigma_{f}(t, u) d u \tag{14.11}
\end{equation*}
$$
[[Andersen_Ch14#^eq14-11|Eq. (14.11)]]

This mapping also explains why discrete-tenor calibration does not uniquely pin down continuous-forward volatility: only interval integrals are observed through Libor dynamics.

Stochastic-volatility state driver is:

$$
\begin{equation*}
d z(t)=\theta\left(z_{0}-z(t)\right) d t+\eta \psi(z(t)) d Z(t), \quad z(0)=z_{0} \tag{14.15}
\end{equation*}
$$
[[Andersen_Ch14#^eq14-15|Eq. (14.15)]]

The positivity-preserving variance specification is essential for keeping simulated volatility levels economically interpretable across stressed and benign rate environments.

The frozen-weight swap-rate projection bridge is:

$$
\begin{align*}
w_{n}(t)=\frac{\varphi\left(L_{n}(t)\right)}{\varphi(S(t))} & \times \frac{\partial S(t)}{\partial L_{n}(t)}=\frac{\varphi\left(L_{n}(t)\right)}{\varphi(S(t))} \times \frac{S(t) \tau_{n}}{1+\tau_{n} L_{n}(t)} \\
& \times\left[\frac{P\left(t, T_{k}\right)}{P\left(t, T_{j}\right)-P\left(t, T_{k}\right)}+\frac{\sum_{i=n}^{k-1} \tau_{i} P\left(t, T_{i+1}\right)}{A(t)}\right] \tag{14.31}
\end{align*}
$$
[[Andersen_Ch14#^eq14-31|Eq. (14.31)]]

These weights provide the exact channel through which individual forward-rate shocks are aggregated into swap-rate diffusion under approximation, preserving tenor-structure intuition.

Martingale-preserving deflated-bond dynamics are:

$$
\begin{equation*}
\frac{d U\left(t, T_{n+1}\right)}{U\left(t, T_{n+1}\right)}=-\sum_{j=q(t)}^{n} \tau_{j} \frac{U\left(t, T_{j+1}\right)}{U\left(t, T_{j}\right)} \sigma_{j}(t)^{\top} d W^{B}(t), \quad n=q(t), \ldots, N-1 \tag{14.75}
\end{equation*}
$$
[[Andersen_Ch14#^eq14-75|Eq. (14.75)]]

#### Key Topical Extensions and Implementation Details

* **Introduction, Setup, and Setting (Sections 14.1 & 14.2)**: The motivation and historical notes for the LMM model arise from the need to overcome HJM's infinite-dimensional complexity and the log-normal rate explosion issue. The LMM is formulated on a fixed tenor structure $0 = T_0 < T_1 < \dots < T_N$. Unlike HJM which uses unobservable instantaneous rates, LMM models discretely compounded, directly observable Libor rates $L_n(t)$ spanning $[T_n, T_{n+1}]$, avoiding the Sandmann-Sondermann explosion issue. The dynamics can be parameterized using separable deterministic volatility functions or stochastic volatility.
* **Correlation (Section 14.3)**:
    * **Empirical PCA**: Principal Component Analysis (PCA) of empirical USD forward rate correlations shows that correlation decays exponentially as $|T_k - T_j|$ increases, with an asymptote that increases with $\min(T_k, T_j)$.
    * **Correlation Estimation and Smoothing**: To smooth noisy empirical correlations, parametric forms are fit. The triple-decaying form $q_2(x, y) = \rho_{\infty}(\min(x, y)) + (1 - \rho_{\infty}(\min(x, y))) \exp(-a(\min(x, y))|y - x|)$ (where the decay rate $a(z)$ and asymptote $\rho_{\infty}(z)$ are decaying exponential functions of tenor maturity) outperforms the four-parameter form $q_1(x, y)$, as demonstrated in the **Example: Fit to USD Data (Section 14.3.2.1)**.
    * **Negative Eigenvalues**: Parametric correlation formulas or rounding errors can generate matrices that are not positive definite. The matrix is repaired by finding the nearest positive semidefinite correlation matrix under a chosen matrix norm.
    * **Correlation PCA & Poor Man's PCA**: To reduce simulation dimension, the correlation matrix is approximated by a low-dimensional factor model ($m \ll N$) using the largest eigenvalues, illustrated in the **Example: USD Data (Section 14.3.4.1)**. A simple alternative is "poor man's PCA" which enforces a simple parametric factor structure.
* **Pricing of European Options (Section 14.4)**:
    * **Caplets & Swaptions**: Caplets are priced via Black's formula. Swaptions are priced using a projected swap-rate volatility approximation.
    * **Spread Options**: Pricing swaptions and CMS spread options requires modeling the term correlation between distinct forward rates. Spread options are priced under the appropriate annuity measure using joint log-normal approximations.
* **Calibration (Section 14.5)**:
    * **Basic Principles**: The basic principles of LMM calibration involve specifying the parametric form of the volatility norms, mapping them to correlation matrices, choosing target benchmark instruments, and optimizing the regularized objective function.
    * **Interpolation on the Whole Grid**: Flat volatility parameters defined on the calendar-tenor grid are interpolated to cover all simulation steps.
    * **Construction of Volatility Vectors**: Volatility vectors $\lambda_k(t)$ are constructed from their norms $\|\lambda_k(t)\|$ and the correlation matrix using Covariance PCA or Correlation PCA. A detailed **Discussion and Recommendation (Section 14.5.4.3)** analyzes the choices, preferring Correlation PCA due to its stability.
    * **Objective Function and Algorithms**: The calibration minimizes a regularized least-squares objective function penalizing deviations from time-homogeneity and smoothness in both calendar time and tenor directions.
    * **Speed-Up Through Sub-Problem Splitting**: Grouping swaption targets by option expiry allows splitting the high-dimensional optimization problem into a sequence of row-by-row sub-problems, significantly reducing computational cost.
    * **Correlation Calibration**: Correlation parameters can be calibrated to market-implied spread option prices to capture correlation smiles.
* **Monte Carlo Simulation & Discretization (Section 14.6)**:
    * **Euler-Type Schemes**: Under the spot measure $\mathrm{Q}^B$, rates are simulated using Euler or log-Euler discretization. Log-Euler prevents rates from becoming negative.
    * **Other Simulation Schemes**: To reduce discretization bias, Euler with Predictor-Corrector (which averages drift at the start and end of the step) or Lagging Predictor-Corrector schemes are preferred. High-order schemes and Brownian-Bridge simulation are also available.
    * **Martingale Discretization**: Simulating deflated bond prices $U(t, T_{n+1}) = P(t, T_{n+1})/B(t)$ instead of Libor rates directly preserves the martingale property of assets under the spot measure, preventing numerical drift. A review of **Comments and Alternatives (Section 14.6.3.2)** covers alternative schemes such as local correlation freezing.
    * **Variance Reduction**:
        * **Antithetic Sampling**: Simulates a second path by reflecting the Gaussian draws $Z \to -Z$.
        * **Control Variates**: Uses portfolios of vanilla zero-coupon bonds or caplets (whose continuous expected values are known analytically) to reduce variance.
        * **Importance Sampling**: Shifts the probability measure to ensure more paths survive (e.g., conditional on a barrier not being breached), which is crucial when pricing knock-out exotics like TARNs.

### Logical Bridge to Chapter 15

Chapter 15 extends the LMM framework with interpolation, advanced skew/smile approximation machinery, and multi-curve plus correlated stochastic-volatility generalizations.

---

### Chapter 15: Libor Market Model II

Chapter 15 takes the LMM framework into production-oriented extensions: stub interpolation, high-order swap-rate approximations, near-Markov reductions, and post-crisis multi-curve modeling. It starts with interpolation formulas that connect discrete-tenor models to off-grid maturities:

The chapter's motivation is implementation realism: front and back stubs, curve segmentation, and multi-curve consistency are not optional market details, they determine whether a model can be used on actual traded products.

$$
\begin{equation*}
P\left(t, T, T_{m}\right)=\frac{T-T_{m-1}}{T_{m}-T_{m-1}}+\frac{T_{m}-T}{T_{m}-T_{m-1}} P\left(t, T_{m-1}, T_{m}\right) \tag{15.1}
\end{equation*}
$$
[[Andersen_Ch15#^eq15-1|Eq. (15.1)]]

Generalized stochastic-volatility LMM dynamics then define the base state pair $(z,L_n)$ used throughout advanced approximation sections:

$$
\begin{align*}
d z(t) & =\theta\left(z_{0}-z(t)\right) d t+\eta(t) \psi(z(t)) d Z(t)  \tag{15.22}\\
d L_{n}(t) & =\sqrt{z(t)} \varphi_{n}\left(t, L_{n}(t)\right) \lambda_{n}(t)^{\top} d W^{T_{n+1}}(t), \quad n=1, \ldots, N-1
\end{align*}
$$
[[Andersen_Ch15#^eq15-23|Eq. (15.22)-(15.23)]]

This setup is the launch point for advanced approximation formulas. It keeps the familiar LMM structure while introducing enough flexibility to fit skew and smile simultaneously.

Projected swap-rate approximations are then upgraded to capture level, skew, and smile effects efficiently in calibration loops:

$$
\begin{equation*}
d S(t) \approx \sqrt{z(t)} p_{S}(t)\left(1+b_{S}(t)(S(t)-S(0))\right) d Y^{A}(t) \tag{15.28}
\end{equation*}
$$
[[Andersen_Ch15#^eq15-28|Eq. (15.28)]]

Calibration is formulated as a weighted multi-objective fitting problem over level/skew/smile target surfaces:

$$
\begin{align*}
\mathcal{I}\left(G_{\lambda}, G_{b}, G_{\eta}\right)= & \frac{w_{S, \lambda}}{N_{S}} \sum_{i=1}^{N_{S}}\left(\bar{\lambda}_{S_{i}}\left(G_{\lambda}, G_{b}, G_{\eta}\right)-\widehat{\lambda}_{S_{i}}\right)^{2}  \tag{15.32}\\
& +\frac{w_{S, b}}{N_{S}} \sum_{i=1}^{N_{S}}\left(\bar{b}_{S_{i}}\left(G_{\lambda}, G_{b}, G_{\eta}\right)-\widehat{b}_{S_{i}}\right)^{2} \\
& +\frac{w_{S, \eta}}{N_{S}} \sum_{i=1}^{N_{S}}\left(\bar{\eta}_{S_{i}}\left(G_{\lambda}, G_{b}, G_{\eta}\right)-\widehat{\eta}_{S_{i}}\right)^{2} \\
& +\ldots
\end{align*}
$$
[[Andersen_Ch15#^eq15-32|Eq. (15.32)]]

The explicit objective decomposition is important because it makes trade-offs visible: improving one region of the volatility surface can degrade another unless weights and parameter regularization are tuned carefully.

Multi-curve extensions are built by separating discounting and forwarding curve objects while retaining HJM consistency at the forward-rate layer:

$$
\begin{equation*}
d f(t, T)=\sigma_{f}(t, T)^{\top} \sigma_{P}(t, T) d t+\sigma_{f}(t, T)^{\top} d W(t) \tag{15.39}
\end{equation*}
$$
[[Andersen_Ch15#^eq15-39|Eq. (15.39)]]

This transition reflects post-crisis market structure directly. The model must now propagate curve basis effects through both pricing identities and stochastic dynamics, rather than treating basis as a static add-on.

$$
\begin{equation*}
\widetilde{P}(t, T)=P(t, T) / Z(t, T), \quad 0 \leq t \leq T \tag{15.45}
\end{equation*}
$$
[[Andersen_Ch15#^eq15-45|Eq. (15.45)]]

The ratio representation clarifies where basis risk lives in the model and how it enters valuation. It also provides a clean algebraic interface for integrating discounting and forwarding conventions.

Applications to forwarding Libors under spot measure then follow:

$$
\begin{equation*}
d \widetilde{L}_{n}(t)=\widetilde{\sigma}_{n}(t)^{\top}\left(\sum_{j=q(t)}^{n} \frac{\tau_{j} \sigma_{j}(t)}{1+\tau_{j} L_{j}(t)} d t+d W^{B}(t)\right) \tag{15.50}
\end{equation*}
$$
[[Andersen_Ch15#^eq15-50|Eq. (15.50)]]

Finally, non-zero-correlation stochastic-volatility dynamics add measure-adjusted mean-reversion effects and close the chapter's modern LMM extension stack:

$$
\begin{equation*}
d z(t)=\theta\left(z_{0}-\alpha(t) z(t)\right) d t+\eta \sqrt{z(t)} d Z^{n+1}(t) \tag{15.58}
\end{equation*}
$$
[[Andersen_Ch15#^eq15-58|Eq. (15.58)]]

Five bridge equations complete the production-oriented extension path in Chapter 15. First, local DVF linearization is made explicit by:

$$
\begin{align*}
\varphi_{n}(t, x) & \approx 1+b_{n}(t)\left(x-L_{n}(0)\right)  \tag{15.24}\\
b_{n}(t) & \triangleq \frac{\partial}{\partial x} \varphi_{n}\left(t, L_{n}(0)\right)
\end{align*}
$$
[[Andersen_Ch15#^eq15-24|Eq. (15.24)]]

This first-order expansion is what turns an abstract local-vol map into a calibration-ready pair of level and slope controls with direct interpretation.

Projected volatility/skew controls are then parameterized by:

$$
\begin{align*}
& p_{S}(t)=\left\|\lambda_{S}\left(t, \mathrm{E}^{A}(\mathbf{L}(t))\right)\right\|  \tag{15.29}\\
& b_{S}(t)=\frac{1}{p_{S}(t)} \sum_{n=j}^{k-1} \frac{\partial\left\|\lambda_{S}\left(t, \mathrm{E}^{A}(\mathbf{L}(t))\right)\right\|}{\partial L_{n}(t)} \frac{\int_{0}^{t}\left(\lambda_{n}(u)^{\top} \lambda_{S}(u, \mathbf{L}(0))\right) d u}{\int_{0}^{t}\left\|\lambda_{S}(u, \mathbf{L}(0))\right\|^{2} d u} \tag{15.30}
\end{align*}
$$
[[Andersen_Ch15#^eq15-29|Eq. (15.29)]]

Together these quantities isolate effective variance scale and local skew sensitivity, allowing optimization routines to target distinct smile features more transparently.

Multi-curve spread construction is anchored by:

$$
\begin{equation*}
\widetilde{P}(t, T)=P(t, T) e^{\int_{t}^{T} \varepsilon(t, u) d u}, \quad T \geq t \tag{15.40}
\end{equation*}
$$
[[Andersen_Ch15#^eq15-40|Eq. (15.40)]]

Writing forwarding discount factors this way makes basis propagation explicit in both valuation identities and downstream dynamic equations for floating-rate products.

The spread-ratio process dynamics are then:

$$
\begin{equation*}
d Z(t, T) / Z(t, T)=(r(t)-\widetilde{r}(t)) d t+\alpha_{Z}(t, T) d t-\sigma_{Z}(t, T)^{\top} d W(t) \tag{15.46}
\end{equation*}
$$
[[Andersen_Ch15#^eq15-46|Eq. (15.46)]]

The decomposition separates pure short-rate differential effects from additional spread-drift and spread-volatility terms, which is crucial when diagnosing multi-curve hedging behavior.

Finally, non-zero-correlation SV closure is given by the measure-adjusted drift multiplier:

$$
\begin{equation*}
\alpha(t)=1+\frac{\eta}{\theta} \rho(t)^{\top} \mu_{n}(t) \approx 1+\frac{\eta}{\theta} \rho(t)^{\top} \sum_{j=q(t)}^{n} \frac{\tau_{j} \varphi\left(L_{j}(0)\right) \lambda_{j}(t)}{1+\tau_{j} L_{j}(0)} . \tag{15.59}
\end{equation*}
$$
[[Andersen_Ch15#^eq15-59|Eq. (15.59)]]

#### Key Topical Extensions and Implementation Details

* **Advanced Approximations via Markovian Projection (Section 15.2)**: To improve on basic swaption pricing approximations (which assume constant swap rate skew), the high-dimensional stochastic-volatility LMM dynamics are projected onto a lower-dimensional Markovian system. This Markovian projection technique maps the swap rate's local volatility function to match the model's actual Libor rate skews, providing a displaced log-normal stochastic volatility approximation that is highly accurate for out-of-the-money options.
* **Swap Market Models (Section 15.5)**:
    * **Basic Ideas**: Swap market models (SMM) specify the SDEs of swap rates directly as martingales under their corresponding swap annuity measures. This contrasts with LMM where Libor rates are modeled. SMM is especially convenient for pricing Bermudan swaptions because the underlying swap rate is a primary state variable with a log-normal or local-volatility distribution, eliminating the need for complex rate-mapping approximations.
* **Multi-Curve QG Models (Section 15.7)**:
    * **Pricing Caplets and Swaptions**: In post-crisis multi-curve QG modeling, the discounting curve is separated from forwarding curves. Under this framework, vanilla options (caplets and swaptions) are priced by adjusting the state covariance matrix to incorporate curve basis spreads, allowing option pricing to remain analytical or highly tractable.

