---
tags: [quant/knowledge, quant/books, interest-rates]
---

[Master Index](../index.md)

# Index for Andersen & Piterbarg (Book 3)

### Chapter 16: Smile Dynamics and Hedging

Chapter 16 deals with the smile dynamics and hedging of interest rate swaptions. In particular, it stress-tests the impact of volatility dynamics on hedging efficiency, comparing the classical Black-Scholes framework with advanced stochastic volatility models.

Option pricing in these models requires accounting for model-consistent sensitivities rather than relying on standard Black model deltas and vegas. The model-consistent delta is defined by the following relation:

$$
\begin{equation*}
\Delta(t)=\frac{\partial c_{\mathrm{B}}}{\partial S}+\frac{\partial c_{\mathrm{B}}}{\partial \sigma_{\mathrm{B}}} \frac{\partial \sigma_{\mathrm{B}}}{\partial S}=\Delta_{\mathrm{B}}(t)+\Upsilon_{\mathrm{B}}(t) \frac{\partial \sigma_{\mathrm{B}}}{\partial S}, \tag{16.1}
\end{equation*}
$$
[[Andersen_Ch16#^eq16-1|Eq. (16.1)]]

This delta accounts for the model backbone, which defines how implied volatility changes when the underlying swap rate moves. Implied volatility curves derived from market data often exhibit a downward sloping smile (skew). If we hedge options assuming a flat backbone, the actual delta we use will deviate from the model-consistent delta, exposing the portfolio to unhedged directional risks. We can calibrate this smile using strike-specific volatilities as shown in:

$$
\begin{equation*}
c(t, S(t) ; T, K)=c_{\mathrm{B}}(t, S(t) ; T, K ; \psi(T, K)) \tag{16.2}
\end{equation*}
$$
[[Andersen_Ch16#^eq16-2|Eq. (16.2)]]

To capture this smile dynamics leakage, the chapter models the evolution of the swaption smile by specifying local volatility structures for swap rates. This allows calibrating the model to the swaption volatility matrix at time 0. The smile dynamics can be written as:

$$
\begin{equation*}
d S(t)=\lambda(b S(t)+(1-b) L) d W(t) \tag{16.3}
\end{equation*}
$$
[[Andersen_Ch16#^eq16-3|Eq. (16.3)]]

The model can be parameterized using a local volatility function of the displaced log-normal type, which is capable of capturing the skew of the volatility surface. Specifically, we define the local volatility of the short rate as:

$$
\begin{equation*}
d S(t)=\lambda(b S(t)+(1-b) S(0)) d W(t) \tag{16.4}
\end{equation*}
$$
[[Andersen_Ch16#^eq16-4|Eq. (16.4)]]

To calibrate the model to market swaptions, the implied volatility is parameterized as a function of the underlying rate and the backbone exponent. This backbone relation is represented by:

$$
\begin{equation*}
d S(t)=\lambda(b S(t)+(m-b) S(0)+(1-m) L) d W(t) \tag{16.5}
\end{equation*}
$$
[[Andersen_Ch16#^eq16-5|Eq. (16.5)]]

Alternatively, a more flexible backbone shape can be achieved by specifying the implied volatility curve with a time-dependent adjustment function, which is given by:

$$
\begin{equation*}
\left.\frac{\partial \sigma_{\mathrm{B}}(0, S(0) ; T, K)}{\partial K}\right|_{K=S(0)} \approx-\frac{1-b}{2} \frac{\lambda}{S(0)} . \tag{16.6}
\end{equation*}
$$
[[Andersen_Ch16#^eq16-6|Eq. (16.6)]]

To calibrate the model to a swaption strip, the local volatility of the swap rate is parameterized under the annuity measure. This swap rate local volatility function is written as:

$$
\begin{equation*}
\left.\frac{\partial \sigma_{\mathrm{B}}(0, S ; T, S)}{\partial S}\right|_{S=S(0)} \approx-(1-m) \frac{\lambda}{S(0)} \tag{16.7}
\end{equation*}
$$
[[Andersen_Ch16#^eq16-7|Eq. (16.7)]]

The SDE for the swap rate under the annuity measure is driven by a multi-dimensional Brownian motion, where the volatility is local. The swap rate SDE is written as:

$$
\begin{equation*}
d z_{n, m}(t)=\theta\left(1-z_{n, m}(t)\right) d t+\eta_{n, m} \sqrt{z_{n . m}(t)} d Z^{n, m}(t) \tag{16.8}
\end{equation*}
$$
[[Andersen_Ch16#^eq16-8|Eq. (16.8)]]

The swaption price is computed using the swap rate volatility function. The swaption price under the annuity measure is given by:

$$
\begin{equation*}
\mathcal{I}_{1}(\lambda, b, \eta)=\sum_{j=1}^{J} w_{j}\left(V\left(K_{j} ; \lambda, b, \eta\right)-\widehat{V}_{j}\right)^{2} \tag{16.10}
\end{equation*}
$$
[[Andersen_Ch16#^eq16-10|Eq. (16.10)]]

To value European swaptions efficiently, we can use the log-normal swaption pricing formula under the annuity measure, which is written as:

$$
\begin{equation*}
\mathcal{I}_{2}(\lambda, b, \eta)=\sum_{j=1}^{J} w_{j}\left(\frac{V\left(K_{j} ; \lambda, b, \eta\right)-\widehat{V}_{j}}{\widehat{\Upsilon}_{j}}\right)^{2} \tag{16.11}
\end{equation*}
$$
[[Andersen_Ch16#^eq16-11|Eq. (16.11)]]

Alternatively, the normal swaption pricing formula is used when swap rates can become negative, which is given by:

$$
\begin{equation*}
\mathcal{I}_{3}(\lambda, b, \eta)=\sum_{j=1}^{J} w_{j}\left(\sigma_{\mathrm{B}}\left(K_{j} ; \lambda, b, \eta\right)-\widehat{\sigma}_{j}\right)^{2} \tag{16.12}
\end{equation*}
$$
[[Andersen_Ch16#^eq16-12|Eq. (16.12)]]

Stochastic volatility models are introduced to capture the dynamics of the smile more accurately, where the variance process follows a square-root diffusion. The swap rate SDE under stochastic volatility is written as:

$$
\begin{equation*}
\sigma_{\mathrm{B}}(0, S(0) ; T, K)=a+b\left(\rho(k-h)+\sqrt{(k-h)^{2}+s^{2}}\right) \tag{16.13}
\end{equation*}
$$
[[Andersen_Ch16#^eq16-13|Eq. (16.13)]]

The variance process itself follows a mean-reverting square-root process, which is written as:

$$
\begin{align*}
& \mathcal{I}_{\text {simooth }}(\mathcal{X})=w^{\lambda} \sum_{n=2}^{N-1}\left(\lambda_{n}(\mathcal{X})-\lambda_{n-1}(\mathcal{X})\right)^{2} \\
& +w^{b} \sum_{n=2}^{N-1}\left(b_{n}(\mathcal{X})-b_{n-1}(\mathcal{X})\right)^{2}+w^{\eta} \sum_{n=2}^{N-1}\left(\eta_{n}(\mathcal{X})-\eta_{n-1}(\mathcal{X})\right)^{2} \tag{16.14}
\end{align*}
$$
[[Andersen_Ch16#^eq16-14|Eq. (16.14)]]

The joint calibration to caps and swaptions requires low-dimensional models that can be solved efficiently. The local model method projects high-dimensional Libor market dynamics onto low-dimensional Markovian models to value path-dependent exotics. This is essential for pricing callable exotics, which are sensitive to both volatility smiles and terminal swap correlations. The caplet price under the local model is written as:

$$
\begin{align*}
& A(t) \triangleq A_{0, N}(t)=\sum_{n=0}^{N-1} \tau_{n} P\left(t, T_{n+1}\right)  \tag{16.15}\\
& S(t) \triangleq S_{0, N}(t)=\frac{P(t, T)-P\left(t, T_{N}\right)}{A(t)} \tag{16.16}
\end{align*}
$$
[[Andersen_Ch16#^eq16-15|Eq. (16.15)]]

The joint calibration to caps and swaptions is performed by minimizing a regularized objective function that penalizes deviations from time-homogeneity and smoothness. The joint objective function is written as:

$$
\begin{equation*}
V(0)=A(0) \mathrm{E}^{A}\left(\frac{X}{A(T)}\right) . \tag{16.17}
\end{equation*}
$$
[[Andersen_Ch16#^eq16-17|Eq. (16.17)]]

The swaption grid calibration is implemented using a row-by-row splitting algorithm, which significantly reduces the computational cost of the optimization. This swaption calibration relation is given by:

$$
\begin{equation*}
P(T, M)=\pi(S(T), M), \quad M \geq T \tag{16.18}
\end{equation*}
$$
[[Andersen_Ch16#^eq16-18|Eq. (16.18)]]

Finally, the chapter reviews alternative interpolation rules for the swaption matrix, which are used to construct a smooth volatility surface from sparse market quotes. The interpolation function is written as:

$$
\begin{equation*}
P(0, M)=A(0) \mathrm{E}^{A}\left(\frac{\pi(S(T), M)}{\sum_{n=0}^{N-1} \tau_{n} \pi\left(S(T), T_{n+1}\right)}\right), \quad M \geq T . \tag{16.19}
\end{equation*}
$$
[[Andersen_Ch16#^eq16-19|Eq. (16.19)]]

Chapter 17 moves from single-rate options to multi-rate vanilla derivatives, introducing joint distributions of multiple rates and copula methods under a single reference pricing measure.

In the next chapter, we will see how these multi-rate options are formulated and how their joint dependencies are modeled.

The shift from a single annuity measure to a joint pricing measure is a key step in this development.

This allows us to handle more complex exotic payouts that depend on multiple swap rates simultaneously.

We will explore how the correlation structure between swap rates affects the pricing of CMS spread options.

---

### Chapter 17: Multi-Rate Vanilla Derivatives

Chapter 17 expands the vanilla option framework to multi-rate vanilla derivatives, including spread options, CMS options, and mid-curve swaptions. In particular, it is important to analyze multi-rate exotics under a single reference measure.

Unlike standard options that depend on a single swap rate, multi-rate options depend on the joint distribution of multiple rates, requiring a rigorous modeling of the terminal correlation structure under a single reference pricing measure. The SDE of the forward Libor rates under the spot measure is written as:

$$
\begin{equation*}
V\left(T_{p}\right)=f\left(S_{1}\left(t_{1}\right), \ldots, S_{d}\left(t_{d}\right)\right) \text { paid at time } T_{p} . \tag{17.1}
\end{equation*}
$$
[[Andersen_Ch17#^eq17-1|Eq. (17.1)]]

The marginal distribution of each rate is assumed to be log-normal under its own forward measure, allowing options to be priced using Black's formula. This marginal pricing relation is written as:

$$
\begin{equation*}
V(0)=\mathrm{E}\left(\beta\left(T_{p}\right)^{-1} f\left(S_{1}\left(t_{1}\right), \ldots, S_{d}\left(t_{d}\right)\right)\right) \tag{17.2}
\end{equation*}
$$
[[Andersen_Ch17#^eq17-2|Eq. (17.2)]]

To describe the joint distribution of rates, the chapter introduces copula methods, which separate marginal distributions from the joint dependency structure. The joint Gaussian distribution of the rates is defined by:

$$
\begin{equation*}
V(0)=P\left(0, T_{p}\right) \mathrm{E}^{T_{r}}\left(f\left(S_{1}\left(t_{1}\right), \ldots, S_{d}\left(t_{d}\right)\right)\right) . \tag{17.3}
\end{equation*}
$$
[[Andersen_Ch17#^eq17-3|Eq. (17.3)]]

The Gaussian copula method projects marginal rates onto Gaussian variables, using a correlation matrix to represent their joint dependency. This state variable SDE is written as:

$$
\begin{equation*}
\mathrm{Q}^{T_{p}}\left(S_{i}\left(t_{i}\right) \in d s\right)=\frac{A_{i}(0)}{P\left(0, T_{p}\right)} \alpha_{i}(s) \mathrm{Q}^{A}\left(S_{i}\left(t_{i}\right) \in d s\right) \tag{17.4}
\end{equation*}
$$
[[Andersen_Ch17#^eq17-4|Eq. (17.4)]]

The discount bond price in the Gaussian framework is recovered as an exponential-quadratic function of the state vector. This bond reconstitution formula is given by:

$$
\begin{equation*}
\alpha_{i}(s)=\mathrm{E}^{A_{\imath}}\left(\left.\frac{P\left(t_{i}, T_{p}\right)}{A_{i}\left(t_{i}\right)} \right\rvert\, S_{i}\left(t_{i}\right)=s\right), \quad i=1, \ldots, d \tag{17.5}
\end{equation*}
$$
[[Andersen_Ch17#^eq17-5|Eq. (17.5)]]

The forward swap rate is approximated as a weighted sum of the underlying forward Libor rates, where the weights are annuity-based. This swap rate approximation is written as:

$$
\begin{equation*}
X_{\imath}=\Psi_{\imath}^{-1}\left(\Phi\left(Z_{\imath}\right)\right), \quad i=1, \ldots, d \tag{17.6}
\end{equation*}
$$
[[Andersen_Ch17#^eq17-6|Eq. (17.6)]]

The valuation of multi-rate derivatives requires evaluating the expected value of the option payoff under the reference pricing measure. This multi-rate pricing equation is written as:

$$
\begin{equation*}
\Psi_{\text {gauss }}\left(x_{1}, \ldots, x_{d}\right)=\Phi_{d}\left(\Phi^{-1}\left(\Psi_{1}\left(x_{1}\right)\right), \ldots, \Phi^{-1}\left(\Psi_{d}\left(x_{d}\right)\right) ; R\right) . \tag{17.7}
\end{equation*}
$$
[[Andersen_Ch17#^eq17-7|Eq. (17.7)]]

CMS spread options, which pay the difference between two swap rates, are priced under the appropriate annuity measure. The CMS spread option price is given by:

$$
\begin{align*}
& \psi_{\text {gauss }}\left(x_{1}, \ldots, x_{d}\right)  \tag{17.8}\\
& \quad=\frac{\partial^{d}}{\partial x_{1} \ldots \partial x_{d}} \Psi_{\text {gauss }}\left(x_{1}, \ldots, x_{d}\right)  \tag{17.9}\\
& \quad=\left.\frac{\partial^{d}}{\partial z_{1} \ldots \partial z_{d}} \Phi_{d}\left(z_{1}, \ldots, z_{d} ; R\right)\right|_{z^{\prime}=\Phi^{-1}\left(\Psi_{i}(x,)\right), \forall i} \times \prod_{\imath=1}^{d} \frac{\Psi_{i}^{\prime}\left(x_{i}\right)}{\Phi^{\prime}\left(\Phi^{-1}\left(\Psi_{i}\left(x_{i}\right)\right)\right)} \\
& \quad=\phi_{d}\left(\Phi^{-1}\left(\Psi_{1}\left(x_{1}\right)\right), \ldots, \Phi^{-1}\left(\Psi_{d}\left(x_{d}\right)\right) ; R\right) \times \prod_{i=1}^{d} \frac{\psi_{\imath}\left(x_{\imath}\right)}{\phi\left(\Phi^{-1}\left(\Psi_{i}\left(x_{i}\right)\right)\right)}
\end{align*}
$$
[[Andersen_Ch17#^eq17-8|Eq. (17.8)]]

The terminal correlation between the swap rates is calibrated to historical correlation data or market-implied spread option prices. This correlation calibration formula is written as:

$$
\begin{equation*}
V=\int \cdots \int f\left(x_{1}, \ldots, x_{d}\right) \psi_{\text {gauss }}\left(x_{1}, \ldots, x_{d}\right) d x_{1} \ldots d x_{d} \tag{17.10}
\end{equation*}
$$
[[Andersen_Ch17#^eq17-10|Eq. (17.10)]]

To capture the tail dependence of the rates, Archimedean copulas are introduced, which are defined through a generator function. This Archimedean copula density is given by:

$$
\begin{equation*}
V=\int \cdots \int f\left(\Psi_{1}^{-1}\left(\Phi\left(z_{1}\right)\right), \ldots, \Psi_{d}^{-1}\left(\Phi\left(z_{d}\right)\right)\right) \phi_{d}\left(z_{1}, \ldots, z_{d} ; R\right) d z_{1} \ldots d z_{d} \tag{17.11}
\end{equation*}
$$
[[Andersen_Ch17#^eq17-11|Eq. (17.11)]]

The joint density can also be constructed using a mixture of copulas, which allows capturing both upper and lower tail dependence. This mixture copula density is written as:

$$
\begin{equation*}
\Psi_{\text {gauss }}\left(x_{1}, \ldots, x_{d}\right)=C_{\text {gauss }}\left(\Psi_{1}\left(x_{1}\right), \ldots, \Psi_{d}\left(x_{d}\right) ; R\right) \tag{17.12}
\end{equation*}
$$
[[Andersen_Ch17#^eq17-12|Eq. (17.12)]]

The valuation of mid-curve swaptions requires modeling the joint distribution of the forward swap rate and the discount bond price. This mid-curve swaption pricing formula is written as:

$$
\begin{gather*}
C\left(u_{1}, u_{2}, \ldots, u_{i-1}, 0, u_{i+1}, \ldots, u_{d}\right)=0,  \tag{17.13}\\
C\left(1,1, \ldots, 1, u_{\imath}, 1,1, \ldots, 1\right)=u_{i} . \tag{17.14}
\end{gather*}
$$
[[Andersen_Ch17#^eq17-13|Eq. (17.13)]]

Finally, the chapter reviews numerical integration methods and copula approximations to calculate the multi-rate option prices efficiently. This multi-rate integration formula is given by:

$$
\begin{equation*}
C_{\mathrm{ID}}=\prod_{i=1}^{d} u_{i} . \tag{17.15}
\end{equation*}
$$
[[Andersen_Ch17#^eq17-15|Eq. (17.15)]]

Chapter 18 establishes the valuation theory and numerical methods for callable Libor exotics, introducing optimal stopping problems and regression-based continuation value approximations.

We will see how continuation values are approximated using regression onto simulated state paths.

This provides a highly tractable framework to estimate lower bounds on exotic option prices.

---

### Chapter 18: Callable Libor Exotics

Chapter 18 establishes the valuation theory and numerical methods for callable Libor exotics (CLEs), such as target redemption notes and Bermudan swaptions. In particular, we formulate the valuation under backward induction.

CLEs are path-dependent exotics that can be cancelled or exercised at the option of the holder or issuer at multiple pre-specified exercise dates. The value of a callable exotic under the spot measure is written as:

$$
\begin{equation*}
X_{n}=\tau_{n} \times\left(C_{n}-L_{n}\left(T_{n}\right)\right) \tag{18.1}
\end{equation*}
$$
[[Andersen_Ch18#^eq18-1|Eq. (18.1)]]

Valuing callable exotics requires solving an optimal stopping problem, which can be formulated as a backward dynamic programming recursion. At each exercise date, the holder compares the immediate exercise value with the continuation value, which represents the expected value of holding the option. This exercise decision rule is given by:

$$
\begin{equation*}
U_{n}(t)=B(t) \sum_{i=n}^{N-1} \mathrm{E}_{t}\left(B\left(T_{i+1}\right)^{-1} X_{i}\right) . \tag{18.2}
\end{equation*}
$$
[[Andersen_Ch18#^eq18-2|Eq. (18.2)]]

The continuation value is defined as the expected value of the next-period option value discounted to the current period. This continuation value recursion is written as:

$$
\begin{equation*}
H_{0}(t) \geq H_{1}(t) \geq \ldots \geq H_{N-2}(t) \tag{18.3}
\end{equation*}
$$
[[Andersen_Ch18#^eq18-3|Eq. (18.3)]]

The optimal stopping time is defined as the first exercise date where the exercise value exceeds the continuation value. This optimal stopping time relation is written as:

$$
\begin{equation*}
H_{n}\left(T_{n}\right)=B\left(T_{n}\right) \sup _{\xi \in \mathcal{T}_{n}} \mathrm{E}_{T_{n}}\left(B\left(T_{\xi}\right)^{-1} U_{\xi}\left(T_{\xi}\right)\right), \quad n=0, \ldots, N-1 \tag{18.4}
\end{equation*}
$$
[[Andersen_Ch18#^eq18-4|Eq. (18.4)]]

To price CLEs in high-dimensional Monte Carlo simulation, the chapter utilizes regression methods, such as the Longstaff-Schwartz algorithm. This method approximates the continuation value by regressing the future option payoffs onto a set of basis functions of the current state variables. This regression approximation is written as:

$$
\begin{equation*}
\eta_{n}(\omega)=\min \left\{k>n: U_{k}\left(T_{k}\right) \geq H_{k}\left(T_{k}\right)\right\} \wedge N, \tag{18.5}
\end{equation*}
$$
[[Andersen_Ch18#^eq18-5|Eq. (18.5)]]

To assess the accuracy of the regression approximation, the chapter develops dual upper bound algorithms, which construct an upper bound by introducing a martingale that penalizes suboptimal exercise decisions. This dual upper bound estimator is given by:

$$
\begin{equation*}
H_{0}(0)=\mathrm{E}\left(B\left(T_{\eta}\right)^{-1} U_{\eta}\left(T_{\eta}\right)\right)=\mathrm{E}\left(\sum_{n=\eta}^{N-1} B\left(T_{n+1}\right)^{-1} X_{n}\right) \tag{18.6}
\end{equation*}
$$
[[Andersen_Ch18#^eq18-6|Eq. (18.6)]]

The martingale is constructed by taking the pathwise derivative of the option value or using the hedging portfolio. This martingale SDE is written as:

$$
\begin{equation*}
H_{n-1}\left(T_{n-1}\right)=B\left(T_{n-1}\right) \mathrm{E}_{T_{n-1}}\left(B\left(T_{n}\right)^{-1} \max \left(H_{n}\left(T_{n}\right), U_{n}\left(T_{n}\right)\right)\right) \tag{18.7}
\end{equation*}
$$
[[Andersen_Ch18#^eq18-7|Eq. (18.7)]]

Evaluating both the regression lower bound and the dual upper bound allows establishing a tight confidence interval for the exotic option price. This dual pricing inequality is written as:

$$
\begin{equation*}
H_{0}(0)=\sum_{n=1}^{N-1} \mathrm{E}\left(B\left(T_{n}\right)^{-1}\left(U_{n}\left(T_{n}\right)-H_{n}\left(T_{n}\right)\right)^{+}\right) \tag{18.8}
\end{equation*}
$$
[[Andersen_Ch18#^eq18-8|Eq. (18.8)]]

The regression variable choice is critical for the convergence of the algorithm, and the chapter recommends using forward swap rates as explanatory variables. This regression variable relation is written as:

$$
\begin{equation*}
\left\|\left(X\left(\omega_{1}\right)-\zeta\left(T, \omega_{1}\right)^{\top} \beta, \ldots, X\left(\omega_{K}\right)-\zeta\left(T, \omega_{K}\right)^{\top} \beta\right)\right\|^{2} \rightarrow \min \tag{18.9}
\end{equation*}
$$
[[Andersen_Ch18#^eq18-9|Eq. (18.9)]]

The regression coefficients are solved using least-squares optimization on the simulated paths. This least-squares coefficient estimator is given by:

$$
\begin{equation*}
\mathcal{R}_{T}(X)=\mathcal{C}\left(\mathcal{R}_{T}(X)\right)^{\top} \zeta(T) . \tag{18.10}
\end{equation*}
$$
[[Andersen_Ch18#^eq18-10|Eq. (18.10)]]

The optimal stopping policy is then used to simulate the option value on a separate set of paths to avoid look-ahead bias. This unbiased lower bound estimator is written as:

$$
\begin{equation*}
\widetilde{H}_{n-1}\left(T_{n-1}\right)=\mathcal{R}_{T_{n-1}}\left(\frac{B\left(T_{n-1}\right)}{B\left(T_{n}\right)} \max \left(\widetilde{H}_{n}\left(T_{n}\right), U_{n}\left(T_{n}\right)\right)\right), \tag{18.11}
\end{equation*}
$$
[[Andersen_Ch18#^eq18-11|Eq. (18.11)]]

The dual upper bound can be further optimized by minimizing the variance of the martingale penalty. This martingale optimization formula is given by:

$$
\begin{equation*}
V_{n}\left(\omega_{k}\right)=\frac{B\left(T_{n-1}, \omega_{k}\right)}{B\left(T_{n}, \omega_{k}\right)} \max \left(\widetilde{H}_{n}\left(T_{n}, \omega_{k}\right), U_{n}\left(T_{n}, \omega_{k}\right)\right) \tag{18.12}
\end{equation*}
$$
[[Andersen_Ch18#^eq18-12|Eq. (18.12)]]

Finally, the chapter reviews numerical examples of pricing TARNs and snowballs using the regression and dual algorithms. This numerical pricing formula is written as:

$$
\begin{equation*}
V_{n}\left(\omega_{k}\right)=\frac{B\left(T_{n-1}, \omega_{k}\right)}{B\left(T_{n}, \omega_{k}\right)} \max \left(\tilde{H}_{n}\left(T_{n}, \omega_{k}\right), \tilde{U}_{n}\left(T_{n}, \omega_{k}\right)\right), \quad k=1, \ldots, K, \tag{18.13}
\end{equation*}
$$
[[Andersen_Ch18#^eq18-13|Eq. (18.13)]]

Chapter 19 focuses on Bermudan swaptions, developing local projection methods to project Libor market dynamics onto low-dimensional Markovian models for fast PDE pricing.

We will analyze the calibration of local volatilities to European swaption smiles.

This enables consistent valuation of standard and non-standard exercise structures.

---

### Chapter 19: Bermudan Swaptions

Chapter 19 focuses on Bermudan swaptions, which are swaptions that can be exercised at multiple pre-specified exercise dates. Bermudan swaptions are the most common interest rate exotics and are used by financial institutions to hedge interest rate risk.

Their valuation is highly sensitive to the term structure of interest rates, forward volatility, and swap rate correlations. The value of a Bermudan swaption under the annuity measure is written as:

$$
\begin{equation*}
U_{n}(t)=\sum_{i=n}^{N-1} \tau_{i} P\left(t, T_{i+1}\right)\left(L_{\imath}(t)-k\right) \tag{19.1}
\end{equation*}
$$
[[Andersen_Ch19#^eq19-1|Eq. (19.1)]]

To value Bermudan swaptions in a Libor market model, the chapter presents local projection methods, which project high-dimensional Libor market dynamics onto a low-dimensional Markovian model. This swap rate mapping is represented by:

$$
\begin{equation*}
V_{\text {Berm }}(0 ; k) \geq \max _{n=1, \ldots, N-1} V_{\text {swaption }, n}(0 ; k), \tag{19.2}
\end{equation*}
$$
[[Andersen_Ch19#^eq19-2|Eq. (19.2)]]

The SDE for the forward swap rate under its own swap measure is driven by a volatility norm and a Brownian motion, which is written as:

$$
\begin{equation*}
\tilde{U}_{n}(t)=\sum_{i=n}^{N-1} R_{i} \tau_{i} P\left(t, T_{i+1}\right)\left(L_{i}(t)-k\right) \tag{19.3}
\end{equation*}
$$
[[Andersen_Ch19#^eq19-3|Eq. (19.3)]]

The terminal correlation between distinct forward swap rates is modeled using a parametric correlation function. This swap rate correlation formula is written as:

$$
\begin{equation*}
v_{n, m}=R_{n+m-1}-R_{n+m}, \quad m=1, \ldots, N-n, \tag{19.4}
\end{equation*}
$$
[[Andersen_Ch19#^eq19-4|Eq. (19.4)]]

The optimal stopping problem for the Bermudan swaption is solved using backward induction on the finite difference grid. This Bermudan dynamic programming recursion is given by:

$$
\begin{equation*}
\widetilde{S}_{n}\left(T_{n}\right)=\sum_{m=1}^{N-n} w_{n, m}\left(T_{n}\right) S_{n, m}\left(T_{n}\right) \tag{19.5}
\end{equation*}
$$
[[Andersen_Ch19#^eq19-5|Eq. (19.5)]]

To calibrate the model to the Bermudan swaption smile, the local volatility parameters are fit to a strip of co-terminal European swaptions. This calibration objective function is written as:

$$
\begin{equation*}
\widetilde{S}_{n}\left(T_{n}\right) \approx \sum_{m=1}^{N-n} w_{n, m}(0) S_{n, m}\left(T_{n}\right), \quad w_{n, m}(0)=\left(R_{n+m-1}-R_{n+m}\right) \frac{A_{n, m}(0)}{\widetilde{A}_{n}(0)} . \tag{19.6}
\end{equation*}
$$
[[Andersen_Ch19#^eq19-6|Eq. (19.6)]]

The exercise boundary of the Bermudan swaption is represented as a function of the underlying swap rate. This exercise boundary relation is given by:

$$
\begin{equation*}
\mu_{n}=\underset{m}{\operatorname{argmin}}\left\{\left|R_{n} A_{n, m}(0)-\widetilde{A}_{n}(0)\right|\right\} \tag{19.7}
\end{equation*}
$$
[[Andersen_Ch19#^eq19-7|Eq. (19.7)]]

The mapping function approach maps the high-dimensional Libor rates to a single Markovian state process. This state process mapping is written as:

$$
\begin{align*}
V\left(x_{0} ; R, q, m\right) & =\widetilde{U}_{n}\left(x_{0}\right)  \tag{19.8}\\
\frac{\partial}{\partial x} V\left(x_{0} ; R, q, m\right) & =\frac{\partial}{\partial x} \widetilde{U}_{n}\left(x_{0}\right)  \tag{19.9}\\
\frac{\partial^{2}}{\partial x^{2}} V\left(x_{0} ; R, q, m\right) & =\frac{\partial^{2}}{\partial x^{2}} \widetilde{U}_{n}\left(x_{0}\right) \tag{19.10}
\end{align*}
$$
[[Andersen_Ch19#^eq19-8|Eq. (19.8)]]

The drift of the Markovian state variable is adjusted to preserve the martingale property of assets under the spot measure. This drift adjustment formula is written as:

$$
\begin{equation*}
R_{i}=\prod_{j=0}^{i-1}\left(1+\tau_{j} k\right), \quad i=1, \ldots, N-1 \tag{19.11}
\end{equation*}
$$
[[Andersen_Ch19#^eq19-11|Eq. (19.11)]]

The local volatility function of the state process is parameterized to match the caplet volatility strip. This local volatility parameterization is given by:

$$
\begin{align*}
\tilde{U}_{n}\left(T_{n}\right) & =\sum_{i=n}^{N-1} R_{i} \tau_{i} P\left(T_{n}, T_{i+1}\right)\left(k-L_{i}\left(T_{n}\right)\right) \\
& =\sum_{i=n}^{N-1}\left(\prod_{j=0}^{i-1}\left(1+\tau_{j} k\right)\right)\left(P\left(T_{n}, T_{i+1}\right)\left(1+\tau_{i} k\right)-P\left(T_{n}, T_{i}\right)\right) \\
& =\sum_{i=n}^{N-1}\left(R_{i+1} P\left(T_{n}, T_{i+1}\right)-R_{i} P\left(T_{n}, T_{i}\right)\right) \\
& =R_{N} P\left(T_{n}, T_{N}\right)-R_{n} \tag{19.12}
\end{align*}
$$
[[Andersen_Ch19#^eq19-12|Eq. (19.12)]]

The continuation value on the grid is solved using Crank-Nicolson finite difference scheme. This finite difference recursion is written as:

$$
\begin{equation*}
U_{n}^{A}\left(T_{n}+\right) \neq U_{n}^{A}\left(T_{n}\right)=U_{n}\left(T_{n}\right) . \tag{19.13}
\end{equation*}
$$
[[Andersen_Ch19#^eq19-13|Eq. (19.13)]]

Non-standard Bermudan swaptions, such as amortizing swaptions, are priced by mapping the non-standard swap rates to standard swap rates. This non-standard swap rate mapping is written as:

$$
\begin{align*}
u_{n}^{B}\left(t_{m}\right)=L\left(t_{m}, t_{m}, T_{n+1}\right) & \left(T_{n+1}-t_{m}\right) P\left(t_{m}, T_{n+1}\right) \\
& \quad-k \sum_{i=m}^{M-1}\left(t_{i+1}-t_{i}\right) P\left(t_{m}, t_{i+1}\right) \tag{19.14}
\end{align*}
$$
[[Andersen_Ch19#^eq19-14|Eq. (19.14)]]

Finally, the chapter reviews Monte Carlo regression methods and parametric boundary methods to price Bermudan swaptions in high dimensions. This parametric boundary estimator is given by:

$$
\begin{align*}
u_{n}^{A}\left(t_{m}\right)=L\left(T_{n}, T_{n}, T_{n+1}\right)\left(T_{n+1}\right. & \left.-t_{m}\right) P\left(t_{m}, T_{n+1}\right) \\
& -k\left(T_{n+1}-t_{m}\right) P\left(t_{m}, T_{n+1}\right) \tag{19.15}
\end{align*}
$$
[[Andersen_Ch19#^eq19-15|Eq. (19.15)]]

Chapter 20 covers TARNs and volatility swaps, analyzing their valuation under global and local projection models using PDE and Monte Carlo schemes.

We will see how target redemption notes accumulate inverse floating coupons.

This triggers early termination when the cumulative sum reaches a barrier.

---

### Chapter 20: TARNs, Volatility Swaps, and Other Derivatives

Chapter 20 covers target redemption notes (TARNs), volatility swaps, and other path-dependent derivatives. In particular, we model volatility structures under local projection.

TARNs are swap contracts that pay a sequence of inverse floating coupons subject to a cap and floor, which terminate automatically once the accumulated coupons reach a target gain. The value of a TARN under the spot measure is written as:

$$
\begin{equation*}
C_{n}=\left(s-g \times L_{n}\left(T_{n}\right)\right)^{+} \tag{20.1}
\end{equation*}
$$
[[Andersen_Ch20#^eq20-1|Eq. (20.1)]]

The coupon payoff is a non-linear function of the underlying swap rate, subject to a cap, floor, and target barrier. This coupon payoff formula is written as:

$$
\begin{align*}
V_{\operatorname{tarn}}(0) & =\mathrm{E}\left(\sum_{n=1}^{N-1} B\left(T_{n+1}\right)^{-1} \tau_{n}\left(C_{n}-L_{n}\left(T_{n}\right)\right) 1_{\left\{Q_{. .}<R\right\}}\right)  \tag{20.2}\\
Q_{n} & =\sum_{i=1}^{n-1} \tau_{i} C_{i}, \quad Q_{1}=0
\end{align*}
$$
[[Andersen_Ch20#^eq20-2|Eq. (20.2)]]

To capture volatility smile effects on TARNs, the chapter uses a local projection method to project the local volatility surface of swap rates onto a low-dimensional model. This swap rate local volatility SDE is written as:

$$
\begin{array}{r}
V_{\operatorname{tarn}}(0)=\mathrm{E}\left(\sum_{n=1}^{N-1} B\left(T_{n+1}\right)^{-1} \tau_{n}\left(\left(s-g \times L_{n}\left(T_{n}\right)\right)^{+}-L_{n}\left(T_{n}\right)\right)\right. \\
\left.\times 1_{\left\{\sum_{i=1}^{n-1} \tau_{i}\left(s-g \times L_{i}\left(T_{i}\right)\right)^{+}<R\right\}}\right) \tag{20.3}
\end{array}
$$
[[Andersen_Ch20#^eq20-3|Eq. (20.3)]]

Volatility swaps, which pay the difference between the realized volatility of a swap rate and a fixed strike, are priced under the spot measure. This volatility swap pricing formula is written as:

$$
\begin{equation*}
\tau_{1}\left(s-g b_{1}\right)^{+}=R \tag{20.4}
\end{equation*}
$$
[[Andersen_Ch20#^eq20-4|Eq. (20.4)]]

The realized volatility is calculated as the square root of the sum of squared log returns of the underlying rate over a discrete grid. This realized volatility formula is given by:

$$
\begin{equation*}
V\left(T_{N}-, I\right)=0 \tag{20.5}
\end{equation*}
$$
[[Andersen_Ch20#^eq20-5|Eq. (20.5)]]

The underlying rate process is modeled using local stochastic volatility, which allows capturing the dynamics of the smile. This local volatility SDE is written as:

$$
\begin{equation*}
V\left(T_{n}-, I\right)=V\left(T_{n}-, I\right)+P\left(T_{n}, T_{n+1}\right) \tau_{n}\left(C_{n}-L_{n}\left(T_{n}\right)\right) 1_{\{I<R\}} \tag{20.6}
\end{equation*}
$$
[[Andersen_Ch20#^eq20-6|Eq. (20.6)]]

Shout options, which allow the holder to lock in a minimum option payoff at a chosen shout date, are priced using PDE methods. This shout option pricing formula is written as:

$$
\begin{equation*}
R+\max _{n, \omega}\left\{\tau_{n} C_{n}(\omega)\right\} \tag{20.7}
\end{equation*}
$$
[[Andersen_Ch20#^eq20-7|Eq. (20.7)]]

Forward swaption straddles are priced under the annuity measure, and are used to trade forward volatility. This swaption straddle pricing formula is written as:

$$
\begin{equation*}
C_{n}=\left|S_{n+1}\left(T_{n+1}\right)-S_{n}\left(T_{n}\right)\right| \tag{20.8}
\end{equation*}
$$
[[Andersen_Ch20#^eq20-8|Eq. (20.8)]]

The joint density of the swap rate and its variance process is solved from the forward Kolmogorov equation. This joint density formula is written as:

$$
\begin{align*}
V_{\text {swaplet }}\left(x, T_{n}\right) \triangleq \mathrm{E}_{T_{n}}^{T_{n+1}} & \left(C_{n} \mid S_{n}\left(T_{n}\right)=x\right) \\
& =\mathrm{E}_{T_{n}}^{T_{n+1}}\left(\left|S_{n+1}\left(T_{n+1}\right)-x\right| \mid S_{n}\left(T_{n}\right)=x\right) \tag{20.9}
\end{align*}
$$
[[Andersen_Ch20#^eq20-9|Eq. (20.9)]]

The drift of the swap rate under the annuity measure contains a correlation term that links swap rate and volatility. This drift formula is given by:

$$
\begin{equation*}
\mathrm{E}\left(\beta\left(T_{n}\right)^{-1} P\left(T_{n}, T_{n+1}\right) V_{\text {swaplet }}\left(S_{n}\left(T_{n}\right), T_{n}\right)\right), \tag{20.10}
\end{equation*}
$$
[[Andersen_Ch20#^eq20-10|Eq. (20.10)]]

The local volatility of the volatility swap is calibrated to the swaption smile. This local volatility calibration formula is written as:

$$
\begin{equation*}
C_{n}^{\prime}=\min \left(\left|S_{n+1}\left(\eta_{n}\right)-S_{n}\left(T_{n}\right)\right|, c\right) \tag{20.11}
\end{equation*}
$$
[[Andersen_Ch20#^eq20-11|Eq. (20.11)]]

The numerical pricing of TARNs is implemented on a two-dimensional grid using finite difference methods. This finite difference PDE recursion is written as:

$$
\begin{align*}
\mathrm{E}_{T_{n}}^{T_{n+1}}\left(\operatorname { m i n } \left(\mid S_{n+1}\left(\eta_{n}^{c}\right)-\right.\right. & \left.\left.S_{n}\left(T_{n}\right) \mid, c\right)\right) \\
& \leq \mathrm{E}_{T_{n}}^{T_{n+1}}\left(\min \left(\left|S_{n+1}\left(\eta_{n}\right)-S_{n}\left(T_{n}\right)\right|, c\right)\right) \tag{20.12}
\end{align*}
$$
[[Andersen_Ch20#^eq20-12|Eq. (20.12)]]

Finally, the chapter analyzes the impact of volatility dynamics and correlations on the pricing of volatility swaps. This volatility dynamics relation is given by:

$$
\begin{align*}
& \mathrm{E}_{T_{n}}^{T_{n+1}}\left(\min \left(\left|S_{n+1}\left(\eta_{n}\right)-S_{n}\left(T_{n}\right)\right|, c\right)\right) \\
& \quad=\mathrm{E}_{T_{n}}^{T_{n+1}}\left(\int_{T_{n}}^{T_{n+1}} \min \left(\left|S_{n+1}(t)-S_{n}\left(T_{n}\right)\right|, c\right) \delta\left(\eta_{n}-t\right) d t\right) \\
& \quad \leq \mathrm{E}_{T_{n}}^{T_{n+1}}\left(\int_{T_{n}}^{T_{n+1}} \mathrm{E}_{t}^{T_{n+1}}\left(\min \left(\left|S_{n+1}\left(\eta_{n}^{c}\right)-S_{n}\left(T_{n}\right)\right|, c\right)\right) \delta\left(\eta_{n}-t\right) d t\right) \\
& \quad=\mathrm{E}_{T_{n}}^{T_{n+1}}\left(\min \left(\left|S_{n+1}\left(\eta_{n}^{c}\right)-S_{n}\left(T_{n}\right)\right|, c\right) \int_{T_{n}}^{T_{n+1}} \delta\left(\eta_{n}-t\right) d t\right) \\
& \quad=\mathrm{E}_{T_{n}}^{T_{n+1}}\left(\min \left(\left|S_{n+1}\left(\eta_{n}^{c}\right)-S_{n}\left(T_{n}\right)\right|, c\right)\right) \tag{20.13}
\end{align*}
$$
[[Andersen_Ch20#^eq20-13|Eq. (20.13)]]

Chapter 21 addresses out-of-model adjustments, presenting systematic techniques like adjusters and path re-weighting to account for risks and cash flows not captured by standard models.

These adjustments allow traders to modify model outputs to match market realities.

This ensures that sparse pricing data and transactions are correctly adjusted.

---

### Chapter 21: Out-of-Model Adjustments

Chapter 21 addresses out-of-model adjustments, which are corrections applied by practitioners to option prices to account for risks and cash flows not captured by standard pricing models. In particular, adjusters are applied to the SDE dynamics.

These adjustments are essential for valuation and risk management of complex exotic derivatives in real-world trading environments. The SDE of the adjusted Libor rate is written as:

$$
\begin{equation*}
V_{\mathrm{mdl}}\left(C_{n} ; \xi_{0}\right)+\frac{\partial V_{\mathrm{mdl}}}{\partial \xi}\left(C_{n} ; \xi_{0}\right)\left(\xi^{*}-\xi_{0}\right) \approx V_{\mathrm{mkt}}\left(C_{n}\right), \tag{21.1}
\end{equation*}
$$
[[Andersen_Ch21#^eq21-1|Eq. (21.1)]]

Model calibration to coupon payoffs and other cash flows is the first step in out-of-model adjustments. The adjusted local volatility function is written as:

$$
\begin{equation*}
\xi^{*} \approx \xi_{0}+\left[\frac{\partial \bar{V}_{\mathrm{mdl}}}{\partial \xi}\left(C ; \xi_{0}\right)\right]^{-1}\left(\bar{V}_{\mathrm{mkt}}(C)-\bar{V}_{\mathrm{mdl}}\left(C ; \xi_{0}\right)\right) \tag{21.2}
\end{equation*}
$$
[[Andersen_Ch21#^eq21-2|Eq. (21.2)]]

The proxy model method is another key adjustment technique. The adjusted option value is computed using path weights in a Monte Carlo simulation as shown in:

$$
\begin{equation*}
\frac{\partial \bar{V}_{\mathrm{mdl}}}{\partial \xi}(C ; \xi) \tag{21.3}
\end{equation*}
$$
[[Andersen_Ch21#^eq21-3|Eq. (21.3)]]

The path weights are solved by minimizing a relative entropy distance subject to matching target benchmark prices. This entropy minimization formula is written as:

$$
\begin{align*}
V_{\mathrm{adj}}\left(H_{0}\right) \approx & V_{\mathrm{mdl}}\left(H_{0} ; \xi_{0}\right)+\frac{\partial V_{\mathrm{mdl}}}{\partial \xi}\left(H_{0}\right)\left(\xi^{*}-\xi_{0}\right) \\
= & V_{\mathrm{mdl}}\left(H_{0} ; \xi_{0}\right) \\
& +\frac{\partial V_{\mathrm{mdl}}}{\partial \xi}\left(H_{0}\right)\left[\frac{\partial \bar{V}_{\mathrm{mdl}}}{\partial \xi}\left(C ; \xi_{0}\right)\right]^{-1}\left(\bar{V}_{\mathrm{mkt}}(C)-\bar{V}_{\mathrm{mdl}}\left(C ; \xi_{0}\right)\right) . \tag{21.4}
\end{align*}
$$
[[Andersen_Ch21#^eq21-4|Eq. (21.4)]]

The adjusted option value is also computed using the difference between the base model and the proxy model. This proxy adjustment formula is given by:

$$
\begin{equation*}
V_{\mathrm{mdl}}\left(C_{n}\right)=\frac{1}{K} \sum_{l=1}^{K} C_{n}^{l} . \tag{21.5}
\end{equation*}
$$
[[Andersen_Ch21#^eq21-5|Eq. (21.5)]]

Asset-based adjustments are used to correct for structural discrepancies, such as basis spreads and multi-curve structures. This adjusted asset pricing relation is written as:

$$
\begin{gather*}
0 \leq p_{\imath} \leq 1 \quad \forall i=1, \ldots, K,  \tag{21.6}\\
\sum_{\imath=1}^{K} p_{i}=1 . \tag{21.7}
\end{gather*}
$$
[[Andersen_Ch21#^eq21-6|Eq. (21.6)]]

Mapping function adjustments are applied to correct for smile and correlation mismatches. This mapping adjustment function is written as:

$$
\begin{equation*}
\sum_{i=1}^{K} C_{n}^{i} p_{i}, \tag{21.8}
\end{equation*}
$$
[[Andersen_Ch21#^eq21-8|Eq. (21.8)]]

Fee adjustments are applied to account for transaction costs and credit valuation adjustments in swap trades. This fee adjustment formula is written as:

$$
\begin{equation*}
\sum_{\imath=1}^{K} C_{n}^{i} p_{i}=V_{\mathrm{mkt}}\left(C_{n}\right), \quad n=1, \ldots, N-1 \tag{21.9}
\end{equation*}
$$
[[Andersen_Ch21#^eq21-9|Eq. (21.9)]]

The fee cash flows are discounted under the appropriate risk-adjusted measure. This fee discounting formula is given by:

$$
\begin{equation*}
I(p) \triangleq \sum_{\imath=1}^{K} p_{\imath} \ln \left(p_{i}\right) \rightarrow \min \tag{21.10}
\end{equation*}
$$
[[Andersen_Ch21#^eq21-10|Eq. (21.10)]]

Strike adjustments are applied to option strikes to correct for smile shifts. This strike adjustment formula is written as:

$$
\begin{equation*}
\min _{\lambda} \max _{p, \mu} J(p, \lambda, \mu), \tag{21.11}
\end{equation*}
$$
[[Andersen_Ch21#^eq21-11|Eq. (21.11)]]

The strike offset is calibrated to the difference between the model and market implied volatilities. This strike calibration formula is given by:

$$
\begin{align*}
J(p, \lambda, \mu) & \triangleq-I(p) \\
& +\sum_{n=1}^{N-1} \lambda_{n}\left(\sum_{i=1}^{K} C_{n}^{i} p_{i}-V_{\mathrm{mkt}}\left(C_{n}\right)\right)+\mu\left(\sum_{i=1}^{K} p_{i}-1\right) \tag{21.12}
\end{align*}
$$
[[Andersen_Ch21#^eq21-12|Eq. (21.12)]]

The adjusted sensitivities (Greeks) are calculated by differentiating the adjusted option pricing formula. This adjusted delta formula is written as:

$$
\begin{equation*}
p_{\imath}^{*}=\frac{1}{Z(\lambda)} \exp \left(\sum_{n=1}^{N-1} \lambda_{n} C_{n}^{i}\right), \quad i=1, \ldots, K, \tag{21.13}
\end{equation*}
$$
[[Andersen_Ch21#^eq21-13|Eq. (21.13)]]

Finally, the chapter reviews the practical implementation of these adjustments in a trading environment. This adjustment implementation relation is written as:

$$
\begin{equation*}
G(\lambda) \triangleq J\left(p^{*}, \lambda, \mu^{*}\right)=\ln (Z(\lambda))-\sum_{n=1}^{N-1} \lambda_{n} V_{\mathrm{mkt}}\left(C_{n}\right) \tag{21.14}
\end{equation*}
$$
[[Andersen_Ch21#^eq21-14|Eq. (21.14)]]

Chapter 22 introduces risk management and sensitivity computations, providing the theoretical and practical foundations for option Greeks, predicting P&L, and VaR.

Traders calculate sensitivities to manage the risk of complex interest rate options.

This ensures that directional and volatility exposures are hedged effectively.

---

### Chapter 22: Introduction to Risk Management

Chapter 22 introduces risk management and sensitivity computations, providing the theoretical and practical foundations for calculating option Greeks. In this chapter, we outline basic sensitivity parameters.

Greeks, such as delta, gamma, and theta, represent the sensitivities of the option value to changes in risk factors. The delta sensitivity is defined as:

$$
\begin{equation*}
\Theta_{\mathrm{mdl}}(t)=C\left(\Theta_{\mathrm{mkt}}(t) ; \Theta_{\mathrm{prml}}(t)\right), \tag{22.1}
\end{equation*}
$$
[[Andersen_Ch22#^eq22-1|Eq. (22.1)]]

The gamma sensitivity represents the second derivative of the option value with respect to the risk factors, which is written as:

$$
\begin{equation*}
V(t)=M\left(\Theta_{\mathrm{mdl}}(t) ; \Theta_{\mathrm{num}}(t)\right), \tag{22.2}
\end{equation*}
$$
[[Andersen_Ch22#^eq22-2|Eq. (22.2)]]

The theta sensitivity represents the derivative of the option value with respect to calendar time, which is given by:

$$
\begin{equation*}
V(t)=H\left(\Theta_{\mathrm{mkt}}(t) ; \Theta_{\mathrm{prm}}(t), \Theta_{\mathrm{num}}(t)\right), \tag{22.3}
\end{equation*}
$$
[[Andersen_Ch22#^eq22-3|Eq. (22.3)]]

The analytical delta of the Black-Scholes model is used as a benchmark for sensitivity calculations. This Black-Scholes delta formula is written as:

$$
\begin{equation*}
V\left(t^{\prime}\right)=L\left(X\left(t^{\prime}\right) ; \Theta_{\mathrm{mdl}}(t)\right), \quad t^{\prime}>t, \tag{22.4}
\end{equation*}
$$
[[Andersen_Ch22#^eq22-4|Eq. (22.4)]]

The analytical gamma of the Black-Scholes model is given by the standard formula:

$$
\begin{equation*}
d X(t) / X(t)=r d t+\sigma d W^{\prime}(t) \tag{22.5}
\end{equation*}
$$
[[Andersen_Ch22#^eq22-5|Eq. (22.5)]]

The analytical vega of the Black-Scholes model is given by the derivative with respect to volatility:

$$
\begin{equation*}
\frac{\partial V}{\partial t}+r X \frac{\partial V}{\partial X}+\frac{1}{2} \sigma^{2} X^{2} \frac{\partial^{2} V}{\partial X^{2}}=r V \tag{22.6}
\end{equation*}
$$
[[Andersen_Ch22#^eq22-6|Eq. (22.6)]]

Calculating Greeks in practice requires finite-difference approximations when analytical formulas are unavailable. The finite-difference delta is computed as:

$$
\begin{equation*}
\frac{\partial V(t)}{\partial \sigma}=(T-t) X(t)^{2} \sigma \frac{\partial^{2} V(t)}{\partial X(t)^{2}}, \quad \text { or } \quad \text { vega }=(T-t) X(t)^{2} \sigma \times \text { gamma. } \tag{22.7}
\end{equation*}
$$
[[Andersen_Ch22#^eq22-7|Eq. (22.7)]]

The finite-difference gamma is computed using a central difference scheme, which is written as:

$$
\begin{equation*}
\frac{\partial V(t)}{\partial X(t)}=e^{-r(T-t)} \mathrm{E}_{t}\left(\frac{\partial g(X(T))}{\partial X(T)} \cdot \frac{X(T)}{X(t)}\right) \tag{22.8}
\end{equation*}
$$
[[Andersen_Ch22#^eq22-8|Eq. (22.8)]]

The finite-difference vega is computed by bumping the volatility surface, which is written as:

$$
\begin{align*}
\frac{\partial V(t)}{\partial t}= & r e^{-r(T-t)} \mathrm{E}_{t}(g(X(T)))+e^{-r(T-t)} \mathrm{E}_{t}\left(\frac{\partial g(X(T))}{\partial X(T)} \cdot \frac{\partial X(T)}{\partial t}\right) \\
= & r V(t)-e^{-r(T-t)} \mathrm{E}_{t}\left(\frac{\partial g(X(T))}{\partial X(T)} X(T)\left(r-\frac{1}{2} \sigma^{2}+\frac{\sigma}{2 \sqrt{T-t}} Z\right)\right) \\
= & r V(t)-X(t) \frac{\partial V(t)}{\partial X(t)}\left(r-\frac{1}{2} \sigma^{2}\right) \\
& -\frac{\sigma}{2 \sqrt{T-t}} e^{-r(T-t)} \mathrm{E}_{t}\left(\frac{\partial g(X(T))}{\partial X(T)} \cdot X(T) Z\right) \tag{22.9}
\end{align*}
$$
[[Andersen_Ch22#^eq22-9|Eq. (22.9)]]

The theta sensitivity can also be computed using a finite-difference approximation in calendar time. This finite-difference theta is given by:

$$
\begin{align*}
\frac{\partial V(t)}{\partial \sigma} & =e^{-r(T-t)} \mathrm{E}_{t}\left(\frac{\partial g(X(T))}{\partial X(T)} \cdot \frac{\partial X(T)}{\partial \sigma}\right) \\
& =e^{-r(T-t)} \mathrm{E}_{t}\left(\frac{\partial g(X(T))}{\partial X(T)} \cdot X(T)(-\sigma(T-t)+\sqrt{T-t} Z)\right) \\
& =-X(t) \frac{\partial V(t)}{\partial X(t)} \sigma(T-t)+e^{-r(T-t)} \sqrt{T-t} \mathrm{E}_{t}\left(\frac{\partial g(X(T))}{\partial X(T)} X(T) Z\right) \tag{22.10}
\end{align*}
$$
[[Andersen_Ch22#^eq22-10|Eq. (22.10)]]

Choosing the optimal bump size is critical to balance truncation error and round-off error. This optimal bump size relation is written as:

$$
\begin{equation*}
\sigma_{\mathrm{ATM}}(t, u)^{2}=(u-t)^{-1} \int_{t}^{u} \sigma(s)^{2} d s \tag{22.11}
\end{equation*}
$$
[[Andersen_Ch22#^eq22-11|Eq. (22.11)]]

P&L explain analysis is used to decompose the daily changes in portfolio value into contributions from the Greeks. This P&L predict formula is written as:

$$
\begin{equation*}
V(\delta)=H\left(\Theta_{\mathrm{mkt}}+\delta\right) \approx H\left(\Theta_{\mathrm{mkt}}\right)+\nabla^{H} \cdot \delta+\frac{1}{2} \delta^{\top} \cdot A^{H} \cdot \delta \tag{22.12}
\end{equation*}
$$
[[Andersen_Ch22#^eq22-12|Eq. (22.12)]]

Traders use waterfall explain and bump-and-reset explain methods to verify that actual P&L matches the predicted risk. This P&L explain formula is written as:

$$
\begin{equation*}
d \Pi(t)=n_{X}(t) d X(t)+N(t) r d t \tag{22.13}
\end{equation*}
$$
[[Andersen_Ch22#^eq22-13|Eq. (22.13)]]

Chapter 23 develops payoff smoothing techniques to improve the stability and convergence of Greeks in numerical integration and Monte Carlo simulation.

We will see how payoff discontinuities introduce substantial noise in sensitivities.

Smoothing the payoff function recovers smooth convergence and stable sensitivities.

---

### Chapter 23: Payoff Smoothing and Related Methods

Chapter 23 develops payoff smoothing and related numerical methods to improve the stability and convergence of Greeks. In particular, discontinuous payoffs cause numerical instabilities.

Numerical pricing schemes, such as finite difference PDE grids and Monte Carlo simulations, often exhibit slow convergence and high variance when pricing options with discontinuous payoffs. The option value under numerical integration is written as:

$$
\begin{equation*}
N=\left\lfloor\frac{10 \sigma \sqrt{T}}{\delta}\right\rfloor, \tag{23.1}
\end{equation*}
$$
[[Andersen_Ch23#^eq23-1|Eq. (23.1)]]

The discretization error is analyzed as a function of the grid step size and payoff smoothness. This discretization error formula is given by:

$$
\begin{align*}
V(S) & =\int_{-\infty}^{\infty} f(x) \pi(x-S) d x  \tag{23.2}\\
& =\int_{-\infty}^{\infty} f(S+x) \pi(x) d x \tag{23.3}
\end{align*}
$$
[[Andersen_Ch23#^eq23-2|Eq. (23.2)]]

To resolve these issues, the chapter introduces payoff smoothing techniques, which replace the discontinuous payoff function with a smooth approximation. This smoothed option value is written as:

$$
\begin{equation*}
\mathrm{E}\left(f(S(T))=\mathrm{E}\left(\mathrm{E}_{T-\delta}(f(S(T)))\right)=\int_{-\infty}^{\infty} \pi(x, T-\delta ; S, 0) V(x, T-\delta) d x\right. \tag{23.4}
\end{equation*}
$$
[[Andersen_Ch23#^eq23-4|Eq. (23.4)]]

One-dimensional box smoothing replaces the step function with a smooth cumulative distribution function. This smoothed digital payoff is given by:

$$
\begin{equation*}
V_{\mathrm{cpn}, n}(0)=\mathrm{E}\left(B\left(T_{n-2}\right)^{-1} V_{\mathrm{cpn}, n}\left(T_{n-2}\right)\right), \tag{23.5}
\end{equation*}
$$
[[Andersen_Ch23#^eq23-5|Eq. (23.5)]]

Alternatively, Gaussian smoothing is applied to smooth digital options, which is written as:

$$
\begin{equation*}
f_{\mathrm{smooth}}(x)=\frac{1}{\epsilon} \int_{x-\epsilon / 2}^{x+\epsilon / 2} f(y) d y \tag{23.6}
\end{equation*}
$$
[[Andersen_Ch23#^eq23-6|Eq. (23.6)]]

Payoff smoothing in multiple dimensions requires specifying a smoothing covariance matrix. This smoothed multi-D option value is written as:

$$
\begin{equation*}
\tilde{V}(S)=\sum_{n=0}^{N} w_{n} f\left(x_{n}+S\right) . \tag{23.7}
\end{equation*}
$$
[[Andersen_Ch23#^eq23-7|Eq. (23.7)]]

For Monte Carlo simulation, the chapter develops Tube Monte Carlo, which restricts the simulation paths to a narrow tube around the barrier. This Tube Monte Carlo estimator is written as:

$$
\begin{equation*}
\tilde{V}(S)=\sum_{n=0}^{N} w_{n} f_{\text {smooth }}\left(x_{n}+S\right) . \tag{23.8}
\end{equation*}
$$
[[Andersen_Ch23#^eq23-8|Eq. (23.8)]]

The tube density is chosen to minimize the variance of the pathwise estimator. This tube density function is given by:

$$
\begin{equation*}
f_{\mathrm{smooth}}(x)=\frac{1}{\epsilon} \int_{x-\epsilon / 2}^{x+\epsilon / 2} f(y) d y \approx f(x) \tag{23.9}
\end{equation*}
$$
[[Andersen_Ch23#^eq23-9|Eq. (23.9)]]

The likelihood ratio weight is applied to option payoffs to restore unbiasedness under the tube measure. This likelihood ratio weight is written as:

$$
\begin{align*}
f_{\mathrm{smooth}}(x)=\frac{1}{\epsilon} \int_{x-\epsilon / 2}^{x+\epsilon / 2} f(y) d y & \approx \frac{x^{*}-x+\epsilon / 2}{\epsilon} f\left(\frac{x-\epsilon / 2+x^{*}}{2}\right) \\
& +\frac{x+\epsilon / 2-x^{*}}{\epsilon} f\left(\frac{x+\epsilon / 2+x^{*}}{2}\right) \tag{23.10}
\end{align*}
$$
[[Andersen_Ch23#^eq23-10|Eq. (23.10)]]

The variance of the Tube Monte Carlo estimator is compared with standard Monte Carlo. This variance reduction ratio is given by:

$$
\begin{equation*}
f_{n}=\int_{-\infty}^{\infty} \kappa_{n}(x) f(x) d x, \quad n=0, \ldots, N \tag{23.11}
\end{equation*}
$$
[[Andersen_Ch23#^eq23-11|Eq. (23.11)]]

The conditional independence approximation is used to simplify the multidimensional expectation. This conditional expectation is written as:

$$
\begin{align*}
\kappa_{n}(x) & =\left(c_{n+1}-c_{n}\right)^{-1} 1_{\left\{x \in\left[c_{n}, c_{n+1}\right]\right\}}  \tag{23.12}\\
c_{n} & =\left(x_{n}+x_{n-1}\right) / 2, \quad n=1, \ldots, N-1
\end{align*}
$$
[[Andersen_Ch23#^eq23-12|Eq. (23.12)]]

The drift of the simulation paths is shifted to keep them inside the tube. This drift shift formula is given by:

$$
\begin{equation*}
f_{n}=\frac{1}{c_{n+1}-c_{n}} \int_{c_{n}}^{c_{n+1}} f(x) d x, \quad c_{n}=\left(x_{n}+x_{n-1}\right) / 2, \quad n=1, \ldots, N-1 \tag{23.13}
\end{equation*}
$$
[[Andersen_Ch23#^eq23-13|Eq. (23.13)]]

Finally, the chapter reviews numerical tests of Tube Monte Carlo for digital and barrier options. This numerical test relation is written as:

$$
\begin{equation*}
V(T-, x)=1_{\{g(x) \leq h(x)\}} p(x)+1_{\{g(x)>h(x)\}} q(x), \tag{23.14}
\end{equation*}
$$
[[Andersen_Ch23#^eq23-14|Eq. (23.14)]]

Chapter 24 focuses on pathwise differentiation, which calculates exact Greeks in a single Monte Carlo simulation run by differentiating state transition equations along each path.

We will explore pathwise sensitivity calculations in the Libor Market Model.

This yields exact sensitivities without running multiple bumped simulations.

---

### Chapter 24: Pathwise Differentiation

Chapter 24 focuses on pathwise differentiation, which is a powerful method for calculating Greeks in Monte Carlo simulation. In this chapter, we outline pathwise derivatives.

Unlike finite difference methods that require running multiple bumped simulations, pathwise differentiation calculates the sensitivities of individual paths. The pathwise derivative is defined as:

$$
\begin{equation*}
H_{n-1}\left(T_{n-1}\right)=B\left(T_{n-1}\right) \mathrm{E}_{T_{n-1}}\left(B\left(T_{n}\right)^{-1} \max \left(H_{n}\left(T_{n}\right), U_{n}\left(T_{n}\right)\right)\right) \tag{24.1}
\end{equation*}
$$
[[Andersen_Ch24#^eq24-1|Eq. (24.1)]]

The pathwise derivative of the state variable is solved from the transition SDE, which is written as:

$$
\begin{equation*}
\Delta_{\alpha} H_{0}(0)=\mathrm{E}\left(\sum_{n=\eta}^{N-1} \Delta_{\alpha}\left(B\left(T_{n}\right)^{-1} X_{n}\right)\right) . \tag{24.2}
\end{equation*}
$$
[[Andersen_Ch24#^eq24-2|Eq. (24.2)]]

For callable exotics, PDM must account for the sensitivity of the optimal exercise time to the underlying parameters. This optimal stopping time derivative is given by:

$$
\begin{align*}
H_{0}(0) & =\mathrm{E}\left(\sum_{n=\eta}^{N-1} B\left(T_{n+1}\right)^{-1} X_{n}\right)  \tag{24.3}\\
\Delta_{\alpha} H_{0}(0) & =\mathrm{E}\left(\sum_{n=\eta}^{N-1} \Delta_{\alpha}\left(B\left(T_{n+1}\right)^{-1} X_{n}\right)\right)
\end{align*}
$$
[[Andersen_Ch24#^eq24-3|Eq. (24.3)]]

The pathwise derivative of the spot measure numeraire is solved by differentiating the short rate integral. This numeraire derivative formula is written as:

$$
\begin{equation*}
\Delta_{\alpha} V_{\mathrm{ki}}(\zeta, X)=\frac{\partial}{\partial \zeta} V_{\mathrm{ki}}(\zeta, X) \Delta_{\alpha} \zeta+\frac{\partial}{\partial X} V_{\mathrm{ki}}(\zeta, X) \Delta_{\alpha} X \tag{24.4}
\end{equation*}
$$
[[Andersen_Ch24#^eq24-4|Eq. (24.4)]]

The pathwise derivatives of forward Libor rates in the LMM are solved recursively. This pathwise Libor derivative recursion is given by:

$$
\begin{equation*}
\left.\frac{\partial}{\partial \zeta} V_{\mathrm{ki}}(\zeta, X)\right|_{\zeta=\eta}=0 \tag{24.5}
\end{equation*}
$$
[[Andersen_Ch24#^eq24-5|Eq. (24.5)]]

The pathwise delta of a European option under the LMM is computed using the derivatives of the rates and annuity. This pathwise delta formula is written as:

$$
\begin{equation*}
X=1_{\{G>h\}} R, \tag{24.6}
\end{equation*}
$$
[[Andersen_Ch24#^eq24-6|Eq. (24.6)]]

The Adjoint Method (AEM) is introduced to calculate Greeks when the number of risk factors is large. This adjoint gradient vector is written as:

$$
\begin{align*}
& \Delta_{\alpha} \mathrm{E}\left(B(T)^{-1} 1_{\{G>h\}} R\right)=\mathrm{E}\left(\Delta_{\alpha}\left(B(T)^{-1}\right) 1_{\{G>h\}} R\right) \\
& \quad+\mathrm{E}\left(B(T)^{-1} 1_{\{G>h\}} \Delta_{\alpha} R\right)+\gamma_{G}(h) \mathrm{E}\left(B(T)^{-1} R \Delta_{\alpha} G \mid G=h\right) \tag{24.7}
\end{align*}
$$
[[Andersen_Ch24#^eq24-7|Eq. (24.7)]]

The pathwise delta of a callable Libor exotic is approximated by keeping the exercise boundary fixed. This pathwise CLE delta formula is given by:

$$
\begin{equation*}
V_{\mathrm{ki}, n}(t)=B(t) \mathrm{E}_{t}\left(B\left(T_{\eta_{n}}\right)^{-1} R_{\eta_{n}}\right), \tag{24.8}
\end{equation*}
$$
[[Andersen_Ch24#^eq24-8|Eq. (24.8)]]

The adjoint state variable is solved recursively backward from the terminal payoff to the initial state. This adjoint recurrence relation is written as:

$$
\begin{align*}
& \Delta_{\alpha} V_{\mathrm{ki}}(0)=\mathrm{E}\left(\left.B\left(T_{\eta}\right)^{-1} \Delta_{\alpha} R_{n}\right|_{n=\eta}\right) \\
& \quad+\mathrm{E}\left(\left.B\left(T_{\eta}\right)^{-1} \gamma_{\eta}\left(h_{\eta}\right)\left(R_{\eta}-V_{\mathrm{ki}, \eta}\left(T_{\eta}\right)\right) \Delta_{\alpha} G_{n}\right|_{n=\eta} \mid G_{\eta}=h_{\eta}\right) \tag{24.9}
\end{align*}
$$
[[Andersen_Ch24#^eq24-9|Eq. (24.9)]]

The adjoint sensitivity of the Libor rates is computed by accumulating the pathwise derivatives. This adjoint Libor sensitivity is written as:

$$
\begin{equation*}
\mathcal{L}=(y(t)-\varkappa(t) x) \frac{\partial}{\partial x}+\frac{1}{2} \sigma_{r}(t)^{2} \frac{\partial^{2}}{\partial x^{2}} \tag{24.10}
\end{equation*}
$$
[[Andersen_Ch24#^eq24-10|Eq. (24.10)]]

The likelihood ratio method (LRM) is introduced as an alternative when the option payoff is discontinuous. This LRM derivative formula is written as:

$$
\begin{align*}
\Delta_{m} H_{0}(0)= & \sum_{n=1}^{N-1} \mathrm{E}\left(1_{\{\eta \geq n\}} \Delta_{m}\left(\beta\left(T_{n+1}\right)^{-1} X_{n}\right)\right) \\
= & \sum_{n=1}^{N-1} \mathrm{E}\left(1_{\{\eta \geq n\}} \beta\left(T_{n+1}\right)^{-1} \frac{\Delta_{m}\left(\beta\left(T_{n+1}\right)^{-1}\right)}{\beta\left(T_{n+1}\right)^{-1}} X_{n}\right) \\
& +\sum_{n=1}^{N-1} \mathrm{E}\left(1_{\{\eta \geq n\}} \beta\left(T_{n+1}\right)^{-1} \Delta_{m}\left(X_{n}\right)\right) \tag{24.11}
\end{align*}
$$
[[Andersen_Ch24#^eq24-11|Eq. (24.11)]]

The hybrid pathwise-likelihood ratio method combines PDM and LRM to price options with discontinuous payoffs. This hybrid derivative formula is written as:

$$
\begin{equation*}
\frac{\partial}{\partial t} \pi_{s, x}(t, y)=\mathcal{L}^{*} \pi_{s, x}(t, y)-r(t) \pi_{s, x}(t, y), \quad(s, x) \text { fixed, } \tag{24.12}
\end{equation*}
$$
[[Andersen_Ch24#^eq24-12|Eq. (24.12)]]

Finally, the chapter reviews numerical examples of pricing Bermudan swaption Greeks using PDM and AEM. This numerical Greeks relation is written as:

$$
\begin{equation*}
\psi\left(y ; T_{n}\right)=\psi\left(y ; T_{n}-\right) \times 1_{\left\{H_{n}\left(T_{n}, y\right)>U_{n}\left(, T_{n}, y\right)\right\}} . \tag{24.13}
\end{equation*}
$$
[[Andersen_Ch24#^eq24-13|Eq. (24.13)]]

Chapter 25 covers importance sampling and control variates, presenting systematic variance reduction techniques to speed up Monte Carlo pricing of complex interest rate portfolios.

We will see how shifted measures direct path simulation to region of interest.

This dramatically improves simulation efficiency and reduces error.

---

### Chapter 25: Importance Sampling and Control Variates

Chapter 25 covers importance sampling (IS) and control variates, which are the two primary variance reduction techniques used in Monte Carlo option pricing. In this chapter, we analyze importance sampling measure shifts.

Importance sampling reduces variance by shifting the probability measure. The option price under the importance sampling measure is written as:

$$
\begin{equation*}
X(0, T)=\mathrm{E}\left(\exp \left(-\int_{0}^{T} x(u) d u\right)\right) \triangleq \mathrm{E}(Y(T)) \tag{25.1}
\end{equation*}
$$
[[Andersen_Ch25#^eq25-1|Eq. (25.1)]]

The likelihood ratio (Radon-Nikodym derivative) is used to correct the expectation under the shifted measure. This likelihood ratio is written as:

$$
\begin{align*}
\varsigma(t) & =\exp \left(-\int_{0}^{t} x(u) d u\right) X(t, T, x(t)) / X(0, T) \\
q(t, x(t)) & =-X(t, T, x(t))^{-1} \sigma_{x}(t, x(t)) \frac{\partial X(t, T, x(t))}{\partial x(t)} \tag{25.2}
\end{align*}
$$
[[Andersen_Ch25#^eq25-2|Eq. (25.2)]]

In short-rate models, the drift of the state process is shifted to generate more paths near the strike. This shifted short-rate SDE is written as:

$$
\begin{equation*}
q(t, x)=\sigma_{x}(t, x) \frac{1-e^{-\varkappa_{G}(T-t)}}{\varkappa_{G}} . \tag{25.3}
\end{equation*}
$$
[[Andersen_Ch25#^eq25-3|Eq. (25.3)]]

For trigger products like TARNs, the trigger event is smoothed using importance sampling. This trigger indicator smoothing formula is written as:

$$
\begin{equation*}
V=\mathrm{E}\left(g(X) 1_{\{X<b\}}\right), \tag{25.4}
\end{equation*}
$$
[[Andersen_Ch25#^eq25-4|Eq. (25.4)]]

Simulating rates under the survival measure ensures that more paths survive the barrier. This survival measure option price is given by:

$$
\begin{equation*}
V=\mathrm{E}(g(X) \mid X<b) \mathrm{P}(X<b) . \tag{25.5}
\end{equation*}
$$
[[Andersen_Ch25#^eq25-5|Eq. (25.5)]]

The drift of the forward rates under the survival measure is adjusted to account for the survival conditioning. This survival-adjusted SDE is written as:

$$
\begin{equation*}
X=\Phi^{-1}(U) \tag{25.6}
\end{equation*}
$$
[[Andersen_Ch25#^eq25-6|Eq. (25.6)]]

Control variates reduce variance by adding a zero-expectation variable to the estimator. The control variates estimator is written as:

$$
\begin{equation*}
X \mid\{X<b\}=\Phi^{-1}\left(U^{\prime}\right), \quad U^{\prime} \sim \mathcal{U}(0, \Phi(b)) . \tag{25.7}
\end{equation*}
$$
[[Andersen_Ch25#^eq25-7|Eq. (25.7)]]

The optimal control variate coefficient minimizes the variance of the estimator. This optimal coefficient is given by:

$$
\begin{equation*}
V=\mathrm{E}\left(g\left(\Phi^{-1}\left(U^{\prime}\right)\right)\right) \Phi\left(\frac{b-\mu}{\sigma}\right), \quad U^{\prime} \sim \mathcal{U}(0, \Phi(b)) . \tag{25.8}
\end{equation*}
$$
[[Andersen_Ch25#^eq25-8|Eq. (25.8)]]

The low-dimensional Markov approximation of the LMM is used as a model-based control variate. This low-dimensional Markov approximation SDE is written as:

$$
\begin{equation*}
V_{\mathrm{tarn}}(0)=\mathrm{E}\left(\sum_{n=1}^{N-1} B\left(T_{n+1}\right)^{-1} X_{n}\left(T_{n}\right) 1_{\left\{Q_{n}<R\right\}}\right), \tag{25.9}
\end{equation*}
$$
[[Andersen_Ch25#^eq25-9|Eq. (25.9)]]

The volatility structure of the proxy model is calibrated to match the base model. This proxy volatility calibration is given by:

$$
\begin{equation*}
C_{n}=\left(s-g \times L_{n}\left(T_{n}\right)\right)^{+} . \tag{25.10}
\end{equation*}
$$
[[Andersen_Ch25#^eq25-10|Eq. (25.10)]]

The analytical option price in the proxy model is used to define the control variate. This proxy pricing relation is written as:

$$
\begin{equation*}
b_{n}=\left(s-\left(R-Q_{n}\right) / \tau_{n}\right) / g \tag{25.11}
\end{equation*}
$$
[[Andersen_Ch25#^eq25-11|Eq. (25.11)]]

The variance reduction factor is analyzed as a function of the correlation between the option and the control variate. This variance reduction formula is given by:

$$
\begin{equation*}
\mathrm{E}(\mathcal{V})=\mathrm{E}\left(\mathcal{V} \mid L_{1}\left(T_{1}\right)>b_{1}\right) \mathrm{Q}^{B}\left(L_{1}\left(T_{1}\right)>b_{1}\right) . \tag{25.12}
\end{equation*}
$$
[[Andersen_Ch25#^eq25-12|Eq. (25.12)]]

Finally, the chapter reviews instrument-based and dynamic control variates. This dynamic control variate relation is written as:

$$
\begin{align*}
V_{\operatorname{tarn}}(0) & =\sum_{n=1}^{N-1} \widetilde{\mathrm{E}}\left(\psi_{n} \mathrm{E}_{T_{n-1}}\left(B\left(T_{n+1}\right)^{-1} X_{n}\left(T_{n}\right)\right)\right)  \tag{25.13}\\
\psi_{n} & =\prod_{k=1}^{n-1} \mathrm{Q}_{T_{k-1}}^{B}\left(L_{k}\left(T_{k}\right)>b_{k}\right) \tag{25.14}
\end{align*}
$$
[[Andersen_Ch25#^eq25-13|Eq. (25.13)]]

Chapter 26 focuses on Vegas in Libor Market Models, comparing direct, indirect, and hybrid Vega calculation methods and discussing Vega projections.

We will analyze the high-dimensional volatility parameter space in LMMs.

This requires advanced projection methods to hedge model sensitivities.

---

### Chapter 26: Vegas in Libor Market Models

Chapter 26 focuses on Vegas in Libor Market Models (LMMs). In this chapter, we define the LMM vega problem.

Vega represents the sensitivity of option prices to changes in the volatility parameters, and is crucial for hedging volatility risk. The Vega sensitivity is defined as:

$$
\begin{equation*}
\sigma_{f}(t, T)=g(t) h(T) \tag{26.1}
\end{equation*}
$$
[[Andersen_Ch26#^eq26-1|Eq. (26.1)]]

Direct Vegas are computed by differentiating the option price with respect to the model parameters, which is written as:

$$
\begin{equation*}
(\Lambda(G))_{i, j}, \quad i=1, \ldots, N_{t}, \quad j=1, \ldots, N_{x} . \tag{26.2}
\end{equation*}
$$
[[Andersen_Ch26#^eq26-2|Eq. (26.2)]]

Cross-Vegas represent the sensitivity of the option price to joint changes in two volatility parameters. This cross-Vega formula is given by:

$$
\begin{equation*}
G^{*}(\widehat{\Lambda})=\underset{G}{\operatorname{argmin}} \mathcal{I}(G ; \widehat{\Lambda}) . \tag{26.3}
\end{equation*}
$$
[[Andersen_Ch26#^eq26-3|Eq. (26.3)]]

Indirect Vegas are computed by differentiating the option price with respect to the market targets, which captures refitting effects. This indirect Vega formula is written as:

$$
\begin{equation*}
\nu_{\mathrm{mkt}}(\delta)=\left.\epsilon^{-1}\left(V\left(G^{*}(\widehat{\Lambda}+\epsilon \delta)\right)-V\left(G^{*}(\widehat{\Lambda})\right)\right) \approx \frac{d}{d u} V\left(G^{*}(\widehat{\Lambda}+u \delta)\right)\right|_{u=0} \tag{26.4}
\end{equation*}
$$
[[Andersen_Ch26#^eq26-4|Eq. (26.4)]]

Hybrid Vegas combine direct and indirect methods to balance speed and accuracy. This hybrid Vega formula is written as:

$$
\begin{equation*}
\left(\delta_{n, m}\right)_{i, j}=1_{\{i=n\}} 1_{\{j=m\}} \tag{26.5}
\end{equation*}
$$
[[Andersen_Ch26#^eq26-5|Eq. (26.5)]]

The model-consistent delta with backbone accounts for the sensitivity of the volatility surface to the underlying rate. This swap rate local volatility SDE is written as:

$$
\left(\delta_{n, m}^{\prime}\right)_{i, j}=\left\{\begin{array}{l}
1, i<n \text { or } i=n, j \leq m  \tag{26.6}\\
0, \text { otherwise }
\end{array}\right.
$$
[[Andersen_Ch26#^eq26-6|Eq. (26.6)]]

The delta with backbone is computed by adding a vega-leakage term to the Black delta. This delta with backbone formula is given by:

$$
\begin{equation*}
\nu_{\mathrm{mkt}}\left(\delta_{n, m}\right)=\nu_{\mathrm{mkt}}\left(\delta_{n, m}^{\prime}\right)-\nu_{\mathrm{mkt}}\left(\delta_{n, m-1}^{\prime}\right) \tag{26.7}
\end{equation*}
$$
[[Andersen_Ch26#^eq26-7|Eq. (26.7)]]

Vega projections are introduced to project the high-dimensional model Vegas onto liquid benchmark swaptions. This projected Vega formula is written as:

$$
\begin{equation*}
\left(\delta_{n, m}^{\prime \prime}\right)_{i, j}=1_{\{i \leq n\}} 1_{\{j \leq m\}} \tag{26.8}
\end{equation*}
$$
[[Andersen_Ch26#^eq26-8|Eq. (26.8)]]

The projection matrix is solved by minimizing the projection error under a chosen norm. This projection matrix formula is written as:

$$
\begin{equation*}
\Lambda\left(G^{*}(\widehat{\Lambda})\right)-\widehat{\Lambda} \approx \Lambda\left(G^{*}(\widehat{\Lambda}+\delta)\right)-(\widehat{\Lambda}+\delta) . \tag{26.9}
\end{equation*}
$$
[[Andersen_Ch26#^eq26-9|Eq. (26.9)]]

The term correlation between distinct rates is kept constant during calibration to capture correlation smile effects. This correlation constraint is written as:

$$
\begin{equation*}
\frac{\partial V}{\partial G}=\frac{\partial V}{\partial \Lambda} \frac{\partial \Lambda}{\partial G} \tag{26.10}
\end{equation*}
$$
[[Andersen_Ch26#^eq26-10|Eq. (26.10)]]

The correlation Vega represents the sensitivity of the option price to changes in correlation parameters. This correlation Vega formula is given by:

$$
\begin{equation*}
\frac{\partial V}{\partial \Lambda}=\left(\frac{\partial \Lambda}{\partial G}\right)^{-1} \frac{\partial V}{\partial G} \tag{26.11}
\end{equation*}
$$
[[Andersen_Ch26#^eq26-11|Eq. (26.11)]]

The numerical efficiency of the Vega algorithms is compared in a test calibration. This Vega efficiency relation is written as:

$$
\begin{equation*}
\sum_{m=1}^{N_{i}} \sum_{n=1}^{N_{i}}\left(\nu_{\mathrm{mdl}}\left(\delta_{n, m}\right)-\sum_{j=1}^{N_{i}} \sum_{i=1}^{N_{i}}\left(\nu_{\mathrm{mkt}}\right)_{2, j} \frac{\partial \Lambda_{i, j}}{\partial \delta_{n, m}}\right)^{2} \rightarrow \min \tag{26.12}
\end{equation*}
$$
[[Andersen_Ch26#^eq26-12|Eq. (26.12)]]

Finally, the chapter reviews practical guidelines for managing Vega risk in LMM trading desks. This Vega management relation is written as:

$$
\begin{align*}
\sum_{n, m} W_{n, m}^{2}\left(\nu_{\mathrm{mdl}}\left(\delta_{n, m}\right)\right. & \left.-\sum_{i, j}\left(\nu_{\mathrm{mkt}}\right)_{i, j} \frac{\partial \Lambda_{l, j}}{\partial \delta_{n, m}}\right)^{2} \\
& +\sum_{i, j} U_{l, j}^{2}\left(\left(\nu_{\mathrm{mkt}}\right)_{i, j}\right)^{2} \rightarrow \mathrm{~min} \tag{26.13}
\end{align*}
$$
[[Andersen_Ch26#^eq26-13|Eq. (26.13)]]

Appendix A presents the mathematical foundations of Markovian projection, defining how multidimensional Ito processes can be mimicked by low-dimensional Markov diffusions.

We will see how Gyongy's theorem forms the core of LSV projection.

This maps complex volatility dynamics to simple proxy diffusions.

---

### Chapter 27: Markovian Projection

Appendix A presents the mathematical foundations of Markovian projection, which is a powerful method to simplify complex multi-dimensional process dynamics. In this section, we analyze Ito process marginals.

This simplifies complex multi-dimensional process dynamics into a low-dimensional Markovian form that preserves the marginal distributions. The SDE of the multi-dimensional Ito process is written as:

$$
\begin{equation*}
d X(t)=\lambda(t) d W(t) \tag{(A. 1 !}
\end{equation*}
$$
[[Andersen_AppA#^eqa-1|Eq. (A.1)]]

The diffusion coefficient of the mimicking process is defined as the conditional expectation of the original volatility process. This conditional expected variance is written as:

$$
\begin{equation*}
b(t, x)^{2}=\mathrm{E}\left(\lambda(t)^{2} \mid X(t)=x\right) \tag{A.2}
\end{equation*}
$$
[[Andersen_AppA#^eqa-2|Eq. (A.2)]]

Gyöngy's theorem states that the mimicking SDE admits a weak solution that has the same marginal distributions. This mimicking SDE is written as:

$$
\begin{equation*}
d Y(t)=b(t, Y(t)) d W(t), \quad Y(0)=X(0) \tag{A.3}
\end{equation*}
$$
[[Andersen_AppA#^eqa-3|Eq. (A.3)]]

The Dupire local volatility function is recovered from option prices using the forward PDE. This local volatility SDE is written as:

$$
\begin{equation*}
d Y(t)=b(t, Y(t)) d W(t) \tag{A.4}
\end{equation*}
$$
[[Andersen_AppA#^eqa-4|Eq. (A.4)]]

The local variance function is evaluated as a function of the derivative of option prices. This Dupire local variance formula is written as:

$$
\begin{equation*}
b(t, K)^{2}=\frac{2 \frac{\partial}{\partial t} c(t, K)}{\frac{\partial^{2}}{\partial K^{2}} c(t, K)}, \tag{A.5}
\end{equation*}
$$
[[Andersen_AppA#^eqa-5|Eq. (A.5)]]

The multi-factor extension of Gyöngy's theorem maps a multi-dimensional volatility loading to a one-dimensional Brownian driver. This multi-factor SDE is written as:

$$
\begin{equation*}
d X(t)=\lambda(t)^{\top} d W(t) \tag{A.6}
\end{equation*}
$$
[[Andersen_AppA#^eqa-6|Eq. (A.6)]]

The mimicking process for the multi-factor SDE is driven by a one-dimensional Brownian motion. This mimicking SDE is written as:

$$
\begin{equation*}
d Y(t)=b(t, Y(t)) d \widetilde{W}(t), \quad Y(0)=X(0) \tag{A.7}
\end{equation*}
$$
[[Andersen_AppA#^eqa-7|Eq. (A.7)]]

Analytical approximations, such as the Gaussian approximation, are used to compute the conditional expected values in closed form. This Gaussian local variance approximation is written as:

$$
\begin{align*}
b(t, x)^{2} & \approx \overline{\lambda(t)^{2}}+s(t)(x-X(0))  \tag{A.8}\\
s(t) & =\frac{\int_{0}^{t} \bar{\varepsilon}(s) \bar{\lambda}(s) \bar{\rho}(s) d s}{\int_{0}^{t} \overline{\lambda(s)^{2}} d s}
\end{align*}
$$
[[Andersen_AppA#^eqa-8|Eq. (A.8)]]

The conditional expectation is evaluated using the standard conditioning formula for Gaussian variables. This Gaussian expectation is written as:

$$
\begin{equation*}
\mathrm{E}(U \mid V)=\mathrm{E}(U)+\frac{\operatorname{Cov}(U, V)}{\operatorname{Var}(V)}(V-\mathrm{E}(V)) \tag{A.9}
\end{equation*}
$$
[[Andersen_AppA#^eqa-9|Eq. (A.9)]]

The least-squares projection method is an alternative approach that projects the variance onto a linear subspace. This least-squares local variance is written as:

$$
\begin{equation*}
b(t, x)^{2} \approx \overline{\lambda(t)^{2}}+s(t)(x-X(0)) \tag{A.10}
\end{equation*}
$$
[[Andersen_AppA#^eqa-10|Eq. (A.10)]]

The slope of the least-squares projection is solved by minimizing the projection error. This least-squares slope is written as:

$$
\begin{equation*}
s(t)=\frac{\mathrm{E}\left(\lambda(t)^{2}(X(t)-X(0))\right)}{\mathrm{E}\left((X(t)-X(0))^{2}\right)}, \quad \overline{\lambda(t)^{2}}=\mathrm{E}\left(\lambda(t)^{2}\right) \tag{A.11}
\end{equation*}
$$
[[Andersen_AppA#^eqa-11|Eq. (A.11)]]

The displaced log-normal approximation replaces the local volatility function with a linear function of spot. This displaced log-normal SDE is written as:

$$
\begin{equation*}
d Y(t)=\sigma(t)(1+b(t)(Y(t)-Y(0))) d W(t) \tag{A.12}
\end{equation*}
$$
[[Andersen_AppA#^eqa-12|Eq. (A.12)]]

The local volatility parameters of the displaced log-normal process are solved from the Gaussian moments. This parameter estimation is given by:

$$
\begin{equation*}
d X_{2}(t)=b_{2}\left(t, X_{2}(t)\right) \sqrt{z_{2}(t)} d W(t) . \tag{A.13}
\end{equation*}
$$
[[Andersen_AppA#^eqa-13|Eq. (A.13)]]

Finally, we conclude the technical appendix on Markovian projection methods.

These projection techniques are widely used to build local stochastic volatility models from high-dimensional Libor market models.

In the next section, we present a detailed comparison of pathwise differentiation and Malliavin calculus for sensitivity analysis.

---

### Chapter 28: Sensitivity Deep Dive: PDM vs. Malliavin Calculus in LMMs

In this deep dive, we compare the Pathwise Differentiation Method (PDM) and Malliavin calculus for calculating Monte Carlo Greeks in Libor Market Models (LMMs). Option Greeks are essential for risk-managing interest rate portfolios, but their calculation in a Monte Carlo framework is often computationally expensive and subject to high numerical noise.

We compare the two methods in terms of variance, speed, and analytical complexity. Let us begin with the formulation of the LMM rate SDE under spot measure.

#### Pathwise Differentiation Method (PDM)

The pathwise differentiation method (PDM) computes Greeks by differentiating the state transition equations and the option payoff along each simulated path. Under the spot measure, the joint SDE of the forward Libor rates is given by Eq. (24.14):

$$
\begin{align*}
d L_{n}(t) & =\varphi\left(L_{n}(t)\right) \lambda_{n}(t)^{\top}\left(\mu_{n}(t) d t+d W(t)\right)  \tag{24.14}\\
\mu_{n}(t) & =\sum_{j=q(t)}^{n} \frac{\tau_{j} \varphi\left(L_{j}(t)\right) \lambda_{j}(t)}{1+\tau_{j} L_{j}(t)}
\end{align*}
$$
[[Andersen_Ch24#^eq24-14|Eq. (24.14)]]

Differentiating this transition SDE with respect to the initial rate $L_m(0)$ yields the pathwise derivative of the Libor rates at future times. This recursion is given by Eq. (24.15):

$$
\begin{align*}
& d\left(\Delta_{m} L_{n}(t)\right)=\varphi^{\prime}\left(L_{n}(t)\right) \lambda_{n}(t)^{\top} \Delta_{m} L_{n}(t)\left(\mu_{n}(t) d t+d W(t)\right) \\
& \quad+\varphi\left(L_{n}(t)\right) \lambda_{n}(t)^{\top} \sum_{j} \frac{\partial \mu_{n}(t)}{\partial L_{j}(t)} \Delta_{m} L_{j}(t) d t \tag{24.15}
\end{align*}
$$
[[Andersen_Ch24#^eq24-15|Eq. (24.15)]]

For options with continuous payoffs, pathwise differentiation is highly efficient and yields exact Greeks in a single simulation run. However, when the option payoff is discontinuous (e.g., digital options or trigger products), PDM fails because the derivative of the indicator function is a Dirac delta function, which has zero value almost everywhere on the path space and cannot be simulated directly.

#### Malliavin Calculus Approach

To resolve the discontinuity issue, the Malliavin calculus approach (also known as the likelihood ratio method in simple settings) applies integration by parts on the path space to shift the derivative from the payoff function to the probability density function. For example, the expected value of the derivative of the payoff is evaluated by multiplying the payoff by a Malliavin weight. This integration by parts relation is written as:

$$
\mathrm{E}\left(\frac{\partial}{\partial x} f(X)\right) = \mathrm{E}\left(f(X) \cdot \text{Weight}\right)
$$

Here, the Malliavin weight is defined as the stochastic integral of the Malliavin derivative of the log-density over the path space, which is given by:

$$
\text{Weight} = \int_0^T h(t) dW(t)
$$

The Malliavin weight is independent of the option payoff, allowing it to price options with discontinuous or even non-analytic payoffs. However, because the Malliavin weight involves the sum of stochastic integrals along each path, it has significantly higher variance than the pathwise estimator, requiring a larger number of simulation paths to achieve the same level of accuracy.

#### Computational Complexity and Adjoint Method (AEM)

When the number of risk factors is large, direct pathwise differentiation requires simulating the derivative of every rate with respect to every initial parameter, which is computationally expensive. The Adjoint Method (AEM) resolves this by solving the recurrence relation backward from the terminal payoff to the initial state, reducing the computational cost of delta and gamma calculations to a fraction of the direct method. This adjoint gradient vector is given by Eq. (24.25):

$$
\begin{equation*}
\Delta U=\frac{\partial U}{\partial \mathbf{L}\left(T_{k}\right)} \mathbf{D}\left(T_{k-1}\right) \ldots \mathbf{D}\left(T_{0}\right) \tag{24.25}
\end{equation*}
$$
[[Andersen_Ch24#^eq24-25|Eq. (24.25)]]

In summary, PDM combined with AEM is the preferred choice for risk-managing portfolios of vanilla and continuous exotic options due to its speed and low variance, while Malliavin calculus is reserved for options with hard discontinuous barriers and digital payoffs where PDM is inapplicable.

