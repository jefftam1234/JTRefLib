---
tags: [quant/knowledge, quant/books, volatility-modeling, stochastic-volatility]
---

# Stochastic Volatility Modeling — Lorenzo Bergomi

[Master Index](../index.md) · [Preface, Epilogue & Bibliography](Front and Back Matter.md)

**Author:** Lorenzo Bergomi
**Scope:** A practitioner-focused guide to the pricing, hedging, and risk management of volatility under local and stochastic volatility frameworks. The book prioritizes the **dynamics and break-even P&L** of option portfolios over abstract measure theory, building up from Black-Scholes carry-P&L accounting to local volatility, general forward-variance stochastic volatility models, and local-stochastic ("mixed") volatility models — always asking whether a model is *usable*: does it supply a genuine, payoff-independent break-even accounting for a hedged position's P&L?

---

### Chapter 1: Introduction

Why would a trader use a stochastic volatility model at all? The chapter's answer reframes the entire question: a pricing model earns the right to be used for risk management not by predicting the underlying's dynamics accurately, but by supplying a *usable accounting identity* for the P&L of a hedged position — one where a break-even condition on realized variance exists and is payoff-independent. Rather than a "wrong" physical assumption about lognormal dynamics, Black-Scholes is best understood as an accounting device.

Consider a trading desk pricing and risk-managing a short position in a European option of maturity $T$ using only the underlying $S$ as hedge. The first sanity check any pricing function $P(t,S)$ must pass is that it reproduces the payoff at maturity:

$$
\begin{equation*}
P(t=T, S)=f(S), \forall S \tag{1.1}
\end{equation*}
$$
[[Ch01#^eq1-1|Eq. (1.1)]]

Delta-hedging with $\Delta=dP/dS$ shares of the underlying and expanding the resulting carry P&L to order $\delta t$ and order two in $\delta S$ (the cross term $\delta S\,\delta t$ and higher orders in $\delta S$ are negligible as $\delta t\to0$, since $\delta S\sim\sqrt{\delta t}$) yields:

$$
\begin{equation*}
P \& L=-\left(\frac{d P}{d t}-r P+(r-q) S \frac{d P}{d S}\right) \delta t-\frac{1}{2} S^{2} \frac{d^{2} P}{d S^{2}}\left(\frac{\delta S}{S}\right)^{2} \tag{1.2}
\end{equation*}
$$
[[Ch01#^eq1-2|Eq. (1.2)]]

which is written compactly as $P\&L=-A(t,S)\delta t - B(t,S)(\delta S/S)^2$:

$$
\begin{equation*}
P \& L=-A(t, S) \delta t-B(t, S)\left(\frac{\delta S}{S}\right)^{2} \tag{1.3}
\end{equation*}
$$
[[Ch01#^eq1-3|Eq. (1.3)]]

If $A$ and $B$ have the same sign, the position systematically loses or makes money regardless of the realized path — a model with this property is unusable, since the price charged must have been mispriced. The model is *usable* only if $A$ and $B$ have opposite signs, $\forall t,\forall S$. Demanding further that the P&L vanish on average — with $\widehat\sigma$ the break-even volatility such that $\langle(\delta S/S)^2\rangle=\widehat\sigma^2\delta t$ — pins down the unique relationship $A=-\widehat\sigma^2 B$, and substituting $A,B$ back in gives the Black-Scholes PDE:

$$
\begin{equation*}
\frac{d P_{\widehat{\sigma}}}{d t}-r P_{\widehat{\sigma}}+(r-q) S \frac{d P_{\widehat{\sigma}}}{d S}=-\frac{\widehat{\sigma}^{2}}{2} S^{2} \frac{d^{2} P_{\widehat{\sigma}}}{d S^{2}} \tag{1.4}
\end{equation*}
$$
[[Ch01#^eq1-4|Eq. (1.4)]]

Substituting (1.4) back into the P&L expansion (1.2) gives the chapter's central accounting identity: the carry P&L is the dollar gamma times the gap between realized and implied (break-even) variance —

$$
\begin{equation*}
P \& L=-\frac{S^{2}}{2} \frac{d^{2} P_{\widehat{\sigma}}}{d S^{2}}\left(\frac{\delta S^{2}}{S^{2}}-\widehat{\sigma}^{2} \delta t\right) \tag{1.5}
\end{equation*}
$$
[[Ch01#^eq1-5|Eq. (1.5)]]

— a derivation that never assumed Brownian motion or a specific process for $S_t$; the probabilistic (diffusion) interpretation of (1.4) is a *consequence* of the break-even condition, not its starting point. This "market model" criterion generalizes cleanly to multiple hedge instruments. With a pricing function $P(t,S_1,\ldots,S_n)$ and cross-gamma matrix $\phi_{ij}=S_iS_j\,d^2P/dS_idS_j$, the multi-asset carry P&L reads:

$$
\begin{equation*}
P \& L=-A(t, S) \delta t-\frac{1}{2} \Sigma_{i j} \phi_{i j}(t, S) \frac{\delta S_{i}}{S_{i}} \frac{\delta S_{j}}{S_{j}} \tag{1.6}
\end{equation*}
$$
[[Ch01#^eq1-6|Eq. (1.6)]]

Diagonalizing $\phi$ into eigen-baskets and requiring $n$ positive break-even numbers $\omega_k$ so that

$$
\begin{equation*}
A=-\frac{1}{2} \Sigma_{k} \varphi_{k} \omega_{k} \tag{1.7}
\end{equation*}
$$
[[Ch01#^eq1-7|Eq. (1.7)]]

leads, after rewriting $A$ via the positive matrix $C=T\omega T^\top$, to the multi-asset market-model identity:

$$
\begin{equation*}
P \& L=-\frac{1}{2} \Sigma_{i j} \phi_{i j}\left(\frac{\delta S_{i}}{S_{i}} \frac{\delta S_{j}}{S_{j}}-C_{i j} \delta t\right) \tag{1.8}
\end{equation*}
$$
[[Ch01#^eq1-8|Eq. (1.8)]]

Any positive covariance matrix $C$ is admissible here — a model is usable, in the multi-asset case, whenever *some* positive break-even covariance matrix relates theta to the cross-gammas. Local volatility (Chapter 2) and forward variance models (Chapter 7) will be shown to satisfy this; most local-stochastic volatility models (Chapter 12) will not.

Having established *when* a model is well-defined, the chapter turns to whether delta-hedging alone is *good enough* in practice — and delivers a central paradox of the book. Even though continuous-time Black-Scholes delta-hedging drives the residual P&L in (1.5) to zero, real desks rehedge discretely (typically daily). Summing the daily P&Ls (1.9) over the option's life, decomposing each day's return as $r_i=\sigma_i\sqrt{\delta t}z_i$ (1.10) to separate the *scale* of the return from its *distribution*, and computing the variance of the resulting sum (1.11) — introducing the excess kurtosis $\kappa$ of daily returns and the variance/variance autocorrelation function $f_{ij}$ — yields, after converting to the vega via the Black-Scholes relation (1.12), a general formula for the standard deviation of the final P&L (1.13). In the idealized Black-Scholes case ($\Omega=0,\kappa=0$) this collapses to $\mathrm{StDev}(P\&L)=|\widehat\sigma\,dP_{\widehat\sigma}/d\widehat\sigma|/\sqrt{2N}$ (1.14) — about 4.5% of a typical option's vega for $N=250$ daily rehedges. In the real case, though, a second contribution survives the continuum limit (1.15), driven by the autocorrelation of daily variance itself. Fitting an exponential decay $f(\tau)=\rho e^{-k\tau}$ to real equity data ($\rho\simeq0.78$, $1/k\simeq45$ days) gives the chapter's headline quantitative result:

$$
\begin{equation*}
\frac{\operatorname{StDev}(P \& L)}{\left|\widehat{\sigma} \frac{d P_{\widehat{\sigma}}}{d \widehat{\sigma}}\right|} \simeq \sqrt{\frac{2+\kappa}{4 N}+\frac{\rho \Omega}{2} \frac{k T-1+e^{-k T}}{(k T)^{2}}} \tag{1.16}
\end{equation*}
$$
[[Ch01#^eq1-16|Eq. (1.16)]]

For a one-year at-the-money option this evaluates to roughly 35% of the option's premium — about eight times larger than the idealized Black-Scholes estimate of 4.5% — and comparison of the two terms shows the dispersion is dominated not by fat tails but by the long-ranged autocorrelation of volatility itself. Delta-hedging is necessary but nowhere near sufficient.

The natural fix — gamma-hedging with a single vanilla option $O$ of implied volatility $\widehat\sigma_O$ — does not eliminate this residual risk; it *relocates* it. The vanilla hedge's own P&L takes the same form as (1.5):

$$
\begin{equation*}
P \& L_{O}=-\frac{S^{2}}{2} \frac{d^{2} O}{d S^{2}}\left(\frac{\delta S^{2}}{S^{2}}-\widehat{\sigma}_{O}^{2} \delta t\right) \tag{1.17}
\end{equation*}
$$
[[Ch01#^eq1-17|Eq. (1.17)]]

and the hedge ratio needed to cancel the gamma of the exotic option $P$ is $\lambda=(d^2O/dS^2)^{-1}d^2P/dS^2$:

$$
\begin{equation*}
\lambda=\frac{1}{\frac{d^{2} O}{d S^{2}}} \frac{d^{2} P}{d S^{2}} \tag{1.18}
\end{equation*}
$$
[[Ch01#^eq1-18|Eq. (1.18)]]

Risk-managing $P$ and $O$ with two different implied volatilities produces an incongruous theta even when gamma vanishes — the same nonsensical situation as $A,B$ sharing a sign in Section 1.1 — so consistency forces $\widehat\sigma=\widehat\sigma_O$. But $\widehat\sigma_O$, unlike a freely-chosen constant $\widehat\sigma$, is a *market* implied volatility: as $\lambda$ is continually rebalanced, $O$ is bought and sold at prevailing prices, so $\widehat\sigma_O$ moves and cannot be held fixed. The vanilla hedge's own P&L, expanded to second order in $\delta\widehat\sigma_O$, acquires new terms:

$$
\begin{equation*}
P \& L_{O}=-\frac{S^{2}}{2} \frac{d^{2} O}{d S^{2}}\left(\frac{\delta S^{2}}{S^{2}}-\widehat{\sigma}_{O}^{2} \delta t\right)-\frac{d O}{d \widehat{\sigma}_{O}} \delta \widehat{\sigma}_{O}-\frac{1}{2} \frac{d^{2} O}{d \widehat{\sigma}_{O}^{2}} \delta \widehat{\sigma}_{O}^{2}-\frac{d^{2} O}{d S d \widehat{\sigma}_{O}} \delta S \delta \widehat{\sigma}_{O} \tag{1.19}
\end{equation*}
$$
[[Ch01#^eq1-19|Eq. (1.19)]]

and these vega/vanna/volga terms do not simply cancel between $P$ and $O$: since $d^2P/dS^2=\lambda\,d^2O/dS^2$ by construction, the gamma-hedged position's residual P&L is

$$
\begin{align*}
P \& L & =-\left(\frac{d P}{d \widehat{\sigma}_{O}}-\lambda \frac{d O}{d \widehat{\sigma}_{O}}\right) \delta \widehat{\sigma}_{O}  \tag{1.20}\\
& -\frac{1}{2}\left(\frac{d^{2} P}{d \widehat{\sigma}_{O}^{2}}-\lambda \frac{d^{2} O}{d \widehat{\sigma}_{O}^{2}}\right) \delta \widehat{\sigma}_{O}^{2}-\left(\frac{d^{2} P}{d S d \widehat{\sigma}_{O}}-\lambda \frac{d^{2} O}{d S d \widehat{\sigma}_{O}}\right) \delta S \delta \widehat{\sigma}_{O}
\end{align*}
$$
[[Ch01#^eq1-20|Eq. (1.20)]]

This is the chapter's central "relocated, not solved" result: unlike the $\delta S^2$ term in (1.5), no deterministic $\delta t$ term offsets these vega/vanna/volga contributions, because the Black-Scholes framework has no concept of implied-volatility *dynamics* at all — no implied volatility-of-volatility, no implied spot/volatility correlation. By using options to gamma-hedge, the desk lowers its exposure to realized volatility but acquires exposure to the dynamics of implied volatilities instead, and nothing in Black-Scholes prices that exposure. Supplying exactly this missing break-even machinery — for the *joint* dynamics of $S$ and its volatility surface — is the task of a stochastic volatility model.

Two closing examples make this abstract "exposure to implied-volatility dynamics" concrete. First, a barrier option (paying 1 at maturity unless $S$ hits a barrier $L$) can be statically hedged, in Black-Scholes, by a European payoff built from Carr-Chou reflection — for vanishing rates, essentially a double European digital struck at $L$, whose Black-Scholes value near $S=L$ is dominated by the local skew:

$$
\begin{equation*}
\mathcal{D}=2\left(\mathcal{D}_{L}^{B S}\left(\widehat{\sigma}_{L}\right)+\left.\frac{d \mathcal{P}_{L}^{B S}}{d \widehat{\sigma}} \frac{d \widehat{\sigma}_{K}}{d K}\right|_{L}\right)-1 \tag{1.24}
\end{equation*}
$$
[[Ch01#^eq1-24|Eq. (1.24)]]

Since the digital's Black-Scholes value at $S=L$ is close to 50% and nearly insensitive to $\widehat\sigma_L$ itself, the entire P&L on unwind is governed by the prevailing ATM skew at the moment the barrier is hit — a single, isolable feature of implied-volatility dynamics a usable stochastic volatility model must get right. Second, a forward-start ("cliquet") call paying $(S_{T_2}/S_{T_1}-k)^+$ at $T_2$:

$$
\begin{equation*}
\left(\frac{S_{T_{2}}}{S_{T_{1}}}-k\right)^{+} \tag{1.25}
\end{equation*}
$$
[[Ch01#^eq1-25|Eq. (1.25)]]

has, for $k=100\%$ and vanishing rates, the approximate Black-Scholes price

$$
\begin{equation*}
P \simeq \frac{1}{\sqrt{2 \pi}} \sigma \sqrt{T_{2}-T_{1}} \tag{1.26}
\end{equation*}
$$
[[Ch01#^eq1-26|Eq. (1.26)]]

which is *independent of spot* — a genuine oddity for the one instrument (the underlying) whose dynamics Black-Scholes actually models. The cliquet's real underlying is the forward implied volatility observed at $T_1$ for maturity $T_2$; it is, in effect, a pure option on volatility.

#### Chapter's Digest and Logical Bridge to Chapter 2

Delta-hedging cancels the order-one $\delta S$ term, but the dispersion of the resulting P&L — driven more by variance autocorrelation than by fat tails — remains large. Gamma-hedging with vanilla options immunizes against realized volatility but substitutes an unpriced exposure to the dynamics of implied volatilities, since Black-Scholes supplies no break-even mechanism for that exposure at all. Chapter 2 begins supplying one: the local volatility model is the simplest example of a genuine *market model* for the spot and vanilla options — one whose implied-volatility dynamics, however restricted, is at least internally consistent and priced.

---

### Chapter 2: Local volatility

Local volatility was introduced historically as a variant of Black-Scholes in which the instantaneous volatility is a deterministic function $\sigma(t,S)$ of time and spot, and it is still routinely described that way. Bergomi insists this framing is backwards: calibration to the observed vanilla smile *fully determines* the model, so the local volatility function only serves an ancillary purpose and has no physical significance — it is a by-product of the model possessing a one-dimensional Markov representation in $(t,S)$, not a modeling choice a trader makes on physical grounds. It is, in this sense, a *peculiar* market model: unlike a stochastic volatility model with free parameters to be estimated, its dynamics are not explicit and must be extracted after the fact, entirely from the shape of the calibrated smile. Traditional presentations calibrate $\sigma(t,S)$ once, at $t=0$, and freeze it — which then makes the common trading practice of *recalibrating* the local volatility function daily look like an improper use of the model. One of this chapter's central payoffs is showing that recalibration is in fact exactly how the model should be used, and that the resulting carry P&L has the standard market-model form of offsetting spot/volatility gamma and theta contributions with well-defined, payoff-independent break-even levels.

Before any of this can be established, the chapter needs the basic machinery relating local volatility to the market it is calibrated to. In the local volatility model, the spot follows

$$
\begin{equation*}
d S_{t}=(r-q) S_{t} d t+\sigma\left(t, S_{t}\right) S_{t} d W_{t} \tag{2.1}
\end{equation*}
$$
[[Ch02#^eq2-1|Eq. (2.1)]]

with $r$ the interest rate and $q$ the dividend-inclusive repo rate, and the corresponding pricing equation is the direct analogue of Black-Scholes with $\sigma(t,S)$ replacing the constant $\widehat\sigma$:

$$
\begin{equation*}
\frac{d P}{d t}+(r-q) S \frac{d P}{d S}+\frac{\sigma(t, S)^{2}}{2} S^{2} \frac{d^{2} P}{d S^{2}}=r P \tag{2.2}
\end{equation*}
$$
[[Ch02#^eq2-2|Eq. (2.2)]]

Given $\sigma(t,S)$, solving (2.2) backwards from maturity generates vanilla option prices; the harder question is the reverse one — given a full configuration of vanilla prices, does a function $\sigma(t,S)$ exist that reproduces them under (2.2)? Bruno Dupire's answer is the chapter's foundational identity, expressing local variance as a ratio of maturity and strike derivatives of call prices:

$$
\begin{equation*}
\sigma(t, S)^{2}=\left.2 \frac{\frac{d C}{d T}+q C+(r-q) K \frac{d C}{d K}}{K^{2} \frac{d^{2} C}{d K^{2}}}\right|_{\substack{K=S \\ T=t}} \tag{2.3}
\end{equation*}
$$
[[Ch02#^eq2-3|Eq. (2.3)]]

The derivation is worth tracing because it clarifies exactly what local volatility can and cannot claim to know. Writing an Itô expansion of the call payoff across an infinitesimal maturity step and taking expectations under an arbitrary — not necessarily deterministic — instantaneous volatility process $\sigma_t$ yields the general Dupire equation, which links the *conditional expectation* of instantaneous variance to option price derivatives:

$$
\begin{equation*}
E\left[\sigma_{T}^{2} \mid S_{T}=K\right]=\frac{E\left[\sigma_{T}^{2} \delta\left(S_{T}-K\right)\right]}{E\left[\delta\left(S_{T}-K\right)\right]}=2 \frac{\frac{d C}{d T}+q C+(r-q) K \frac{d C}{d K}}{K^{2} \frac{d^{2} C}{d K^{2}}} \tag{2.6}
\end{equation*}
$$
[[Ch02#^eq2-6|Eq. (2.6)]]

This identity holds in any diffusive model for $S_t$; knowledge of vanilla option prices is not sufficient to pin down the process $\sigma_t$ itself, but it characterizes the entire class of diffusive processes that yield the same vanilla smile. Local volatility is simply the simplest member of that class — the one obtained by setting $\sigma_t \equiv \sigma(t,S)$, a deterministic function, so that the conditional expectation on the left of (2.6) collapses to $\sigma(t,S)^2$ and (2.3) follows directly. For the Dupire formula to be usable, though, both the numerator and denominator of (2.3) must be positive, and this is where the chapter ties local volatility to no-arbitrage. The denominator is, up to discounting, the risk-neutral density of $S_T$:

$$
\begin{equation*}
\frac{d^{2} C(K, T)}{d K^{2}}=e^{-r T} E\left[\delta\left(S_{T}-K\right)\right] \tag{2.8}
\end{equation*}
$$
[[Ch02#^eq2-8|Eq. (2.8)]]

which is positive as soon as butterfly spreads — a static, model-free portfolio — are not arbitrageable; violation of this is called a strike arbitrage. The numerator is positive as soon as, for any moneyness $k$ and $T_1\le T_2$, forward-scaled call prices increase with maturity:

$$
\begin{equation*}
e^{q T_{1}} C\left(k F_{T_{1}}, T_{1}\right) \leq e^{q T_{2}} C\left(k F_{T_{2}}, T_{2}\right) \tag{2.9}
\end{equation*}
$$
[[Ch02#^eq2-9|Eq. (2.9)]]

Bergomi proves this via an explicit static-hedge argument — buying the $T_2$ option, selling a scaled quantity of the $T_1$ option, and entering a barrier-triggered delta position at $T_1$ — and shows the resulting P&L is strictly positive at $T_2$ by convexity of the payoff; a violation of (2.9) is thus arbitrageable and is called a maturity arbitrage. The conclusion carries a genuinely counter-intuitive payoff, stated bluntly in the text: contrary to a frequently heard assertion, steep skews observed for short-maturity equity smiles are not evidence that jumps are needed to generate them — as long as the smile is arbitrage-free, local volatility "will be happy to oblige."

Since (2.9) is a statement about prices, the next step translates it into a statement about implied volatilities directly. Writing the price of a fixed-moneyness option as a function of $\tau=\widehat\sigma_{kT}^2T$ through the standardized martingale $U_\tau$, and exploiting that its Black-Scholes forward PDE has a strictly positive risk-neutral density (hence a strictly increasing price in $\tau$), yields the equivalence between (2.9) and the convex order condition on integrated variance:

$$
\begin{equation*}
e^{q T_{1}} C_{B S}\left(k F_{T_{1}}, T_{1}, \widehat{\sigma}_{k T_{1}}\right) \leq e^{q T_{2}} C_{B S}\left(k F_{T_{2}}, T_{2}, \widehat{\sigma}_{k T_{2}}\right) \Leftrightarrow T_{1} \widehat{\sigma}_{k T_{1}}^{2} \leq T_{2} \widehat{\sigma}_{k T_{2}}^{2} \tag{2.14}
\end{equation*}
$$
[[Ch02#^eq2-14|Eq. (2.14)]]

— the pivotal technical step being that the one-directional monotonicity established first has to be strengthened into a genuine equivalence, "the reverse implication," before it can be used as a design constraint on interpolation schemes. The same convex-order logic is not tied to hockey-stick payoffs: for any family of European payoffs $h(S_T/F_T)$ with $h$ convex, absence of arbitrage forces the analogous condition on the corresponding lognormal implied volatilities,

$$
\begin{equation*}
T_{2} \widehat{\sigma}_{T_{2}}^{2} \geq T_{1} \widehat{\sigma}_{T_{1}}^{2} \tag{2.17}
\end{equation*}
$$
[[Ch02#^eq2-17|Eq. (2.17)]]

a generalization later reused in Chapter 4's power-payoff framework. To turn all of this into a computational recipe, the chapter adopts a $(y,t)$ coordinate system tailored to how vanilla smiles are quoted,

$$
\begin{align*}
y & =\ln \left(\frac{K}{F_{t}}\right)  \tag{2.18a}\\
f(t, y) & =\left(t-t_{0}\right) \widehat{\sigma}_{K t}^{2} \tag{2.18b}
\end{align*}
$$
[[Ch02#^eq2-18a|Eq. (2.18a)]] [[Ch02#^eq2-18b|Eq. (2.18b)]]

and substitutes the Black-Scholes formula for $C$ into (2.3), computing all derivatives analytically, to arrive at the Dupire formula expressed directly in terms of the smile:

$$
\begin{equation*}
\sigma(t, S)^{2}=\left.\frac{\frac{d f}{d t}}{\left(\frac{y}{2 f} \frac{d f}{d y}-1\right)^{2}+\frac{1}{2} \frac{d^{2} f}{d y^{2}}-\frac{1}{4}\left(\frac{1}{4}+\frac{1}{f}\right)\left(\frac{d f}{d y}\right)^{2}}\right|_{y=\ln \left(\frac{S}{F_{t}}\right)} \tag{2.19}
\end{equation*}
$$
[[Ch02#^eq2-19|Eq. (2.19)]]

Because market quotes only exist at discrete strikes and maturities, using (2.19) in practice requires an interpolation/extrapolation of the smile that itself respects the convex-order condition — the simple rule that successive $f_i(y)$ profiles must not cross, enforced here by a piecewise-affine interpolation in maturity, guarantees the positivity (2.19) needs without introducing spurious maturity-crossing sensitivities. Cash-amount dividends complicate this picture, since option prices then depart from the pure Black-Scholes formula; the chapter's exact solution maps the dividend-paying asset onto a driftless proxy $X$,

$$
\begin{equation*}
S_{t}=\alpha(t) X_{t}-\delta(t) \tag{2.21}
\end{equation*}
$$
[[Ch02#^eq2-21|Eq. (2.21)]]

whose own smile can be Dupire-inverted directly via (2.19), with the local volatility of $S$ recovered as

$$
\begin{equation*}
\sigma(t, S)=\frac{S+\delta(t)}{S} \sigma_{X}(t, X(S, t)) \tag{2.22}
\end{equation*}
$$
[[Ch02#^eq2-22|Eq. (2.22)]]

Across dividend dates $\sigma_X$ stays continuous while $\sigma$ jumps, because $\delta(t)$ jumps — an artifact some readers taking local volatility "seriously" as a physical volatility might object to, but which the text defends on the grounds that just before a dividend date the cash portion of $S$ carries no volatility, so the lognormal volatility of $S$ has no reason to be continuous either. Because this exact mapping is numerically heavier than desks typically want, Section 2.3.1.2 develops a Bos-Vandermark-style approximate scheme that achieves the same cross-dividend matching condition using only a modified moneyness coordinate, without ever needing the exact $X$-mapping machinery.

With local volatility now expressible from the smile, the chapter turns the question around: given a local volatility function, how do implied volatilities respond as spot moves? The answer starts from a completely general delta-hedging identity. Pricing a vanilla option in a base "model I" but letting the true dynamics of $S_t$ be governed by an arbitrary "model II" instantaneous volatility $\sigma_{2t}$, and tracking the discounted gamma/theta mismatch accumulated while delta-hedging with model I, gives

$$
\begin{equation*}
P_{2}\left(0, S_{0}, \bullet\right)=P_{1}\left(0, S_{0}\right)+E_{2}\left[\int_{0}^{T} e^{-r t} \frac{S_{t}^{2}}{2} \frac{d^{2} P_{1}}{d S^{2}}\left(\sigma_{2 t}^{2}-\sigma_{1}\left(t, S_{t}\right)^{2}\right) d t\right] \tag{2.30}
\end{equation*}
$$
[[Ch02#^eq2-30|Eq. (2.30)]]

Specializing model I to Black-Scholes at the (unknown) implied volatility $\widehat\sigma_{KT}$ and model II to the local volatility model itself turns this generic identity into the chapter's central "weighted-average" representation: implied variance is the dollar-gamma-weighted average of local variance, over paths generated by the local volatility model, with the weighting gamma itself computed at the implied volatility being solved for:

$$
\begin{equation*}
\widehat{\sigma}_{K T}^{2}=\frac{E_{\sigma(t, S)}\left[\int_{0}^{T} e^{-r t} S_{t}^{2} \frac{d^{2} P_{\hat{\sigma}_{K T}}}{d S^{2}} \sigma(t, S)^{2} d t\right]}{E_{\sigma(t, S)}\left[\int_{0}^{T} e^{-r t} S_{t}^{2} \frac{d^{2} P_{\hat{\sigma}_{K T}}}{d S^{2}} d t\right]} \tag{2.32}
\end{equation*}
$$
[[Ch02#^eq2-32|Eq. (2.32)]]

This formula is implicit — $\widehat\sigma_{KT}$ appears on both sides — so Section 2.4.2 linearizes it around a Black-Scholes base case, expanding at order one in the local variance perturbation $\delta u(t,S)=\sigma(t,S)^2-\sigma_0^2(t)$ and using the martingale property of the discounted dollar gamma to collapse the resulting double integral into a clean Gaussian-kernel average. Specializing further to an expansion around a *constant* base volatility $\sigma_0$ produces the chapter's most-used working formula:

$$
\begin{equation*}
\widehat{\sigma}_{K T}=\frac{1}{T} \int_{0}^{T} d t \int_{-\infty}^{+\infty} d y \frac{e^{-\frac{y^{2}}{2}}}{\sqrt{2 \pi}} \sigma\left(t, F_{t} e^{\frac{t}{T} x_{K}+\sigma_{0} \sqrt{\frac{(T-t) t}{T}} y}\right) \tag{2.42}
\end{equation*}
$$
[[Ch02#^eq2-42|Eq. (2.42)]]

Bergomi is unusually blunt about this formula's limits: however appealing (2.42) may be as an order-one approximation, market smiles on equity underlyings are strong enough, and bid-offer spreads tight enough, that its numerical accuracy is not adequate for real trading applications — the naive fix of using a lognormal density with a different implied volatility per time slice ("most-likely-path" tricks) does not solve the problem, because the dependence of the averaging density itself on $\sigma(t,S)$ is not a second-order effect that can be neglected. What survives is *relative* use of (2.42): computing the skew and curvature of the smile near the forward. Assuming a locally quadratic local volatility function in log-moneyness $x=\ln(S/F_t)$,

$$
\begin{equation*}
\sigma(t, S)=\bar{\sigma}(t)+\alpha(t) x+\frac{\beta(t)}{2} x^{2} \tag{2.44}
\end{equation*}
$$
[[Ch02#^eq2-44|Eq. (2.44)]]

and carrying out the Gaussian integral in (2.42) at order one in $\alpha,\beta$ produces the ATMF (at-the-money-forward) skew and curvature as time-weighted averages of the local volatility function's own skew and curvature:

$$
\begin{align*}
\mathcal{S}_{T}=\left.\frac{d \widehat{\sigma}_{K T}}{d \ln K}\right|_{F_{T}} & =\frac{1}{T} \int_{0}^{T} \frac{t}{T} \alpha(t) d t  \tag{2.48}\\
\left.\frac{d^{2} \widehat{\sigma}_{K T}}{d \ln K^{2}}\right|_{F_{T}} & =\frac{1}{T} \int_{0}^{T}\left(\frac{t}{T}\right)^{2} \beta(t) d t \tag{2.49}
\end{align*}
$$
[[Ch02#^eq2-48|Eq. (2.48)]] [[Ch02#^eq2-49|Eq. (2.49)]]

For constant $\alpha,\beta$ these collapse to the well-known "half" and "third" scaling rules ($\mathcal{S}_T=\alpha/2$, curvature $=\beta/3$); for a power-law-decaying local skew $\alpha(t)=\alpha_0(\tau_0/t)^\gamma$ (with a short-maturity cutoff $\tau_0$), the long-maturity ATMF skew inherits the same decay exponent, rescaled by $1/(2-\gamma)$:

$$
\begin{equation*}
\mathcal{S}_{T} \simeq \frac{1}{2-\gamma} \alpha_{0}\left(\frac{\tau_{0}}{T}\right)^{\gamma} \tag{2.53}
\end{equation*}
$$
[[Ch02#^eq2-53|Eq. (2.53)]]

matching the $\gamma\simeq1/2$ decay typically observed in equity index smiles. All of this, though, is a perturbative story that breaks down as $T\to0$, since the whole Gaussian-kernel machinery is built on a first-order expansion. For the true short-maturity limit, Section 2.4.6 derives an *exact*, non-perturbative result directly from (2.19): taking $T\to0$ and keeping only the leading powers of $T$ in numerator and denominator collapses the Dupire relation to a simple ODE in $y$, whose solution shows that the inverse of implied volatility is the *harmonic* average — not the usual quadratic-variation-style average — of the inverse of local volatility along the path from spot to strike:

$$
\begin{equation*}
\frac{1}{\widehat{\sigma}\left(T=0, S e^{y}\right)}=\frac{1}{y} \int_{0}^{y} \frac{d u}{\sigma\left(T=0, S e^{u}\right)} \tag{2.54}
\end{equation*}
$$
[[Ch02#^eq2-54|Eq. (2.54)]]

The text flags this as genuinely surprising at first sight, since squared-volatility averages are the norm everywhere else in the book; the harmonic form appears here specifically because there is no time-averaging left once maturity vanishes, and it has the reassuring property that a region of zero local volatility between spot and strike forces the implied volatility for that strike to vanish too, since the spot could never cross it in the $T\to0$ limit anyway.

Having pinned down the *static* map between local and implied volatility, the chapter turns to *dynamics*: since implied volatility in this model is a deterministic function of $(t,S)$ alone, its response to a spot move is entirely dictated by the calibrated smile — no separate assumption about the volatility of volatility is available to the user. This is, in fact, the source of one of the chapter's deepest points, flagged explicitly in the text and only resolved in Section 2.7: the volatilities of volatilities and spot/volatility covariances that fall out of this analysis are exactly the break-even levels of the P&L of a delta- and vega-hedged position, an aspect that looks counter-intuitive at first since nobody chose those levels — they are a byproduct of calibration. Differentiating (2.42) with respect to $\ln S_0$ and specializing to the strike equal to the forward gives the rate at which the ATMF implied volatility responds to a spot move, which for a maturity-independent skew reduces to twice the skew itself — the first appearance of a "factor of two" that recurs throughout the section. Normalizing this sensitivity by the ATMF skew of the same maturity defines the Skew Stickiness Ratio,

$$
\begin{equation*}
\mathcal{R}_{T}=\frac{1}{\mathcal{S}_{T}} \frac{d \widehat{\sigma}_{F_{T} T}}{d \ln S_{0}} \tag{2.61}
\end{equation*}
$$
[[Ch02#^eq2-61|Eq. (2.61)]]

which quantifies, in units of the local skew, how much of the spot/implied-volatility covariance the model actually generates: $\mathcal{R}_T=1$ is the "sticky-strike" regime where fixed-strike implied volatilities stay frozen as spot moves, while $\mathcal{R}_T=0$ is "sticky-delta," where the whole smile translates with the spot. Substituting the ATMF sensitivity formula into this definition gives the order-one SSR formula for the local volatility model,

$$
\begin{equation*}
\mathcal{R}_{T}=1+\frac{1}{T} \int_{0}^{T} \frac{\mathcal{S}_{t}}{\mathcal{S}_{T}} d t \tag{2.64}
\end{equation*}
$$
[[Ch02#^eq2-64|Eq. (2.64)]]

which for strong real-world skews tends to overestimate the true SSR, since it omits higher-order contributions from the local skew. The formula's most quoted consequence, though, is exact in a limiting sense that goes beyond the order-one expansion used to derive it. When the ATMF skew does not depend on maturity, (2.64) collapses to

$$
\begin{equation*}
\mathcal{R}_{T}=2 \tag{2.66}
\end{equation*}
$$
[[Ch02#^eq2-66|Eq. (2.66)]]

and the same value is reached in the $T\to0$ limit for *any* smooth skew term structure, since as $T\to0$ the ATMF sensitivity and the ATMF skew both converge to the same short-maturity slope of the smile at the forward:

$$
\begin{equation*}
\lim _{T \rightarrow 0} \frac{d \widehat{\sigma}_{F_{T} T}}{d \ln S_{0}}=\left.2 \lim _{T \rightarrow 0} \frac{d \widehat{\sigma}_{K T}}{d \ln S_{0}}\right|_{K=F_{T}}=\left.2 \lim _{T \rightarrow 0} \frac{d \widehat{\sigma}_{K T}}{d \ln K}\right|_{K=F_{T}} \tag{2.67}
\end{equation*}
$$
[[Ch02#^eq2-67|Eq. (2.67)]]

giving the short-maturity "R=2 rule":

$$
\begin{equation*}
\lim _{T \rightarrow 0} \mathcal{R}_{T}=2 \tag{2.68}
\end{equation*}
$$
[[Ch02#^eq2-68|Eq. (2.68)]]

Section 2.5.3 then shows both statements are not artifacts of the order-one expansion at all but are *exactly* true for any time-independent local volatility function $\sigma(S/F_t)$. The proof exploits a backward/forward pricing symmetry: rewriting the pricing PDE in rescaled moneyness coordinates shows call and put prices at time zero obey mirror-image forward equations in $(s,k)$ versus $(k,s)$, which — after invoking the Black-Scholes identity $f_{BS}(\tau s;k;\Sigma)=(s-k)+f_{BS}(\tau k;s;\Sigma)$ and call/put parity — forces the two implied-volatility surfaces indexed by spot and strike to coincide exactly. Differentiating that symmetry with respect to $\ln K$ at $K=F_T$, $t=0$ yields the exact identity

$$
\begin{equation*}
\left.\frac{d \widehat{\sigma}_{K T}}{d \ln K}\right|_{K=F_{T}}=\left.\frac{d \widehat{\sigma}_{K T}}{d \ln S_{0}}\right|_{K=F_{T}} \tag{2.79}
\end{equation*}
$$
[[Ch02#^eq2-79|Eq. (2.79)]]

from which $\mathcal{R}_T=2$ follows identically, with no perturbative approximation anywhere in the argument. Section 2.5.4 completes the picture by showing that once the skew's power-law decay is not idealized down to a vanishing cutoff, the long-maturity SSR instead settles at $(2-\gamma)/(1-\gamma)$ (2.80) — so $\mathcal{R}=2$ is specifically a short-maturity, model-independent statement about local (and, as Chapter 9 will show, stochastic) volatility models, not a universal constant true at every maturity.

The next question is how these skews evolve at *future* dates, since Section 1.3.1 already showed barrier option prices are dominated by the local skew at the barrier for the residual maturity as generated by whichever model prices them. Re-expressing the quadratic local volatility parametrization around a future reference spot $S_\tau$ at date $\tau$ and using (2.48) gives the ATMF skew for residual maturity $\theta$ observed at $\tau$:

$$
\begin{equation*}
\mathcal{S}_{\theta}(\tau)=\left.\frac{d \widehat{\sigma}_{K \tau+\theta}\left(S_{\tau}, \tau\right)}{d \ln K}\right|_{K=F_{\tau+\theta}\left(S_{\tau}\right)}=\frac{1}{\theta} \int_{\tau}^{\tau+\theta} \frac{t-\tau}{\theta} \alpha(t) d t \tag{2.90}
\end{equation*}
$$
[[Ch02#^eq2-90|Eq. (2.90)]]

which, rewritten purely in terms of the initial vanilla smile's term structure of ATMF skews, gives the chapter's forward-skew formula:

$$
\begin{equation*}
\mathcal{S}_{\theta}(\tau)=\mathcal{S}_{\tau+\theta}-\frac{\tau}{\theta}\left(\frac{1}{\theta} \int_{\tau}^{\tau+\theta} \mathcal{S}_{t} d t-\mathcal{S}_{\tau+\theta}\right) \tag{2.91}
\end{equation*}
$$
[[Ch02#^eq2-91|Eq. (2.91)]]

This is a second deliberately counter-intuitive result: for a typical, decreasing term structure of ATMF skew, the forward-starting skew generated by the model is much weaker than the spot-starting skew of the same residual maturity, $|\mathcal{S}_\theta(\tau)|\ll|\mathcal{S}_\theta|$ — and for a power-law skew this weakening scales explicitly as $(\theta/\tau)^\gamma$ (2.92). Since local volatility is not time-homogeneous, its own calibration-day smile silently encodes assumptions about future skews (and hence future volatility of volatility, since that too is pinned by the skew) that will change unpredictably every time the model is recalibrated — a mispricing of forward-smile risk that Chapter 3's treatment of forward-start options and Chapter 12's comparison of local, stochastic, and local-stochastic smiles both return to.

This sets up the chapter's reckoning with whether local volatility is a *legitimate* model to trade at all. The naive candidate delta is the "local volatility delta," computed with the local volatility function held fixed:

$$
\begin{equation*}
\Delta^{\mathrm{LV}}=\left.\frac{d P^{\mathrm{LV}}}{d S}\right|_{\sigma} \tag{2.93}
\end{equation*}
$$
[[Ch02#^eq2-93|Eq. (2.93)]]

but the P&L this delta generates only has the clean market-model form of (2.94) if market implied volatilities move exactly as $\sigma(t,S)$ prescribes — and in reality, implied volatilities will move however they wish, so the local volatility function recalibrated at $t+\delta t$ will differ from the one calibrated at $t$, and the naive P&L formula (2.94) is not the desk's actual P&L. Sections 2.7.2–2.7.3 resolve this by changing state variables: representing the pricing function directly in terms of the traded implied volatilities $\widehat\sigma_{KT}$ (or, equivalently, vanilla option prices $O_{KT}$) rather than the frozen local volatility function. Differentiating the local volatility price with respect to $S$ while holding the smile fixed defines the sticky-strike delta,

$$
\begin{equation*}
\Delta^{\mathrm{SS}}=\left.\frac{d P}{d S}\right|_{\widehat{\sigma}_{K T}} \tag{2.98}
\end{equation*}
$$
[[Ch02#^eq2-98|Eq. (2.98)]]

and a careful accounting of the P&L over $\delta t$ — expanding the option's total change to second order in $\delta S$ and $\delta\widehat\sigma_{KT}$, then using the local volatility pricing PDE to eliminate the time-derivative term and isolate a well-defined drift $\mu_{KT}$ for implied volatility — produces the full carry P&L decomposition (2.97)–(2.106), with $\Delta^{\mathrm{SS}}$ as its linear spot sensitivity. Replacing implied volatilities with vanilla option prices $O_{KT}$ as state variables removes the model-dependence of the drift entirely (an asset's drift is simply its financing cost) and yields the chapter's central proof of legitimacy: the carry P&L of a delta- and vega-hedged local volatility position, expressed in option prices, is

$$
\begin{align*}
& P \& L= \\
& -\frac{d \mathcal{P}}{d S}(\delta S-(r-q) S \delta t)-\frac{d \mathcal{P}}{d O_{K T}} \bullet\left(\delta O_{K T}-r O_{K T} \delta t\right)  \tag{2.107a}\\
& -\frac{1}{2} \frac{d^{2} \mathcal{P}}{d S^{2}}\left[\delta S^{2}-\sigma^{2}(t, S) S^{2} \delta t\right]  \tag{2.107b}\\
& -\frac{d^{2} \mathcal{P}}{d S d O_{K T}} \bullet\left[\delta S \delta O_{K T}-\sigma^{2}(t, S) S^{2} \frac{d \Omega_{K T}^{L V}}{d S} \delta t\right]  \tag{2.107c}\\
& -\frac{1}{2} \frac{d^{2} \mathcal{P}}{d O_{K T} d O_{K^{\prime} T^{\prime}}} \bullet\left[\delta O_{K T} \delta O_{K^{\prime} T^{\prime}}-\sigma^{2}(t, S) S^{2} \frac{d \Omega_{K T}^{L V}}{d S} \frac{d \Omega_{K^{\prime} T^{\prime}}^{L V}}{d S} \delta t\right] \tag{2.107d}
\end{align*}
$$
[[Ch02#^eq2-107a|Eq. (2.107a)]] [[Ch02#^eq2-107b|Eq. (2.107b)]] [[Ch02#^eq2-107c|Eq. (2.107c)]] [[Ch02#^eq2-107d|Eq. (2.107d)]]

This is exactly the market-model form of Section 1.1: every gamma and cross-gamma — spot/spot, spot/option, option/option — is paired with an offsetting theta at a payoff-independent break-even level, vindicating local volatility, in the chapter's own words, as "the simplest market model," despite $\sigma(t,S)$ itself carrying no physical meaning and despite the model being used, unapologetically, with daily recalibration. Equation (2.107) also settles which delta actually matters: the natural market-model delta, obtained by moving $S$ while holding all vanilla option prices fixed, is

$$
\begin{equation*}
\Delta^{\mathrm{MM}}=\left.\frac{d \mathcal{P}}{d S}\right|_{O_{K T}} \tag{2.108}
\end{equation*}
$$
[[Ch02#^eq2-108|Eq. (2.108)]]

and $\Delta^{\mathrm{SS}}$ turns out to be nothing more than $\Delta^{\mathrm{MM}}$ augmented by the Black-Scholes deltas of the vanilla options used as vega hedges — it has no special status of its own. That local volatility genuinely is a diffusive market model, and not merely coincidentally reproducing (2.107)'s form, is confirmed directly in Section 2.7.5 by exhibiting joint SDEs for $S_t$ and the vanilla prices $O_{KT,t}$ themselves (2.109)–(2.111) whose solution, when the local volatility function is held frozen inside them, is exactly the local volatility model — the carry P&L (2.107) then follows automatically from the model being diffusive, without repeating the derivation of Section 2.7.2. Section 2.7.6's two-underlying metaphor makes the earlier point about $\Delta^{\mathrm{LV}}$ vivid: an ATM basket option on two perfectly correlated assets $S_1,S_2$ admits a one-dimensional Markov representation in $S_1$ alone, exactly as local volatility does in $(t,S)$, and differentiating that reduced representation gives a "delta" with respect to $S_1$ of roughly 50% and a delta with respect to $S_2$ of exactly zero — obviously nonsensical, since a real move in $S_2$ is entirely unhedged. The local volatility delta $\Delta^{\mathrm{LV}}$ is the same kind of artifact: any move of $\widehat\sigma_{KT}$ not already implied by the frozen $\sigma(t,S)$ passes straight through it unhedged. Trading $\Delta^{\mathrm{LV}}$ instead of $\Delta^{\mathrm{MM}}$ or $\Delta^{\mathrm{SS}}$ does make the residual vega P&L (2.112)–(2.113) vanish identically with probability one, since $S$ and $\widehat\sigma_{KT}$ are perfectly correlated inside the model — but the text is explicit that this has no practical relevance, because in real markets $S$ and implied volatility are separate instruments that move independently of any model's internal correlation assumption. The chapter's remaining qualification on legitimacy is that the volatilities of volatilities and spot/volatility correlations the model generates are neither chosen by the user nor even constant — they are dictated, unpredictably, by whatever the market smile happens to look like at each future recalibration date, an issue with no analogue in the exogenously parameterized stochastic volatility models of Chapter 7.

A short digression clarifies why payoff-independence of break-even levels matters operationally. Risk-managing a short 90/110 call spread by delta-hedging each leg with its *own* Black-Scholes implied volatility gives a P&L

$$
\begin{equation*}
P \& L=-\frac{1}{2} S^{2} \frac{d^{2} P_{B S}^{1}}{d S^{2}}\left(\frac{\delta S^{2}}{S^{2}}-\widehat{\sigma}_{K_{1} T}^{2} \delta t\right)+\frac{1}{2} S^{2} \frac{d^{2} P_{B S}^{2}}{d S^{2}}\left(\frac{\delta S^{2}}{S^{2}}-\widehat{\sigma}_{K_{2} T}^{2} \delta t\right) \tag{2.114}
\end{equation*}
$$
[[Ch02#^eq2-114|Eq. (2.114)]]

which, at the particular spot level where the two legs' gammas cancel exactly, still pays out "free" theta from the skew $(\widehat\sigma_{K_1T}>\widehat\sigma_{K_2T})$ in a region with no gamma risk whatsoever — a payoff-dependent break-even level with no economic justification. Risk-managing both legs inside the *same* local volatility model instead gives

$$
\begin{align*}
P \& L & =-\frac{1}{2} S^{2} \frac{d^{2} P_{\sigma}^{1}}{d S^{2}}\left(\frac{\delta S^{2}}{S^{2}}-\sigma(t, S)^{2} \delta t\right)+\frac{1}{2} S^{2} \frac{d^{2} P_{\sigma}^{2}}{d S^{2}}\left(\frac{\delta S^{2}}{S^{2}}-\sigma(t, S)^{2} \delta t\right) \\
& =-\frac{1}{2} S^{2}\left(\frac{d^{2} P_{\sigma}^{1}}{d S^{2}}-\frac{d^{2} P_{\sigma}^{2}}{d S^{2}}\right)\left(\frac{\delta S^{2}}{S^{2}}-\sigma(t, S)^{2} \delta t\right) \tag{2.115}
\end{align*}
$$
[[Ch02#^eq2-115|Eq. (2.115)]]

in which cancellation of the net gamma forces cancellation of the associated theta as well — a much more defensible outcome, and a direct illustration of why a single, model-consistent break-even level beats a payoff-dependent patchwork of separately-quoted implied volatilities.

The chapter closes with the machinery that turns all of this into a tradeable hedge. Since the vanilla hedge $d\mathcal{P}/dO_{KT}$ immunizes the position against any perturbation $\delta\sigma^2(t,S)$ of the local volatility surface, the practical question — first solved by Dupire and later by Dupire and Henry-Labordère — is which portfolio of vanilla options reproduces a given exotic option's dollar-gamma exposure to that perturbation. Perturbing the pricing PDE and applying the Feynman-Kac theorem expresses the resulting price change as an average of $\delta\sigma^2$ weighted by the density and by $\phi(t,S)$, the conditional expectation of the dollar gamma. Requiring a continuous portfolio $\Pi=\int d\tau\,dK\,\mu(\tau,K)\,C_{K\tau}$ of vanilla options to reproduce this exposure exactly reduces to matching dollar gammas through the linear operator

$$
\begin{equation*}
\mathcal{L} f=\frac{d f}{d t}+(r-q) S \frac{d f}{d S}+\frac{1}{2} S^{2} \frac{d^{2}}{d S^{2}}\left(\sigma^{2}(t, S) f\right)-r f \tag{2.120}
\end{equation*}
$$
[[Ch02#^eq2-120|Eq. (2.120)]]

which annihilates the dollar gamma of any European option away from its strike, and whose action on a discrete vanilla portfolio concentrates entirely at the Dirac-delta-like jumps in $\phi_\Pi$ at each option's maturity. Matching this discontinuity to a continuous density of hedging options yields the chapter's practical vega-hedge formula:

$$
\begin{equation*}
\mu(\tau, K)=-\frac{1}{K^{2}} \mathcal{L} \phi(\tau, K) \tag{2.121}
\end{equation*}
$$
[[Ch02#^eq2-121|Eq. (2.121)]]

Computing $\phi$ itself for a path-dependent exotic requires Malliavin-style transition-density weights (2.122)–(2.123), which are numerically delicate — noisy whenever spot observation dates are closely spaced — but present no conceptual difficulty. Once the hedge portfolio $\mu(\tau,K)$ is known, quoting the exotic option's price at a chosen flat volatility $\widehat\sigma_0$ and then re-marking the hedge at *market* vanilla prices rather than model prices produces the calibration/mark-to-market interpretation of the final quoted price:

$$
\begin{equation*}
P=P^{0}+\int_{0}^{T} d \tau \int_{0}^{\infty} d K \mu(\tau, K)\left(C_{K \tau}-C_{K \tau}^{0}\right) \tag{2.124}
\end{equation*}
$$
[[Ch02#^eq2-124|Eq. (2.124)]]

so that the mismatch between market and model vanilla prices is passed on to the client as an explicit hedging cost — the price produced by a calibrated model is only as credible as the hedge it implies. Section 2.10 briefly identifies Markov-functional models as a special case of local volatility, admitting a representation $S_t=f(t,W_t)$ calibrated to only a single maturity's smile (2.125)–(2.128), which is enough for exotics whose payoff only depends on a single observation date but not for anything genuinely path-dependent, since the model then offers no separate handle on the dynamics between that one calibrated date and any other; Appendix A develops the Uncertain Volatility Model as a related but distinct construction, driven by a trading criterion (a bound on realized volatility) rather than by smile calibration, and is a self-contained supplement outside this chapter's core Dupire/no-arbitrage/dynamics/market-model arc.

A useful way to summarize the chapter's arc is that local volatility answers two separate questions with two separate levels of rigor: the *static* question (what local variance function reproduces a given vanilla smile?) is answered exactly by the Dupire formula, subject only to the standard no-arbitrage conditions on the input smile; the *dynamic* question (what joint dynamics of spot and implied volatility does that same local variance function generate?) is answered only approximately, via the Gaussian-kernel expansion of Section 2.4 and its various specializations, except in the two genuinely exact corners of short maturity (the harmonic-average and $\mathcal R=2$ results) and the market-model carry P&L itself. This split between an exact static calibration and an approximate, entirely calibration-derived dynamics is precisely what later chapters use as the point of contrast: forward variance models (Chapter 7) invert the split, sacrificing exact calibration for direct, exogenous control of the dynamics.

#### Chapter's Digest and Logical Bridge to Chapter 3

Local volatility earns its keep as "the simplest market model" precisely because its one-dimensional Markov representation in $(t,S)$, while stripping $\sigma(t,S)$ of any physical meaning, also guarantees that the daily-recalibrated carry P&L of a delta- and vega-hedged position has the clean, payoff-independent gamma/theta form that Chapter 1 identified as the hallmark of a usable model. The price of this frugality is that every dynamical property the model generates — the ATMF skew, the SSR, the volatility of volatility, the shape of future forward-starting skews — is dictated entirely by the calibrated smile and cannot be controlled independently, and because local volatility is not time-homogeneous, these implied dynamical assumptions silently drift every time the model is recalibrated to a new smile. Chapter 3 picks up exactly this thread, examining forward-start options — instruments whose value depends directly on the future, forward-starting smile that Section 2.6 showed local volatility systematically underprices in variability — and develops the machinery of forward-start option pricing and the power-payoff/log-contract framework needed to analyze them, extending the convex-order results of Section 2.2.2.2 to a new class of payoffs.

---

### Chapter 3: Forward-start options

Forward-start options — cliquets — pay a function of the ratio of two future spot observations, $g(S_{T_2}/S_{T_1})$, and the chapter's opening move is to introduce a tempting shortcut and then discard it: a "forward smile" $\widehat\sigma_k^{T_1T_2}$, obtained by pricing the cliquet in some model and back-implying a Black-Scholes volatility from it. This concept is rejected almost as soon as it is introduced — it has no historical counterpart, conflates two genuinely distinct risks, and turns out to be of limited practical use — setting a skeptical tone that runs through the whole chapter.

Because vanilla smiles of maturity $T_1$ and $T_2$ only price payoffs that depend on $S_{T_1}$ or $S_{T_2}$ separately, it is not even obvious a priori what information, if any, in the two vanilla smiles is relevant to a payoff that genuinely couples both dates together. This is the chapter's real governing question, and the two-part structure of what follows — first a model-independent replication analysis (Section 3.1), then a model-specific critique of local volatility (Section 3.2) — is organized around answering it from two different directions: what can be said about hedging a cliquet without committing to any particular model at all, and what goes wrong once a specific, familiar model is asked to do the job.

Working first in a Black-Scholes setting with deterministic time-dependent volatility, homogeneity in $S$ means the cliquet's price depends on nothing but the integrated variance between $T_1$ and $T_2$ — a single number, the forward volatility $\widehat\sigma_{T_1T_2}$:

$$
\begin{equation*}
P=e^{-r\left(T_{1}-t\right)} G\left(\widehat{\sigma}_{T_{1} T_{2}}\right) \tag{3.1}
\end{equation*}
$$
[[Ch03#^eq3-1|Eq. (3.1)]]

defined as the residual integrated variance between the two maturities:

$$
\begin{equation*}
\widehat{\sigma}_{T_{1} T_{2}}^{2}=\frac{1}{T_{2}-T_{1}} \int_{T_{1}}^{T_{2}} \sigma(t)^{2} d t=\frac{\left(T_{2}-t\right) \widehat{\sigma}_{T_{2}}^{2}-\left(T_{1}-t\right) \widehat{\sigma}_{T_{1}}^{2}}{T_{2}-T_{1}} \tag{3.2}
\end{equation*}
$$
[[Ch03#^eq3-2|Eq. (3.2)]]

Two questions follow immediately: how to hedge this option, and which vanilla payoffs $\widehat\sigma_{T_1},\widehat\sigma_{T_2}$ should even be calibrated to. This second question is not a formality — implied volatilities can be defined for any convex payoff, but only if the corresponding volatilities satisfy the same convex order condition established in Chapter 2:

$$
\begin{equation*}
\left(T_{2}-t\right) \widehat{\sigma}_{T_{2}}^{2} \geq\left(T_{1}-t\right) \widehat{\sigma}_{T_{1}}^{2} \tag{3.3}
\end{equation*}
$$
[[Ch03#^eq3-3|Eq. (3.3)]]

so (3.2)'s well-definedness is itself a genuine question, one the chapter revisits explicitly a second time later on. Section 3.1.2 poses a deliberately counter-intuitive construction problem to answer the hedging question: a single vanilla option has spot-dependent vega, so does a portfolio of vanilla options exist whose *aggregate* vega is independent of $S$? Writing the vega of a portfolio with density $\rho(K)$ and changing variables to $u=S/K$,

$$
\begin{equation*}
\operatorname{Vega}_{\Pi}(S)=\int \frac{d u}{u^{3}} S^{2} \rho\left(\frac{S}{u}\right) \varphi(u) \tag{3.4}
\end{equation*}
$$
[[Ch03#^eq3-4|Eq. (3.4)]]

shows that $S$-independence requires precisely

$$
\begin{equation*}
\rho(K) \propto \frac{1}{K^{2}} \tag{3.5}
\end{equation*}
$$
[[Ch03#^eq3-5|Eq. (3.5)]]

— a $1/K^2$-weighted vanilla portfolio, which Section 3.1.3's digression on payoff replication shows synthesizes exactly the *log contract*. Replicating a general payoff $f$ via cash, forwards, and a continuum of puts/calls struck at the forward,

$$
\begin{equation*}
P_{f}=f\left(F_{T}\right) e^{-r(T-t)}+\int_{0}^{F_{T}} \frac{d^{2} f}{d K^{2}} P_{K} d K+\int_{F_{T}}^{\infty} \frac{d^{2} f}{d K^{2}} C_{K} d K \tag{3.7}
\end{equation*}
$$
[[Ch03#^eq3-7|Eq. (3.7)]]

recovers, upon integrating the $1/K^2$ density twice, the payoff $-2\ln S$ (the factor of 2 is a convention anticipating variance swaps in Chapter 5), whose Black-Scholes price is

$$
\begin{equation*}
Q^{T}(t, S)=-2 e^{-r(T-t)}\left(\ln S+(r-q)(T-t)-\frac{\widehat{\sigma}_{T}^{2}}{2}(T-t)\right) \tag{3.8}
\end{equation*}
$$
[[Ch03#^eq3-8|Eq. (3.8)]]

Combining $T_1$- and $T_2$-maturity log contracts in the ratio that matches $dG/d\widehat\sigma_{T_1T_2}$ gives the static hedge:

$$
\begin{equation*}
\Pi=-P+\frac{\mathcal{N}}{2}\left(e^{r\left(T_{2}-T_{1}\right)} Q^{T_{2}}-Q^{T_{1}}\right) \tag{3.9}
\end{equation*}
$$
[[Ch03#^eq3-9|Eq. (3.9)]]

whose hedge ratios depend on neither $t$ nor $S$ — genuinely stable, needing rebalancing only when $\widehat\sigma_{T_1T_2}$ itself moves — since the log contract's dollar gamma is a constant,

$$
\begin{equation*}
S^{2} \frac{d^{2} Q^{T}}{d S^{2}}=2 e^{-r(T-t)} \tag{3.10}
\end{equation*}
$$
[[Ch03#^eq3-10|Eq. (3.10)]]

Even this perfect first-order hedge is not risk-free, though. The exact P&L of $\Pi$, expanded before any approximation,

$$
\begin{align*}
& P \& L_{\Pi}=-e^{-r\left(T_{1}-t\right)}\left(G\left(\widehat{\sigma}_{T_{1} T_{2}}+\delta \widehat{\sigma}_{T_{1} T_{2}}\right)-G\left(\widehat{\sigma}_{T_{1} T_{2}}\right)\right) \\
&+\frac{1}{2 \widehat{\sigma}_{T_{1} T_{2}}} \frac{d G}{d \widehat{\sigma}_{T_{1} T_{2}}} e^{-r\left(T_{1}-t\right)}\left(\left(\widehat{\sigma}_{T_{1} T_{2}}+\delta \widehat{\sigma}_{T_{1} T_{2}}\right)^{2}-\widehat{\sigma}_{T_{1} T_{2}}^{2}\right) \\
&=e^{-r\left(T_{1}-t\right)}\left[-\left(G\left(\widehat{\sigma}_{T_{1} T_{2}}+\delta \widehat{\sigma}_{T_{1} T_{2}}\right)-G\left(\widehat{\sigma}_{T_{1} T_{2}}\right)\right)+\frac{d G}{d\left(\widehat{\sigma}_{T_{1} T_{2}}^{2}\right)} \delta\left(\widehat{\sigma}_{T_{1} T_{2}}^{2}\right)\right] \tag{3.11}
\end{align*}
$$
[[Ch03#^eq3-11|Eq. (3.11)]]

vanishes exactly only if $G$ happens to be affine in $\widehat\sigma_{T_1T_2}^2$ (a case of no commercial interest); generically the leading residual term is second order in $\delta(\widehat\sigma_{T_1T_2}^2)$:

$$
\begin{equation*}
P \& L_{\Pi}=-\frac{e^{-r\left(T_{1}-t\right)}}{2} \frac{d^{2} G}{d\left(\widehat{\sigma}_{T_{1} T_{2}}^{2}\right)^{2}}\left(\delta\left(\widehat{\sigma}_{T_{1} T_{2}}^{2}\right)\right)^{2} \tag{3.12}
\end{equation*}
$$
[[Ch03#^eq3-12|Eq. (3.12)]]

This is "volatility-of-volatility risk" before $T_1$: a real, unhedgeable exposure to the *variation* of the forward volatility itself, distinct from anything a vanilla vega hedge addresses. A second, independent risk appears exactly at $T_1$: the cliquet turns into a plain vanilla option of maturity $T_2$, whose market price uses the smile actually observed at $T_1$, while the hedge portfolio was built assuming the log-contract implied volatility $\widehat\sigma_{T_2}(T_1)$ — a gap called forward-smile risk. Section 3.1.6's conclusion prices both honestly rather than pretending they cancel:

$$
\begin{equation*}
P=e^{-r T_{1}} G\left(\widehat{\sigma}_{T_{1} T_{2}}(t=0)\right)+\left(\delta P_{1}+\delta P_{2}\right) \tag{3.17}
\end{equation*}
$$
[[Ch03#^eq3-17|Eq. (3.17)]]

with $\delta P_1$ an explicit charge for volatility-of-volatility risk and $\delta P_2$ for forward-smile risk — the quoted price is the toy-model price *plus* two named, separately-justified risk premiums, not a single clean number. Section 3.1.7 goes further and shows, via a model-independent super/sub-replication argument using only the market's $T_1,T_2$ vanilla smiles, that these smiles alone bound the cliquet's price only very loosely (an ATM forward call's implied vol can range roughly 9%–25% against a flat 20% input smile) — calibrating a stochastic volatility model to the vanilla smile "on the grounds that... vanilla option prices provide information that should, somehow, be used" is explicitly called fallacious, since only instruments with genuinely congruent risk narrow that range.

Section 3.1.9 repeats the whole construction for payoffs linear, rather than independent, in spot — forward volatility agreements common in FX, $S_{T_1}g(S_{T_2}/S_{T_1})$, whose Black-Scholes value is

$$
\begin{equation*}
P(t, S)=S e^{-q\left(T_{1}-t\right)} G\left(\widehat{\sigma}_{T_{1} T_{2}}\right) \tag{3.18}
\end{equation*}
$$
[[Ch03#^eq3-18|Eq. (3.18)]]

The matching hedge instrument this time is the $S\ln S$ contract,

$$
\begin{equation*}
R^{T}(t, S)=S e^{-q(T-t)}\left(\ln S+(r-q)(T-t)+\frac{(T-t) \widehat{\sigma}_{T}^{2}}{2}\right) \tag{3.19}
\end{equation*}
$$
[[Ch03#^eq3-19|Eq. (3.19)]]

whose vega, $dR^T/d\widehat\sigma_T=Se^{-q(T-t)}(T-t)\widehat\sigma_T$, matches the $S$-dependence and repo-rate sensitivity of the FVA's own $\widehat\sigma_{T_1},\widehat\sigma_{T_2}$ vegas exactly, giving an analogous static hedge $\Pi=-P+\mathcal N(e^{q(T_2-T_1)}R^{T_2}-R^{T_1})$:

$$
\begin{equation*}
\Pi=-P+\mathcal{N}\left(e^{q\left(T_{2}-T_{1}\right)} R^{T_{2}}-R^{T_{1}}\right) \tag{3.23}
\end{equation*}
$$
[[Ch03#^eq3-23|Eq. (3.23)]]

and a residual P&L of exactly the same form as (3.12):

$$
\begin{equation*}
P \& L_{\Pi}=-\frac{S e^{-q\left(T_{1}-t\right)}}{2} \frac{d^{2} G}{d\left(\widehat{\sigma}_{T_{1} T_{2}}^{2}\right)^{2}}\left(\delta\left(\widehat{\sigma}_{T_{1} T_{2}}^{2}\right)\right)^{2} \tag{3.24}
\end{equation*}
$$
[[Ch03#^eq3-24|Eq. (3.24)]]

A hasty assessment of the FVA payoff might have singled out the cross spot/forward-volatility covariance as a main risk — but with the correct hedge instrument, that covariance simply does not appear in the residual P&L at all; the reapplication of the well-definedness condition (3.3) to the $S\ln S$ payoff shows it, too, holds here, closing the chapter's recurring "is $\widehat\sigma_{T_1T_2}$ well-defined?" question a second time.

Section 3.2 turns the same lens on the local volatility model specifically, and the news is not good. Expanding around a constant volatility $\sigma_0$ with a locally quadratic local-vol parametrization, dropping curvature, gives the log-contract implied vol as a function of the ATMF skew:

$$
\begin{equation*}
\widehat{\sigma}_{T}=\widehat{\sigma}_{F_{T} T}-\frac{\widehat{\sigma}_{F_{T} T}^{2} T}{2} \mathcal{S}_{T} \tag{3.30}
\end{equation*}
$$
[[Ch03#^eq3-30|Eq. (3.30)]]

and the model's forward-volatility dynamics, at leading order, collapses to a simple average of the *current* local skew $\alpha(u)$ over the residual interval:

$$
\begin{align*}
d \widehat{\sigma}_{T_{1} T_{2}} & =\frac{T_{2}-t}{T_{2}-T_{1}} \frac{\widehat{\sigma}_{T_{2}}}{\widehat{\sigma}_{T_{1} T_{2}}} d \widehat{\sigma}_{T_{2}}-\frac{T_{1}-t}{T_{2}-T_{1}} \frac{\widehat{\sigma}_{T_{1}}}{\widehat{\sigma}_{T_{1} T_{2}}} d \widehat{\sigma}_{T_{1}} \\
& \simeq\left(\frac{1}{T_{2}-T_{1}} \int_{t}^{T_{2}} \alpha(u) d u-\frac{1}{T_{2}-T_{1}} \int_{t}^{T_{1}} \alpha(u) d u\right) d \ln S_{t} \\
& \simeq\left(\frac{1}{T_{2}-T_{1}} \int_{T_{1}}^{T_{2}} \alpha(u) d u\right) d \ln S_{t} \tag{3.33}
\end{align*}
$$
[[Ch03#^eq3-33|Eq. (3.33)]]

Both the model's implied vol-of-vol *and* its forward skew are thus dictated purely by the calibration-day skew of the local volatility surface — the model is not time-homogeneous, and both risks are systematically mispriced as a result. The corresponding forward ATMF skew, for residual maturity $\theta=T_2-T_1$ observed at $T_1$,

$$
\begin{equation*}
\mathcal{S}_{\theta}\left(S_{T_{1}}, T_{1}\right)=\int_{0}^{1} \alpha\left(T_{1}+u \theta\right) u d u \tag{3.34}
\end{equation*}
$$
[[Ch03#^eq3-34|Eq. (3.34)]]

is systematically weaker than the spot-starting skew of the same residual maturity — the same "local volatility misprices forward-smile risk" finding as Chapter 2's Section 2.6, now demonstrated concretely for a cliquet. Section 3.2.4 makes this failure vivid by deriving the model's own *implied* vega hedge for a forward-start call explicitly. Perturbing the pricing PDE for the payoff $(S_{T_2}/S_{T_1}-k)^+$ and tracking the resulting dollar-gamma density $\phi$ shows it has a genuine discontinuity at $T_1$ (a static bucket of vanilla options of maturity $T_1$),

$$
\begin{align*}
\Psi_{1}(K) & =-\frac{1}{K^{2}} \phi\left(T_{1}^{+}, S\right) \\
& =-\frac{1}{K^{2}} \frac{1}{\sqrt{2 \pi \sigma_{0}^{2}\left(T_{2}-T_{1}\right)}} e^{-\frac{\left(-\ln k+\frac{\sigma_{0}^{2}\left(T_{2}-T_{1}\right)}{2}\right)^{2}}{2 \sigma_{0}^{2}\left(T_{2}-T_{1}\right)}} \tag{3.36}
\end{align*}
$$
[[Ch03#^eq3-36|Eq. (3.36)]]

but also a second discontinuity at $T_2$, requiring a *second*, differently-weighted static bucket of $T_2$-maturity vanilla options,

$$
\begin{equation*}
\Psi_{2}(K)=\frac{1}{K^{2}} \frac{k}{\sqrt{2 \pi \sigma_{0}^{2} T_{1}\left(1-\frac{T_{1}}{T_{2}}\right)}} e^{-\frac{\left(\left(1-\frac{T_{1}}{T_{2}}\right) \ln \frac{K}{S_{0}}-\ln k\right)^{2}}{2 \sigma_{0}^{2} T_{1}\left(1-\frac{T_{1}}{T_{2}}\right)}} \tag{3.38}
\end{equation*}
$$
[[Ch03#^eq3-38|Eq. (3.38)]]

plus a continuous band of intermediate-maturity vanilla options between $T_1$ and $T_2$ that must be perpetually rebalanced as time elapses and spot moves. Comparing the relative sizes of these three pieces,

$$
\begin{equation*}
\frac{\Psi_{2}(K)}{\left|\Psi_{1}(K)\right|}, \frac{\left(T_{2}-T_{1}\right) \mu(T, K)}{\left|\Psi_{1}(K)\right|} \tag{3.39}
\end{equation*}
$$
[[Ch03#^eq3-39|Eq. (3.39)]]

confirms local volatility's own hedge for a cliquet is nothing like the clean, static log-contract portfolio of Section 3.1.4 — evidence, the chapter argues, that both the model's price and its implied hedge for this instrument are suspicious.

It is worth stepping back to see why this chapter's two main threads — the model-independent replication analysis of Section 3.1 and the local-volatility-specific critique of Section 3.2 — are, in a sense, telling the same story from two different vantage points. Section 3.1 shows that *no* model, however it is parametrized, can make volatility-of-volatility risk and forward-smile risk disappear through clever hedging alone: they are structural features of any payoff that couples two future dates, not an artifact of a poorly chosen model. Section 3.2 then shows what happens when a *specific*, non-time-homogeneous model is asked to price and hedge that same structural risk: it does not fail outright, but it systematically understates both risks and generates an implausibly complicated hedge in the process, because its own dynamics were never actually chosen — they are a byproduct of whatever smile happened to be quoted on the calibration day. The natural response, deferred to Chapter 12's comparison of local, stochastic, and local-stochastic volatility models, is to ask whether a model that *is* time-homogeneous, and whose forward-smile behavior can be controlled directly rather than inherited passively from calibration, would price cliquets more honestly.

#### Chapter's Digest and Logical Bridge to Chapter 4

Cliquets need a dynamic hedge (log contracts, or $S\ln S$ contracts for FVA-type payoffs) that leaves two irreducible, separately-priced charges — volatility-of-volatility risk and forward-smile risk — even after the static hedge is in place; and the local volatility model, precisely because it is not time-homogeneous, systematically misprices both, producing an implausible, non-static implied hedge of its own. Chapter 4 turns from this diagnostic (what goes wrong, and why) to a genuinely constructive question: can a model be built that lets the modeler *freely specify* the joint dynamics of spot and implied volatility, rather than having that dynamics forced on it by a single day's calibrated smile?

---

### Chapter 4: Stochastic volatility — introduction

Chapter 2's local volatility model, in this book's framing, is an extremely constrained special case of stochastic volatility: every implied volatility moves in 100%-correlated lockstep with $S$ and with every other implied volatility, and whatever volatility-of-volatility the model happens to generate is a byproduct of calibration, not something the modeler ever chose. This chapter's governing question is whether a model can instead let the modeler *freely specify* the joint dynamics of $S$ and its implied volatilities — and it proceeds by testing three successively less ambitious targets, in order, to find out which one is actually tractable.

The first, most direct attempt is to model vanilla option prices themselves. Assuming zero rates and repo, write joint SDEs for $S$ and a call price $C_{KT}$:

$$
\left\{\begin{align*}
d S & =\bar{\sigma} S d W^{S}  \tag{4.1}\\
d C_{K T} & =\Lambda_{K T} d W^{K T}
\end{align*}\right.
$$
[[Ch04#^eq4-1|Eq. (4.1)]]

where $W^S,W^{KT}$ are Brownian motions, $\Lambda_{KT}$ is the (as yet unspecified) volatility of $C_{KT}$, and $\bar\sigma$ is the instantaneous volatility of $S$. Since $r=0$ by assumption, $C_{KT}$'s pricing drift must vanish, and the option's payoff must be recovered exactly at maturity — the terminal condition:

$$
\begin{equation*}
C_{K T}(t=T)=\left(S_{T}-K\right)^{+} \tag{4.2}
\end{equation*}
$$
[[Ch04#^eq4-2|Eq. (4.2)]]

Working instead with implied volatilities $\widehat\sigma_{KT}$ directly — since a directional position on $\widehat\sigma_{KT}$ is just a delta-hedged vanilla position, and the zero-drift condition on $C_{KT}$ translates into a specific risk-neutral drift for $\widehat\sigma_{KT}$ — gives the equivalent joint SDEs

$$
\left\{\begin{align*}
d S & =\bar{\sigma} S d W^{S}  \tag{4.3}\\
d \widehat{\sigma}_{K T} & =\mu_{K T} d t+\lambda_{K T} d W^{K T}
\end{align*}\right.
$$
[[Ch04#^eq4-3|Eq. (4.3)]]

with terminal condition on the implied volatility itself,

$$
\begin{equation*}
\lim _{t \rightarrow T}(T-t) \widehat{\sigma}_{K T, t}^{2}=0 \tag{4.4}
\end{equation*}
$$
[[Ch04#^eq4-4|Eq. (4.4)]]

A necessary consistency constraint falls out immediately: as maturity shrinks to zero, the at-the-money implied volatility must converge to the instantaneous volatility of $S$ itself,

$$
\begin{equation*}
\lim _{T \rightarrow t} \widehat{\sigma}_{S T, t}=\bar{\sigma}_{t} \tag{4.5}
\end{equation*}
$$
[[Ch04#^eq4-5|Eq. (4.5)]]

but the deeper problem is structural, not just a boundary condition to satisfy: any solution of (4.3) over $[0,T]$ automatically determines every implied volatility $\widehat\sigma_{K\tau}$ for $\tau<T$ as well, so the initial smile the modeler wanted to specify freely turns out to be a *consequence* buried inside the process for $\bar\sigma$, not an independent input. Direct modeling of vanilla implied volatilities is declared impractical.

Section 4.2 retreats to a related but distinct target: rather than the smile itself, model the dynamics of the *local volatility function* that the smile implies via Dupire's formula at each instant. In a model where the local volatility function is itself dynamic, its short end is

$$
\begin{equation*}
\bar{\sigma}_{t}=\sigma_{t}\left(t, S_{t}\right) \tag{4.6}
\end{equation*}
$$
[[Ch04#^eq4-6|Eq. (4.6)]]

— note the subscript $t$ rather than a fixed calibration date $t_0$: the local volatility function itself is now allowed to move. Writing joint SDEs for $S$ and the whole surface $\sigma_{\tau\mathcal S}^2$, and imposing that every vanilla price has vanishing drift, requires $\sigma_{\tau\mathcal S}^2$ to still satisfy the Dupire ratio at each instant,

$$
\begin{equation*}
\sigma_{\tau \mathcal{S}}^{2}=\left.2 \frac{\frac{d C_{K T}}{d T}}{K^{2} \frac{d^{2} C_{K T}}{d K^{2}}}\right|_{\substack{K=\mathcal{S} \\ T=\tau}} \tag{4.8}
\end{equation*}
$$
[[Ch04#^eq4-8|Eq. (4.8)]]

Propagating the zero-drift requirement through the transition density $\rho_{\tau\mathcal S}=E_t[\delta(S_\tau-\mathcal S)]$ — which itself solves the backward equation $d\rho_{\tau\mathcal S}/dt+\tfrac12\sigma_{tS}^2S^2\,d^2\rho_{\tau\mathcal S}/dS^2=0$ — and using the resulting functional derivative of $\rho_{\tau\mathcal S}$ with respect to a perturbation of the local variance surface produces an explicit, but forbidding, drift for $\sigma_{\tau\mathcal S}^2$:

$$
\begin{align*}
\mu_{\tau \mathcal{S}}=-\frac{\lambda_{\tau \mathcal{S}}}{\rho_{\tau \mathcal{S}}} & \left(S \frac{d \rho_{\tau \mathcal{S}}\left(t S, \sigma^{2}\right)}{d S} \sigma_{t S} \frac{\left\langle d W^{S} d W^{\tau \mathcal{S}}\right\rangle}{d t}\right.  \tag{4.14}\\
& \left.+\frac{1}{2} \int_{t}^{\tau} d u \int_{0}^{\infty} d x \rho_{u x} x^{2} \frac{d^{2} \rho_{\tau \mathcal{S}}\left(u x, \sigma^{2}\right)}{d x^{2}} \lambda_{u x} \frac{\left\langle d W^{u x} d W^{\tau \mathcal{S}}\right\rangle}{d t}\right)
\end{align*}
$$
[[Ch04#^eq4-14|Eq. (4.14)]]

This drift is non-local — it depends on the *whole* local volatility function through a double integral, not just its value at $(\tau,\mathcal S)$ — and requires forward transition densities for every intermediate date, making it computationally far too expensive to use as a practical model. Worse, inspecting (4.14) for special cases where the drift vanishes turns up a paradox rather than a solution: freezing every shock $\lambda_{\tau\mathcal S}\equiv0$ trivially recovers the static local volatility model, as expected, but so does the seemingly much richer choice of giving *every point* of the local volatility surface its own independent Brownian shock, uncorrelated with $S$ and with every other point on the surface. Because the expectation $E[\bar\sigma_t^2]$ is preserved by construction in that case, the randomness averages out over any finite horizon and vanilla prices come out numerically identical to the static local volatility model's — a seemingly non-trivial stochastic construction collapses to no new content at all. Modeling the dynamic local volatility function directly is, once again, an impasse.

Section 4.3 finally succeeds, by abandoning vanilla options as the modeling target altogether in favor of a deliberately simple payoff family: *power payoffs* $S_T^p$, whose Black-Scholes price (with zero cash dividends) is closed-form in a single implied volatility $\widehat\sigma_{pT}$:

$$
\begin{equation*}
Q^{p T}=e^{-r(T-t)} F_{T}^{p} e^{\frac{p(p-1)}{2}(T-t) \widehat{\sigma}^{2}} \tag{4.15}
\end{equation*}
$$
[[Ch04#^eq4-15|Eq. (4.15)]]

Using Matytsin's moneyness change of variable $z(K)=\ln(F_T/K)/(\widehat\sigma_{KT}\sqrt T)-\widehat\sigma_{KT}\sqrt T/2$ applied to the characteristic-function representation of $Q^{pT}$, and re-expressing everything in terms of a $p$-dependent moneyness $y=z+p\omega(z)$,

$$
\begin{equation*}
y(K)=\frac{\ln \left(\frac{F_{T}}{K}\right)}{\widehat{\sigma}_{K T} \sqrt{T}}+\left(p-\frac{1}{2}\right) \widehat{\sigma}_{K T} \sqrt{T} \tag{4.19}
\end{equation*}
$$
[[Ch04#^eq4-19|Eq. (4.19)]]

yields a genuinely model-free relationship — no dynamics assumed at all — expressing a power payoff's implied volatility as a Gaussian-moneyness-weighted average of the vanilla smile:

$$
\begin{equation*}
e^{\frac{p(p-1)}{2} \widehat{\sigma}_{p}^{2} T}=\int_{-\infty}^{+\infty} d y \frac{e^{-\frac{y^{2}}{2}}}{\sqrt{2 \pi}} e^{\frac{p(p-1)}{2} \widehat{\sigma}_{K(y, p) T}^{2}} \tag{4.20}
\end{equation*}
$$
[[Ch04#^eq4-20|Eq. (4.20)]]

Taking the limits $p\to0$ and $p\to1$ in this averaging formula isolates two special, practically important degenerate payoffs. As $p\to0$, expanding (4.20) to first order gives the implied volatility of the *log contract*:

$$
\begin{align*}
\widehat{\sigma}_{\ln S}^{2} & =\int_{-\infty}^{+\infty} d y \frac{e^{-\frac{y^{2}}{2}}}{\sqrt{2 \pi}} \widehat{\sigma}_{K(y) T}^{2}  \tag{4.21a}\\
y(K) & =\frac{\ln \left(\frac{F_{T}}{K}\right)}{\widehat{\sigma}_{K T} \sqrt{T}}-\frac{\widehat{\sigma}_{K T} \sqrt{T}}{2} \tag{4.21b}
\end{align*}
$$
[[Ch04#^eq4-21a|Eq. (4.21a)]] [[Ch04#^eq4-21b|Eq. (4.21b)]]

whose Black-Scholes dollar gamma and vega are both independent of $S$, and which is the direct theoretical ancestor of the variance swap studied in Chapter 5. As $p\to1$, the analogous limit gives the $S\ln S$ contract instead, whose dollar gamma is proportional to $S$ rather than constant:

$$
\begin{align*}
\widehat{\sigma}_{S \ln S}^{2} & =\int_{-\infty}^{+\infty} d y \frac{e^{-\frac{y^{2}}{2}}}{\sqrt{2 \pi}} \widehat{\sigma}_{K(y) T}^{2}  \tag{4.22a}\\
y(K) & =\frac{\ln \left(\frac{F_{T}}{K}\right)}{\widehat{\sigma}_{K T} \sqrt{T}}+\frac{\widehat{\sigma}_{K T} \sqrt{T}}{2} \tag{4.22b}
\end{align*}
$$
[[Ch04#^eq4-22a|Eq. (4.22a)]] [[Ch04#^eq4-22b|Eq. (4.22b)]]

Because power payoffs are convex (or concave) functions of $S_T/F_T$, their implied volatilities obey the same convex-order (calendar-monotonicity) condition established in Chapter 2:

$$
\begin{equation*}
\left(T_{2}-t\right) \widehat{\sigma}_{p T_{2}}^{2} \geq\left(T_{1}-t\right) \widehat{\sigma}_{p T_{1}}^{2} \tag{4.23}
\end{equation*}
$$
[[Ch04#^eq4-23|Eq. (4.23)]]

which licenses defining a genuinely *positive* forward variance — discretely,

$$
\begin{equation*}
\xi^{p T_{1} T_{2}}=\frac{\left(T_{2}-t\right) \widehat{\sigma}_{p T_{2}}^{2}-\left(T_{1}-t\right) \widehat{\sigma}_{p T_{1}}^{2}}{T_{2}-T_{1}} \tag{4.24}
\end{equation*}
$$
[[Ch04#^eq4-24|Eq. (4.24)]]

or, in the continuous limit used throughout the rest of the book,

$$
\begin{equation*}
\xi^{p T}=\frac{d}{d T}\left((T-t) \hat{\sigma}_{p T}^{2}\right) \tag{4.25}
\end{equation*}
$$
[[Ch04#^eq4-25|Eq. (4.25)]]

This is the chapter's central conceptual pivot: forward variance, not spot (instantaneous) variance, becomes the fundamental state variable whose dynamics the model gets to specify — precisely the freedom that eluded Sections 4.1 and 4.2. Writing joint SDEs for $(S,\xi^{pT})$ and imposing that $Q^{pT}$ have vanishing drift pins down the forward-variance drift explicitly:

$$
\begin{equation*}
\mu^{T}=-\left(p \bar{\sigma} \lambda^{T} \rho^{S T}+\frac{p(p-1)}{2} \int_{t}^{T} \rho^{T u} \lambda^{T} \lambda^{u} d u\right) \tag{4.29}
\end{equation*}
$$
[[Ch04#^eq4-29|Eq. (4.29)]]

and taking $T\to t$ in the price SDE leaves only the drift's first term, forcing the short end of every such curve, for *every* value of $p$, to equal the instantaneous variance of $S$:

$$
\begin{equation*}
\xi_{t}^{t}=\bar{\sigma}_{t}^{2} \tag{4.30}
\end{equation*}
$$
[[Ch04#^eq4-30|Eq. (4.30)]]

Substituting (4.29) back in gives the final, fully explicit joint SDEs — a genuinely tractable stochastic-volatility structure, in sharp contrast to the impasses of Sections 4.1 and 4.2:

$$
\begin{align*}
d S_{t} & =\sqrt{\xi_{t}^{t}} S_{t} d W_{t}^{S}  \tag{4.31a}\\
d \xi_{t}^{p T} & =-\left(p \sqrt{\xi_{t}^{t}} \lambda_{t}^{T} \rho^{S T}+\frac{p(p-1)}{2} \int_{t}^{T} \rho^{T u} \lambda_{t}^{T} \lambda_{t}^{u} d u\right) d t+\lambda_{t}^{T} d W_{t}^{T} \tag{4.31b}
\end{align*}
$$
[[Ch04#^eq4-31a|Eq. (4.31a)]] [[Ch04#^eq4-31b|Eq. (4.31b)]]

Section 4.3.4 then asks whether this whole infinite-dimensional curve $\xi_t^T$ can be collapsed to a low-dimensional Markov representation, borrowing an HJM-style trick from Cheyette. Assuming a specific separable covariance structure,

$$
\begin{align*}
d S_{t} & =\sqrt{\xi_{t}^{t}} S_{t} d W_{t}^{S}  \tag{4.32a}\\
d \xi_{t}^{T} & =\mu_{t}^{T} d t+\sum_{i=0}^{n} \alpha_{i}(T) \beta_{i t} d W_{t}^{i} \tag{4.32b}
\end{align*}
$$
[[Ch04#^eq4-32a|Eq. (4.32a)]] [[Ch04#^eq4-32b|Eq. (4.32b)]]

gives an explicit representation of the whole curve in terms of a handful of state variables $x_{it}$, whose short end (setting $T=t$) is

$$
\begin{equation*}
\xi_{t}^{t}=\xi_{0}^{t}+\Sigma_{i} \alpha_{i}(t) x_{i t} \tag{4.34}
\end{equation*}
$$
[[Ch04#^eq4-34|Eq. (4.34)]]

This is a fruitful failure rather than a clean success: the Cheyette-style construction is shown to be structurally incompatible with the non-negativity constraint on instantaneous variance — unlike a forward *rate*, which can go negative without breaking a model, a forward *variance* cannot, and the ansatz (4.32b) provides no mechanism to enforce this. The problem is left explicitly open here, resolved only in Chapter 7's more careful treatment of finite-factor forward-variance models.

A related practical question, taken up in Section 4.3.5, is whether the power-payoff family for different values of $p$ could each be calibrated to their own, independently-observed market data, giving the modeler several separate variance curves $\xi^{pT}$ to work with simultaneously. The boundary condition (4.30) rules this out: since the short end of *every* such curve must coincide with the same instantaneous variance of $S$ regardless of $p$, the curves cannot be freely and independently specified — they are all tied together at $T=t$ by the same underlying dynamics of $S$ itself, so in practice only one power-payoff family (almost always $p\to0$, the log contract) is used as the model's actual state variable.

The chapter closes by revisiting the $p\to0$ case once more, since the log contract is the cleanest possible instance of the whole state-variable idea. Taking the limit $p\to0$ in (4.31) directly gives:

$$
\left\{\begin{align*}
d S_{t} & =\sqrt{\zeta_{t}^{t}} S_{t} d W_{t}^{S}  \tag{4.36}\\
d \zeta_{t}^{T} & =\lambda_{t}^{T} d W_{t}^{T}
\end{align*}\right.
$$
[[Ch04#^eq4-36|Eq. (4.36)]]

so the log contract's forward variance curve $\zeta_t^T$ has *no drift at all* — the cleanest possible realization of "forward variance is the fundamental state variable." Log contracts themselves are not directly traded (they would require an unbounded range of strikes to replicate exactly), but variance swaps are, and Chapter 5 shows that (4.36) carries over essentially unchanged once $\zeta_t^T$ is replaced by variance-swap forward variances.

#### Chapter's Digest and Logical Bridge to Chapter 5

Three successively narrower modeling targets were tried: vanilla prices/implied volatilities directly (impractical — the smile becomes a derived consequence, not a free input), the dynamic local volatility function (an impasse — the drift is non-local and a seemingly rich construction collapses to the static model), and power payoffs (successful — a model-free averaging relationship to the vanilla smile, and forward variance as a state variable with an explicit, tractable joint SDE). Chapter 5 picks up exactly at the $p\to0$ log-contract limit and shows it is, for all practical purposes, identical to the market's actual variance-swap contract.

---

### Chapter 5: Variance swaps

This chapter is explicitly billed as a prerequisite for everything that follows: variance swaps (VS) and their connection to delta-hedged log contracts are the concrete bridge between Chapter 4's abstract "forward variance as fundamental state variable" pivot and the general stochastic volatility models of Chapter 7. The chapter's overall shape is a controlled demolition — it first shows VS and log contracts are, in an idealized diffusive world, literally the same instrument, then spends the rest of its length carefully cataloguing every real-world effect (jumps, dividends, discreteness, stochastic rates) that breaks this equivalence, insisting each time that the break be priced as an explicit, separately-justified charge rather than smuggled in through a convenient model choice.

A variance swap (VS) pays at maturity the realized variance of the underlying, computed as the sum of squared daily log-returns, less a fixed strike $\widehat\sigma_{\mathrm{VS},T}^2$ set so the contract has zero initial value. The market convention,

$$
\begin{equation*}
\frac{252}{N} \sum_{i=0}^{N-1} \ln ^{2}\left(\frac{S_{i+1}}{S_{i}}\right)-\widehat{\sigma}_{\mathrm{VS}, T}^{2}(t) \tag{5.1}
\end{equation*}
$$
[[Ch05#^eq5-1|Eq. (5.1)]]

is biased for short maturities because trading days aren't uniformly distributed through the year, so the chapter works instead with the cleaner convention

$$
\begin{equation*}
\frac{1}{T-t} \sum_{i=0}^{N-1} \ln ^{2}\left(\frac{S_{i+1}}{S_{i}}\right)-\widehat{\sigma}_{\mathrm{VS}, T}^{2}(t) \tag{5.2}
\end{equation*}
$$
[[Ch05#^eq5-2|Eq. (5.2)]]

Before anything else, the chapter establishes that forward VS variance is genuinely costless to trade: a calendar-spread position — long $(T_2-t)$ VSs of maturity $T_2$, short $(T_1-t)e^{-r(T_2-T_1)}$ VSs of maturity $T_1$ — has a payoff, capitalized at $T_2$, of

$$
\begin{align*}
& \sum_{T_{1}}^{T_{2}} \ln ^{2}\left(\frac{S_{i+1}}{S_{i}}\right)-\left(\left(T_{2}-t\right) \widehat{\sigma}_{\mathrm{VS}, T_{2}}^{2}(t)-\left(T_{1}-t\right) \widehat{\sigma}_{\mathrm{VS}, T_{1}}^{2}(t)\right) \\
& =\sum_{T_{1}}^{T_{2}} \ln ^{2}\left(\frac{S_{i+1}}{S_{i}}\right)-\left(T_{2}-T_{1}\right) \widehat{\sigma}_{\mathrm{VS}, T_{1} T_{2}}^{2}(t) \tag{5.4}
\end{align*}
$$
[[Ch05#^eq5-4|Eq. (5.4)]]

and unwinding this position later at $t'<T_1$ cancels the realized-variance term entirely, leaving a pure P&L in the change of forward variance,

$$
\begin{equation*}
\left(T_{2}-T_{1}\right)\left(\widehat{\sigma}_{\mathrm{VS}, T_{1} T_{2}}^{2}\left(t^{\prime}\right)-\widehat{\sigma}_{\mathrm{VS}, T_{1} T_{2}}^{2}(t)\right) \tag{5.5}
\end{equation*}
$$
[[Ch05#^eq5-5|Eq. (5.5)]]

generated at no upfront cost — establishing that the continuous forward variance $\xi_t^T$ is driftless,

$$
\begin{equation*}
d \xi_{t}^{T}=\bullet d W_{t}^{T} \tag{5.6}
\end{equation*}
$$
[[Ch05#^eq5-6|Eq. (5.6)]]

Section 5.2 then asks whether a European payoff exists whose Black-Scholes dollar gamma, risk-managed at zero implied volatility, exactly reproduces the VS payoff's gamma:

$$
\begin{equation*}
\frac{1}{2} S^{2} \frac{d^{2} P_{\widehat{\sigma}=0}}{d S^{2}}=e^{-r(T-t)} \tag{5.9}
\end{equation*}
$$
[[Ch05#^eq5-9|Eq. (5.9)]]

The answer is yes, up to a factor: the log contract, whose Black-Scholes price is

$$
\begin{equation*}
Q^{T}(t, S)=-2 e^{-r(T-t)}\left(\ln S+(r-q)(T-t)-\frac{\widehat{\sigma}^{2}}{2}(T-t)\right) \tag{5.10}
\end{equation*}
$$
[[Ch05#^eq5-10|Eq. (5.10)]]

whose delta and dollar gamma at zero implied vol, $dQ_{\widehat\sigma=0}^T/dS=-e^{-r(T-t)}\cdot2/S$ and $\tfrac12S^2d^2Q_{\widehat\sigma=0}^T/dS^2=e^{-r(T-t)}$, satisfy (5.9) exactly. This gives the chapter's first landmark identity: in a diffusive world, a VS and the corresponding log contract are literally the same object viewed two ways —

$$
\begin{align*}
\widehat{\sigma}_{\mathrm{VS}, T} & =\widehat{\sigma}_{T}  \tag{5.13}\\
\xi_{t}^{T} & =\zeta_{t}^{T} \tag{5.14}
\end{align*}
$$
[[Ch05#^eq5-13|Eq. (5.13)]] [[Ch05#^eq5-14|Eq. (5.14)]]

replicable, since the log contract's own replication uses a $2/K^2$-weighted vanilla portfolio, purely from market vanilla prices:

$$
\begin{equation*}
\widehat{\sigma}_{\mathrm{VS}, T}^{2}=\frac{e^{r T}}{T} \int_{0}^{\infty} \frac{2}{K^{2}}\left(P_{\text {market }}^{K T}-P_{\widehat{\sigma}=0}^{K T}\right) d K \tag{5.16}
\end{equation*}
$$
[[Ch05#^eq5-16|Eq. (5.16)]]

Section 5.2.1 supplies a second, independent route to the same volatility, quoting Chapter 4's power-payoff result directly: a Gaussian-moneyness-weighted average of the vanilla smile,

$$
\begin{align*}
\widehat{\sigma}_{\mathrm{VS}, T}^{2} & =\int_{-\infty}^{+\infty} d y \frac{e^{-\frac{y^{2}}{2}}}{\sqrt{2 \pi}} \widehat{\sigma}_{K(y) T}^{2}  \tag{5.17a}\\
y(K) & =\frac{\ln \left(\frac{F_{T}}{K}\right)}{\widehat{\sigma}_{K T} \sqrt{T}}-\frac{\widehat{\sigma}_{K T} \sqrt{T}}{2} \tag{5.17b}
\end{align*}
$$
[[Ch05#^eq5-17a|Eq. (5.17a)]] [[Ch05#^eq5-17b|Eq. (5.17b)]]

Section 5.3 is the chapter's central reversal: how robust is this VS/log-contract equivalence once returns are not infinitesimally small? This matters because the whole derivation of the equivalence relied on an expansion to second order in the daily return $\delta S/S$ — exactly the order at which delta-hedging works perfectly in the continuous-time limit, and exactly the order that Chapter 1 already flagged as insufficient once hedging happens at a discrete, daily frequency rather than continuously. Section 5.3.1 confirms the equivalence nonetheless holds exactly for *any* diffusive model calibrated to the smile, in the limit of frequent observations — a genuinely clean, model-independent result, since diffusive returns scale like $\sqrt{\Delta t}$ and their third moment is subleading. Section 5.3.2 then breaks it by introducing jumps, which do not obey that same scaling. Under a jump-diffusion SDE, the log contract's implied volatility becomes

$$
\begin{equation*}
\widehat{\sigma}_{T}^{2}=E\left[\frac{1}{T} \int_{0}^{T} \bar{\sigma}_{t}^{2} d t\right]-2 \lambda \overline{\ln (1+J)-J} \tag{5.25}
\end{equation*}
$$
[[Ch05#^eq5-25|Eq. (5.25)]]

while the VS's implied volatility, computed directly from the quadratic variation of $\ln S_t$, gets a different jump contribution,

$$
\begin{equation*}
\widehat{\sigma}_{\mathrm{VS}, T}^{2}=E\left[\frac{1}{T} \int_{0}^{T} \bar{\sigma}_{t}^{2} d t\right]+\lambda \overline{\ln ^{2}(1+J)} \tag{5.27}
\end{equation*}
$$
[[Ch05#^eq5-27|Eq. (5.27)]]

so the two are no longer equal — the clean equivalence of Main Node 4 genuinely breaks, with an explicit, non-zero gap:

$$
\begin{equation*}
\widehat{\sigma}_{\mathrm{VS}, T}^{2}-\widehat{\sigma}_{T}^{2}=\lambda \overline{\ln ^{2}(1+J)+2 \ln (1+J)-2 J} \tag{5.28}
\end{equation*}
$$
[[Ch05#^eq5-28|Eq. (5.28)]]

The natural next question is why this gap only appears once jumps are added: (5.26)'s derivation shows that, for small $dt$, at most one jump can occur with probability proportional to $\lambda\,dt$, so its contribution to instantaneous quadratic variation is genuinely linear in the jump intensity — but this also means the same identity requires *both* $\bar\sigma\sqrt{\Delta t}\ll1$ and $\lambda\Delta t\ll1$ to hold; if observations are not frequent enough, more than one jump can land within a single interval and the whole linear-in-$\lambda$ structure of (5.27)-(5.28) breaks down.

Section 5.3.3-5.3.4 refuses to let a jump model quietly "explain away" this gap, though. Expanding (5.28) in powers of jump size $J$, the leading contribution is cubic,

$$
\begin{equation*}
\widehat{\sigma}_{\mathrm{VS}, T}^{2}-\widehat{\sigma}_{T}^{2} \simeq-\frac{1}{3} \lambda \overline{J^{3}} \tag{5.29}
\end{equation*}
$$
[[Ch05#^eq5-29|Eq. (5.29)]]

which, combined with the small-jump ATMF skew formula quoted from Chapter 10,

$$
\begin{equation*}
\mathcal{S}_{T} \simeq \frac{\lambda \overline{J^{3}}}{6 \widehat{\sigma}_{T}^{3} T} \tag{5.30}
\end{equation*}
$$
[[Ch05#^eq5-30|Eq. (5.30)]]

implies the gap is approximately proportional to the ATMF skew. But Section 5.3.4 then derives the *same* gap in a completely model-free way, directly from the P&L of a static VS-vs-log-contract position, landing on

$$
\begin{equation*}
\widehat{\sigma}_{\mathrm{VS}, T}^{2}-\widehat{\sigma}_{T}^{2} \simeq-\frac{\left\langle r^{3}\right\rangle}{3 \Delta t} \simeq-\frac{s_{\Delta t}}{3} \widehat{\sigma}_{T}^{3} \sqrt{\Delta t} \tag{5.37}
\end{equation*}
$$
[[Ch05#^eq5-37|Eq. (5.37)]]

or, as a relative adjustment,

$$
\begin{equation*}
\frac{\widehat{\sigma}_{\mathrm{VS}, T}}{\widehat{\sigma}_{T}}-1 \simeq-\frac{s_{\Delta t}}{6} \widehat{\sigma}_{T} \sqrt{\Delta t} \tag{5.38}
\end{equation*}
$$
[[Ch05#^eq5-38|Eq. (5.38)]]

with $s_{\Delta t}$ the actual skewness of daily returns — no jump model required at all. That this model-free route and the earlier jump-model route land on structurally the same order-three-in-$J$ (equivalently order-one-in-skewness) result is itself a useful consistency check, but it does not mean the two routes are equally trustworthy as a basis for pricing: only the model-free route ties the adjustment directly to something actually measurable in historical data.

Section 5.3.5 (posed, tellingly, as a question) confronts the two routes to the same gap directly: inverting the jump-model route to infer the implied daily-return skewness a given ATMF skew and maturity would demand,

$$
\begin{equation*}
\widehat{s}_{\Delta t} \simeq 6 \frac{\mathcal{S}_{T} T}{\sqrt{\Delta t}} \tag{5.39}
\end{equation*}
$$
[[Ch05#^eq5-39|Eq. (5.39)]]

gives an implausibly large number (≈−19 for typical index parameters) — and the answer depends on which maturity $T$ happened to be used for the calibration, unless the ATMF skew scales exactly like $1/T$, which — via the identity $\mathcal S_T=s_T/(6\sqrt T)$ derived by perturbing the lognormal density (Appendix B) —

$$
\begin{equation*}
\mathcal{S}_{T}=\frac{s_{T}}{6 \sqrt{T}} \tag{5.40}
\end{equation*}
$$
[[Ch05#^eq5-40|Eq. (5.40)]]

real markets' skew (closer to $1/\sqrt T$) simply does not do. The chapter concludes this model-mediated inference is "unreasonably model-dependent," and backs this with data in Section 5.3.7: the historically *realized* mismatch, measured directly on 20 years of S&P 500 daily returns,

$$
\begin{equation*}
\frac{1}{2}\left(\frac{\left\langle r^{2}\right\rangle}{\left\langle 2\left(e^{r}-1\right)-2 r\right\rangle}-1\right) \tag{5.41}
\end{equation*}
$$
[[Ch05#^eq5-41|Eq. (5.41)]]

is tiny — about 0.2% for typical index parameters, even though the *implied* skewness a naive jump-model calibration would back out is wildly implausible. The chapter's practical resolution is to price any large-return adjustment as an explicit, separately-calibrated stress-test reserve, never as a byproduct of fitting a jump model to the vanilla smile:

$$
\begin{equation*}
\widehat{\sigma}_{\mathrm{VS}, T}^{2}=\widehat{\sigma}_{T}^{2}+\varepsilon\left(\ln ^{2}(1+J)+2 \ln (1+J)-2 J\right) \tag{5.42}
\end{equation*}
$$
[[Ch05#^eq5-42|Eq. (5.42)]]

with $\varepsilon$ an explicit, desk-chosen stress-test frequency rather than anything backed out of the vanilla smile.

Before turning to dividends, Section 5.4 adds one more practical caveat to the model-free story of Section 5.3.4: in reality the log contract is never traded directly, only approximated by a finite strip of vanilla options at discrete strikes. Comparing the VS payoff against the P&L of delta-hedging such a discretized replicating portfolio shows that coarser strike spacing (5% steps rather than 1%) introduces its own mismatch, of a similar order of magnitude to the skewness-driven effect of Section 5.3 — but, unlike that effect, this one is pure discretization noise with no directional bias, and shrinks steadily as the strike grid is refined.

The remainder of the chapter turns from principle to practice. Cash-amount dividends have a negligible effect on the VS payoff itself (Section 5.6.1), since the extra term they contribute is second order and vanishes for large baskets of dividend dates:

$$
\begin{equation*}
E\left[\ln ^{2}\left(\frac{S_{i+1}}{S_{i}}\right)\right]=\bar{\sigma}^{2} \Delta t+\ln ^{2}\left(\frac{S_{i}-d}{S_{i}}\right)+2\left(r-q-\frac{\bar{\sigma}^{2}}{2}\right) \Delta t \ln \left(\frac{S_{i}-d}{S_{i}}\right) \tag{5.44}
\end{equation*}
$$
[[Ch05#^eq5-44|Eq. (5.44)]]

but they *do* break the log contract's constant dollar gamma,

$$
\begin{equation*}
\frac{S^{2}}{2} \frac{d^{2} Q^{T}}{d S^{2}}=\frac{S^{2}}{\left(S-\sum_{t<T_{j}<T} d_{j}\right)^{2}} \tag{5.45}
\end{equation*}
$$
[[Ch05#^eq5-45|Eq. (5.45)]]

which Section 5.6.2 repairs by supplementing the log contract with a strip of extra European payoffs at each dividend date, working backward from the last dividend date $T_N$,

$$
\begin{equation*}
E^{N}(S)=2 \ln \left(\frac{S-d_{N}}{S}\right) \tag{5.46}
\end{equation*}
$$
[[Ch05#^eq5-46|Eq. (5.46)]]

giving the fully dividend-adjusted VS replication formula:

$$
\begin{equation*}
\widehat{\sigma}_{\mathrm{VS}, T}^{2}=\frac{e^{r T}}{T}\left(Q_{\text {market }}^{T}-Q_{\widehat{\sigma}=0}^{T}+\sum_{T_{j}<T}\left(E_{\text {market }}^{j}-E_{\widehat{\sigma}=0}^{j}\right)\right) \tag{5.47}
\end{equation*}
$$
[[Ch05#^eq5-47|Eq. (5.47)]]

For index baskets with many dividend dates, replicating each $E^j$ payoff by hand is impractical; Section 5.7's backward PDE offers a computationally cheaper route to the same $\widehat\sigma_{\mathrm{VS},T}$, since any diffusive model calibrated to the smile gives the identical answer. Finally, Section 5.8 addresses a subtlety often overlooked entirely: with stochastic interest rates, the Black-Scholes implied volatility one backs out of a market price is really the volatility of the *forward*, not of the spot itself, giving a genuine gap between the log-contract volatility $\widehat\sigma_T$ and the VS volatility $\widehat\sigma_{\mathrm{VS},T}$ of the spot,

$$
\begin{equation*}
\widehat{\sigma}_{T}^{2}=\widehat{\sigma}_{\mathrm{VS}, T}^{2}+\rho \widehat{\sigma}_{\mathrm{VS}, T} \sigma_{r} T+\frac{\sigma_{r}^{2} T^{2}}{3} \tag{5.54}
\end{equation*}
$$
[[Ch05#^eq5-54|Eq. (5.54)]]

or, at order one in the interest-rate volatility $\sigma_r$:

$$
\begin{equation*}
\widehat{\sigma}_{\mathrm{VS}, T}=\widehat{\sigma}_{T}-\frac{\rho}{2} \sigma_{r} T \tag{5.55}
\end{equation*}
$$
[[Ch05#^eq5-55|Eq. (5.55)]]

For a 5-year VS with 5bps/day interest-rate volatility and 50% correlation, this is about one point of volatility — "not a small effect," and one easily missed if the VS/log-contract equivalence of Main Node 4 is trusted uncritically at long maturities. This is a genuinely separate source of the VS/log-contract gap from anything discussed so far: it has nothing to do with jumps, dividends, or discreteness, and survives even in a perfectly continuous, dividend-free, pure-diffusion setting, simply because the VS payoff is defined on the spot while the vanilla smile is quoted on the forward.

Section 5.9's weighted variance swaps extend the standard VS by weighting realized variance with an arbitrary function of spot — corridor variance swaps that only accrue variance while spot stays in a band being one popular example — and can be priced with the same backward PDE machinery of Section 5.7, simply replacing the local variance $\sigma^2(t,S)$ with its weighted counterpart $w(S)\sigma^2(t,S)$. The two appendices (timer options; the lognormal-density perturbation underlying (5.40)) extend the machinery further but are supplementary to this chapter's core VS/log-contract arc.

It is worth noting how much of this chapter's methodology recurs almost verbatim later in the book. The insistence on pricing jump/large-return risk as an explicit, separately-controlled reserve rather than folding it into a recalibrated model reappears in Chapter 10's treatment of fat-tailed returns and jump-diffusion models generally; the model-free, replication-based derivation of (5.37)-(5.38) is the template Chapter 9 follows when it links the static ATMF skew to the dynamic Skew Stickiness Ratio; and the dividend-adjustment machinery of Section 5.6 is reused directly whenever a later chapter's variance-swap-like instrument must be priced on a dividend-paying single stock rather than an index.

#### Chapter's Digest and Logical Bridge to Chapter 6

Variance swaps and log contracts are, in a diffusive world, the same object — but jumps, dividends, and stochastic rates each break that equivalence in a different, quantifiable way, and the chapter's consistent methodological stance is to price each break as an explicit, separately-justified adjustment rather than let a convenient model quietly launder it into the vanilla-smile calibration. Chapter 6 turns to the book's first genuine stochastic volatility model, Heston, using its tractability to build the suitability criteria — term-structure flexibility, skew flexibility, vol-of-vol flexibility — that a general class of models will need to satisfy.

---

### Chapter 6: An example of one-factor dynamics: the Heston model

Before building the general forward-variance framework of Chapter 7, this chapter takes a deliberate detour through the Heston model — an analytically tractable, "first-generation" stochastic volatility model written in terms of the instantaneous variance $V_t$ rather than a forward-variance curve. The point is diagnostic: use Heston's tractability to build the suitability criteria (term-structure flexibility, skew flexibility, vol-of-vol flexibility) that Chapter 7's more general models will need to satisfy.

Heston owes its enduring popularity to tractability rather than realism: being an affine model, the Laplace transform of the moment-generating function of $\ln S$ is known analytically, so vanilla option prices follow from a single numerical inversion rather than a full Monte Carlo simulation or PDE solve. The chapter deliberately sets this analytic machinery aside, however, and instead concentrates on the *joint spot/volatility dynamics* the model generates — since it is that dynamics, not the pricing formula, that will determine whether Heston is a suitable building block for the general theory to come. Assuming zero interest rate and repo for simplicity, the defining SDEs are:

$$
\left\{\begin{align*}
d S_{t} & =\sqrt{V_{t}} S_{t} d W_{t}  \tag{6.1}\\
d V_{t} & =-k\left(V_{t}-V^{0}\right) d t+\sigma \sqrt{V_{t}} d Z_{t}
\end{align*}\right.
$$
[[Ch06#^eq6-1|Eq. (6.1)]]

where $\sigma$ ("volatility of volatility") has dimension time$^{-1}$, and $W_t,Z_t$ are correlated Brownian motions with correlation $\rho$. Taking the expectation of the variance SDE gives an ODE for the forward variance $\bar V_T=E_t[V_T]$:

$$
\begin{equation*}
d \bar{V}_{T}=-k\left(\bar{V}_{T}-V^{0}\right) d T \tag{6.2}
\end{equation*}
$$
[[Ch06#^eq6-2|Eq. (6.2)]]

whose solution gives the forward variance in closed form:

$$
\begin{equation*}
\xi_{t}^{T}=V^{0}+e^{-k(T-t)}\left(\xi_{t}^{t}-V^{0}\right) \tag{6.3}
\end{equation*}
$$
[[Ch06#^eq6-3|Eq. (6.3)]]

and the corresponding VS volatility follows by integrating over the residual maturity:

$$
\begin{equation*}
\widehat{\sigma}_{T}^{2}(t)=\frac{1}{T-t} \int_{t}^{T} \xi_{t}^{\tau} d \tau=V^{0}+\frac{1-e^{-k(T-t)}}{k(T-t)}\left(V_{t}-V^{0}\right) \tag{6.4}
\end{equation*}
$$
[[Ch06#^eq6-4|Eq. (6.4)]]

Rewriting the Heston SDEs in terms of this forward-variance curve gives a *Markov-functional* representation, since $\xi_t^T$ is a deterministic function of $\xi_t^t$ via (6.3):

$$
\left\{\begin{align*}
d S_{t} & =\sqrt{\xi_{t}^{t}} S_{t} d W_{t}  \tag{6.5}\\
d \xi_{t}^{T} & =\sigma e^{-k(T-t)} \sqrt{\xi_{t}^{t}} d Z_{t}
\end{align*}\right.
$$
[[Ch06#^eq6-5|Eq. (6.5)]]

This surfaces a subtle but easily-missed point: Heston should *not* be regarded as simply a special case of the forward variance models Chapter 7 will build. A one-dimensional Markov representation like (6.5) exists only because the initial forward-variance curve happens to satisfy the specific ODE constraint (6.2) — a structural restriction, not a free modeling choice — so Heston cannot be recalibrated to an arbitrary term structure of VS volatilities the way a genuine forward-variance model can.

Section 6.3 makes an unusually blunt methodological point about first-generation models generally: the traditional "market price of risk" story invoked to justify $V_t$'s risk-neutral drift is dismissed outright, since $V_t$ is not a physical quantity with its own risk premium — it is simply the short end of the forward-variance curve, and its drift is nothing more than that curve's initial slope, $d\xi_t^T/dT|_{T=t}$. This is a direct consequence of treating forward variance, not instantaneous variance, as fundamental (the pivot made in Chapter 4).

Differentiating (6.4) gives the SDE for the VS volatility itself:

$$
\begin{equation*}
d \widehat{\sigma}_{T}=\bullet d t+\frac{\sigma}{2} \frac{1-e^{-k(T-t)}}{k(T-t)} \frac{\widehat{\sigma}_{t}}{\widehat{\sigma}_{T}} d Z_{t} \tag{6.6}
\end{equation*}
$$
[[Ch06#^eq6-6|Eq. (6.6)]]

with two limiting regimes for short and long residual maturities:

$$
\begin{array}{ll}
T-t \ll \frac{1}{k} & d \widehat{\sigma}_{T} \simeq \bullet d t+\frac{\sigma}{2}\left(1-\frac{k(T-t)}{2}\right) \frac{\widehat{\sigma}_{t}}{\widehat{\sigma}_{T}} d Z_{t} \\
T-t \gg \frac{1}{k} & d \widehat{\sigma}_{T} \simeq \bullet d t+\frac{\sigma}{2} \frac{1}{k(T-t)} \frac{\widehat{\sigma}_{t}}{\widehat{\sigma}_{T}} d Z_{t} \tag{6.8}
\end{array}
$$
[[Ch06#^eq6-8|Eq. (6.8)]]

so for a flat term structure of VS volatilities, the instantaneous (lognormal) volatility of $\widehat\sigma_T$ itself decays like $\frac{1-e^{-k(T-t)}}{k(T-t)}$:

$$
\begin{equation*}
\operatorname{vol}\left(\widehat{\sigma}_{T}\right) \propto \frac{1-e^{-k(T-t)}}{k(T-t)} \tag{6.9}
\end{equation*}
$$
[[Ch06#^eq6-9|Eq. (6.9)]]

Compared against real Euro Stoxx 50 vol-of-vol data, this decays noticeably faster than what is observed empirically (which fits a power law instead) — the model's single embedded time scale $1/k$ is already showing signs of strain.

The chapter's main technical payload derives the ATMF skew perturbatively. Starting from the pricing PDE for a European option $P(t,S,V)$,

$$
\begin{equation*}
\frac{d P}{d t}-k\left(V-V^{0}\right) \frac{d P}{d V}+\frac{V}{2} S^{2} \frac{d^{2} P}{d S^{2}}+\frac{\sigma^{2}}{2} V \frac{d^{2} P}{d V^{2}}+\rho \sigma V S \frac{d^{2} P}{d S d V}=0 \tag{6.10}
\end{equation*}
$$
[[Ch06#^eq6-10|Eq. (6.10)]]

the zero-vol-of-vol ($\sigma=0$) solution $P^0$ solves the simpler PDE

$$
\begin{equation*}
\frac{d P^{0}}{d t}-k\left(V-V^{0}\right) \frac{d P^{0}}{d V}+\frac{V}{2} S^{2} \frac{d^{2} P^{0}}{d S^{2}}=0 \tag{6.11}
\end{equation*}
$$
[[Ch06#^eq6-11|Eq. (6.11)]]

and, since $V$ is then deterministic, reduces exactly to a Black-Scholes price with the deterministic time-dependent volatility $V_\tau(V_t)=V^0+(V_t-V^0)e^{-k(\tau-t)}$:

$$
\begin{equation*}
V_{\tau}\left(V_{t}\right)=V^{0}+\left(V_{t}-V^{0}\right) e^{-k(\tau-t)} \tag{6.12}
\end{equation*}
$$
[[Ch06#^eq6-12|Eq. (6.12)]]

giving the order-zero implied volatility

$$
\begin{equation*}
\widehat{\sigma}_{T}^{2}(t, V)=\frac{1}{T-t} \int_{t}^{T} V_{\tau}(V) d \tau=V^{0}+\left(V-V^{0}\right) \frac{1-e^{-k(T-t)}}{k(T-t)} \tag{6.13}
\end{equation*}
$$
[[Ch06#^eq6-13|Eq. (6.13)]]

Expanding $P=P^0+\delta P$ to order one in $\sigma$ and inserting into (6.10) leaves $\delta P$ solving a forced version of (6.11), sourced entirely by the correlation term:

$$
\begin{equation*}
\frac{d \delta P}{d t}-k\left(V-V^{0}\right) \frac{d \delta P}{d V}+\frac{V}{2} S^{2} \frac{d^{2} \delta P}{d S^{2}}=-\rho \sigma V S \frac{d^{2} P^{0}}{d S d V} \tag{6.14}
\end{equation*}
$$
[[Ch06#^eq6-14|Eq. (6.14)]]

Solving this via the Feynman-Kac representation and using the fact that the discounted dollar gamma is a martingale under the $\sigma=0$ dynamics yields $\delta P$ purely in terms of the second- and third-order log-$S$ derivatives of $P^0$:

$$
\begin{equation*}
\delta P=\frac{\rho \sigma}{2}\left[\int_{0}^{T} V_{\tau}(V) \frac{1-e^{-k(T-\tau)}}{k} d \tau\right]\left(\frac{d^{3} P^{0}}{d \ln S^{3}}-\frac{d^{2} P^{0}}{d \ln S^{2}}\right)_{t=0, S, V} \tag{6.15}
\end{equation*}
$$
[[Ch06#^eq6-15|Eq. (6.15)]]

which is not a coincidence: since the ODE for $E[V_t]$ in (6.2) does not involve $\sigma$ at all, perturbing in powers of $\sigma$ leaves forward variances unchanged at every order — perturbing at order one in volatility of volatility is equivalent to perturbing at order one in the third cumulant with forward variances held fixed. Converting $\delta P$ into an implied-volatility correction via the vega-gamma relation and reinstating rates/repo gives:

$$
\begin{equation*}
\delta \widehat{\sigma}_{K T}=\frac{1}{\widehat{\sigma}_{T}^{3} T^{2}}\left[\frac{\rho \sigma}{2} \int_{0}^{T} V_{\tau}(V) \frac{1-e^{-k(T-\tau)}}{k} d \tau\right]\left(\frac{\widehat{\sigma}_{T}^{2} T}{2}+\ln \frac{K}{F_{T}}\right) \tag{6.16}
\end{equation*}
$$
[[Ch06#^eq6-16|Eq. (6.16)]]

Specializing to the ATMF strike and skew gives the chapter's central quantitative pair:

$$
\begin{align*}
\widehat{\sigma}_{F_{T} T} & =\widehat{\sigma}_{T}\left(1+\frac{\widehat{\sigma}_{T} T}{2} \mathcal{S}_{T}\right)  \tag{6.17a}\\
\mathcal{S}_{T} & =\left.\frac{d \widehat{\sigma}_{K T}}{d \ln K}\right|_{F_{T}}=\frac{1}{\widehat{\sigma}_{T}^{3} T^{2}} \frac{\rho \sigma}{2} \int_{0}^{T} V_{\tau} \frac{1-e^{-k(T-\tau)}}{k} d \tau \tag{6.17b}
\end{align*}
$$
[[Ch06#^eq6-17a|Eq. (6.17a)]] [[Ch06#^eq6-17b|Eq. (6.17b)]]

Taking the short- and long-maturity limits of (6.17b) exposes the model's second empirical weakness. For $kT\ll1$:

$$
\begin{align*}
\widehat{\sigma}_{F_{T} T} & =\sqrt{V}\left(1+\frac{\rho \sigma T}{8}\right)  \tag{6.18a}\\
\mathcal{S}_{T} & =\frac{\rho \sigma}{4 \sqrt{V}}=\frac{\rho \sigma}{4 \widehat{\sigma}_{F_{T} T}} \tag{6.18b}
\end{align*}
$$
[[Ch06#^eq6-18a|Eq. (6.18a)]] [[Ch06#^eq6-18b|Eq. (6.18b)]]

so at short maturities the ATMF skew is predicted to be *inversely proportional* to the level of ATMF volatility. This is a sharp, falsifiable claim, and Section 6.6.4 tests it directly against Euro Stoxx 50 data: rather than the predicted inverse relationship, skew and volatility are observed to be, if anything, *positively* correlated — it is "unreasonable to hard-wire" the inverse dependence the model insists on. For $kT\gg1$, by contrast:

$$
\begin{align*}
\widehat{\sigma}_{F_{T} T} & =\sqrt{V^{0}}\left(1+\frac{\rho \sigma}{4 k}\right)+\frac{1}{2 k T}\left(\frac{V-V^{0}}{\sqrt{V^{0}}}+\frac{\rho \sigma}{4 k} \frac{V-3 V^{0}}{\sqrt{V^{0}}}\right)  \tag{6.19a}\\
\mathcal{S}_{T} & =\frac{\rho \sigma}{2 \sqrt{V^{0}}} \frac{1}{k T} \tag{6.19b}
\end{align*}
$$
[[Ch06#^eq6-19a|Eq. (6.19a)]] [[Ch06#^eq6-19b|Eq. (6.19b)]]

and in the special case where the term structure of VS volatilities is exactly flat ($V=V^0$), (6.16) collapses to a single closed-form skew term structure spanning all maturities:

$$
\begin{equation*}
\left.\frac{d \widehat{\sigma}_{K T}}{d \ln K}\right|_{F_{T}}=\frac{\rho \sigma}{2 \sqrt{V^{0}}} \frac{k T+e^{-k T}-1}{(k T)^{2}} \tag{6.20}
\end{equation*}
$$
[[Ch06#^eq6-20|Eq. (6.20)]]

Fitted against real index skew term structures, this single-time-scale form cannot match both the short and long end simultaneously — real market skews decay like a power law, not the shape (6.20) forces. This sets up the chapter's real conclusion, deliberately reframed from what one might expect: what disqualifies Heston for trading exotics is not primarily its imprecise fit to historical dynamics, but its *lack of flexibility* — a one-factor model simply has too few free parameters to be independently recalibrated to whatever forward-skew or vol-of-vol view a desk actually needs to hold, especially since forward-skew risk (Section 3.1.7) cannot itself be hedged by trading vanilla options. The chapter closes on a pointed rhetorical question that motivates everything that follows: which practical pricing or hedging issue actually calls for the specific SDE (6.5)?

#### Chapter's Digest and Logical Bridge to Chapter 7

Heston is tractable and pedagogically useful, but its single time scale $1/k$ structurally caps its flexibility — in vol-of-vol term structure, in ATMF skew term structure, and in the (empirically false) skew/volatility-level relationship it forces. Chapter 7 answers the chapter's closing question directly: it builds the general class of forward-variance models — with $N$ independent factors and independently-tunable time scales — of which Heston is only the most restrictive, one-factor member.

---
### Chapter 7: Forward variance models

Chapter 6's diagnosis of Heston left a specific, actionable complaint: a single mean-reversion speed $k$ forces vol-of-vol term structure, ATMF skew term structure, and the skew/volatility-level relationship all to share one time scale, and no amount of recalibration can decouple them. Chapter 7 answers that complaint by building the general class of models in which the state variable is not an instantaneous variance $V_t$ but the whole forward variance curve $\xi_t^T$ — fulfilling the promise, made at the end of Chapter 4, that this curve is the natural object to model directly, since it is what a trader can actually observe and hedge with, maturity by maturity.

The chapter opens exactly as Chapter 1 did: by writing down the carry P&L of a hedged position and reading off what break-even levels make it driftless. With a book of vanilla and exotic options hedged in delta and vega against every point of the variance curve, the P&L over $\delta t$, expanded to order one in $\delta t$ and order two in $\delta S,\delta\xi$, reads

$$
\begin{align*}
P \& L= & -\frac{d P}{d t} \delta t+r P \delta t-(r-q) S \frac{d P}{d S} \delta t \\
& -\frac{S^{2}}{2} \frac{d^{2} P}{d S^{2}} \frac{\delta S^{2}}{S^{2}}-\frac{1}{2} \int_{t}^{T} d u \int_{t}^{T} d u^{\prime} \frac{d^{2} P}{\delta \xi^{u} \delta \xi^{u^{\prime}}} \delta \xi^{u} \delta \xi^{u^{\prime}}-\int_{t}^{T} d u S \frac{d^{2} P}{d S \delta \xi^{u}} \frac{\delta S}{S} \delta \xi^{u} \tag{7.1}
\end{align*}
$$
[[Ch07#^eq7-1|Eq. (7.1)]]

Because $\xi$ is a whole curve rather than a single number, the position now carries one gamma-type term per pair of maturities $(u,u')$ and one cross term per maturity $u$. Break-even levels for these random quadratic terms are given names — a spot/variance covariance $\mu$ and a variance/variance covariance $\nu$ — defined by

$$
\begin{align*}
\mu(t, u) \delta t & =\left\langle\frac{\delta S}{S} \delta \xi^{u}\right\rangle_{t}  \tag{7.2a}\\
\nu\left(t, u, u^{\prime}\right) \delta t & =\left\langle\delta \xi^{u} \delta \xi^{u^{\prime}}\right\rangle_{t} \tag{7.2b}
\end{align*}
$$
[[Ch07#^eq7-2a|Eq. (7.2a)]] [[Ch07#^eq7-2b|Eq. (7.2b)]]

exactly as Chapter 1's single break-even variance $\sigma^2$ was defined by demanding the gamma term be driftless. Demanding that the whole random part of the P&L vanish at these break-even levels pins down the desired form

$$
\begin{align*}
P \& L= & -\frac{S^{2}}{2} \frac{d^{2} P}{d S^{2}}\left(\frac{\delta S^{2}}{S^{2}}-\sigma(t, S, \xi)^{2} \delta t\right)  \tag{7.3a}\\
& -\frac{1}{2} \int_{t}^{T} d u \int_{t}^{T} d u^{\prime} \frac{d^{2} P}{\delta \xi^{u} \delta \xi^{u^{\prime}}}\left(\delta \xi^{u} \delta \xi^{u^{\prime}}-\nu\left(t, u, u^{\prime}, \xi\right) \delta t\right)  \tag{7.3b}\\
& -\int_{t}^{T} d u S \frac{d^{2} P}{d S \delta \xi^{u}}\left(\frac{\delta S}{S} \delta \xi^{u}-\mu(t, u, \xi) \delta t\right) \tag{7.3c}
\end{align*}
$$
[[Ch07#^eq7-3a|Eq. (7.3a)]] [[Ch07#^eq7-3b|Eq. (7.3b)]] [[Ch07#^eq7-3c|Eq. (7.3c)]]

with the gamma break-even $\sigma(t,S,\xi)^2$ pinned to the instantaneous forward variance $\xi^t$ itself, since a short-maturity variance swap is the natural hedge for the spot gamma term. Reading off the coefficient of $\delta t$ in (7.3) gives the chapter's founding object: an infinite-dimensional pricing PDE in state variables $(t,S,\xi)$,

$$
\begin{align*}
& \frac{d P}{d t}+(r-q) S \frac{d P}{d S}+\frac{\xi^{t}}{2} S^{2} \frac{d^{2} P}{d S^{2}} \\
& +\frac{1}{2} \int_{t}^{T} d u \int_{t}^{T} d u^{\prime} \nu\left(t, u, u^{\prime}, \xi\right) \frac{d^{2} P}{\delta \xi^{u} \delta \xi^{u^{\prime}}}+\int_{t}^{T} d u \mu(t, u, \xi) S \frac{d^{2} P}{d S \delta \xi^{u}}=r P \tag{7.4}
\end{align*}
$$
[[Ch07#^eq7-4|Eq. (7.4)]]

with functional derivatives taken with respect to the whole curve $\xi^u$ rather than a finite list of variables. This PDE is exact, but it is also, as stated, useless: nobody can solve a PDE with functional derivatives, or trade a continuum of hedge instruments. Its probabilistic dual is no better off — the corresponding SDE system for $(S_t,\xi_t^u)$ has one equation per maturity $u$, with $\mu,\nu$ reappearing as instantaneous covariations,

$$
\begin{align*}
\lim _{d t \rightarrow 0} \frac{1}{d t} E_{t}\left[d \ln S_{t} d \xi_{t}^{u}\right] & =\sqrt{\xi_{t}^{t}} \lambda_{t}^{u} \frac{1}{d t} E_{t}\left[d W_{t}^{S} d W_{t}^{u}\right]=\mu(t, u, \xi)  \tag{7.5}\\
\lim _{d t \rightarrow 0} \frac{1}{d t} E_{t}\left[d \xi_{t}^{u} d \xi_{t}^{u^{\prime}}\right] & =\lambda_{t}^{u} \lambda_{t}^{u^{\prime}} \frac{1}{d t} E_{t}\left[d W_{t}^{u} d W_{t}^{u^{\prime}}\right]=\nu\left(t, u, u^{\prime}, \xi\right) \tag{7.6}
\end{align*}
$$
[[Ch07#^eq7-5|Eq. (7.5)]] [[Ch07#^eq7-6|Eq. (7.6)]]

Section 7.2 confronts this uselessness head-on by asking when the family $\{\xi_t^u\}$ collapses to a small, finite set of Markov state variables. The natural starting point, motivated by the historical evidence for lognormal VS-volatility behavior surveyed in Chapter 6, is a lognormal, time-homogeneous ansatz for each forward variance separately,

$$
\begin{equation*}
d \xi_{t}^{T}=\omega(T-t) \xi_{t}^{T} d W_{t}^{T} \tag{7.7}
\end{equation*}
$$
[[Ch07#^eq7-7|Eq. (7.7)]]

whose solution is

$$
\begin{equation*}
\ln \left(\xi_{t}^{T}\right)=\ln \left(\xi_{0}^{T}\right)-\frac{1}{2} \int_{0}^{t} \omega^{2}(T-\tau) d \tau+\int_{0}^{t} \omega(T-\tau) d W_{\tau}^{T} \tag{7.8}
\end{equation*}
$$
[[Ch07#^eq7-8|Eq. (7.8)]]

Even driving every $\xi^T$ with the *same* Brownian motion $W_t$ does not, in general, collapse the curve to a finite state: (7.8) expresses $\ln\xi^T$ as a weighted average of the increments of $W_t$ with a $T$-dependent weight $\omega(T-t)$, so simulating the curve at time $t$ requires the full path of $W_t$ up to $t$ — every $\xi^T$ has to be tracked individually, which is exactly the intractability (7.4) already displayed. The chapter's founding technical insight is that this collapses only for one specific parametric family of vol-of-vol term structures: an exponentially decaying $\omega$,

$$
\begin{align*}
\omega(u) & =\omega e^{-k u}  \tag{7.9}\\
\int_{0}^{t} \omega(T-\tau) d W_{\tau} & =e^{-k T} \int_{0}^{t} e^{k \tau} d W_{\tau}
\end{align*}
$$
[[Ch07#^eq7-9|Eq. (7.9)]]

which factors the $T$-dependence out of the integral entirely, so that a single scalar quantity — the value of the integral at time $t$ — suffices to reconstruct the whole curve. Markov tractability, in other words, is not free: it is bought by imposing this specific exponential form, equivalent to driving $\xi_t^T$ with one Ornstein-Uhlenbeck (OU) process $X_t$, giving the explicit one-factor solution

$$
\begin{equation*}
\xi_{t}^{T}=\xi_{0}^{T} \exp \left(\omega e^{-k(T-t)} X_{t}-\frac{\omega^{2}}{2} e^{-2 k(T-t)} E\left[X_{t}^{2}\right]\right) \tag{7.10}
\end{equation*}
$$
[[Ch07#^eq7-10|Eq. (7.10)]]

with $\omega$ the lognormal volatility of the vanishing-maturity forward variance $\xi_t^{T=t}$. Section 7.3 generalizes this to $N$ independent OU factors, each with its own mean-reversion speed $k_i$ and weight $w_i$, giving

$$
\begin{equation*}
d X_{t}^{i}=-k_{i} X_{t}^{i} d t+d W_{t}^{i}, X_{t=0}^{i}=0 \tag{7.14}
\end{equation*}
$$
[[Ch07#^eq7-14|Eq. (7.14)]]

Being driven by OU processes, this $N$-factor model is exactly and efficiently simulable in closed form at any time step — a decisive practical advantage over Heston's non-central chi-square variance process, which admits no such shortcut. Section 7.3.2 works out the resulting closed-form instantaneous volatility and correlation structure of the forward variances themselves, and, for the practically important case of a flat term structure of VS volatilities, the instantaneous volatility of a VS volatility $\widehat\sigma_T$ collapses to

$$
\begin{equation*}
\nu_{T}(t)=\nu \alpha_{w} \sqrt{\sum_{i j} w_{i} w_{j} \rho_{i j} I\left(k_{i}(T-t)\right) I\left(k_{j}(T-t)\right)} \tag{7.24}
\end{equation*}
$$
[[Ch07#^eq7-24|Eq. (7.24)]]

where

$$
\begin{equation*}
I(x)=\frac{1-e^{-x}}{x} \tag{7.25}
\end{equation*}
$$
[[Ch07#^eq7-25|Eq. (7.25)]]

This closed form — a weighted combination of the individual factors' $I(k_i(T-t))$ profiles — is what will be fitted to a target vol-of-vol term structure further on. Section 7.3.3 then makes a point easy to get backwards: even though an $N$-factor model has only $N$ independent sources of randomness, a vega hedge must still immunize the position, at order one, against *all* deformations $\delta\xi^T$ of the curve — not merely those spanned by the model's own $N$ factors — because the deltas $dP/d\xi^T$ that the pricing PDE (7.4) generates are computed independently of whatever covariance structure a particular model happens to impose. Only if the whole set of deltas is traded does the position actually realize, during $\delta t$, the gamma/theta P&L with break-even levels $\mu,\nu$ that the pricing equation assumes. The number of factors sets the *rank* of the break-even covariance matrix the model generates — not which instruments must be hedged against. This is spelled out explicitly by writing the $N$-factor model's covariance functions,

$$
\begin{align*}
\mu(t, u, \xi) & =\omega \alpha_{w} \sqrt{\xi_{t}^{t}} \xi_{t}^{u} \sum_{i} \rho_{S X^{i}} w_{i} e^{-k_{i}(u-t)}  \tag{7.27a}\\
\nu\left(t, u, u^{\prime}, \xi\right) & =\omega^{2} \alpha_{w}^{2} \xi_{t}^{u} \xi_{t}^{u^{\prime}} \sum_{i j} \rho_{i j} w_{i} w_{j} e^{-k_{i}(u-t)} e^{-k_{j}\left(u^{\prime}-t\right)} \tag{7.27b}
\end{align*}
$$
[[Ch07#^eq7-27a|Eq. (7.27a)]] [[Ch07#^eq7-27b|Eq. (7.27b)]]

which are perfectly well-defined for every pair of maturities $(u,u')$, whether or not $u,u'$ coincide with any of the model's own characteristic time scales $1/k_i$.

Section 7.4 now asks the practical question Chapter 6 left hanging: how many factors, and what time scales, does a genuinely usable model need? A one-factor model, it turns out, reproduces Heston's own vol-of-vol term-structure formula (6.9) exactly — inheriting its single-timescale flaw wholesale, since $\nu_T(t)=\nu I(k(T-t))$ is literally the same functional form. Two OU factors, with mean-reversion speeds $k_1,k_2$, a correlation $\rho_{12}$, and a mixing parameter $\theta$, resolve this:

$$
\begin{align*}
d \xi_{t}^{T} & =(2 \nu) \xi_{t}^{T} \alpha_{\theta}\left((1-\theta) e^{-k_{1}(T-t)} d W_{t}^{1}+\theta e^{-k_{2}(T-t)} d W_{t}^{2}\right)  \tag{7.28}\\
\alpha_{\theta} & =1 / \sqrt{(1-\theta)^{2}+\theta^{2}+2 \rho_{12} \theta(1-\theta)} \tag{7.29}
\end{align*}
$$
[[Ch07#^eq7-28|Eq. (7.28)]] [[Ch07#^eq7-29|Eq. (7.29)]]

Introducing the driftless Gaussian mixing process that the two factors combine into,

$$
\begin{equation*}
x_{t}^{T}=\alpha_{\theta}\left[(1-\theta) e^{-k_{1}(T-t)} X_{t}^{1}+\theta e^{-k_{2}(T-t)} X_{t}^{2}\right] \tag{7.30}
\end{equation*}
$$
[[Ch07#^eq7-30|Eq. (7.30)]]

with quadratic variation

$$
\begin{align*}
\left\langle\left(d x_{t}^{T}\right)^{2}\right\rangle & =\eta^{2}(T-t) d t  \tag{7.31a}\\
\eta(u) & =\alpha_{\theta} \sqrt{(1-\theta)^{2} e^{-2 k_{1} u}+\theta^{2} e^{-2 k_{2} u}+2 \rho_{12} \theta(1-\theta) e^{-\left(k_{1}+k_{2}\right) u}} \tag{7.31b}
\end{align*}
$$
[[Ch07#^eq7-31a|Eq. (7.31a)]] [[Ch07#^eq7-31b|Eq. (7.31b)]]

turns (7.28) into the compact SDE $d\xi_t^T=(2\nu)\xi_t^T dx_t^T$, whose Markov-functional solution is

$$
\begin{align*}
\xi_{t}^{T} & =\xi_{0}^{T} f^{T}\left(t, x_{t}^{T}\right)  \tag{7.33}\\
f^{T}(t, x) & =e^{\omega x-\frac{\omega^{2}}{2} \chi(t, T)} \tag{7.34}
\end{align*}
$$
[[Ch07#^eq7-33|Eq. (7.33)]] [[Ch07#^eq7-34|Eq. (7.34)]]

expressing $\xi_t^T$ as a function of the single Gaussian process $x_t^T$ — a genuine Markov-functional model, with the mapping function $f^T$ presently just an exponential (so forward variances remain lognormal), though Section 7.7.1 will later relax this. Taking $k_1>k_2$ and calling $X^1,X^2$ the short and long factors respectively, the resulting instantaneous volatility of a VS volatility works out, for a flat term structure, to the closed form

$$
\begin{equation*}
\nu_{T}(t)=\nu \alpha_{\theta} \sqrt{(1-\theta)^{2} A_{1}^{2}+\theta^{2} A_{2}^{2}+2 \rho_{12} \theta(1-\theta) A_{1} A_{2}} \tag{7.39}
\end{equation*}
$$
[[Ch07#^eq7-39|Eq. (7.39)]]

with $A_i=I(k_i(T-t))$. Because $k_1,k_2$ and $\theta$ are three genuinely independent knobs, this two-time-scale $\nu_T(t)$ can be tuned to match a realistic market benchmark — a power-law vol-of-vol term structure

$$
\begin{equation*}
\nu_{T}^{\mathrm{B}}(t)=\sigma_{0}\left(\frac{\tau_{0}}{T-t}\right)^{\alpha} \tag{7.40}
\end{equation*}
$$
[[Ch07#^eq7-40|Eq. (7.40)]]

over a wide range of maturities, from one month to five years, with an exponent $\alpha$ typically between 0.3 and 0.6 for equity indexes — something a single time scale $1/k$ is structurally incapable of doing, since $I(k(T-t))$ has only one inflection region, not the smoothly decaying power-law shape (7.40) requires across decades of maturity. This is the chapter's central worked resolution of Chapter 6's diagnosed weakness: the flexibility that was missing is not a vague qualitative improvement but an explicit, parametrized, curve-fittable one.

Numerically, the three sets are calibrated with widely separated mean-reversion speeds — a short factor with $k_1$ ranging from roughly 2.5 to 7.5 per year and a long factor with $k_2$ well under 0.5 — combined with a mixing weight $\theta$ and a correlation $\rho_{12}$ that can be positive, negative, or zero without materially degrading the fit to (7.40). This is itself worth pausing on: the same target vol-of-vol curve is reproducible from qualitatively different underlying factor correlations, which foreshadows the chapter's later point that matching spot-starting vol-of-vol says very little about the model's forward vol-of-vol or VIX-future behavior.

Sections 7.4.2–7.5 are candid about what this added flexibility does not buy. Correlations of forward variances in the two-factor model still collapse to a single time scale $1/(k_1-k_2)$ — the relevant scale for correlation is a *difference* of mean-reversion speeds, not either speed on its own — which is exactly what would motivate a third factor in a still more elaborate model. And "calibrating" the model to the VS or ATMF volatility term structure is explicitly distinguished from calibrating to the vanilla smile: the former merely inputs the model's own state variables $\xi_0^T$ and is barely a calibration at all, while the latter, addressed at length in Chapter 8, is where the model's genuine explanatory content is tested.

This distinction matters practically because it determines what a trading desk is actually agreeing to when it says a model is "calibrated." Inputting the VS curve as $\xi_0^T$ commits a desk only to today's term structure, leaving the vol-of-vol and correlation parameters $\nu,\theta,k_1,k_2,\rho_{12}$ free to be set from a genuinely independent view — typically the desk's own read of the vol-of-vol and skew markets — rather than forced by whatever the vanilla smile alone would imply.

Section 7.6 turns to a full worked example — pricing and hedging options on realized variance — that puts the machinery to concrete use. The key observation is that although the realized quadratic variation $Q_t$ itself is not directly hedgeable, a specific combination of $Q_t$ and a VS mark-to-market is:

$$
\begin{equation*}
U_{t}=\frac{Q_{t}+(T-t) \widehat{\sigma}_{T}^{2}(t)}{T} \tag{7.43}
\end{equation*}
$$
[[Ch07#^eq7-43|Eq. (7.43)]]

is exactly replicated, over any interval, by holding $\tfrac{T-t}{T}$ variance-swap contracts of maturity $T$ — a position with no cash outlay and hence no pricing drift, making $U_t$ a legitimate driftless underlying in its own right. In terms of forward variances, this reads

$$
\begin{equation*}
U_{t}=\frac{1}{T}\left(\int_{0}^{t} \xi_{\tau}^{\tau} d \tau+\int_{t}^{T} \xi_{t}^{\tau} d \tau\right) \tag{7.44}
\end{equation*}
$$
[[Ch07#^eq7-44|Eq. (7.44)]]

with the payoff of the call on realized variance reducing simply to $(U_T-\widehat{\boldsymbol\sigma}^2)^+$. Differentiating (7.43) and keeping only the driftless part gives

$$
\begin{equation*}
d U_{t}=\frac{T-t}{T} \underbrace{d\left(\widehat{\sigma}_{T}^{2}(t)\right)}_{\text {Diffusive portion }} \tag{7.45}
\end{equation*}
$$
[[Ch07#^eq7-45|Eq. (7.45)]]

and, assuming $\widehat\sigma_T(t)$ is lognormal with instantaneous volatility $\nu_T(t)$, this becomes an SDE for $U_t$ itself,

$$
\begin{equation*}
\frac{d U_{t}}{U_{t}}=2 R_{t} \frac{T-t}{T} \nu_{T}(t) d W_{t} \tag{7.46}
\end{equation*}
$$
[[Ch07#^eq7-46|Eq. (7.46)]]

with the ratio

$$
\begin{equation*}
R_{t}=\frac{\widehat{\sigma}_{T}^{2}(t)}{U_{t}}=\frac{T \widehat{\sigma}_{T}^{2}(t)}{Q_{t}+(T-t) \widehat{\sigma}_{T}^{2}(t)}=\frac{T \widehat{\sigma}_{T}^{2}(t)}{t \sigma_{r}^{2}(t)+(T-t) \widehat{\sigma}_{T}^{2}(t)} \tag{7.47}
\end{equation*}
$$
[[Ch07#^eq7-47|Eq. (7.47)]]

The presence of $R_t$ — which mixes past realized volatility with future implied volatility — is precisely what prevents (7.46) from being simply lognormal, and its quadratic variation over $[t,T]$,

$$
\begin{equation*}
\int_{t}^{T} 4 R_{\tau}^{2}\left(\frac{T-\tau}{T}\right)^{2} \nu_{T}^{2}(\tau) d \tau \tag{7.48}
\end{equation*}
$$
[[Ch07#^eq7-48|Eq. (7.48)]]

is not autonomous in $U_t$ alone. But $R_\tau$ is a ratio of two integrated variances and so is a number of order one; freezing it at its value $R_t$ on the pricing date turns $U_t$ into a genuinely lognormal process and yields the "simple model" (SM)'s closed-form Black-Scholes-type price,

$$
\begin{align*}
P(t, U) & =P_{\mathrm{BS}}\left(t, U, \sigma_{\mathrm{eff}}, T\right)  \tag{7.49a}\\
\sigma_{\mathrm{eff}}^{2} & =\frac{1}{T-t} \int_{t}^{T} 4 R_{t}^{2}\left(\frac{T-\tau}{T}\right)^{2} \nu_{T}^{2}(\tau) d \tau \tag{7.49b}
\end{align*}
$$
[[Ch07#^eq7-49a|Eq. (7.49a)]] [[Ch07#^eq7-49b|Eq. (7.49b)]]

The remarkable practical conclusion, checked numerically against a full Monte Carlo simulation of the two-factor model, is that options on realized variance are simple instruments after all: their only volatility-of-volatility risk is exposure to the realized volatility of the single VS volatility $\widehat\sigma_T$ for the option's own maturity — the detailed correlation structure of the underlying forward variance curve is almost irrelevant to the price, as long as $\nu_T(\tau)$ itself is matched. This is, in a real sense, Chapter 7's answer to "does all this factor-model machinery actually matter for pricing?" — mostly not, for spot-starting options on realized variance, provided the curve $\nu_T(\tau)$ is right.

This is checked concretely by pricing an at-the-money call on realized variance in each of the three parameter sets of Table 7.1 and comparing to the SM: prices computed with Set I, Set II, Set III — which differ substantially in correlation $\rho_{12}$ and hence in the correlation structure of forward variances — agree with each other and with the SM to within a few basis points, despite the underlying dynamics being visibly different in each case.

Two refinements matter more for risk management than for price. First, the baseline approximation $R_\tau\equiv R_t$ implicitly assumes a flat VS term structure; for a genuinely non-flat one, $R_\tau$ needs the amended approximation

$$
\begin{equation*}
R_{\tau} \simeq \frac{E_{t}\left[\widehat{\sigma}_{T}^{2}(\tau)\right]}{E_{t}\left[U_{\tau}\right]}=\frac{T \widehat{\sigma}_{\tau T}^{2}(t)}{Q_{t}+(T-t) \widehat{\sigma}_{T}^{2}(t)} \tag{7.54}
\end{equation*}
$$
[[Ch07#^eq7-54|Eq. (7.54)]]

which feeds into an amended effective volatility

$$
\begin{equation*}
\sigma_{\mathrm{eff}}^{2}=\frac{4}{T-t} \int_{t}^{T}\left(\frac{T-\tau}{T}\right)^{2}\left(\frac{T \widehat{\sigma}_{\tau T}^{2}(t)}{Q_{t}+(T-t) \widehat{\sigma}_{T}^{2}(t)}\right)^{2} \nu_{T}^{2}(\tau) d \tau \tag{7.55}
\end{equation*}
$$
[[Ch07#^eq7-55|Eq. (7.55)]]

which now depends on the full VS term structure up to $T$, not just $\widehat\sigma_T$ alone — and this dependence must itself be hedged with intermediate VSs, since the option's mark-to-market P&L is now exposed to $\widehat\sigma_\tau$ for every $\tau\in[0,T]$, not merely the terminal VS. Rewritten in a suggestive form,

$$
\begin{equation*}
\sigma_{\mathrm{eff}}^{2}=\frac{4}{T-t} \int_{t}^{T}\left(\frac{\left(Q_{t}+(T-t) \widehat{\sigma}_{T}^{2}(t)\right)-\left(Q_{t}+(\tau-t) \widehat{\sigma}_{\tau}^{2}(t)\right)}{Q_{t}+(T-t) \widehat{\sigma}_{T}^{2}(t)}\right)^{2} \nu_{T}^{2}(\tau) d \tau \tag{7.59}
\end{equation*}
$$
[[Ch07#^eq7-59|Eq. (7.59)]]

reveals that the intermediate VS hedge that offsets the exposure of $\sigma_{\text{eff}}^2$ to the term structure is exactly the same position that offsets the exposure of $R_\tau$ to the accrued quadratic variation $Q_t$: the vega hedge doubles, automatically, as the gamma/theta hedge. For a flat term structure this collapses to a clean statement about aggregate vega,

$$
\begin{equation*}
\int_{0}^{T} \lambda(\tau) d \tau+\mu_{T}=0 \tag{7.60}
\end{equation*}
$$
[[Ch07#^eq7-60|Eq. (7.60)]]

— the total vega of the intermediate-VS hedge exactly cancels, leaving only the vega of the terminal-maturity VS that hedges $U_t$ itself. Second, quadratic variation does not, in reality, accrue continuously: it is a discrete sum of squared daily log-returns,

$$
\begin{align*}
r & =\sigma \sqrt{\Delta t} Z-\frac{\sigma^{2} \Delta t}{2}  \tag{7.64a}\\
& \simeq \sigma \sqrt{\Delta t} Z \tag{7.64b}
\end{align*}
$$
[[Ch07#^eq7-64a|Eq. (7.64a)]] [[Ch07#^eq7-64b|Eq. (7.64b)]]

and this discreteness contributes an extra, purely statistical piece to $\sigma_{\text{eff}}^2$,

$$
\begin{equation*}
\sigma_{\mathrm{eff}}^{2}=\frac{1}{T} \int_{0}^{T} 4\left(\frac{T-\tau}{T}\right)^{2}\left(\frac{\widehat{\sigma}_{\tau T}^{2}(0)}{\widehat{\sigma}_{T}^{2}(0)}\right)^{2} \nu_{T}^{2}(\tau) d \tau+\frac{2+\kappa}{N T} \tag{7.70}
\end{equation*}
$$
[[Ch07#^eq7-70|Eq. (7.70)]]

driven by the conditional excess kurtosis $\kappa$ of daily returns and the number of return observations $N=T/\Delta t$ — largest, relative to the diffusive piece, for the shortest maturities, where there are fewest returns to average over.

A close cousin of the spot-starting option is the option on *forward* realized variance, priced and hedged the same way but starting its clock at a future date $T_1$ rather than today,

$$
\begin{align*}
\sigma_{\text {eff }}^{2} & =\frac{4}{T_{2}-T_{1}} \int_{T_{1}}^{T_{2}}\left(\frac{T_{2}-\tau}{T_{2}-T_{1}}\right)^{2}\left(\frac{\widehat{\sigma}_{\tau T_{2}}^{2}\left(T_{1}\right)}{\widehat{\sigma}_{T_{2}}^{2}\left(T_{1}\right)}\right)^{2} \nu_{T_{2}}^{2}(\tau) d \tau \\
P\left(t=T_{1}\right) & =P_{\mathrm{BS}}\left(t=T_{1}, \widehat{\sigma}_{T_{2}}^{2}\left(T_{1}\right), \sigma_{\text {eff }}, T_{2}\right) \tag{7.82}
\end{align*}
$$
[[Ch07#^eq7-82|Eq. (7.82)]]

Unlike the spot-starting case, this price genuinely depends on the volatility of *forward* VS volatilities over $[0,T_1]$ — information that cannot be backed out of the curve $\nu_T(\tau)$ alone, and that only a full multi-factor model, not the SM, can supply. This is precisely the distinction Section 7.4 drew between parameter sets matched on spot-starting vol-of-vol but differing on forward vol-of-vol: it is exactly forward-starting instruments like this one where that distinction becomes economically visible, not merely academic.

Sections 7.7–7.8 close the chapter with two remaining classes of instrument: VIX futures and options, and discrete-time forward variance models. The VIX future's settlement value is a 30-day log-contract implied volatility computed from listed S&P 500 options,

$$
\begin{equation*}
\widehat{\sigma}_{\mathrm{VIX}, T_{i}}^{2}\left(T_{i}\right)=\frac{2 e^{r \Delta}}{\Delta}\left(\int_{0}^{F_{T_{i}}^{T_{i}+\Delta}} P_{\mathrm{mkt}, T_{i}}^{K, T_{i}+\Delta} \frac{d K}{K^{2}}+\int_{F_{T_{i}}^{T_{i}+\Delta}}^{\infty} C_{\mathrm{mkt}, T_{i}}^{K, T_{i}+\Delta} \frac{d K}{K^{2}}\right) \tag{7.83}
\end{equation*}
$$
[[Ch07#^eq7-83|Eq. (7.83)]]

and, being the value of a future, it must be driftless prior to expiry:

$$
\begin{equation*}
F_{t}^{i}=E\left[F_{T_{i}}^{i}\right] \tag{7.86}
\end{equation*}
$$
[[Ch07#^eq7-86|Eq. (7.86)]]

But VIX futures display pronounced, strongly positively-sloped smiles — visibly at odds with the near-lognormality of variance swaptions established earlier in the chapter — so the lognormal ansatz (7.34) needs to be relaxed while preserving the Markov-functional structure. Driftlessness of $f^T(t,x)$ translates into a PDE,

$$
\begin{equation*}
\frac{d f^{T}}{d t}+\frac{\eta^{2}(T-t)}{2} \frac{d^{2} f^{T}}{d x^{2}}=0 \tag{7.87}
\end{equation*}
$$
[[Ch07#^eq7-87|Eq. (7.87)]]

whose solution, for any admissible terminal profile, produces a general local-volatility-type instantaneous volatility function for $\xi_t^T$,

$$
\begin{equation*}
\operatorname{vol}\left(\xi_{t}^{T}\right)=\eta(T-t)\left|\frac{d \ln f^{T}}{d x}\right|_{x=f^{T^{-1}}\left(\xi_{t}^{T}, t\right)} \tag{7.88}
\end{equation*}
$$
[[Ch07#^eq7-88|Eq. (7.88)]]

A simple two-exponential mixture — positively sloping because a linear combination of positive exponentials is dominated, at large $x$, by whichever term has the larger effective vol-of-vol — supplies exactly the flexibility needed:

$$
\begin{equation*}
f^{T}(t, x)=\left(1-\gamma_{T}\right) e^{\omega_{T} x-\frac{\omega_{T}^{2}}{2} \chi(t, T)}+\gamma_{T} e^{\beta_{T} \omega_{T} x-\frac{\left(\beta_{T} \omega_{T}\right)^{2}}{2} \chi(t, T)} \tag{7.89}
\end{equation*}
$$
[[Ch07#^eq7-89|Eq. (7.89)]]

with mixing parameter $\gamma_T\in[0,1]$ and displacement parameter $\beta_T\in[0,1]$, and with $\omega_T$ conveniently reparametrized through a dimensionless scale factor $\zeta_T$,

$$
\begin{equation*}
\omega_{T}=\frac{2 \nu}{\left(1-\gamma_{T}\right)+\beta_{T} \gamma_{T}} \zeta_{T} \tag{7.90}
\end{equation*}
$$
[[Ch07#^eq7-90|Eq. (7.90)]]

Calibrated against real VIX futures and options, this parametrization fits observed smiles well within the bid/offer spread. Different parameter sets, all matched on spot-starting vol-of-vol, still leave a genuine choice: the expected instantaneous variance of the VIX future itself,

$$
\begin{equation*}
\sigma_{t}^{2}=E_{0}\left[\frac{\left\langle\left(d \ln F_{t}^{i}\right)^{2}\right\rangle}{d t}\right] \tag{7.92}
\end{equation*}
$$
[[Ch07#^eq7-92|Eq. (7.92)]]

is a curve in $t$ whose *shape* — not merely its integral — differs across otherwise-equivalent parameter sets, and this shape is the criterion for choosing between them. For path-dependent payoffs on VIX-linked instruments (VIX ETFs and ETNs, for instance) that require the future's value at intermediate dates, not merely at expiry,

$$
\begin{equation*}
F_{t}^{i}=E_{t}\left[F_{T_{i}}^{i}\left(X_{T_{i}}^{1}, X_{T_{i}}^{2}\right)\right] \tag{7.93}
\end{equation*}
$$
[[Ch07#^eq7-93|Eq. (7.93)]]

must be evaluated by two-dimensional quadrature at every observation date — computationally awkward enough that, for payoffs needing frequent observation, Section 7.8.2's discrete models, which parametrize VIX futures directly rather than deriving them from a continuous variance curve, become the practical choice.

Section 7.7.4 closes the VIX discussion with a genuine model-independent result. Since any European payoff on the VIX future $F_{T_i}^i$ can be statically replicated by cash, forwards, and vanilla options on $F_{T_i}^i$, applying this replication to the payoff $f(S)=S^2$,

$$
\begin{equation*}
S^{2}=S_{*}^{2}+2 S_{*}\left(S-S_{*}\right)+2 \int_{0}^{S_{*}}(K-S)^{+} d K+2 \int_{S_{*}}^{\infty}(S-K)^{+} d K \tag{7.97}
\end{equation*}
$$
[[Ch07#^eq7-97|Eq. (7.97)]]

and translating the resulting equality of payoffs into an equality of prices yields a consistency condition linking the S&P 500's own forward log-contract volatility directly to market prices of VIX futures and options,

$$
\begin{equation*}
\widehat{\sigma}_{T_{i}, T_{i+\Delta}}^{2}(t)=\left(F_{t}^{i}\right)^{2}+2 \int_{0}^{F_{t}^{i}} \mathcal{P}_{t}^{K i} d K+2 \int_{F_{t}^{i}}^{\infty} \mathcal{C}_{t}^{K i} d K \tag{7.98}
\end{equation*}
$$
[[Ch07#^eq7-98|Eq. (7.98)]]

This condition is model-independent — it holds regardless of which forward variance model, if any, is used — yet volatilities derived from the S&P 500 VS market and from the VIX market via (7.98) do not, in practice, coincide, suggesting a genuine arbitrage strategy rather than a modeling artifact: two markets pricing what is, by (7.98), the same quantity, disagree.

In practice, executing on this discrepancy is not trivial: the replication behind (7.98) requires a static position in a continuum of VIX options plus a position in the VIX future, both facing bid/offer spreads and both, unlike the underlying S&P 500 vanilla book, comparatively illiquid at extreme strikes — so the observed gap between the two markets' implied volatilities is a genuine trading opportunity only to the extent that it exceeds these frictions.

Section 7.8 closes the chapter by building discrete-time forward variance models, in which VS volatilities for a fixed schedule of maturities $T_i$, rather than a continuum $\xi_t^\tau$, are the state variables — mirroring the continuous two-factor SDE directly:

$$
\begin{align*}
d \xi_{t}^{i} & =\left(2 \nu_{i}\right) \xi_{t}^{i} \alpha_{\theta_{i}}\left(\left(1-\theta_{i}\right) e^{-k_{1}\left(T_{i}-t\right)} d W_{t}^{1}+\theta_{i} e^{-k_{2}\left(T_{i}-t\right)} d W_{t}^{2}\right)  \tag{7.99}\\
\alpha_{\theta_{i}} & =1 / \sqrt{\left(1-\theta_{i}\right)^{2}+\theta_{i}^{2}+2 \rho_{12} \theta_{i}\left(1-\theta_{i}\right)}
\end{align*}
$$
[[Ch07#^eq7-99|Eq. (7.99)]]

with solution

$$
\begin{align*}
\xi_{t}^{i} & =\xi_{0}^{i} e^{\omega_{i} x_{t}^{T_{i}}-\frac{\omega_{i}^{2}}{2}} \chi\left(t, T_{i}\right)  \tag{7.100}\\
x_{t}^{T_{i}} & =\alpha_{\theta_{i}}\left[\left(1-\theta_{i}\right) e^{-k_{1}\left(T_{i}-t\right)} X_{t}^{1}+\theta_{i} e^{-k_{2}\left(T_{i}-t\right)} X_{t}^{2}\right]
\end{align*}
$$
[[Ch07#^eq7-100|Eq. (7.100)]]

Demanding that the discrete and continuous models generate identical dynamics for forward volatilities, at least for a curve equal to its own initial values, pins down the mapping between the discrete parameters $\theta_i,\nu_i$ and the continuous model's $\theta,\nu$:

$$
\left\{\begin{align*}
\theta_{i} & =\frac{\theta A_{i}^{2}(0)}{\theta A_{i}^{2}(0)+(1-\theta) A_{i}^{1}(0)}  \tag{7.102}\\
\nu_{i} & =\nu \frac{\alpha_{\theta} \theta}{\alpha_{\theta_{i}} \theta_{i}} A_{i}^{2}(0)
\end{align*}\right.
$$
[[Ch07#^eq7-102|Eq. (7.102)]]

This mapping is what lets a book priced with a discrete model for one exotic and a continuous model for another remain internally consistent — both are, in the sense of (7.102), the same underlying dynamics viewed at different time granularities. VIX futures, rather than variances, can also be modeled directly and driven by the identical two-factor SDE, since they too are driftless:

$$
\begin{align*}
d F_{t}^{i} & =\nu_{i} F_{t}^{i} \alpha_{\theta_{i}}\left(\left(1-\theta_{i}\right) e^{-k_{1}\left(T_{i}-t\right)} d W_{t}^{1}+\theta_{i} e^{-k_{2}\left(T_{i}-t\right)} d W_{t}^{2}\right)  \tag{7.103}\\
\alpha_{\theta_{i}} & =1 / \sqrt{\left(1-\theta_{i}\right)^{2}+\theta_{i}^{2}+2 \rho_{12} \theta_{i}\left(1-\theta_{i}\right)}
\end{align*}
$$
[[Ch07#^eq7-103|Eq. (7.103)]]

and, more generally, as a mapping function of a single Gaussian driver,

$$
\begin{equation*}
F_{t}^{i}=F_{0}^{i} f^{i}\left(t, x_{t}^{T_{i}}\right) \tag{7.105}
\end{equation*}
$$
[[Ch07#^eq7-105|Eq. (7.105)]]

subject to the same driftlessness PDE encountered before,

$$
\begin{equation*}
\frac{d f^{i}}{d t}+\frac{\eta_{i}^{2}\left(T_{i}-t\right)}{2} \frac{d^{2} f^{i}}{d x^{2}}=0 \tag{7.106}
\end{equation*}
$$
[[Ch07#^eq7-106|Eq. (7.106)]]

Squaring $f^i$ and solving the same PDE for the squared quantity $g^i=(f^i)^2$ produces a model calibrated to arbitrary, non-arbitrageable VIX smiles in which the entire dynamics reduces to simulating just the two underlying Ornstein-Uhlenbeck processes $X_t^1,X_t^2$ — the discrete-model analogue of the earlier claim that a two-factor model is "just" two OU processes with a nonlinear read-out. The resulting local-volatility function for the VIX future itself,

$$
\begin{equation*}
\operatorname{vol}\left(F^{i}\right)=\eta\left(T_{i}-t\right)\left|\frac{d \ln f^{i}}{d x}\right|_{x=f^{-1}\left(t, F_{t}^{i}\right)} \tag{7.107}
\end{equation*}
$$
[[Ch07#^eq7-107|Eq. (7.107)]]

is the direct discrete-model analogue of (7.88), read off the same driftlessness PDE with $F^i$ in place of $\xi^T$.

Having specified a dynamics for the forward variance curve, Section 7.8.3 asks the remaining question: what dynamics should the underlying $S_t$ itself follow, consistent with these forward volatilities? The requirements are that the VS volatility observed at $T_i$ be recovered exactly, that the distribution of returns over $[T_i,T_{i+1}]$ not depend on the level of $S_{T_i}$ (so that future and spot-starting smiles stay decoupled), and that the shape of the future smile itself be controllable. A path-dependent local volatility dynamics,

$$
\begin{equation*}
d S_{t}=(r-q) S_{t} d t+\sigma^{i}\left(\frac{S_{t}}{S_{T_{i}}}\right) S_{t} d W_{t}^{S} \tag{7.108}
\end{equation*}
$$
[[Ch07#^eq7-108|Eq. (7.108)]]

meets the first two requirements automatically, since $s_t=S_t/S_{T_i}$ satisfies an autonomous SDE with initial condition $s_{T_i}=1$ independent of $S_{T_i}$'s own level. A specific three-parameter functional form for $\sigma^i$,

$$
\begin{equation*}
\sigma^{i}(s)=\sigma_{0}^{i} \frac{n^{i}}{n^{i}-1} \frac{\left(n^{i}-1\right) s^{\beta^{i}-1}}{\left(n^{i}-1\right)+s^{\beta^{i}-1}} \tag{7.109}
\end{equation*}
$$
[[Ch07#^eq7-109|Eq. (7.109)]]

decouples the model's three degrees of freedom cleanly: $\sigma_0^i$ sets the at-the-money level, $n^i$ caps the maximum attainable local volatility, and $\beta^i$ controls the future skew of maturity $T_{i+1}$ directly, since at order one in $(\beta^i-1)$ the ATMF skew reduces to

$$
\begin{equation*}
\left.\left.\frac{d \widehat{\sigma}_{K, T_{i+1}}}{d \ln K}\right|_{S_{T_{i}}} \simeq \frac{1}{2} \frac{d \sigma^{i}}{d \ln S}\right|_{S_{T_{i}}}=\frac{\sigma_{0}^{i}}{2} \frac{n^{i}-1}{n^{i}}\left(\beta^{i}-1\right) \tag{7.110}
\end{equation*}
$$
[[Ch07#^eq7-110|Eq. (7.110)]]

This is the payoff the whole chapter has been building toward: a future-skew scenario, for every interval $[T_i,T_{i+1}]$, that can be dialed in independently, at will, as a function of the VS volatility prevailing at $T_i$ — fixed regardless of volatility level, or tied to it by a power law with an arbitrarily chosen exponent.

Section 7.8.4 makes this last point rigorous by deriving, from Chapter 8's general order-one-in-vol-of-vol ATMF skew formula, the discrete model's own skew decomposition. Writing the skew for maturity $T=N\Delta$ built from $N$ intervals of tenor $\Delta$, and separating the contribution of the intra-interval future skew $\mathcal S_\Delta$ from the covariance of spot with the forward-variance factors, the calculation hinges on a function

$$
\begin{equation*}
\zeta(x, N)=\frac{1}{N^{2}} \sum_{i=0, j>i}^{N-1} \int_{i}^{i+1}\left(e^{-x}\right)^{j-u} d u=\frac{e^{x}-1}{x} \sum_{n=1}^{N-1} \frac{N-n}{N^{2}}\left(e^{-x}\right)^{n} \tag{7.118}
\end{equation*}
$$
[[Ch07#^eq7-118|Eq. (7.118)]]

which packages the accumulated covariance between today's spot return and every later interval's forward variance factor. The final result, at order one in $\nu$ and $(\beta-1)$, splits the ATMF skew for maturity $N\Delta$ into exactly two additive pieces:

$$
\begin{align*}
\mathcal{S}_{N \Delta} & =\frac{1}{N} \mathcal{S}_{\Delta}+\frac{1}{2 \widehat{\sigma}^{3} T^{2}} \Delta \sum_{i=0, j>i}^{N-1} \int_{T_{i}}^{T_{i+1}}\left\langle d \ln S_{t} d \xi_{t}^{j}\right\rangle_{0}  \tag{7.119a}\\
& =\frac{1}{N} \mathcal{S}_{\Delta}+\nu_{\Delta} \alpha_{\theta_{\Delta}}\left(\left(1-\theta_{\Delta}\right) \rho_{S X^{1}} \zeta\left(k_{1} \Delta, N\right)+\theta_{\Delta} \rho_{S X^{2}} \zeta\left(k_{2} \Delta, N\right)\right) \tag{7.119b}
\end{align*}
$$
[[Ch07#^eq7-119a|Eq. (7.119a)]] [[Ch07#^eq7-119b|Eq. (7.119b)]]

a forward-smile term, contributed by the intra-interval skew $\mathcal S_\Delta$ and decaying like $1/N$ — the familiar $1/\sqrt{T}$-type decay of skewness under aggregation of independent returns, echoing both Chapter 3's forward-skew discussion and Chapter 5's Appendix B — and a vol-of-vol covariance term, controlled entirely by $\nu_\Delta$ and the spot/factor correlations $\rho_{SX^1},\rho_{SX^2}$, that survives even when the future skew $\mathcal S_\Delta$ is set to exactly zero.

Taking the continuous-time limit $\Delta\to0$ eliminates the first piece entirely — a continuous model has no intra-interval future skew left over, by construction — and leaves a single order-one-in-vol-of-vol expression,

$$
\begin{equation*}
\mathcal{S}_{T}=\nu \alpha_{\theta}\left((1-\theta) \rho_{S X^{1}} \frac{k_{1} T-\left(1-e^{-k_{1} T}\right)}{\left(k_{1} T\right)^{2}}+\theta \rho_{S X^{2}} \frac{k_{2} T-\left(1-e^{-k_{2} T}\right)}{\left(k_{2} T\right)^{2}}\right) \tag{7.120}
\end{equation*}
$$
[[Ch07#^eq7-120|Eq. (7.120)]]

This single formula is the chapter's explicit, parametrized answer to Chapter 6's flexibility critique. By choosing how the future ATMF skew (7.111) depends on the future ATMF volatility level — a fixed skew, a skew proportional to volatility, or a skew inversely proportional to volatility — the discrete forward variance model can be made to mimic, for short maturities, the two-factor continuous model's own skew/volatility independence, the $3/2$ model's behavior, or Heston's own short-maturity skew/volatility relationship (6.18a) — all three, and anything between them, from the same parametric family. What Chapter 6 found structurally impossible in a one-factor model — decoupling the ATMF skew's dependence on volatility level from the vol-of-vol term structure itself — Chapter 7 shows is not merely possible but explicitly dialable, one parameter at a time.

#### Chapter's Digest and Logical Bridge to Chapter 8

Chapter 7 built the general machinery — the pricing PDE (7.4), its Markov-functional resolution (7.9)/(7.14), the two-factor model (7.28)/(7.29) and its VIX and discrete-time extensions — needed to give forward variance models genuine, independently-tunable control over vol-of-vol term structure, VIX smiles, and ATMF skew/volatility relationships. What the chapter has *not* done is characterize, in general, what determines the shape of the vanilla smile these models produce, beyond the specific order-one skew formula (7.119)–(7.120) worked out for the discrete case. Chapter 8 takes up that question directly: it develops the general perturbative machinery for the smile of stochastic volatility models, of which forward variance models are one instance, and asks what features of a model's covariance structure — not just its vol-of-vol level — actually shape the smile that comes out.

---
### Chapter 8: The smile of stochastic volatility models

Chapter 7 built the machinery — forward variance models with independently tunable vol-of-vol term structure, VIX smiles, and ATMF skew term structure — but left one question conspicuously unanswered: given a model's covariance structure, what actually determines the *shape* of the vanilla smile it produces? Chapter 6 answered this only for Heston's specific instantaneous-variance dynamics. Chapter 8 builds the same machinery once, generally, for any forward-variance model, and shows precisely which functionals of the spot/variance and variance/variance covariance functions $\mu,\nu$ fix the near-the-money smile.

The chapter's very first move is to note that this generality costs nothing: any stochastic volatility model, including instantaneous-variance-based ones like Heston, can be rewritten as a forward variance model, so the pricing PDE already derived in Chapter 7 is the correct universal starting point,

$$
\begin{align*}
& \frac{d P}{d t}+(r-q) S \frac{d P}{d S}+\frac{\xi^{t}}{2} S^{2} \frac{d^{2} P}{d S^{2}} \\
& +\frac{1}{2} \int_{t}^{T} d u \int_{t}^{T} d u^{\prime} \nu\left(t, u, u^{\prime}, \xi\right) \frac{d^{2} P}{\delta \xi^{u} \delta \xi^{u^{\prime}}}+\int_{t}^{T} d u \mu(t, u, \xi) S \frac{d^{2} P}{d S \delta \xi^{u}}=r P \tag{8.1}
\end{align*}
$$
[[Ch08#^eq8-1|Eq. (8.1)]]

not something Heston-specific. Section 8.2 develops an expansion of this PDE in powers of a bookkeeping parameter $\varepsilon$ that scales volatility of volatility, using time-dependent perturbation theory rather than naive term-by-term substitution — the same technique used across quantum mechanics for a Hamiltonian split into an exactly-solvable piece and a small perturbation. Splitting the pricing operator into its Black-Scholes piece and two perturbation operators,

$$
\begin{align*}
H_{t}^{0} & =\frac{\xi^{t}}{2}\left(\partial_{x}^{2}-\partial_{x}\right)  \tag{8.3a}\\
\mathcal{W}_{t}^{1} & =\int_{t}^{T} d u \mu(t, u, \xi) \partial_{x \xi^{u}}^{2}  \tag{8.3b}\\
\mathcal{W}_{t}^{2} & =\frac{1}{2} \int_{t}^{T} d u \int_{t}^{T} d u^{\prime} \nu\left(t, u, u^{\prime}, \xi\right) \partial_{\xi^{u} \xi^{u^{\prime}}}^{2} \tag{8.3c}
\end{align*}
$$
[[Ch08#^eq8-3a|Eq. (8.3a)]] [[Ch08#^eq8-3b|Eq. (8.3b)]] [[Ch08#^eq8-3c|Eq. (8.3c)]]

with $x=\ln S$, and writing the price itself as a power series

$$
\begin{equation*}
P=P_{0}+\varepsilon P_{1}+\varepsilon^{2} P_{2}+\cdots \tag{8.4}
\end{equation*}
$$
[[Ch08#^eq8-4|Eq. (8.4)]]

lets the propagator $U_{st}$ be expanded in a Dyson series in the perturbation, exactly as in time-dependent perturbation theory, with the zeroth-order piece given by the explicit free propagator

$$
\begin{equation*}
U_{s t}^{0}=: \exp \left(\int_{s}^{t} H_{\tau}^{0} d \tau\right):=e^{\frac{1}{2}\left(\int_{s}^{t} \xi^{\tau} d \tau\right)\left(\partial_{x}^{2}-\partial_{x}\right)} \tag{8.11}
\end{equation*}
$$
[[Ch08#^eq8-11|Eq. (8.11)]]

Two structural properties of this free propagator — that $\partial_x$ commutes with it (so $\partial_x^nP_0$ is a martingale, recovering the familiar fact that Black-Scholes deltas and dollar gammas are martingales) and that differentiating with respect to $\xi^u$ picks up an extra gamma-type term exactly when $u$ falls inside the propagation window — do all the work of turning the formal Dyson series into closed-form integral operators acting on the plain Black-Scholes price $P_0$. Carried through order by order, the order-one price correction collapses to a single dimensionless number multiplying a fixed differential operator on $P_0$,

$$
\begin{align*}
C_{t}^{x \xi}(\xi) & =\int_{t}^{T} d \tau \int_{\tau}^{T} d u \mu(\tau, u, \xi)  \tag{8.13a}\\
& =\int_{t}^{T}(T-\tau)\left\langle d \ln S_{\tau} d \widehat{\sigma}_{T}^{2}(\tau)\right\rangle \tag{8.13b}
\end{align*}
$$
[[Ch08#^eq8-13a|Eq. (8.13a)]] [[Ch08#^eq8-13b|Eq. (8.13b)]]

a doubly-integrated spot/variance covariance functional — and the order-two correction requires two further such numbers, a doubly-integrated variance/variance covariance,

$$
\begin{align*}
C_{t}^{\xi \xi}(\xi) & =\int_{t}^{T} d \tau \int_{\tau}^{T} d u \int_{\tau}^{T} d u^{\prime} \nu\left(\tau, u, u^{\prime}, \xi\right)  \tag{8.14a}\\
& =\int_{t}^{T}(T-\tau)^{2}\left\langle d \widehat{\sigma}_{T}^{2}(\tau) d \widehat{\sigma}_{T}^{2}(\tau)\right\rangle \tag{8.14b}
\end{align*}
$$
[[Ch08#^eq8-14a|Eq. (8.14a)]] [[Ch08#^eq8-14b|Eq. (8.14b)]]

and a third, subtler quantity that captures how the spot/variance covariance itself changes as the variance curve moves,

$$
\begin{align*}
D_{t}(\xi) & =\int_{t}^{T} d \tau \int_{\tau}^{T} d u \mu_{\tau u} \frac{\delta C_{\tau}^{x \xi}(\xi)}{\delta \xi^{u}}  \tag{8.17a}\\
& =\int_{t}^{T} d \tau \lim _{d \tau \rightarrow 0} \frac{1}{d \tau} E_{\tau}\left[d \ln S_{\tau} d C_{\tau}^{x \xi}\right]  \tag{8.17b}\\
& =\int_{t}^{T} d \tau \int_{\tau}^{T} d u(T-u) \frac{1}{d \tau}\left\langle d \ln S_{\tau} d\left[\frac{\left\langle d \ln S_{u} d \widehat{\sigma}_{T}^{2}(u)\right\rangle}{d u}\right]\right\rangle \tag{8.17c}
\end{align*}
$$
[[Ch08#^eq8-17a|Eq. (8.17a)]] [[Ch08#^eq8-17b|Eq. (8.17b)]]

Assembling these three numbers gives the chapter's founding result: the full order-two-in-volatility-of-volatility price expansion, valid for *any* forward variance model,

$$
\begin{align*}
& P=\left[1+\varepsilon \frac{C_{0}^{x \xi}\left(\xi_{0}\right)}{2} \partial_{x}\left(\partial_{x}^{2}-\partial_{x}\right)\right. \\
& \left.+\varepsilon^{2}\left(\frac{C_{0}^{\xi \xi}\left(\xi_{0}\right)}{8}\left(\partial_{x}^{2}-\partial_{x}\right)^{2}+\frac{C_{0}^{x \xi}\left(\xi_{0}\right)^{2}}{8} \partial_{x}^{2}\left(\partial_{x}^{2}-\partial_{x}\right)^{2}+\frac{D_{0}\left(\xi_{0}\right)}{2} \partial_{x}^{2}\left(\partial_{x}^{2}-\partial_{x}\right)\right)\right] P_{0} \tag{8.18}
\end{align*}
$$
[[Ch08#^eq8-18|Eq. (8.18)]]

with the subscript $0$ signaling that $C_0^{x\xi}, C_0^{\xi\xi}, D_0$ are evaluated on the variance curve observed at $t=0$. Every detail of the model's joint spot/variance dynamics that matters for the near-the-money price, at this order, is compressed into exactly these three numbers — a genuinely universal reduction, since it holds regardless of which particular forward variance model (one-factor, two-factor, or otherwise) generated $\mu$ and $\nu$ in the first place.

The structure of (8.18) is itself worth pausing on before converting it to implied volatilities: the correction terms take the Hermite-polynomial form

$$
\begin{equation*}
\sum_{n=3}^{\infty} \alpha_{n}\left(\partial_{x}^{n}-\partial_{x}^{2}\right) P_{0} \tag{8.19}
\end{equation*}
$$
[[Ch08#^eq8-19|Eq. (8.19)]]

the same form Chapter 5's cumulant expansion assigns to any perturbation of the lognormal density that leaves the log-contract's implied volatility unchanged — consistent with the fact that, order by order in $\varepsilon$, forward variances remain driftless, so VS implied volatilities never shift as volatility of volatility is dialed up or down. Section 8.3 converts (8.18) into an expansion of implied volatilities themselves, and the result is unexpectedly clean: at order $\varepsilon^2$, implied volatility is *exactly* quadratic in log-moneyness,

$$
\begin{equation*}
\widehat{\sigma}(K, T)=\widehat{\sigma}_{F_{T} T}+\mathcal{S}_{T} \ln \left(\frac{K}{F_{T}}\right)+\frac{\mathcal{C}_{T}}{2} \ln ^{2}\left(\frac{K}{F_{T}}\right)+O\left(\varepsilon^{3}\right) \tag{8.20}
\end{equation*}
$$
[[Ch08#^eq8-20|Eq. (8.20)]]

with the ATMF volatility, skew, and curvature all given in closed form directly in terms of the same three numbers $C^{x\xi},C^{\xi\xi},D$ and $Q=\widehat\sigma_T^2T$,

$$
\begin{align*}
\widehat{\sigma}_{F_{T} T}= & \hat{\sigma}_{T}\left[1+\frac{\varepsilon}{4 Q} C^{x \xi}\right.  \tag{8.21a}\\
& \left.+\frac{\varepsilon^{2}}{32 Q^{3}}\left(12\left(C^{x \xi}\right)^{2}-Q(Q+4) C^{\xi \xi}+4 Q(Q-4) D\right)\right] \\
\mathcal{S}_{T}= & \hat{\sigma}_{T}\left[\frac{\varepsilon}{2 Q^{2}} C^{x \xi}+\frac{\varepsilon^{2}}{8 Q^{3}}\left(4 Q D-3\left(C^{x \xi}\right)^{2}\right)\right]  \tag{8.21b}\\
\mathcal{C}_{T}= & \hat{\sigma}_{T} \frac{\varepsilon^{2}}{4 Q^{4}}\left(4 Q D+Q C^{\xi \xi}-6\left(C^{x \xi}\right)^{2}\right) \tag{8.21c}
\end{align*}
$$
[[Ch08#^eq8-21a|Eq. (8.21a)]] [[Ch08#^eq8-21b|Eq. (8.21b)]] [[Ch08#^eq8-21c|Eq. (8.21c)]]

Every one of the three coefficients in the quadratic expansion (8.20) is thus, in principle, computable from nothing more than the model's spot/variance and variance/variance covariance structure evaluated on today's curve.

At order one alone, the ATMF skew collapses further, to a strikingly compact formula built from $C^{x\xi}$ directly,

$$
\begin{equation*}
\mathcal{S}_{T}=\hat{\sigma}_{T} \frac{C^{x \xi}}{2\left(\hat{\sigma}_{T}^{2} T\right)^{2}} \tag{8.22}
\end{equation*}
$$
[[Ch08#^eq8-22|Eq. (8.22)]]

recovered independently by relating the ATMF skew to the skewness $s_T$ of $\ln S_T$ via the identity

$$
\begin{equation*}
\mathcal{S}_{T}=\frac{s_{T}}{6 \sqrt{T}} \tag{8.23}
\end{equation*}
$$
[[Ch08#^eq8-23|Eq. (8.23)]]

— an internal consistency check between the order-one skewness expansion and the cumulant-expansion machinery of Chapter 5. Rewriting (8.22) using the definition (8.13) gives a second, more physically transparent form,

$$
\begin{equation*}
\mathcal{S}_{T}=\frac{1}{2 \widehat{\sigma}_{T}^{3} T} \int_{0}^{T} \frac{T-\tau}{T} \frac{\left\langle d \ln S_{\tau} d \widehat{\sigma}_{T}^{2}(\tau)\right\rangle_{0}}{d \tau} d \tau \tag{8.24}
\end{equation*}
$$
[[Ch08#^eq8-24|Eq. (8.24)]]

the ATMF skew as a weighted average of the instantaneous covariance of the spot with the VS volatility for the residual maturity, evaluated over the option's life on the initial variance curve. This is precisely the identity Chapter 2 derived for the local volatility model (formula (2.89)) and Chapter 6 derived for Heston specifically — the general machinery derived here is shown to reduce correctly to both earlier special cases, closing a loop that had been left open across three prior chapters. Specializing (8.24) to a flat term structure gives a clean working formula,

$$
\begin{equation*}
\mathcal{S}_{T}=\frac{1}{\widehat{\sigma}_{T}^{2} T} \int_{0}^{T} \frac{T-\tau}{T} \frac{\left\langle d \ln S_{\tau} d \widehat{\sigma}_{T}(\tau)\right\rangle_{0}}{d \tau} d \tau \tag{8.25}
\end{equation*}
$$
[[Ch08#^eq8-25|Eq. (8.25)]]

and the corresponding order-one ATMF volatility is simply

$$
\begin{equation*}
\widehat{\sigma}_{F_{T} T}=\hat{\sigma}_{T}+\frac{Q}{2} \mathcal{S}_{T} \tag{8.26}
\end{equation*}
$$
[[Ch08#^eq8-26|Eq. (8.26)]]

matching the relationship between ATMF and VS volatilities already derived in Chapter 3's forward-start-option analysis, this time recovered as a byproduct of the general perturbative machinery rather than assumed.

Section 8.4 is the chapter's methodological centerpiece: it supplies an entirely independent derivation of the order-one result — not via PDE perturbation theory, but via an explicit hedging argument. Consider the auxiliary undiscounted price process $Q_t$ evaluated at the VS variance rather than a constant Black-Scholes volatility; because it solves the ordinary Black-Scholes PDE at every instant,

$$
\begin{equation*}
\frac{d Q}{d t}+(r-q) S \frac{d Q}{d S}+\frac{\omega}{2} S^{2} \frac{d^{2} Q}{d S^{2}}=0 \tag{8.27}
\end{equation*}
$$
[[Ch08#^eq8-27|Eq. (8.27)]]

with the usual Black-Scholes vega/gamma relationship applying at every instant as well,

$$
\begin{equation*}
\frac{d Q}{d \omega}=\frac{T-t}{2} S^{2} \frac{d^{2} Q}{d S^{2}} \tag{8.28}
\end{equation*}
$$
[[Ch08#^eq8-28|Eq. (8.28)]]

the spot-gamma contribution to $Q_t$'s daily P&L cancels automatically against the drift of $\omega_t$ — precisely because the VS used as a vega hedge doubles as a gamma hedge, the exact mechanism already seen at work in Chapter 7's variance-swap hedging discussion. Taking expectations and integrating over $[0,T]$ leaves only the spot/volatility and volatility/volatility cross-gamma contributions, giving a genuinely new representation of a European option's price in *any* diffusive model:

$$
\begin{align*}
P= & P_{B S}\left(0, S_{0}, \widehat{\sigma}_{T}^{2}(0)\right) \\
& +E\left[\int_{0}^{T} e^{-r t}\left(\frac{d^{2} P_{\mathrm{BS}}}{d S d\left(\widehat{\sigma}_{T}^{2}\right)}\left\langle d S_{t} d \widehat{\sigma}_{T}^{2}(t)\right\rangle+\frac{1}{2} \frac{d^{2} P_{\mathrm{BS}}}{\left(d\left(\widehat{\sigma}_{T}^{2}\right)\right)^{2}}\left\langle d \widehat{\sigma}_{T}^{2}(t) d \widehat{\sigma}_{T}^{2}(t)\right\rangle\right)\right] \tag{8.29}
\end{align*}
$$
[[Ch08#^eq8-29|Eq. (8.29)]]

the Black-Scholes price plus the expectation of the spot/volatility and volatility/volatility gamma P&Ls accrued by a position that is delta-hedged and VS-vega-hedged throughout its life. Specializing this to order one in volatility of volatility — dropping the pure volatility/volatility term, which only enters at order two — gives

$$
\begin{equation*}
P=P_{0}+\frac{1}{2} E\left[\int_{0}^{T}(T-t)\left\langle d \ln S_{t} d \widehat{\sigma}_{T}^{2}\right\rangle\right]\left(\frac{d^{3} P_{\mathrm{BS}}}{d x^{3}}-\frac{d^{2} P_{\mathrm{BS}}}{d x^{2}}\right) \tag{8.31}
\end{equation*}
$$
[[Ch08#^eq8-31|Eq. (8.31)]]

which recovers, term for term, the order-one piece of (8.18) already derived from perturbation theory — and reading off the ATMF skew from (8.31) reproduces (8.24) exactly,

$$
\begin{equation*}
\mathcal{S}_{T}=\frac{1}{2 \widehat{\sigma}_{T}^{3}(0) T} \int_{0}^{T} \frac{T-t}{T} \frac{\left\langle d \ln S_{t} d \widehat{\sigma}_{T}^{2}(t)\right\rangle_{0}}{d t} d t \tag{8.32}
\end{equation*}
$$
[[Ch08#^eq8-32|Eq. (8.32)]]

Two entirely independent routes — one a formal PDE expansion, the other an explicit hedging argument — land on the same prefactor. This is not a coincidence to be shrugged off; it is exactly the kind of cross-check the project's counter-intuitive-discovery process exists to flag, and it sets up the chapter's sharpest reversal.

Section 8.4.2 asks the obvious next question: the prefactor in (8.31) is, on its face, an abstract weighted spot/volatility covariance sitting inside a perturbation-theory formula — is it merely a mathematical artifact of the expansion, or can it be *materialized* as an actual, tradeable P&L? The answer is yes, and the payoff that does it is $\ln^2(S/S_0)$: delta-hedging and VS-vega-hedging a short position in this specific payoff, and summing the resulting spot/volatility cross-gamma P&Ls discounted back to $t=0$, produces exactly

$$
\begin{equation*}
P \& L=e^{-r T} \sum_{i}\left(T-t_{i}\right) \delta \ln S_{i} \delta\left(\widehat{\sigma}_{T}^{2}\left(t_{i}\right)\right) \tag{8.33}
\end{equation*}
$$
[[Ch08#^eq8-33|Eq. (8.33)]]

An abstract number buried in a perturbative expansion turns out to correspond to a concrete, hedgeable, market-observable exposure — the integrated spot/volatility covariance implied by option prices is not merely a modeling convenience but something a desk could, in principle, trade directly by running exactly this delta- and vega-hedged position in a log-squared payoff.

This is precisely the kind of reversal the project's counter-intuitive-discovery discipline exists to surface: a quantity that first appears purely as a coefficient inside a perturbation series turns out, on closer inspection, to be exactly what a specific, nameable trading strategy realizes as its P&L. Nothing about the derivation of (8.31) forced this outcome — it is a nontrivial fact about the payoff $\ln^2(S/S_0)$ specifically, not a generic property of every term that appears in a Taylor expansion of an option price.

Section 8.5 takes the opposite limit — vanishing maturity — where the whole variance curve collapses to a single instantaneous variance, and the general expressions (8.21) for ATM volatility, skew, and curvature turn out to be *exact*, not merely order-two approximations. Because the skewness of $\ln S_T$ vanishes only like $\sqrt T$ as $T\to0$,

$$
\begin{equation*}
s_{T}=\frac{3}{2} \frac{\mu_{0}}{\hat{\sigma}_{T}^{3}} \sqrt{T} \tag{8.34}
\end{equation*}
$$
[[Ch08#^eq8-34|Eq. (8.34)]]

the ATMF skew tends to a finite, nonzero constant rather than vanishing — a genuinely non-Gaussian short-maturity limit generated purely by volatility of volatility, since the drivers of $S_t$ remain ordinary Brownian motions throughout. Keeping only leading terms in $T$ gives the short ATM volatility, skew, and curvature entirely in terms of the instantaneous covariance functions $\mu_0,\nu_0$ and their sensitivity to the variance level,

$$
\begin{align*}
\widehat{\sigma}_{S, T=0} & =\hat{\sigma}_{0}  \tag{8.35a}\\
\mathcal{S}_{0} & =\hat{\sigma}_{0} \frac{1}{4\left(\xi_{0}^{0}\right)^{2}} \mu_{0}  \tag{8.35b}\\
\mathcal{C}_{0} & =\hat{\sigma}_{0} \frac{1}{4\left(\xi_{0}^{0}\right)^{4}}\left(\frac{2}{3} \xi_{0}^{0} \mu_{0} \frac{d \mu_{0}}{d \xi_{0}^{0}}+\frac{1}{3} \xi_{0}^{0} \nu_{0}-\frac{3}{2} \mu_{0}^{2}\right) \tag{8.35c}
\end{align*}
$$
[[Ch08#^eq8-35a|Eq. (8.35a)]] [[Ch08#^eq8-35b|Eq. (8.35b)]] [[Ch08#^eq8-35c|Eq. (8.35c)]]

Rewriting (8.35b) in terms of the observable short ATM volatility rather than the instantaneous variance produces the chapter's most quotable model-independent result — the short ATM skew is a *direct, model-independent* measure of the instantaneous covariance of $\ln S$ and the short ATM volatility,

$$
\begin{equation*}
\mathcal{S}_{0}=\frac{1}{2 \hat{\sigma}_{0}^{2}} \frac{\left\langle d \ln S d \hat{\sigma}_{0}\right\rangle}{d t} \tag{8.36}
\end{equation*}
$$
[[Ch08#^eq8-36|Eq. (8.36)]]

equivalent, in the SSR notation Chapter 9 introduces, to the property $\mathcal R_0=2$ — while the short curvature involves an additional, genuinely model-dependent piece, the covariance of $\ln S$ with the short ATMF skew itself,

$$
\begin{equation*}
\mathcal{C}_{0}=\frac{1}{4 \hat{\sigma}_{0}}\left(\frac{8}{3} \frac{\left\langle d \ln S d \mathcal{S}_{0}\right\rangle}{\hat{\sigma}_{0} d t}+\frac{4}{3} \frac{\left\langle d \hat{\sigma}_{0} d \hat{\sigma}_{0}\right\rangle}{\hat{\sigma}_{0}^{2} d t}-8 \mathcal{S}_{0}^{2}\right) \tag{8.37}
\end{equation*}
$$
[[Ch08#^eq8-37|Eq. (8.37)]]

so that, unlike the skew, the short curvature cannot be inverted to yield a model-independent measure of vol-of-vol without first assuming something about how the short skew itself moves with the spot.

Sections 8.5.1–8.5.2 apply (8.36)–(8.37) to two contrasting short-ATM-volatility dynamics. A lognormal short ATM volatility,

$$
\begin{align*}
d S & =\widehat{\sigma}_{0} S d W^{S}  \tag{8.38a}\\
d \widehat{\sigma}_{0} & =\bullet d t+\nu \widehat{\sigma}_{0} d W^{\widehat{\sigma}_{0}} \tag{8.38b}
\end{align*}
$$
[[Ch08#^eq8-38a|Eq. (8.38a)]] [[Ch08#^eq8-38b|Eq. (8.38b)]]

— the short-maturity limit of the two-factor lognormal forward-variance model of Section 7.4 — gives a short skew that is *constant*, entirely independent of the level of ATM volatility, with the curvature inversely proportional to it,

$$
\begin{align*}
\mathcal{S}_{0} & =\frac{\rho}{2} \nu  \tag{8.39a}\\
\mathcal{C}_{0} & =\frac{1}{6 \hat{\sigma}_{0}}\left(2-3 \rho^{2}\right) \nu^{2} \tag{8.39b}
\end{align*}
$$
[[Ch08#^eq8-39a|Eq. (8.39a)]] [[Ch08#^eq8-39b|Eq. (8.39b)]]

with skew, curvature and vol-of-vol tied together by

$$
\begin{equation*}
\nu^{2}=3 \hat{\sigma}_{0} \mathcal{C}_{0}+6 \mathcal{S}_{0}^{2} \tag{8.40}
\end{equation*}
$$
[[Ch08#^eq8-40|Eq. (8.40)]]

A normal short ATM volatility instead,

$$
\begin{align*}
d S & =\widehat{\sigma}_{0} S d W^{S}  \tag{8.41a}\\
d \widehat{\sigma}_{0} & =\bullet d t+\sigma d W^{\widehat{\sigma}_{0}} \tag{8.41b}
\end{align*}
$$
[[Ch08#^eq8-41a|Eq. (8.41a)]] [[Ch08#^eq8-41b|Eq. (8.41b)]]

— the short-maturity limit of Heston — gives a skew that is *inversely proportional* to the ATM volatility instead of constant,

$$
\begin{align*}
\mathcal{S}_{0} & =\frac{\rho}{2} \frac{\sigma}{\hat{\sigma}_{0}}  \tag{8.42a}\\
\mathcal{C}_{0} & =\frac{1}{6 \hat{\sigma}_{0}}\left(2-5 \rho^{2}\right)\left(\frac{\sigma}{\hat{\sigma}_{0}}\right)^{2} \tag{8.42b}
\end{align*}
$$
[[Ch08#^eq8-42a|Eq. (8.42a)]] [[Ch08#^eq8-42b|Eq. (8.42b)]]

The source explicitly confirms that (8.42a) reproduces Chapter 6's Heston-specific short-maturity skew formula (6.18b) verbatim — closing, for the first time within this chapter, the loop between the fully general machinery developed here and the Heston-specific derivation carried out three chapters earlier. The comparison also carries a warning: choosing $\nu=\sigma/\hat\sigma_0$ so that the two models share the same instantaneous vol-of-vol at $t=0$ leaves the skew identical in both, via

$$
\begin{equation*}
\left(\frac{\sigma}{\hat{\sigma}_{0}}\right)^{2}=3 \hat{\sigma}_{0} \mathcal{C}_{0}+10 \mathcal{S}_{0}^{2} \tag{8.43}
\end{equation*}
$$
[[Ch08#^eq8-43|Eq. (8.43)]]

but the *curvature* differs between (8.40) and (8.43) — confirming that curvature genuinely depends on how the covariance structure changes with the level of volatility, a piece of information the skew alone cannot supply, and that backing a vol-of-vol level out of a market smile is never truly model-independent. For the special case of vanishing spot/volatility correlation, the skew vanishes identically and the curvature becomes a clean, direct read of the vol-of-vol level itself,

$$
\begin{equation*}
\mathcal{C}_{0}=\frac{1}{3 \hat{\sigma}_{0}} \frac{\left\langle d \hat{\sigma}_{0} d \hat{\sigma}_{0}\right\rangle}{\hat{\sigma}_{0}^{2} d t} \tag{8.44}
\end{equation*}
$$
[[Ch08#^eq8-44|Eq. (8.44)]]

Section 8.6 closes the general/Heston loop a second time, more thoroughly, by embedding Heston in a general one-factor mean-reverting family,

$$
\left\{\begin{align*}
d S_{t} & =\sqrt{V_{t}} S_{t} d W_{t}^{S}  \tag{8.45}\\
d V_{t} & =-k\left(V_{t}-V^{0}\right) d t+\sigma V_{t}^{\varphi} d W_{t}^{V}
\end{align*}\right.
$$
[[Ch08#^eq8-45|Eq. (8.45)]]

with Heston corresponding to the specific exponent $\varphi=\tfrac12$. Taking conditional expectations pins down the forward variance curve this family generates,

$$
\begin{equation*}
\xi_{t}^{T}=E_{t}\left[V_{T}\right]=V^{0}+\left(V_{t}-V^{0}\right) e^{-k(T-t)} \tag{8.46}
\end{equation*}
$$
[[Ch08#^eq8-46|Eq. (8.46)]]

with driftless dynamics

$$
\begin{equation*}
d \xi_{t}^{T}=e^{-k(T-t)} \sigma \xi_{t}^{t \varphi} d W_{t}^{V} \tag{8.47}
\end{equation*}
$$
[[Ch08#^eq8-47|Eq. (8.47)]]

— and here is exactly the structural flaw Chapter 6 diagnosed, now visible in a single line: the same mean-reversion speed $k$ governs both the term structure of forward variances (8.46) and their volatilities (8.47), a first-generation model's hallmark limitation. The resulting spot/variance and variance/variance covariance functions,

$$
\begin{align*}
\mu(t, u, \xi) & =\rho \sigma e^{-k(u-t)} \xi_{t}^{t \varphi+\frac{1}{2}}  \tag{8.48a}\\
\nu\left(t, u, u^{\prime}, \xi\right) & =\sigma^{2} e^{-k(u-t)} e^{-k\left(u^{\prime}-t\right)} \xi_{t}^{t^{2 \varphi}} \tag{8.48b}
\end{align*}
$$
[[Ch08#^eq8-48a|Eq. (8.48a)]] [[Ch08#^eq8-48b|Eq. (8.48b)]]

feed directly into (8.13)-(8.17) and (8.21) to reproduce, for general $\varphi$, the short ATMF skew

$$
\begin{equation*}
\mathcal{S}_{0}=\frac{\rho \sigma}{4} V_{t=0}^{(\varphi-1)} \tag{8.49}
\end{equation*}
$$
[[Ch08#^eq8-49|Eq. (8.49)]]

which at $\varphi=\tfrac12$ recovers Chapter 6's Heston-specific result exactly, and whose order-one price correction recovers Chapter 6's equation (6.15) verbatim — the abstract prefactor there is now identified, with no ambiguity, as the doubly-integrated spot/variance covariance $C_0^{x\xi}$.

Section 8.7 puts the two-factor model of Chapter 7 through this same machinery, and this is where the chapter's practical payoff lands. The two-factor model's own spot/variance and variance/variance covariance functions,

$$
\begin{align*}
\mu(t, u, \xi)= & 2 \nu \xi^{u} \sqrt{\xi^{t}} \alpha_{\theta}\left[\rho_{S X^{1}}(1-\theta) e^{-k_{1}(u-t)}+\rho_{S X^{2}} \theta e^{-k_{2}(u-t)}\right]  \tag{8.50}\\
\nu\left(t, u, u^{\prime}, \xi\right)= & 4 \nu^{2} \xi^{u} \xi^{u^{\prime}} \alpha_{\theta}^{2}\left[(1-\theta)^{2} e^{-k_{1}\left(u+u^{\prime}-2 t\right)}+\theta^{2} e^{-k_{2}\left(u+u^{\prime}-2 t\right)}\right. \\
& \left.+\rho_{12} \theta(1-\theta)\left(e^{-k_{1}(u-t)} e^{-k_{2}\left(u^{\prime}-t\right)}+e^{-k_{2}(u-t)} e^{-k_{1}\left(u^{\prime}-t\right)}\right)\right] \tag{8.51}
\end{align*}
$$
[[Ch08#^eq8-50|Eq. (8.50)]] [[Ch08#^eq8-51|Eq. (8.51)]]

feed into $C^{x\xi}, C^{\xi\xi}, D$ by Gaussian quadrature, and the resulting order-two smile is checked against exact Monte Carlo prices: excellent near the money, degrading — as expected — for far-away strikes and long maturities, but the ATMF skew alone is already remarkably accurate at order *one* in $\varepsilon$,

$$
\begin{equation*}
\mathcal{S}_{T}^{\text {order 1 }}=\hat{\sigma}_{T} \frac{\varepsilon}{2 Q^{2}} C^{x \xi} \tag{8.52}
\end{equation*}
$$
[[Ch08#^eq8-52|Eq. (8.52)]]

Specializing $\mu$ to a flat forward-variance curve,

$$
\begin{equation*}
\mu(t, u, \xi)=2 \nu \xi_{0}^{\frac{3}{2}} \alpha_{\theta}\left[(1-\theta) \rho_{S X^{1}} e^{-k_{1}(u-t)}+\theta \rho_{S X^{2}} e^{-k_{2}(u-t)}\right] \tag{8.53}
\end{equation*}
$$
[[Ch08#^eq8-53|Eq. (8.53)]]

and inserting into (8.13) gives the order-one ATMF skew as an explicit double integral,

$$
\begin{equation*}
\mathcal{S}_{T}^{\text {order 1 }}=\frac{\nu \alpha_{\theta}}{\widehat{\sigma}_{T}^{3} T^{2}} \int_{0}^{T} d t \sqrt{\xi_{0}^{t}} \int_{t}^{T} d u \xi_{0}^{u}\left[(1-\theta) \rho_{S X^{1}} e^{-k_{1}(u-t)}+\theta \rho_{S X^{2}} e^{-k_{2}(u-t)}\right] \tag{8.54}
\end{equation*}
$$
[[Ch08#^eq8-54|Eq. (8.54)]]

which, for a flat term structure, integrates out to the same closed form already derived directly from the pricing PDE in Chapter 7,

$$
\begin{equation*}
\mathcal{S}_{T}^{\text {order 1 }}=\nu \alpha_{\theta}\left[(1-\theta) \rho_{S X^{1}} \frac{k_{1} T-\left(1-e^{-k_{1} T}\right)}{\left(k_{1} T\right)^{2}}+\theta \rho_{S X^{2}} \frac{k_{2} T-\left(1-e^{-k_{2} T}\right)}{\left(k_{2} T\right)^{2}}\right] \tag{8.55}
\end{equation*}
$$
[[Ch08#^eq8-55|Eq. (8.55)]]

— the same equation (7.120) Chapter 7 arrived at from its own discrete-model skew decomposition, now re-derived from the fully general smile-expansion machinery, a second independent confirmation that Chapter 7's construction is not an isolated trick but a genuine instance of this chapter's general theory. Because forward variances in the two-factor model are lognormal, (8.55) does not depend on the level of VS volatility at all — checked numerically against a full Monte Carlo simulation, ATMF skew levels are practically unchanged as VS volatility is dialed up or down, a striking contrast with Heston, where (8.42a) hard-wires the skew to scale inversely with the ATM volatility level whether a desk wants that coupling or not.

This flexibility, however, is not unlimited. Once $k_1,k_2,\rho_{12}$ are fixed to match a target vol-of-vol term structure, the remaining correlations $\rho_{SX^1},\rho_{SX^2}$ cannot be chosen freely — together with $\rho_{12}$ they must form a valid correlation matrix, which pins $\rho_{SX^2}$ to $\rho_{SX^1}$ up to one free parameter $\chi\in[-1,1]$:

$$
\begin{equation*}
\rho_{S X^{2}}=\rho_{12} \rho_{S X^{1}}+\chi \sqrt{1-\rho_{12}^{2}} \sqrt{1-\rho_{S X^{1}}^{2}} \tag{8.56}
\end{equation*}
$$
[[Ch08#^eq8-56|Eq. (8.56)]]

Within this constraint, though, the model still comfortably spans the range of power-law ATMF skew decay exponents observed in real equity index markets — a concrete, checkable resolution of exactly the flexibility gap Chapter 6 identified in Heston.

Section 8.10 closes the loop back to Chapter 7's VIX-smile machinery by asking how the *smile* of volatility of volatility — not just its overall level — feeds back into the vanilla smile. Since $D_0(\xi_0)$, unlike $C_0^{x\xi}$ and $C_0^{\xi\xi}$, depends on the derivative of $\mu$ with respect to the variance level, it alone is sensitive to whether vol-of-vol itself has a smile; using VIX-calibrated vol-of-vol-smile parameters from Section 7.7.1,

$$
\begin{equation*}
\gamma_{T}=50 \%, \beta_{T}=15 \%, \zeta_{T}=100 \% \tag{8.57}
\end{equation*}
$$
[[Ch08#^eq8-57|Eq. (8.57)]]

confirms exactly the prediction (8.21) makes: the upward-sloping vol-of-vol smile observed in real VIX markets weakens the ATMF skew and steepens the curvature, relative to an unsmiled lognormal benchmark calibrated to the same instantaneous vol-of-vol level — a genuinely new, checkable prediction that only the full order-two machinery, not the order-one skew formula alone, is capable of making.

This is worth restating plainly: at order one in volatility of volatility the vanilla smile cannot distinguish a model with a smiled vol-of-vol structure from one without, since $D_0$ only enters at order two. A trader relying solely on the order-one ATMF skew formula (8.55) would be blind to exactly this effect — a second concrete instance, alongside Section 8.4.2, of an order-two refinement carrying genuine, non-cosmetic financial content rather than merely improving numerical accuracy.

#### Chapter's Digest and Logical Bridge to Chapter 9

Chapter 8 reduced the near-the-money smile of any forward variance model to three numbers — $C_0^{x\xi}$, $C_0^{\xi\xi}$, $D_0$ — built entirely from the model's spot/variance and variance/variance covariance functions, and showed these numbers correctly reproduce both the local volatility model's skew identity (Chapter 2) and Heston's short-maturity skew (Chapter 6) as special cases, while resolving Chapter 6's flexibility complaint explicitly for the two-factor model. What remains unaddressed is the *dynamic* side of the same coin: Chapter 8 characterizes the smile a model produces today, holding the variance curve fixed, but says nothing yet about how that smile itself *evolves* as the spot and the curve move together — the joint static/dynamic consistency that determines, for instance, whether a model's smile dynamics matches what its own vanilla-hedging P&L requires. Chapter 9 takes up exactly this link between the static smile shape characterized here and the dynamic properties — the SSR briefly previewed in Section 8.5 among them — that a genuinely consistent stochastic volatility model must also satisfy.

---
### Chapter 9: Linking static and dynamic properties of stochastic volatility models

A stochastic volatility model can be assessed two ways: by its dynamic properties — the volatilities and correlations of volatilities, spot-starting or forward-starting, and their joint dynamics with the spot — or by its static properties, the shape of the smile it produces today. Both are consequences of the same underlying joint dynamics of spot and forward variances. But how tight is the link between them? Could a model be reparametrized so that today's vanilla smile stays essentially unchanged while its dynamics is materially different? Chapter 9's answer is a qualified no: it exhibits a genuine model-independent constraint tying the two together, and shows that violating it in practice creates a real, extractable trading P&L rather than a purely theoretical inconsistency.

The chapter's starting point is Chapter 8's order-one ATMF skew formula, specialized using the definition of $C^{x\xi}$,

$$
\begin{equation*}
\mathcal{S}_{T}=\frac{1}{2 \sqrt{T}} \frac{1}{\left(\int_{0}^{T} \xi_{0}^{\tau} d \tau\right)^{3 / 2}} \int_{0}^{T} d \tau \int_{\tau}^{T} \mu\left(\tau, u, \xi_{0}\right) d u \tag{9.1}
\end{equation*}
$$
[[Ch09#^eq9-1|Eq. (9.1)]]

rewritten, exactly as in Chapter 8's (8.24), as a weighted average of the instantaneous spot/VS-volatility covariance,

$$
\begin{equation*}
\mathcal{S}_{T}=\frac{1}{2 \widehat{\sigma}_{T}^{3} T} \int_{0}^{T} \frac{T-\tau}{T} \frac{\left\langle d \ln S_{\tau} d \widehat{\sigma}_{T}^{2}(\tau)\right\rangle_{0}}{d \tau} d \tau \tag{9.2}
\end{equation*}
$$
[[Ch09#^eq9-2|Eq. (9.2)]]

This is purely a *static* object — the smile's slope at the money, today. The dynamic counterpart is built by asking how the ATMF volatility itself co-moves with the spot: normalizing the regression coefficient of $\delta\widehat\sigma_{F_TT}$ on $\delta\ln S$ by the ATMF skew defines the chapter's central object, the Skew Stickiness Ratio,

$$
\begin{equation*}
\mathcal{R}_{T}=\frac{1}{\mathcal{S}_{T}} \frac{E\left[d \ln S d \widehat{\sigma}_{F_{T}(S) T}\right]}{E\left[(d \ln S)^{2}\right]} \tag{9.3}
\end{equation*}
$$
[[Ch09#^eq9-3|Eq. (9.3)]]

extending to any stochastic volatility model an SSR concept the book already used for the local volatility model. Two limiting values are known immediately: $\mathcal R_T$ vanishes identically in jump-diffusion or Lévy models with independent stationary increments, since implied volatility there is a pure function of moneyness and never moves when spot moves alone; and it reduces to the local-volatility-model SSR of Section 2.5.2 when implied volatilities are deterministic functions of $(t,S)$. Computing $\mathcal R_T$ at lowest order in volatility of volatility for a general stochastic volatility model — replacing the ATMF volatility with the VS volatility, which differ only at order one — gives the working formula that ties the static skew to the dynamic SSR through the same covariance function $\mu$ that built (9.1),

$$
\begin{equation*}
\mathcal{R}_{T}=\frac{\int_{0}^{T} \xi_{0}^{\tau} d \tau}{T \xi_{0}^{0}} \frac{T \int_{0}^{T} \mu\left(0, u, \xi_{0}\right) d u}{\int_{0}^{T} d \tau \int_{\tau}^{T} \mu\left(\tau, u, \xi_{0}\right) d u} \tag{9.4}
\end{equation*}
$$
[[Ch09#^eq9-4|Eq. (9.4)]]

Static and dynamic quantities share a common microscopic origin: both are functionals of the same $\mu$.

Taking the short-maturity limit of (9.1) and (9.4) recovers, respectively, Chapter 8's short ATMF skew result and a clean universal value for the short SSR: $\mathcal S_0=\mu(0,0,\xi_0)/4(\xi_0^0)^{3/2}$, and $\mathcal R_0=2$ — the same "2" that appears in the denominator of Chapter 8's short-skew formula (8.36), now shown to be a universal short-maturity limit for *any* diffusive stochastic volatility model, not merely a local-volatility-model peculiarity.

Section 9.4 delivers the chapter's central model-independent result. Assuming a flat term structure of VS volatilities and time-homogeneity — $\mu(\tau,u,\xi)$ depends only on $u-\tau$, written $\mu(u-\tau)$ — one of the integrals in (9.1) can be done analytically, giving simple closed forms for both quantities purely in terms of $\mu$,

$$
\begin{align*}
\mathcal{S}_{T} & =\frac{1}{2 \xi_{0}^{3 / 2} T^{2}} \int_{0}^{T}(T-t) \mu(t) d t  \tag{9.5}\\
\mathcal{R}_{T} & =\frac{\int_{0}^{T} \mu(t) d t}{\int_{0}^{T}\left(1-\frac{t}{T}\right) \mu(t) d t} \tag{9.6}
\end{align*}
$$
[[Ch09#^eq9-5|Eq. (9.5)]] [[Ch09#^eq9-6|Eq. (9.6)]]

Defining the cumulative function

$$
\begin{equation*}
g(\tau)=\int_{0}^{\tau} \mu(t) d t \tag{9.7}
\end{equation*}
$$
[[Ch09#^eq9-7|Eq. (9.7)]]

lets both $\mathcal S_T$ and $\mathcal R_T$ be rewritten purely in terms of $g$, and a purely convexity-based argument — assuming only that $\mu(t)$ decays monotonically to zero, $g$ is either increasing-concave or decreasing-convex, so $\mathcal R_T$, being the ratio of $g(T)$ to its own average over $[0,T]$, is bounded above by 2 and below by 1 — pins down a genuinely striking model-independent range,

$$
\begin{equation*}
\mathcal{R}_{T} \in[1,2] \tag{9.9}
\end{equation*}
$$
[[Ch09#^eq9-9|Eq. (9.9)]]

a constraint that holds regardless of which particular diffusive stochastic volatility model generated $\mu$. To go below 1, one presumably needs a jump or Lévy component in $\ln S_t$'s dynamics. The result comes with an immediate caveat, though: because the SSR's denominator is the instantaneous variance of $S$, the ratio can be made artificially small or large simply by shifting the short end of the VS-volatility term structure — a warning anyone reading an SSR off a real market smile needs to keep in mind.

Section 9.5 sharpens this into a genuine classification. Assuming $\mu(t)$ decays for large $t$ with a characteristic power-law exponent,

$$
\begin{equation*}
\mu(t) \propto \frac{1}{t^{\gamma}} \tag{9.10}
\end{equation*}
$$
[[Ch09#^eq9-10|Eq. (9.10)]]

models split cleanly into two classes at long maturity. **Type I** models ($\gamma>1$, fast decay) give

$$
\begin{equation*}
\mathcal{S}_{T} \propto \frac{1}{T} \text { and } \lim _{T \rightarrow \infty} \mathcal{R}_{T}=1 \tag{9.11}
\end{equation*}
$$
[[Ch09#^eq9-11|Eq. (9.11)]]

while **Type II** models ($\gamma<1$, slow decay) give

$$
\begin{equation*}
\mathcal{S}_{T} \propto \frac{1}{T^{\gamma}} \text { and } \lim _{T \rightarrow \infty} \mathcal{R}_{T}=2-\gamma \tag{9.12}
\end{equation*}
$$
[[Ch09#^eq9-12|Eq. (9.12)]]

Both are compactly unified by a single long-maturity relationship linking a purely static exponent to a purely dynamic limit,

$$
\begin{equation*}
\mathcal{S}_{T} \propto \frac{1}{T^{2-\mathcal{R}_{\infty}}} \tag{9.13}
\end{equation*}
$$
[[Ch09#^eq9-13|Eq. (9.13)]]

In Type I models the long-maturity behavior is universal, bearing no trace of the underlying dynamics beyond the fact that it decays fast; Type II models are richer, since both the ATMF-skew exponent and the SSR limit directly encode $\gamma$. Section 9.6 confirms Heston is a textbook Type I model by direct computation. Its exponentially decaying $\mu(t)$ (Chapter 8's (8.48)) forces the long-maturity ATMF skew, already derived in Chapter 6,

$$
\begin{equation*}
\mathcal{S}_{T}=\frac{\rho \sigma}{2 k T} \frac{1}{\sqrt{V^{0}}} \tag{9.14}
\end{equation*}
$$
[[Ch09#^eq9-14|Eq. (9.14)]]

to display Type I scaling directly, and the corresponding spot/ATMF-volatility covariance,

$$
\begin{align*}
\frac{E\left[d \widehat{\sigma}_{F_{T} T} d \ln S_{t}\right]}{E\left[\left(d \ln S_{t}\right)^{2}\right]} & =\frac{\partial \widehat{\sigma}_{F_{T} T}}{\partial V} \frac{E[d V d \ln S]}{V d t} \\
& =\frac{1}{2 k T \sqrt{V^{0}}} \frac{E[d V d \ln S]}{V d t} \\
& =\frac{1}{2 k T \sqrt{V^{0}}} \rho \sigma \tag{9.15}
\end{align*}
$$
[[Ch09#^eq9-15|Eq. (9.15)]]

divides out to give exactly $\mathcal R_T=1$ — confirming by hand what the classification predicted, and closing a loop back to Chapter 6's own Heston formulas.

Section 9.7 is the chapter's most instructive "can this even be built?" episode. Real equity index skews decay like $T^{-1/2}$ — squarely Type II, with $\gamma\approx1/2$ — so a genuine Type II model is called for. The obvious approach, giving forward variances a power-law-decaying vol-of-vol directly,

$$
\begin{equation*}
d \xi_{t}^{T}=\omega \xi_{t}^{T} \frac{1}{(T-t+\theta)^{\gamma}} d W \tag{9.17}
\end{equation*}
$$
[[Ch09#^eq9-17|Eq. (9.17)]]

does produce the desired power-law $\mu(\tau)\propto1/\tau^\gamma$ — but at the cost of sacrificing the finite-dimensional Markov representation that makes a model usable at all: forward variances can no longer be written as a function of finitely many state variables. This is a genuinely "fruitful failure," left as effectively closed rather than pursued further. The two-factor model of Chapters 7-8 instead solves the *practical* version of the problem: although its own $\mu(\tau)$, a linear combination of exponentials, is strictly Type I as $\tau\to\infty$, over the entire range of maturities relevant to real exotics risk management it exhibits genuine Type II scaling, an SSR plateau sitting around $2-\gamma$ rather than converging all the way to 1. Its ATMF skew and SSR read

$$
\begin{align*}
\mathcal{S}_{T} & =\frac{\omega}{2} \sum_{i} w_{i} \rho_{i S} \frac{k_{i} T-\left(1-e^{-k_{i} T}\right)}{\left(k_{i} T\right)^{2}}  \tag{9.16a}\\
\mathcal{R}_{T} & =\frac{\sum_{i} w_{i} \rho_{i S} \frac{1-e^{-k_{i} T}}{k_{i} T}}{\sum_{i} w_{i} \rho_{i S} \frac{k_{i} T-\left(1-e^{-k_{i} T}\right)}{\left(k_{i} T\right)^{2}}} \tag{9.16b}
\end{align*}
$$
[[Ch09#^eq9-16a|Eq. (9.16a)]] [[Ch09#^eq9-16b|Eq. (9.16b)]]

with the general (non-flat term structure) version requiring numerical integration,

$$
\begin{align*}
& \mathcal{S}_{T}=\frac{\nu \alpha_{\theta}}{\widehat{\sigma}_{T}^{3} T^{2}} \int_{0}^{T} d t \sqrt{\xi_{0}^{t}} \int_{t}^{T} \xi_{0}^{u}\left[(1-\theta) \rho_{S X^{1}} e^{-k_{1}(u-t)}+\theta \rho_{S X^{2}} e^{-k_{2}(u-t)}\right] d u \\
& \mathcal{R}_{T}=\frac{1}{\sqrt{\xi_{0}^{0}}} \frac{\widehat{\sigma}_{T}^{2} T \int_{0}^{T} \xi_{0}^{t}\left[(1-\theta) \rho_{S X^{1}} e^{-k_{1} t}+\theta \rho_{S X^{2}} e^{-k_{2} t}\right] d t}{\xi_{0}^{t}\left(\int_{t}^{T} \xi_{0}^{u}\left[(1-\theta) \rho_{S X^{1}} e^{-k_{1}(u-t)}+\theta \rho_{S X^{2}} e^{-k_{2}(u-t)}\right] d u\right)} \tag{9.18}
\end{align*}
$$
[[Ch09#^eq9-18|Eq. (9.18)]]

and the flat-term-structure closed form used for the actual numerical illustration,

$$
\begin{align*}
\mathcal{S}_{T} & =\nu \alpha_{\theta}\left[(1-\theta) \rho_{S X^{1}} \frac{k_{1} T-\left(1-e^{-k_{1} T}\right)}{\left(k_{1} T\right)^{2}}+\theta \rho_{S X^{2}} \frac{k_{2} T-\left(1-e^{-k_{2} T}\right)}{\left(k_{2} T\right)^{2}}\right] \\
\mathcal{R}_{T} & =\frac{(1-\theta) \rho_{S X^{1}} \frac{1-e^{-k_{1} T}}{k_{1} T}+\theta \rho_{S X^{2}} \frac{1-e^{-k_{2} T}}{k_{2} T}}{(1-\theta) \rho_{S X^{1}} \frac{k_{1} T-\left(1-e^{-k_{1} T}\right)}{\left(k_{1} T\right)^{2}}+\theta \rho_{S X^{2}} \frac{k_{2} T-\left(1-e^{-k_{2} T}\right)}{\left(k_{2} T\right)^{2}}} \tag{9.20}
\end{align*}
$$
[[Ch09#^eq9-20|Eq. (9.20)]]

Plotted against maturity, this SSR starts at 2 for short maturities exactly as the universal limit demands, tends to 1 only as $T\to\infty$, but sits on a stable shoulder around 1.5 for maturities from one to five years — precisely $2-\gamma$ for the $\gamma\approx1/2$ this parameter set was calibrated to reproduce in Chapter 8. A model that is "Type I on paper" behaves like Type II exactly where practical risk management lives; the asymptotic classification and the practically relevant behavior can diverge, and it is the latter that matters for pricing real exotics. A separate empirical check directly on Euro Stoxx 50 and S&P 500 historical data confirms real market dynamics is likewise consistent with Type II, not Type I, measured via a realized SSR estimator

$$
\begin{equation*}
\mathcal{R}_{T}^{r}=\frac{\sum_{i} \ln \frac{S_{i+1}}{S_{i}}\left(\widehat{\sigma}_{T, i+1}-\widehat{\sigma}_{T, i}\right)}{\sum_{i} \mathcal{S}_{T, i}\left(\ln \frac{S_{i+1}}{S_{i}}\right)^{2}} \tag{9.22}
\end{equation*}
$$
[[Ch09#^eq9-22|Eq. (9.22)]]

which hovers around 1.5 for both indexes at a two-year maturity — matching the two-factor model's own plateau almost exactly.

Sections 9.8–9.10 turn from analytic classification to the chapter's most consequential empirical result. First, a caveat on which variance correctly normalizes the SSR approximation: it is the *average* level of spot/volatility covariance generated by the model up to the option's maturity, not its particular value at $t=0$, that is practically relevant — using the short-dated instantaneous variance in (9.3)'s denominator, rather than the VS volatility for the actual maturity of interest, silently distorts the comparison. For *short* maturities specifically, the realized SSR is defined using the implied — not realized — variance in the denominator,

$$
\begin{equation*}
\mathcal{R}_{T}^{r, \text { short }}=\frac{\sum_{i} \ln \frac{S_{i+1}}{S_{i}}\left(\widehat{\sigma}_{T, i+1}-\widehat{\sigma}_{T, i}\right)}{\Delta t \sum_{i} \mathcal{S}_{T, i} \widehat{\sigma}_{T, i}^{2}} \tag{9.23}
\end{equation*}
$$
[[Ch09#^eq9-23|Eq. (9.23)]]

and measured on real index data, this realized short SSR sits systematically well below the model-independent implied value of 2 — posing the chapter's sharpest rhetorical question: does this point to a genuine discrepancy between the dynamics of stochastic models and the dynamics of actual market smiles? The resolution is a direct analogy to realized-versus-implied volatility: VS implied volatilities are systematically higher than realized volatilities for equity indexes, and this alone signals no modeling inconsistency, since the gap is exactly what a VS's carry P&L captures. The same logic applies here — no inconsistency exists as long as the gap $2-\mathcal R_T^{r,\text{short}}$ can itself be materialized as a genuine, tradeable option P&L.

Section 9.10 builds exactly that strategy. Risk-managing the position with the short-maturity lognormal model of Chapter 8's Section 8.5.1 — chosen because any model calibrated to the market smile gives the same implied spot/volatility covariance, and the lognormal choice specifically keeps the short ATM skew independent of the ATM volatility level, minimizing spurious mark-to-market P&L from remarking the skew — the relevant pricing PDE is

$$
\begin{equation*}
\frac{d P}{d t}+\frac{\widehat{\sigma}_{0}^{2}}{2} S^{2} \frac{d^{2} P}{d S^{2}}+\frac{\nu^{2}}{2} \widehat{\sigma}_{0}^{2} \frac{d^{2} P}{d \widehat{\sigma}_{0}^{2}}+\rho \nu \widehat{\sigma}_{0}^{2} S \frac{d^{2} P}{d S d \widehat{\sigma}_{0}}=0 \tag{9.24}
\end{equation*}
$$
[[Ch09#^eq9-24|Eq. (9.24)]]

with the vanilla price parametrized directly by ATM level, skew, and curvature,

$$
\begin{align*}
P & =P_{B S}(t, S, \widehat{\sigma}(x), K, T)  \tag{9.25a}\\
\widehat{\sigma}(x) & =\widehat{\sigma}_{0}+\mathcal{S} x+\frac{\mathcal{C}}{2} x^{2} \tag{9.25b}
\end{align*}
$$
[[Ch09#^eq9-25a|Eq. (9.25a)]] [[Ch09#^eq9-25b|Eq. (9.25b)]]

and $\rho,\nu$ recoverable directly from a market smile's observed skew and curvature,

$$
\begin{equation*}
\rho=\frac{2 \mathcal{S}}{\sqrt{3 \hat{\sigma}_{0} \mathcal{C}+6 \mathcal{S}^{2}}} \quad \nu=\sqrt{3 \hat{\sigma}_{0} \mathcal{C}+6 \mathcal{S}^{2}} \tag{9.26}
\end{equation*}
$$
[[Ch09#^eq9-26|Eq. (9.26)]]

A long, delta-hedged position in this option, calibrated daily to the near-the-money smile, has a carry P&L over $\delta t$ at order one in $\delta t$ and order two in $\delta S,\delta\widehat\sigma_0$ given by

$$
\begin{equation*}
P \& L=\frac{d \Pi}{d t} \delta t+\frac{d \Pi}{d \widehat{\sigma}_{0}} \delta \widehat{\sigma}_{0}+\frac{1}{2} \frac{d^{2} \Pi}{d S^{2}} \delta S^{2}+\frac{1}{2} \frac{d^{2} \Pi}{d \widehat{\sigma}_{0}^{2}}\left(\delta \widehat{\sigma}_{0}\right)^{2}+\frac{d^{2} \Pi}{d S d \widehat{\sigma}_{0}} \delta S \delta \widehat{\sigma}_{0} \tag{9.27}
\end{equation*}
$$
[[Ch09#^eq9-27|Eq. (9.27)]]

with the pure-vega term left over since the position is delta- but not vega-hedged. Using (9.24) to eliminate the theta, this splits cleanly into three break-even pieces — a pure vega term, a spot-gamma-plus-vol-gamma term at their respective break-even levels, and a cross-gamma/theta term at its own break-even level $\rho\nu\widehat\sigma_0$,

$$
\begin{align*}
P \& L= & \frac{d \Pi}{d \widehat{\sigma}_{0}} \delta \widehat{\sigma}_{0}  \tag{9.28a}\\
& +\frac{1}{2} S^{2} \frac{d^{2} \Pi}{d S^{2}}\left(\left(\frac{\delta S}{S}\right)^{2}-\widehat{\sigma}_{0}^{2} \delta t\right)+\frac{1}{2} \widehat{\sigma}_{0}^{2} \frac{d^{2} \Pi}{d \widehat{\sigma}_{0}^{2}}\left(\left(\frac{\delta \widehat{\sigma}_{0}}{\widehat{\sigma}_{0}}\right)^{2}-\nu^{2} \delta t\right)  \tag{9.28b}\\
& +S \widehat{\sigma}_{0} \frac{d^{2} \Pi}{d S d \widehat{\sigma}_{0}}\left(\frac{\delta S}{S} \frac{\delta \widehat{\sigma}_{0}}{\widehat{\sigma}_{0}}-\rho \nu \widehat{\sigma}_{0} \delta t\right) \tag{9.28c}
\end{align*}
$$
[[Ch09#^eq9-28a|Eq. (9.28a)]] [[Ch09#^eq9-28b|Eq. (9.28b)]] [[Ch09#^eq9-28c|Eq. (9.28c)]]

Since $\mathcal S=\rho\nu/2$ is exactly the statement $\mathcal R_{T=0}=2$, substituting this into the cross-gamma/theta piece (9.28c) rewrites it directly in terms of how far the realized cross-gamma/theta ratio departs from the universal value 2,

$$
\begin{align*}
S \widehat{\sigma}_{0} \frac{d^{2} \Pi}{d S d \widehat{\sigma}_{0}}\left(\frac{\delta S}{S} \frac{\delta \widehat{\sigma}_{0}}{\widehat{\sigma}_{0}}-\rho \nu \widehat{\sigma}_{0} \delta t\right) & =S \widehat{\sigma}_{0} \frac{d^{2} \Pi}{d S d \widehat{\sigma}_{0}}\left(\frac{\delta S}{S} \frac{\delta \widehat{\sigma}_{0}}{\widehat{\sigma}_{0}}-2 \mathcal{S} \widehat{\sigma}_{0} \delta t\right) \\
& =S \widehat{\sigma}_{0}^{2} \mathcal{S} \frac{d^{2} \Pi}{d S d \widehat{\sigma}_{0}}\left(\frac{\frac{\delta S}{S} \delta \widehat{\sigma}_{0}}{\mathcal{S} \widehat{\sigma}_{0}^{2} \delta t}-2\right) \delta t  \tag{9.29a}\\
& =S \widehat{\sigma}_{0}^{2} \mathcal{S} \frac{d^{2} \Pi}{d S d \widehat{\sigma}_{0}}\left(\mathcal{R}_{T=0}^{r, \text { short }}-2\right) \delta t \tag{9.29b}
\end{align*}
$$
[[Ch09#^eq9-29a|Eq. (9.29a)]] [[Ch09#^eq9-29b|Eq. (9.29b)]]

This is the central identity of the entire arbitraging strategy: the cross-gamma/theta P&L of a delta-hedged, gamma-neutral cross-strike position is directly proportional to $(\mathcal R_{T=0}^{r,\text{short}}-2)$ — exactly the deviation the empirical data showed is systematically negative. Equivalently, defining a realized instantaneous skew

$$
\begin{equation*}
\mathcal{S}^{r}=\frac{1}{2 \widehat{\sigma}_{0} \delta t} \frac{\delta S}{S} \frac{\delta \widehat{\sigma}_{0}}{\widehat{\sigma}_{0}}=\left(\frac{\mathcal{R}_{T=0}^{r, \text { short }}}{2}\right) \mathcal{S} \tag{9.30}
\end{equation*}
$$
[[Ch09#^eq9-30|Eq. (9.30)]]

the same P&L can be written directly as a gap between realized and implied skew,

$$
\begin{equation*}
2 S \widehat{\sigma}_{0}^{2} \frac{d^{2} \Pi}{d S d \widehat{\sigma}_{0}}\left(\mathcal{S}^{r}-\mathcal{S}\right) \delta t \tag{9.31}
\end{equation*}
$$
[[Ch09#^eq9-31|Eq. (9.31)]]

the final, directly backtestable trading result: sell the position whenever the market's implied skew is priced above what the spot/volatility co-movements actually realize.

To make this P&L attribution reliable in practice, the Black-Scholes theta of a strike-dependent smile must first be shown to split cleanly, with strike-independent break-even levels, into the three gamma/theta contributions of (9.28) — checked via

$$
\begin{equation*}
\frac{\widehat{\sigma}(x)^{2}}{2} S^{2} \frac{d^{2} P_{B S}}{d S^{2}}=\frac{\widehat{\sigma}_{0}^{2}}{2} S^{2} \frac{d^{2} P}{d S^{2}}+\frac{\nu^{2}}{2} \widehat{\sigma}_{0}^{2} \frac{d^{2} P}{d \widehat{\sigma}_{0}^{2}}+\rho \nu \widehat{\sigma}_{0}^{2} S \frac{d^{2} P}{d S d \widehat{\sigma}_{0}} \tag{9.32}
\end{equation*}
$$
[[Ch09#^eq9-32|Eq. (9.32)]]

which holds to good accuracy for strikes within about 5% of the money — precisely the 95%/105% range the actual backtest trades. Backtested on nine years of Euro Stoxx 50 one-month implied volatilities — selling 95%-strike options, buying the ratio of 105%-strike options that zeroes out the spot gamma, delta-hedging daily and unwinding the next day — the realized ATM skew averages about 20% below the implied skew, translating directly, via (9.31), into a materialized, positive P&L consistent with exactly the gap the realized-SSR data predicted.

Section 9.11 closes by folding the SSR back into practical model choice. Autocalls supply a concrete real-world stake in this: their upper knockout threshold generates a large spot/volatility cross-gamma, and unwinding the resulting vega hedge as spot rallies produces a cross-gamma P&L that is offset by a theta term quantified directly by the pricing model's own SSR,

$$
\begin{align*}
P \& L & =-S \widehat{\sigma}_{T} \frac{d^{2} \Pi}{d S d \widehat{\sigma}_{T}}\left(\frac{\delta S}{S} \frac{\delta \widehat{\sigma}_{T}}{\widehat{\sigma}_{T}}-\mathcal{R}_{T} \mathcal{S}_{T} \frac{\sigma^{2}}{\widehat{\sigma}_{T}} \delta t\right)  \tag{9.33}\\
& =-S \mathcal{S}_{T} \frac{d^{2} \Pi}{d S d \widehat{\sigma}_{T}}\left(\mathcal{R}_{T}^{r} \sigma_{r}^{2}-\mathcal{R}_{T} \sigma^{2}\right) \delta t
\end{align*}
$$
[[Ch09#^eq9-33|Eq. (9.33)]]

Stochastic volatility and local volatility models produce structurally *opposite* SSR term structures — both start at the universal value 2, but a stochastic volatility model's SSR then decreases toward 1 or a Type-II plateau, while the local volatility model's SSR instead *increases* for longer maturities, typically toward a value near 3 for a $\gamma=1/2$ skew decay. This is not a minor calibration detail: whichever direction a book's net cross-gamma exposure points, one of the two model families is the conservative choice, and the other is not. A final practical caveat closes the loop: because vega-hedging in practice happens far less often than delta-hedging, isolating exactly which piece of the theoretical cross-gamma/theta P&L is actually realized requires re-deriving the spot/volatility cross-gamma P&L under periodic, rather than continuous, vega rehedging,

$$
\begin{equation*}
\frac{d^{2} P}{d S d \widehat{\sigma}} \Sigma_{i}\left(\Sigma_{j \geq i} \delta \widehat{\sigma}_{j}\right) \delta S_{i} \tag{9.34}
\end{equation*}
$$
[[Ch09#^eq9-34|Eq. (9.34)]]

which shows that only about half of the theoretical cross-gamma/theta P&L computed under continuous vega-hedging is actually captured once realistic, weekly-or-so rehedging frequencies are used — a sobering coda to the arbitrage strategy just constructed, and a reminder that the model-independent constraint this chapter derives is a bound on what is *theoretically* extractable, not a guarantee of what a real trading desk will actually collect.

#### Chapter's Digest and Logical Bridge to Chapter 10

Chapter 9 established that a model's static smile and its dynamic SSR are not independently choosable: they are tied together by the same spot/variance covariance function $\mu$, bounded to $\mathcal R_T\in[1,2]$ under mild assumptions, and classified into Type I/Type II behavior by a single decay exponent $\gamma$ — with the gap between a model's implied SSR and the market's realized SSR directly convertible into a tradeable P&L. This closes the book's long-running thread, begun in Chapter 6 and carried through Chapters 7-8, of asking what forward variance models can and cannot decouple: vol-of-vol term structure, VIX smiles, ATMF skew term structure, and now the static/dynamic SSR link are all shown to be governed by the same handful of covariance functionals, with genuine but bounded flexibility in each. What none of the last four chapters has addressed, however, is the more basic empirical question of *why* equity indexes exhibit the pronounced negative skew and fat left tail that motivate stochastic volatility modeling in the first place. Chapter 10 turns to exactly that question — surveying the actual economic and structural causes of equity smiles, as a bridge from the model-building machinery of Chapters 6-9 to the practical modeling choices of the book's remaining chapters.

---
### Chapter 10: What causes equity smiles?

Chapters 6-9 built and characterized an entire apparatus of stochastic volatility models — their vol-of-vol term structure, their smile expansion, their static/dynamic SSR link — all under the working assumption that genuine covariance between spot and implied volatility is what generates the equity smile. Chapter 10 steps back and asks the foundational question the rest of the book has taken for granted: is that actually true? Large historical drawdowns in equity indexes are often blamed for the strong negative skew of implied smiles. Are vanilla smiles in any way related to statistical properties of historical returns? Which payoffs are actually impacted by fat-tailed returns, and which are not?

Section 10.1 starts from an empirical fact. Fitting a century of Dow Jones daily closes to a Student distribution,

$$
\begin{equation*}
\rho_{\mu}(x)=\frac{\Gamma\left(\frac{1+\mu}{2}\right)}{\sqrt{\mu \pi} \Gamma\left(\frac{\mu}{2}\right)} \frac{1}{\left(1+\frac{x^{2}}{\mu}\right)^{\frac{1+\mu}{2}}} \tag{10.1}
\end{equation*}
$$
[[Ch10#^eq10-1|Eq. (10.1)]]

gives a good fit for $\mu$ anywhere in $[3,4]$ — with tails so fat that the lognormal distribution underestimates the probability of a 4-standard-deviation move by a factor of 100. Just as strikingly, and contrary to the common belief that crashes make negative returns fatter-tailed than positive ones, the fit shows negative and positive returns have very similar tail thickness — a pattern that holds not just for the Dow Jones but broadly across equity indexes.

Section 10.1.1 asks whether this century-long unconditional fit is even the right question, since volatility regimes clearly varied enormously over a hundred years. The natural ansatz separates the scale of a daily return from its residual shape,

$$
\begin{equation*}
r_{i}=\sigma_{i} \sqrt{\delta t} z_{i} \tag{10.2}
\end{equation*}
$$
[[Ch10#^eq10-2|Eq. (10.2)]]

with $\sigma_i$ the (stochastic) instantaneous volatility and $z_i$ identically distributed. Even when $z_i$ is itself Gaussian, a stochastic $\sigma_i$ alone generates non-Gaussian, fat-tailed unconditional returns — so the real question is how much of the empirical tail thickness in Section 10.1 is just this scale effect, and how much is genuine residual non-Gaussianity in $z_i$ itself. Conditioning each return on a trailing realized-volatility estimate only modestly thins the tails — $\mu$ rises from 3.2 to 3.8 for negative returns — which means a substantial, *residual* non-Gaussianity survives beyond anything stochastic volatility alone explains. This residual, purely distributional effect is exactly what the rest of the chapter sets out to isolate and price.

Section 10.2.1 builds a deliberately surgical tool for isolating it: take the two-factor forward-variance model of Chapter 7 exactly as is, and replace only its daily Gaussian spot shock with a two-sided Student-distributed one,

$$
\begin{equation*}
S_{t+\Delta}=S_{t}\left[1+(r-q) \Delta+\sigma_{t} \delta Z\right] \tag{10.3}
\end{equation*}
$$
[[Ch10#^eq10-3|Eq. (10.3)]]

with separate degrees-of-freedom parameters $\mu_+,\mu_-$ and probabilities $p_+,p_-$ for the positive and negative tails, chosen so that the shock's first two moments match those of a standard Gaussian increment exactly,

$$
\begin{equation*}
E[\delta Z]=0 \quad E\left[(\delta Z)^{2}\right]=\Delta \tag{10.4}
\end{equation*}
$$
[[Ch10#^eq10-4|Eq. (10.4)]]

To keep this fat-tailed shock correlated with the model's own forward-variance factors, and to let it degenerate smoothly back to the ordinary Gaussian two-factor model as a check, it is built as an explicit deterministic mapping of the Monte Carlo engine's own Brownian increment,

$$
\begin{equation*}
\delta Z=\sqrt{\Delta} f\left(\frac{\delta W^{S}}{\sqrt{\Delta}}\right) \tag{10.5}
\end{equation*}
$$
[[Ch10#^eq10-5|Eq. (10.5)]]

with $f$ defined piecewise via the two Student inverse-CDFs on either side of the split point,

$$
\begin{cases}x \leq \mathcal{N}_{G}^{-1}\left(p_{-}\right) & f(x)=\zeta_{-} \sqrt{\frac{\mu_{-}-2}{\mu_{-}}} \mathcal{N}_{\mu_{-}}^{-1}\left(\frac{\mathcal{N}_{G}(x)}{2 p_{-}}\right)  \tag{10.6}\\ x \geq \mathcal{N}_{G}^{-1}\left(p_{-}\right) & f(x)=\zeta_{+} \sqrt{\frac{\mu_{+}-2}{\mu_{+}}} \mathcal{N}_{\mu+}^{-1}\left(\frac{1}{2}+\frac{\mathcal{N}_{G}(x)-p_{-}}{2 p_{+}}\right)\end{cases}
$$
[[Ch10#^eq10-6|Eq. (10.6)]]

so that setting $p_+=p_-=\tfrac12$ and $\mu_+,\mu_-\to\infty$ collapses $\delta Z$ back to the plain Gaussian $\delta W^S$, exactly recovering the standard two-factor model as a sanity check.

Simply swapping in this fat-tailed shock, however, silently changes the model's spot/volatility covariance structure along with its return distribution — a confound that would make it impossible to attribute any change in the resulting smile to the tail thickness specifically, rather than to an accidental change in the genuine stochastic-volatility channel already characterized in Chapters 7-9. The fix is to rescale the spot/volatility correlations uniformly,

$$
\begin{equation*}
\frac{\rho_{i S}^{*}}{\rho_{i S}}=\frac{\Delta}{E\left[\delta Z \delta W^{S}\right]}=\frac{1}{\int_{-\infty}^{+\infty} \phi(x) x f(x) d x} \tag{10.7}
\end{equation*}
$$
[[Ch10#^eq10-7|Eq. (10.7)]]

so that the instantaneous covariance of $\delta Z$ with each forward-variance factor exactly matches what the standard Gaussian model would have generated. With this rescaling in place, any difference between the fat-tailed model's smile and the standard model's smile can be attributed purely to the shape of the daily return distribution — a clean separation of channels that is the whole point of the construction.

Section 10.2.3 explains why this distinction between *scale* and *distributional shape* is even meaningful for pricing, rather than a purely statistical nicety. The P&L of a short delta-hedged option over one rehedging interval $\Delta$ is

$$
\begin{equation*}
P \& L=-(P(t+\Delta, S+\delta S)-P(t, S))+\frac{d P}{d S} \delta S \tag{10.8}
\end{equation*}
$$
[[Ch10#^eq10-8|Eq. (10.8)]]

and writing the pricing function as a conditional expectation over the next day's fat-tailed shock,

$$
\begin{equation*}
P(t, S)=E[P(t+\Delta, S(1+\sigma \delta Z))] \tag{10.9}
\end{equation*}
$$
[[Ch10#^eq10-9|Eq. (10.9)]]

and expanding in powers of $\delta Z$ gives the full P&L decomposition, order by order,

$$
\begin{align*}
P \& L & =-\left(P(t+\Delta, S+\delta S)-E_{t, S}[P(t+\Delta, S(1+\sigma \delta Z))]\right)+\frac{d P}{d S} \delta S \\
& =-\frac{S^{2}}{2} \frac{d^{2} P}{d S^{2}}\left(\frac{\delta S^{2}}{S^{2}}-\sigma^{2} \Delta\right)-\sum_{k>2} \frac{S^{k}}{k!} \frac{d^{k} P}{d S^{k}}\left(\left(\frac{\delta S}{S}\right)^{k}-\sigma^{k} E\left[\delta Z^{k}\right]\right) \tag{10.10}
\end{align*}
$$
[[Ch10#^eq10-10|Eq. (10.10)]]

This is the chapter's technical crux: in the continuum limit $\Delta\to0$ that every previous chapter's diffusive models implicitly took, the order-$k$ term for $k>2$ scales like $\Delta^{k/2}$, and summed over the $T/\Delta$ rehedges in an option's life the cumulative contribution scales like $\Delta^{k/2-1}$ — vanishing for any $k>2$. This is precisely *why* none of Chapters 6-9's continuous-time machinery ever needed to worry about anything beyond the second moment of the daily shock. But at the actual, finite daily time scale relevant to real trading and to genuinely path-dependent daily-return payoffs, these higher-order terms do not vanish — they survive, and they matter, exactly where the continuum-limit intuition breaks down.

Running the fat-tailed model with this apparatus produces the chapter's central verdict, and it is a genuine debunking. The near-ATMF vanilla smile for standard maturities is barely affected by how fat the daily tails are — only far out-of-the-money strikes move at all, and mostly on the high-strike side, not the low-strike side one might naively expect crash risk to distort. The ATMF skew, exactly as Chapters 8-9 already argued, remains overwhelmingly generated by the genuine covariance of spot with the implied volatility of the residual maturity — not by historical drawdowns, and not by fat daily tails. Payoffs that directly accumulate daily returns are the systematic exception: variance swaps show only a modest gap between VS-implied and log-contract-implied volatility, consistent with the small effect already measured empirically in Chapter 5, while daily cliquets are strongly sensitive to tail thickness — matching their market prices requires an implied $\mu_-\simeq2.2$, well below even the fattest end of the historical $[3,4]$ range. This gap is not read as evidence that markets secretly believe in far fatter tails than history shows; it is read as pricing compensation for a distinct, genuinely unhedgeable risk — daily-gap risk — that a delta-hedging strategy simply cannot neutralize between rehedges.

Appendix A turns from this empirical question to jump-diffusion and Lévy models, and asks a pointed methodological question: when are these models actually needed, and what do they supply that stochastic volatility alone does not? The chapter's answer is that they are best understood not as physical descriptions of returns at all, but as an explicit stress-test reserve mechanism. Starting from the ordinary Black-Scholes gamma/theta P&L,

$$
\begin{equation*}
P \& L=-\frac{1}{2} S^{2} \frac{d^{2} P}{d S^{2}}\left(\frac{\delta S^{2}}{S^{2}}-\sigma^{2} \delta t\right) \tag{10.11}
\end{equation*}
$$
[[Ch10#^eq10-11|Eq. (10.11)]]

imagine a shock of fixed relative magnitude $J$ applied to the underlying; the P&L this generates on a delta-hedged position is

$$
\begin{equation*}
P \& L=-(P(t, S(1+J))-P(t, S))+J S \frac{d P}{d S} \tag{10.12}
\end{equation*}
$$
[[Ch10#^eq10-12|Eq. (10.12)]]

Assigning this shock an annualized intensity $\lambda$, the theta needed to offset its expected impact on average is

$$
\begin{equation*}
\lambda \delta t\left[(P(t, S(1+J))-P(t, S))-J S \frac{d P}{d S}\right] \tag{10.13}
\end{equation*}
$$
[[Ch10#^eq10-13|Eq. (10.13)]]

which, integrated to today's price adjustment for an option of maturity $T$, gives the reserve

$$
\begin{equation*}
\Delta P(t, S)=\lambda(T-t)\left[(P(t, S(1+J))-P(t, S))-J S \frac{d P}{d S}\right] \tag{10.14}
\end{equation*}
$$
[[Ch10#^eq10-14|Eq. (10.14)]]

directly parallel to the variance-swap reserve adjustment already derived in Chapter 5. This ad-hoc reserve, though, has a genuine limitation: it includes no provision for the financing cost of its own delta hedge, and no theta to offset its own gamma — it is a patch bolted onto a Black-Scholes price, not a fully self-consistent pricing model.

Appendix A.2 fixes this by building the reserve directly into the pricing PDE from the start, demanding the P&L read

$$
\begin{equation*}
P \& L=-\frac{1}{2} S^{2} \frac{d^{2} P}{d S^{2}}\left(\frac{\delta S^{2}}{S^{2}}-\sigma^{2} \delta t\right)+\lambda\left[(P(t, S(1+J))-P(t, S))-J S \frac{d P}{d S}\right] \delta t \tag{10.15}
\end{equation*}
$$
[[Ch10#^eq10-15|Eq. (10.15)]]

which yields the jump-diffusion pricing PDE directly,

$$
\begin{equation*}
\frac{d P}{d t}+(r-q) S \frac{d P}{d S}+\frac{\sigma^{2}}{2} S^{2} \frac{d^{2} P}{d S^{2}}+\lambda\left[(P(t, S(1+J))-P(t, S))-J S \frac{d P}{d S}\right]=r P \tag{10.16}
\end{equation*}
$$
[[Ch10#^eq10-16|Eq. (10.16)]]

whose probabilistic dual is an ordinary diffusion with Poisson jumps added,

$$
\begin{equation*}
d S_{t}=(r-q) S_{t} d t+\sigma S_{t} d W_{t}+J S_{t}\left(d N_{t}-\lambda d t\right) \tag{10.17}
\end{equation*}
$$
[[Ch10#^eq10-17|Eq. (10.17)]]

Making the jump amplitude itself random, with density $\rho(J)$, generalizes the PDE to its full non-local form,

$$
\begin{align*}
\frac{d P}{d t} & +(r-q) S \frac{d P}{d S}+\frac{\sigma^{2}}{2} S^{2} \frac{d^{2} P}{d S^{2}} \\
& +\lambda \int_{-1}^{\infty} \rho(J)\left[(P(t, S(1+J))-P(t, S))-J S \frac{d P}{d S}\right] d J=r P \tag{10.18}
\end{align*}
$$
[[Ch10#^eq10-18|Eq. (10.18)]]

For constant volatility and vanilla payoffs, this is homogeneous in $\ln S$ and solvable by Laplace transform: substituting $P=Se^{-q\tau}f(\tau,x)$ and transforming reduces the PDE to a simple first-order ODE,

$$
\begin{equation*}
-\frac{d F}{d \tau}+\left[\frac{\sigma^{2}}{2} p(1+p)+\lambda(\psi(p)-(1+p) \psi(0))\right] F=0 \tag{10.19}
\end{equation*}
$$
[[Ch10#^eq10-19|Eq. (10.19)]]

whose explicit solution is

$$
\begin{align*}
F(\tau, p) & =\frac{1}{p(1+p)} e^{\tau H(p)}  \tag{10.20}\\
H(p) & =\frac{\sigma^{2}}{2} p(1+p)+\lambda(\psi(p)-(1+p) \psi(0))
\end{align*}
$$
[[Ch10#^eq10-20|Eq. (10.20)]]

with the ordinary Black-Scholes vanilla payoff transform recovered exactly when $\lambda=0$. Pricing directly off the solution of (10.18), rather than off a Black-Scholes price plus an ad-hoc reserve bolt-on, ensures the additional theta really is (10.13), that the resulting delta already incorporates the reserve's own hedge, and that the jump's P&L impact on the reserve itself is properly accounted for — exactly the three gaps the ad-hoc policy left open.

Appendix A.3 backs the implied ATMF skew out of this solution. Using the cumulant-generating function of the terminal log-return density,

$$
\begin{align*}
e^{L(\tau, q)} & =\int_{-\infty}^{+\infty} e^{-q z} \rho_{z}(z) d z  \tag{10.21}\\
& =\int_{-\infty}^{+\infty} e^{-q z}\left(e^{x}\left(\partial_{x}+\partial_{x}^{2}\right) f(\tau, x)\right)_{x=-z} d z \\
& =q(1+q) F(\tau,-(1+q))
\end{align*}
$$
[[Ch10#^eq10-21|Eq. (10.21)]]

which at $t=0$ reads

$$
\begin{align*}
L(T, q) & =T\left(\frac{\sigma^{2}}{2} q(1+q)+\lambda(\psi(-(1+q))+q \psi(0))\right) \\
& =T(\phi(q)+q \phi(-1)) \tag{10.22}
\end{align*}
$$
[[Ch10#^eq10-22|Eq. (10.22)]]

with

$$
\begin{align*}
\phi(q) & =\frac{\sigma^{2}}{2} q^{2}+\psi(-(1+q)) \\
& =\frac{\sigma^{2}}{2} q^{2}+\int_{-\infty}^{+\infty} \lambda \rho^{*}(u)\left(e^{-q u}+q u-1\right) d u \tag{10.23}
\end{align*}
$$
[[Ch10#^eq10-23|Eq. (10.23)]]

or, directly in terms of the jump size $J$ rather than the log-jump $u$,

$$
\begin{equation*}
L(T, q)=T\left(\frac{\sigma^{2}}{2} q(1+q)+\int_{-1}^{\infty} \lambda \rho(J)\left((1+J)^{-q}-1+q J\right) d J\right) \tag{10.24}
\end{equation*}
$$
[[Ch10#^eq10-24|Eq. (10.24)]]

Expanding for small jump sizes and stopping at the lowest nontrivial order gives

$$
\begin{equation*}
L(T, q)=T\left(\frac{\sigma^{2}+\lambda \overline{J^{2}}}{2} q(1+q)-\frac{\lambda \overline{J^{3}}}{6} q(1+q)(2+q)\right) \tag{10.25}
\end{equation*}
$$
[[Ch10#^eq10-25|Eq. (10.25)]]

The order-two-in-$J$ piece merely rescales the volatility level, $\sigma^2\to\sigma^2+\lambda\overline{J^2}$, and generates no skew at all — the genuinely new content only arrives at order three, exactly the same cubic-in-jump-size structure the book already used, in reverse, to critique jump-model-mediated skewness inference back in Chapter 5. Converting this cumulant expansion into an ATMF skew via the same skew/skewness relationship used throughout Chapters 8-9 gives the universal small-jump result,

$$
\begin{equation*}
\mathcal{S}_{T}=\frac{\lambda \overline{J^{3}}}{6\left(\sigma^{2}+\lambda \overline{J^{2}}\right)^{3 / 2} T}=\frac{\lambda \overline{J^{3}}}{6 \widehat{\sigma}_{T}^{3} T} \tag{10.26}
\end{equation*}
$$
[[Ch10#^eq10-26|Eq. (10.26)]]

so that, for small jump sizes, jump-diffusion models generate an ATMF skew that decays like $1/T$ — the identical Type I scaling Chapter 9 derived for Heston, since a jump-diffusion process with independent stationary increments is, in exactly the sense Chapter 9's classification cares about, a fast-decaying-covariance model.

Appendix A.4 asks a sharper question: once jumps are layered onto an *existing*, already-calibrated market model — say local volatility — does the resulting P&L still take the clean gamma/theta form, with the jump contribution simply added on as an extra theta? The answer is more subtle than "yes." Adding the jump term to the local volatility model's own pricing equation forces the local volatility function itself to absorb part of the jump's contribution, since it must still be chosen to match market vanilla prices exactly,

$$
\begin{equation*}
\sigma^{2}(T, K)=\sigma_{\text {loc }}^{2}(K, T)-2 \lambda \frac{\left[(1+J) C_{K /(1+J), T}-C_{K, T}\right]-J C+\lambda J K \frac{d C}{d K}}{K^{2} \frac{d^{2} C}{d K^{2}}} \tag{10.27}
\end{equation*}
$$
[[Ch10#^eq10-27|Eq. (10.27)]]

with the jump-free Dupire local volatility $\sigma_{\text{loc}}$ shifted downward, at order two in $J$, by exactly $\lambda J^2$ — the local volatility function is not an independent free parameter once jumps are present, but is itself contaminated by the jump specification. Re-deriving the delta-and-vega-hedged carry P&L for this jump-adjusted local volatility model shows that the extra theta contribution is not simply proportional to a jump in spot alone,

$$
\begin{align*}
+\lambda & {\left[\left(P\left(t, S(1+J), \widehat{\sigma}_{K T}+\Delta \Sigma_{K T}^{\mathrm{LV}}(t, S)\right)-P\left(t, S, \widehat{\sigma}_{K T}\right)\right)\right.} \\
& \left.-J S\left(\frac{d P}{d S}+\frac{d P}{d \widehat{\sigma}_{K T}} \bullet \frac{d \Sigma_{K T}^{\mathrm{LV}}}{d S}\right)\right] \delta t \tag{10.28}
\end{align*}
$$
[[Ch10#^eq10-28|Eq. (10.28)]]

but is instead proportional to the P&L generated by *joint* jumps in spot *and* implied volatility — a jump on $S$ generically forces a simultaneous, calibration-consistent jump in implied volatilities too, not merely in the spot alone. Depending on which market model one starts from, adding jumps is thus never as simple as bolting on an extra theta term; the jump scenario for implied volatility is dictated jointly by the jump scenario for spot and by whatever smile the base model was calibrated to.

Appendix A.5 generalizes from a single Poisson process to a superposition of arbitrarily many, each with its own intensity and jump density. If the total intensity remains finite, this collapses back to an ordinary, effective Poisson process — nothing new. But letting the total intensity diverge, so that jumps occur infinitely frequently, produces the genuinely interesting Lévy-Khinchine case,

$$
\begin{equation*}
\phi(q)=\frac{\sigma^{2}}{2} q^{2}+\int_{-\infty}^{+\infty}\left(e^{-q u}+q u-1\right) k(u) d u \tag{10.29}
\end{equation*}
$$
[[Ch10#^eq10-29|Eq. (10.29)]]

with $k(u)$ the Lévy measure governing infinitesimally frequent jumps of every size simultaneously. The chapter's closing methodological verdict is unambiguous: jump and Lévy processes are useful, deliberately crude reserve-policy and stress-test tools — a clean way to charge a theta for a specified adverse scenario, and nothing more — not a genuine, "bona fide characterization of the actual or implied dynamics of real securities." Section 10.2.1's fat-tailed stochastic volatility construction, built specifically to isolate and price the residual, non-Gaussian shape of daily returns while leaving the genuine covariance structure of spot and implied volatility untouched, is the more honest tool whenever a realistic model of tail risk — rather than a stress-test charge — is actually what is needed.

#### Chapter's Digest and Logical Bridge to Chapter 11

Chapter 10 answered the question the rest of the book had quietly assumed: equity smiles are overwhelmingly generated by genuine spot/implied-volatility covariance, not by historical drawdowns or fat-tailed daily returns — with the sharp exception of directly path-dependent, daily-return-accumulating payoffs like cliquets, where fat tails matter a great deal. Jump and Lévy models, meanwhile, are best understood as a stress-test reserve mechanism, not a physically realistic description of dynamics. Having now settled what actually causes the smile of a single underlying, the book's remaining chapters turn outward. Chapter 11 extends the entire forward-variance and covariance-function apparatus built across Chapters 6-9 from a single equity underlying to the genuinely multi-asset case — index/single-stock dispersion, basket smiles, and the joint dynamics of correlated implied volatility surfaces — asking which of this chapter's single-name conclusions survive, and what genuinely new phenomena appear, once more than one underlying is in play.

---
### Chapter 11: Multi-asset stochastic volatility

Every model built in Chapters 6-10 has priced a single underlying. Chapter 11 extends the machinery to baskets and indexes, and opens with a claim worth taking seriously: multi-asset stochastic volatility models are not mere juxtapositions of single-asset models. Beyond the ordinary asset/asset correlation matrix a multi-asset Black-Scholes model would need, a genuine multi-asset stochastic volatility model requires cross spot/volatility and volatility/volatility correlations that have no single-asset analogue at all. The chapter's task is to parametrize these new degrees of freedom parsimoniously, and connect them to two observable, tradeable quantities: the basket's own ATMF skew, and the correlation swap.

Section 11.1 starts at the short-maturity limit, where a basket $B=\sum_i\alpha_iS_i$ with relative weights $w_i=\alpha_iS_i/B$ has return

$$
\begin{equation*}
\frac{d B}{B}=\sum_{i} w_{i} \frac{d S_{i}}{S_{i}} \tag{11.1}
\end{equation*}
$$
[[Ch11#^eq11-1|Eq. (11.1)]]

and, applying Chapter 8's short-maturity skew formula directly to the basket as if it were a single asset,

$$
\begin{align*}
\widehat{\sigma}_{B} & =\sigma_{B}  \tag{11.2a}\\
\mathcal{S}_{B} & =\frac{1}{2 \widehat{\sigma}_{B}^{2}} \frac{\left\langle d \ln B d \widehat{\sigma}_{B}\right\rangle}{d t} \tag{11.2b}
\end{align*}
$$
[[Ch11#^eq11-2a|Eq. (11.2a)]] [[Ch11#^eq11-2b|Eq. (11.2b)]]

with the basket's own instantaneous volatility built from the components' volatilities and pairwise correlations,

$$
\begin{equation*}
\widehat{\sigma}_{B}=\sqrt{\sum_{i j} w_{i} w_{j} \rho_{i j} \sigma_{i} \sigma_{j}}=\sqrt{\sum_{i j} w_{i} w_{j} \rho_{i j} \widehat{\sigma}_{i} \widehat{\sigma}_{j}} \tag{11.3}
\end{equation*}
$$
[[Ch11#^eq11-3|Eq. (11.3)]]

Differentiating this and substituting into (11.2b) gives the chapter's founding decomposition of the basket skew,

$$
\begin{equation*}
\mathcal{S}_{B}=\frac{1}{\widehat{\sigma}_{B}^{3}} \sum_{i j} w_{i} w_{j}^{2} \rho_{i j} \widehat{\sigma}_{i} \widehat{\sigma}_{j}^{2} \mathcal{S}_{j}+\frac{1}{2 \widehat{\sigma}_{B}^{3}} \sum_{i j \neq k} w_{i} w_{j} w_{k} \rho_{i j} \widehat{\sigma}_{i} \frac{\left\langle d \ln S_{k} d \widehat{\sigma}_{j}\right\rangle}{d t} \tag{11.4}
\end{equation*}
$$
[[Ch11#^eq11-4|Eq. (11.4)]]

split cleanly into a piece generated by each component's *own* ATMF skew, and a piece generated by the cross-covariance of one component's spot with a *different* component's implied volatility — a distinction with no single-asset counterpart at all.

Section 11.1.1 specializes to a homogeneous, equally-weighted basket of $n$ assets, where these formulas simplify to

$$
\begin{align*}
\widehat{\sigma}_{B} & =\sqrt{\frac{1+(n-1) \rho_{S S}}{n}} \widehat{\sigma}  \tag{11.5a}\\
\mathcal{S}_{B} & =\frac{1+(n-1) \rho_{S S}}{n} \frac{\widehat{\sigma}^{3}}{\widehat{\sigma}_{B}^{3}}\left[\frac{\mathcal{S}}{n}+\frac{n-1}{n} \frac{1}{2 \widehat{\sigma}^{2}} \frac{\langle d \ln S d \widehat{\sigma}\rangle_{\text {cross }}}{d t}\right] \tag{11.5b}
\end{align*}
$$
[[Ch11#^eq11-5a|Eq. (11.5a)]] [[Ch11#^eq11-5b|Eq. (11.5b)]]

Immediately visible in (11.5b) is a striking asymptotic fact: the components'-own-skew contribution carries an explicit $1/n$ prefactor and vanishes as the basket grows large, while the cross-covariance term does not. Taking $n\gg1$ — the regime of a real equity index — leaves

$$
\begin{align*}
& \widehat{\sigma}_{B} \simeq \sqrt{\rho_{S S}} \widehat{\sigma}  \tag{11.6a}\\
& \mathcal{S}_{B} \simeq \frac{1}{\sqrt{\rho_{S S}}} \frac{1}{2 \widehat{\sigma}^{2}} \frac{\langle d \ln S d \widehat{\sigma}\rangle_{\text {cross }}}{d t} \tag{11.6b}
\end{align*}
$$
[[Ch11#^eq11-6a|Eq. (11.6a)]] [[Ch11#^eq11-6b|Eq. (11.6b)]]

A large index's ATMF skew is thus governed *entirely* by cross-asset spot/volatility covariance, not by anything intrinsic to the individual stocks' own skew — a genuinely new, purely multi-asset phenomenon. Two limiting cases bracket what this covariance can do: if a component's volatility is uncorrelated with the *other* components, $\mathcal S_B$ vanishes identically, regardless of how skewed each stock's own smile is. At the opposite extreme, if every component's volatility is driven by the same single Brownian motion — the cross covariance equal to the diagonal (own-asset) covariance defining each stock's own skew,

$$
\begin{equation*}
\mathcal{S}=\frac{1}{2 \widehat{\sigma}^{2}} \frac{\langle d \ln S d \widehat{\sigma}\rangle_{\mathrm{diag}}}{d t} \tag{11.7}
\end{equation*}
$$
[[Ch11#^eq11-7|Eq. (11.7)]]

— the basket skew is amplified relative to the component skew,

$$
\begin{equation*}
\mathcal{S}_{B} \simeq \frac{1}{\sqrt{\rho_{S S}}} \mathcal{S} \tag{11.8}
\end{equation*}
$$
[[Ch11#^eq11-8|Eq. (11.8)]]

Section 11.1.2 checks what the ubiquitous multi-asset local volatility model predicts for this ratio. Because local-volatility implied volatilities are driven by the same Brownian motion as the spot itself, cross spot/volatility covariance there is forced to equal $\rho_{SS}$ times the diagonal covariance, and (11.5) collapses to

$$
\begin{align*}
\widehat{\sigma}_{B} & =\sqrt{\frac{1+(n-1) \rho_{S S}}{n}} \widehat{\sigma} \simeq \sqrt{\rho_{S S}} \widehat{\sigma}  \tag{11.9a}\\
\mathcal{S}_{B} & =\sqrt{\frac{1+(n-1) \rho_{S S}}{n}} \mathcal{S} \simeq \sqrt{\rho_{S S}} \mathcal{S} \tag{11.9b}
\end{align*}
$$
[[Ch11#^eq11-9a|Eq. (11.9a)]] [[Ch11#^eq11-9b|Eq. (11.9b)]]

giving the falsifiable, clean short-maturity prediction that the basket-to-component volatility ratio and skew ratio should be *identical*, both equal to $\sqrt{\rho_{SS}}$,

$$
\begin{equation*}
\frac{\mathcal{S}_{B}}{\mathcal{S}} \simeq \frac{\widehat{\sigma}_{B}}{\widehat{\sigma}} \simeq \sqrt{\rho_{S S}} \tag{11.10}
\end{equation*}
$$
[[Ch11#^eq11-10|Eq. (11.10)]]

Section 11.1.3 confronts this prediction with real data, and the result is a genuine, sharp reversal. Typical implied index correlations of $\rho_{SS}\simeq60\%$ would, via (11.10), predict an index skew roughly 20% *lower* than the component skew — but real market index skews are typically about 25% *larger* than component skews, which is in fact close to what (11.8)'s 100%-correlated-volatilities case predicts. The local volatility model gets not merely the wrong number here, but the wrong *sign* of the effect entirely.

Section 11.1.4 addresses a subtlety needed before even applying these large-basket formulas to a real index: real indexes are not equally weighted. Assuming for simplicity that all component volatilities equal $\sigma$, the basket variance is

$$
\begin{equation*}
\sigma_{B}^{2}=\sigma^{2} \sum_{i j} w_{i} w_{j} \rho_{i j}=\sigma^{2}\left(\sum_{i} w_{i}^{2}+\rho_{S S} \sum_{i \neq j} w_{i} w_{j}\right) \tag{11.11}
\end{equation*}
$$
[[Ch11#^eq11-11|Eq. (11.11)]]

which, matched against (11.5a), defines the *effective* number of index components,

$$
\begin{equation*}
n^{\star}=\frac{1}{\Sigma_{i} w_{i}^{2}} \tag{11.12}
\end{equation*}
$$
[[Ch11#^eq11-12|Eq. (11.12)]]

which can be far smaller than the actual component count $n$: the Euro Stoxx 50's $n^\star=37$ against $n=50$ is a modest gap, but the KOSPI's $n^\star=17$ against $n=200$ shows just how concentrated a cap-weighted index's effective diversification can be — a practical correction that must be applied before comparing any real index against the large-$n$ asymptotic formulas above.

Section 11.2 builds the chapter's central modeling contribution: a parsimonious, factor-structure-independent parametrization of the new cross-correlations. For a homogeneous basket where each asset is driven by the two-factor model with identical parameters, the diagonal (own-asset) spot/variance and variance/variance covariance functions are

$$
\begin{equation*}
\mu^{\operatorname{diag}}(t, u, \xi)=2 \nu \alpha_{\theta} \sqrt{\xi^{t}} \xi^{u}\left((1-\theta) \rho_{S X_{1}}^{\operatorname{diag}} e^{-k_{1}(u-t)}+\theta \rho_{S X_{2}}^{\operatorname{diag}} e^{-k_{2}(u-t)}\right) \tag{11.13}
\end{equation*}
$$
[[Ch11#^eq11-13|Eq. (11.13)]]

$$
\begin{align*}
\nu^{\operatorname{diag}}\left(t, u, u^{\prime}, \xi\right) & =4 \nu^{2} \xi^{u} \xi^{u^{\prime}} \alpha_{\theta}^{2}\left[(1-\theta)^{2} \rho_{X_{1} X_{1}}^{\operatorname{diag}} e^{-k_{1}\left(u+u^{\prime}-2 t\right)}\right.  \tag{11.14}\\
& +\theta^{2} \rho_{X_{2} X_{2}}^{\operatorname{diag}} e^{-k_{2}\left(u+u^{\prime}-2 t\right)} \\
& \left.+\rho_{X_{1} X_{2}}^{\operatorname{diag}} \theta(1-\theta)\left(e^{-k_{1}(u-t)} e^{-k_{2}\left(u^{\prime}-t\right)}+e^{-k_{2}(u-t)} e^{-k_{1}\left(u^{\prime}-t\right)}\right)\right]
\end{align*}
$$
[[Ch11#^eq11-14|Eq. (11.14)]]

Rather than parametrizing the cross factor/factor and spot/factor correlations directly — quantities whose meaning depends on the model's specific factor structure — the chapter introduces exactly two dimensionless multipliers, $\chi_{S\sigma}$ and $\chi_{\sigma\sigma}$, that rescale the cross covariance functions relative to their diagonal counterparts uniformly, independent of factor structure entirely. Whether this parametrization is even admissible is not a free choice: since each asset carries three Brownian motions, the global correlation matrix is $3n\times3n$, and demanding it be positive semi-definite is a real constraint. Diagonalizing in the basis of the discrete translation operator (since the correlation structure is invariant under permuting assets) reduces the problem to two eigenvalue conditions on a per-mode quantity,

$$
A_{S S}(k)=\left(1+\rho_{S S} \Sigma_{j=1}^{j=n-1} e^{i j \theta_{k}}\right)= \begin{cases}1+(n-1) \rho_{S S} & k=0  \tag{11.15}\\ 1-\rho_{S S} & k \neq 0\end{cases}
$$
[[Ch11#^eq11-15|Eq. (11.15)]]

with $A_{\sigma\sigma}(k),A_{S\sigma}(k)$ defined analogously in terms of $\chi_{\sigma\sigma},\chi_{S\sigma}$, giving a rescaled effective correlation matrix

$$
\rho(k)=\left(\begin{array}{ccc}
1 & \zeta \rho_{S X_{1}} & \zeta \rho_{S X_{2}}  \tag{11.16}\\
\zeta \rho_{S X_{1}} & 1 & \rho_{X_{1} X_{2}} \\
\zeta \rho_{S X_{2}} & \rho_{X_{1} X_{2}} & 1
\end{array}\right)
$$
[[Ch11#^eq11-16|Eq. (11.16)]]

whose positivity, evaluated at exactly the two modes $k=0$ and $k\neq0$, gives the necessary and sufficient conditions,

$$
\begin{align*}
\zeta & =\frac{1+(n-1) \chi_{S \sigma}}{\sqrt{\left(1+(n-1) \rho_{S S}\right)\left(1+(n-1) \chi_{\sigma \sigma}\right)}}  \tag{11.17a}\\
\zeta & =\frac{1-\chi_{S \sigma}}{\sqrt{\left(1-\rho_{S S}\right)\left(1-\chi_{\sigma \sigma}\right)}} \tag{11.17b}
\end{align*}
$$
[[Ch11#^eq11-17a|Eq. (11.17a)]] [[Ch11#^eq11-17b|Eq. (11.17b)]]

together with the ordinary single-asset correlation-matrix positivity condition applied at each of these two $\zeta$ values. Notably, $\chi_{S\sigma}$ — unlike $\rho_{SS}$ and $\chi_{\sigma\sigma}$, both confined to $[-1/(n-1),1]$ — is allowed to exceed 1, which turns out to be exactly what is needed whenever the index skew is much steeper than the component skew, precisely the empirical situation Section 11.1.3 uncovered. The special case $\chi_{S\sigma}=\chi_{\sigma\sigma}=\rho_{SS}$ collapses $\zeta$ to exactly 1, reproducing the single-asset correlation matrix directly and automatically guaranteeing positivity — this, it turns out, is exactly the multi-asset local volatility model's own implicit correlation structure, a fact Section 11.3.3 returns to.

Sections 11.3-11.3.3 extend the basket skew analysis beyond the short-maturity limit, using Chapter 8's general order-one skew formula applied to the basket,

$$
\begin{equation*}
\mathcal{S}_{T}^{B}=\hat{\sigma}_{T}^{B} \frac{C_{B}^{x \xi}(T)}{2\left(\left(\hat{\sigma}_{T}^{B}\right)^{2} T\right)^{2}} \tag{11.18}
\end{equation*}
$$
[[Ch11#^eq11-18|Eq. (11.18)]]

with the basket's own doubly-integrated spot/variance covariance built from the basket's forward variance curve,

$$
\begin{align*}
C_{B}^{x \xi}(T) & =\int_{0}^{T} d t \int_{t}^{T} d u \mu_{B}(t, u)  \tag{11.19}\\
\mu_{B}(t, u) & =\lim _{d t \rightarrow 0} \frac{1}{d t} E_{t}\left[d \ln B_{t} d \zeta_{t}^{u}\right] \tag{11.20}
\end{align*}
$$
[[Ch11#^eq11-19|Eq. (11.19)]] [[Ch11#^eq11-20|Eq. (11.20)]]

For a homogeneous, equally-weighted, frozen-weight basket using the $\chi_{S\sigma},\chi_{\sigma\sigma}$ parametrization, this reduces to a closed form,

$$
\begin{equation*}
\mathcal{S}_{T}^{B}=\frac{\hat{\sigma}_{0 T}^{B}}{2\left(\left(\hat{\sigma}_{0 T}^{B}\right)^{2} T\right)^{2}}\left[\frac{1+(n-1) \chi_{S \sigma}}{n}\right]\left[\frac{1+(n-1) \rho_{S S}}{n}\right] \int_{0}^{T} d t \int_{t}^{T} d u \mu(u-t) \tag{11.23}
\end{equation*}
$$
[[Ch11#^eq11-23|Eq. (11.23)]]

which, rewritten as a ratio to the single-component skew, gives the chapter's key general working formula,

$$
\begin{equation*}
\mathcal{S}_{T}^{B}=\mathcal{S}_{T} \sqrt{\frac{n}{1+(n-1) \rho_{S S}}}\left[\frac{1+(n-1) \chi_{S \sigma}}{n}\right] \tag{11.24}
\end{equation*}
$$
[[Ch11#^eq11-24|Eq. (11.24)]]

directly invertible to back an implied value of $\chi_{S\sigma}$ out of an observed ratio of basket skew to component skew — turning the chapter's parsimonious parametrization into a calibratable, market-observable quantity. Taking $T\to0$ recovers the earlier short-maturity result exactly,

$$
\begin{align*}
\mathcal{S}_{T=0}^{B} & =\frac{1}{4\left(\hat{\sigma}_{T=0}^{B}\right)^{3}}\left[\frac{1+(n-1) \chi_{S \sigma}}{n}\right]\left[\frac{1+(n-1) \rho_{S S}}{n}\right] \mu(0) \\
& =\frac{1}{4\left(\hat{\sigma}_{T=0}^{B}\right)^{3}}\left[\frac{1+(n-1) \chi_{S \sigma}}{n}\right]\left[\frac{1+(n-1) \rho_{S S}}{n}\right] 4 \widehat{\sigma}_{T=0}^{3} \mathcal{S}_{T=0} \\
& =\mathcal{S}_{T=0}\left(\frac{\widehat{\sigma}_{T=0}}{\hat{\sigma}_{T=0}^{B}}\right)^{3}\left[\frac{1+(n-1) \chi_{S \sigma}}{n}\right]\left[\frac{1+(n-1) \rho_{S S}}{n}\right] \tag{11.25}
\end{align*}
$$
[[Ch11#^eq11-25|Eq. (11.25)]]

a consistency check that the general formula correctly degenerates to the special case already derived, exactly as Chapter 8's general smile-expansion machinery was checked against Chapter 6's Heston-specific results.

Section 11.3.1 works this out explicitly for the two-factor model. With a flat term structure the spot/variance covariance function is

$$
\begin{equation*}
\mu(\tau)=(2 \nu) \hat{\sigma}^{3} \alpha_{\theta}\left((1-\theta) \rho_{S X_{1}} e^{-k_{1} \tau}+\theta \rho_{S X_{2}} e^{-k_{2} \tau}\right) \tag{11.26}
\end{equation*}
$$
[[Ch11#^eq11-26|Eq. (11.26)]]

giving the component's own ATMF skew via Chapter 8's formula,

$$
\begin{equation*}
\mathcal{S}_{T}=\nu \alpha_{\theta}\left[(1-\theta) \rho_{S X^{1}} \frac{k_{1} T-\left(1-e^{-k_{1} T}\right)}{\left(k_{1} T\right)^{2}}+\theta \rho_{S X^{2}} \frac{k_{2} T-\left(1-e^{-k_{2} T}\right)}{\left(k_{2} T\right)^{2}}\right] \tag{11.27}
\end{equation*}
$$
[[Ch11#^eq11-27|Eq. (11.27)]]

and, substituting into (11.24), the basket skew directly,

$$
\begin{align*}
\mathcal{S}_{T}^{B} & =\nu \alpha_{\theta} \sqrt{\frac{n}{1+(n-1) \rho_{S S}}}\left[\frac{1+(n-1) \chi_{S \sigma}}{n}\right] \\
& \times\left[(1-\theta) \rho_{S X^{1}} \frac{k_{1} T-\left(1-e^{-k_{1} T}\right)}{\left(k_{1} T\right)^{2}}+\theta \rho_{S X^{2}} \frac{k_{2} T-\left(1-e^{-k_{2} T}\right)}{\left(k_{2} T\right)^{2}}\right] \tag{11.28}
\end{align*}
$$
[[Ch11#^eq11-28|Eq. (11.28)]]

The basket VS volatility needed alongside this, unlike the skew, can be computed *exactly* rather than perturbatively, since it only requires the volatility of volatility's effect on the zeroth-order basket variance, governed by a cross-asset volatility-of-volatility decay function

$$
\begin{equation*}
h(\tau)=e^{-\nu^{2} \alpha_{\theta}^{2}\left(1-\chi_{\sigma \sigma}\right)\left((1-\theta)^{2} \frac{1-e^{-2 k_{1} \tau}}{2 k_{1}}+\theta^{2} \frac{1-e^{-2 k_{2} \tau}}{2 k_{2}}+2 \theta(1-\theta) \rho_{X_{1} X_{2}} \frac{1-e-\left(k_{1}+k_{2}\right) \tau}{k_{1}+k_{2}}\right)} \tag{11.29}
\end{equation*}
$$
[[Ch11#^eq11-29|Eq. (11.29)]]

giving, for identical flat term structures,

$$
\begin{equation*}
\left(\hat{\sigma}_{T}^{B}\right)^{2}=\widehat{\sigma}^{2} \frac{1+(n-1) \rho_{S S} \frac{1}{T} \int_{0}^{T} h(\tau) d \tau}{n} \tag{11.30}
\end{equation*}
$$
[[Ch11#^eq11-30|Eq. (11.30)]]

Numerical checks (Section 11.3.2) confirm (11.24) and (11.30) are accurate against a full Monte Carlo of the two-factor model, validating the perturbative shortcut for the skew alongside the exact formula for volatility.

Section 11.3.3 asks whether this two-factor machinery can *mimic* the multi-asset local volatility model exactly. Since the local volatility model forces $\rho^{\text{cross}}(S,\widehat\sigma_T)=\rho_{SS}\,\rho^{\text{diag}}(S,\widehat\sigma_T)$ and likewise for the volatility/volatility cross-correlation, it corresponds precisely to the parameter choice

$$
\begin{align*}
& \chi_{S \sigma}=\rho_{S S}  \tag{11.31a}\\
& \chi_{\sigma \sigma}=\rho_{S S} \tag{11.31b}
\end{align*}
$$
[[Ch11#^eq11-31a|Eq. (11.31a)]] [[Ch11#^eq11-31b|Eq. (11.31b)]]

exactly the $\zeta\equiv1$ case identified back in Section 11.2.1. Tested numerically, though, this mapping is only accurate for *short* maturities; for longer maturities the two-factor model's implied volatilities run systematically lower than the local volatility model's — because forward volatilities are themselves genuinely more volatile in the stochastic volatility model, the effective level of spot/spot correlation the SV model delivers ends up lower than the nominal $\rho_{SS}$ used to calibrate it. This same effect, the chapter notes, will resurface directly in the correlation swap.

Section 11.4 turns to that instrument and delivers the chapter's most genuinely counter-intuitive result. A correlation swap of maturity $T$ pays the average realized pairwise correlation minus a fixed strike,

$$
\begin{equation*}
\frac{1}{n(n-1)} \sum_{i \neq j} \rho_{i j}-\widehat{\rho} \tag{11.32}
\end{equation*}
$$
[[Ch11#^eq11-32|Eq. (11.32)]]

with each pairwise correlation estimated from daily log-returns in the standard way,

$$
\begin{equation*}
\rho_{i j}=\frac{\Sigma r_{k}^{i} r_{k}^{j}}{\sqrt{\Sigma r_{k}^{i 2}} \sqrt{\Sigma r_{k}^{j^{2}}}} \tag{11.33}
\end{equation*}
$$
[[Ch11#^eq11-33|Eq. (11.33)]]

For centered, constant-volatility log-returns, this single-pair estimator's bias and standard deviation, at leading order in $1/N$, are

$$
\begin{align*}
E\left[\rho_{i j}\right] & =\rho_{S S}\left(1-\frac{1-\rho_{S S}^{2}}{2 N}\right)  \tag{11.34}\\
\operatorname{Stdev}\left(\rho_{i j}\right) & =\frac{1}{\sqrt{N}}\left(1-\rho_{S S}^{2}\right) \tag{11.35}
\end{align*}
$$
[[Ch11#^eq11-34|Eq. (11.34)]] [[Ch11#^eq11-35|Eq. (11.35)]]

Naively, one might expect that *averaging* many such pairwise estimators over a large basket — letting $n\to\infty$ — drives this standard deviation to zero, exactly as averaging independent measurements ordinarily does. It does not. In the $n\to\infty$ limit, the standard deviation of the averaged correlation-swap estimator is

$$
\begin{align*}
\lim _{n \rightarrow \infty} E\left[\frac{1}{n(n-1)} \sum_{i \neq j} \rho_{i j}\right] & =\rho_{S S}\left(1-\frac{1-\rho_{S S}^{2}}{2 N}\right) \\
\lim _{n \rightarrow \infty} \operatorname{Stdev}\left(\frac{1}{n(n-1)} \sum_{i \neq j} \rho_{i j}\right) & =\frac{1}{\sqrt{N}}\left(1-\rho_{S S}^{2}\right) \frac{\sqrt{2}\left|\rho_{S S}\right|}{1+\rho_{S S}} \tag{11.36}
\end{align*}
$$
[[Ch11#^eq11-36|Eq. (11.36)]]

only *marginally* smaller than the single-pair standard deviation (11.35) at $n=2$ — the averaging that would ordinarily wash out estimation noise almost entirely fails to do so here, because pairwise correlation estimators computed from overlapping return series are themselves strongly correlated with one another. As a measure of correlation, averaging pairwise correlations is thus, despite being the market-standard payoff, a poorly defined estimator — a genuinely surprising fact about an instrument whose entire purpose is supposedly to measure correlation cleanly.

Beyond this measurement pathology, the correlation swap's *fair strike* is itself model-dependent, not simply equal to the underlying $\rho_{SS}$. In the limit of short returns, for a homogeneous basket with instantaneous (lognormal) volatilities $\sigma_t^i$,

$$
\begin{equation*}
\widehat{\rho}(T)=E\left[\rho_{i j}\right]=\rho_{S S} E\left[\frac{\int_{0}^{T} \sigma_{t}^{i} \sigma_{t}^{j} d t}{\sqrt{\int_{0}^{T} \sigma_{t}^{i 2} d t} \sqrt{\int_{0}^{T} \sigma_{t}^{j^{2}} d t}}\right] \tag{11.37}
\end{equation*}
$$
[[Ch11#^eq11-37|Eq. (11.37)]]

which equals $\rho_{SS}$ exactly only in the special case that the two underlyings' instantaneous volatilities are collinear processes — otherwise $\widehat\rho\le\rho_{SS}$ strictly, with the gap governed by exactly the sort of cross-volatility covariance structure Section 11.2 was built to parametrize.

Risk-managing a correlation swap on two underlyings with variance swaps as the natural hedge reveals why this gap is not a mere academic curiosity. Because the swap's value has vanishing sensitivity to implied volatilities at inception, its residual exposure shows up entirely as a cross-gamma/theta P&L; for the special case where realized volatilities and correlation match their implied values, this P&L reduces to a clean quadratic form in the mismatch of the two underlyings' relative VS-variance moves,

$$
\begin{align*}
P \& L & =\frac{\widehat{\rho}_{0}}{8} \frac{t(T-t)}{T^{2}}\left(\frac{\left(\delta\left(\widehat{\sigma}_{1 t}^{2}\right)\right)^{2}}{\left(\widehat{\sigma}_{1 t}^{2}\right)^{2}}+\frac{\left(\delta\left(\widehat{\sigma}_{2 t}^{2}\right)\right)^{2}}{\left(\widehat{\sigma}_{2 t}^{2}\right)^{2}}-2 \frac{\delta\left(\widehat{\sigma}_{1 t}^{2}\right) \delta\left(\widehat{\sigma}_{2 t}^{2}\right)}{\widehat{\sigma}_{1 t}^{2} \widehat{\sigma}_{2 t}^{2}}\right) \\
& =\frac{\widehat{\rho}_{0}}{8} \frac{t(T-t)}{T^{2}}\left(\frac{\delta\left(\widehat{\sigma}_{1 t}^{2}\right)}{\widehat{\sigma}_{1 t}^{2}}-\frac{\delta\left(\widehat{\sigma}_{2 t}^{2}\right)}{\widehat{\sigma}_{2 t}^{2}}\right)^{2} \tag{11.38}
\end{align*}
$$
[[Ch11#^eq11-38|Eq. (11.38)]]

vanishing exactly when the two VS volatilities move in lockstep proportionally, and nonzero otherwise — confirming that a correlation swap is, in fact, a genuine exotic volatility-of-volatility instrument, not merely a static bet on today's correlation level.

Section 11.4.1 works this out quantitatively in the two-factor model. Writing the mixing-process representation of the instantaneous variance at $T=t$,

$$
\begin{equation*}
x_{t}^{t}=\alpha_{\theta}\left[(1-\theta) X_{t}^{1}+\theta X_{t}^{2}\right] \tag{11.39}
\end{equation*}
$$
[[Ch11#^eq11-39|Eq. (11.39)]]

and expanding the lognormal instantaneous variance to first order in the vol-of-vol parameter $\nu$,

$$
\begin{equation*}
\xi_{t}^{t}=\xi_{0}^{t}\left(1+2 \nu x_{t}^{t}+2 \nu^{2}\left[\left(x_{t}^{t}\right)^{2}-E\left[\left(x_{t}^{t}\right)^{2}\right]\right]\right) \tag{11.40}
\end{equation*}
$$
[[Ch11#^eq11-40|Eq. (11.40)]]

then computing $\widehat\rho(T)$ order by order in $\nu$ and re-exponentiating (to guarantee $\widehat\rho$ stays non-negative for large $\nu$ rather than crossing zero) gives

$$
\begin{equation*}
\widehat{\rho}(T)=\rho_{S S} e^{-\frac{\nu^{2}}{2} E\left[\overline{(x-y)^{2}}-\overline{x-y}^{2}\right]} \tag{11.41}
\end{equation*}
$$
[[Ch11#^eq11-41|Eq. (11.41)]]

which, by a pathwise Cauchy-Schwarz-type inequality on the bracketed term, automatically guarantees $\widehat\rho\le\rho_{SS}$ — consistent with the general property already derived in (11.37) — and collapses back to exact equality precisely when the two underlyings' mixing processes coincide. Specialized to the two-factor model's own $x_t^t$ from (11.39), this produces the chapter's final explicit correlation-swap formula,

$$
\begin{equation*}
\hat{\rho}(T)=\rho_{S S} e^{-\left(1-\chi_{\sigma \sigma}\right)\left(\nu^{2} T\right) \frac{1}{T} \int_{0}^{T} f(t) d t} \tag{11.42}
\end{equation*}
$$
[[Ch11#^eq11-42|Eq. (11.42)]]

Section 11.4.2's numerical examples show the gap $\rho_{SS}-\widehat\rho(T)$ implied by this formula *widens* with maturity as vol-of-vol accumulates its effect over a longer horizon — the opposite of the term structure ordinarily observed for market-implied correlation-swap strikes, and sharply different from the local volatility model's near-flat prediction for the identical calibrated smile. A model built to be internally consistent with the vanilla and VIX markets can, on this one instrument, point in a direction opposite to what the correlation-swap market itself typically shows — a genuine, quantifiable tension a desk pricing correlation products needs to be aware of, not a modeling error to be patched away.

Appendix A derives the single-pair bias and variance results (11.34)-(11.35) used above from first principles, expanding the correlation estimator in powers of centered second-moment fluctuations $\varepsilon_{ij}=x_i\cdot x_j-\rho$,

$$
\begin{equation*}
\widehat{\rho}_{i j}=\rho-\left(\rho \frac{\varepsilon_{i i}+\varepsilon_{j j}}{2}-\varepsilon_{i j}\right)-\left(\varepsilon_{i j} \frac{\varepsilon_{i i}+\varepsilon_{j j}}{2}-\rho\left(\frac{3}{8}\left(\varepsilon_{i i}^{2}+\varepsilon_{j j}^{2}\right)+\frac{1}{4} \varepsilon_{i i} \varepsilon_{j j}\right)\right) \tag{11.43}
\end{equation*}
$$
[[Ch11#^eq11-43|Eq. (11.43)]]

and averaging this expansion over all pairs to isolate the leading-order behavior of the multi-asset estimator,

$$
\begin{equation*}
\widehat{\rho}=\rho-\frac{1}{n(n-1)} \sum_{i \neq j}\left(\rho \frac{\varepsilon_{i i}+\varepsilon_{j j}}{2}-\varepsilon_{i j}\right) \tag{11.44}
\end{equation*}
$$
[[Ch11#^eq11-44|Eq. (11.44)]]

whose variance, taken in the $n\to\infty$ limit using Wick's theorem on the Gaussian fourth-order moments of the $\varepsilon_{ij}$, is exactly the surprising result (11.36) derived above — the "poorly defined estimator" paradox is thus not a numerical artifact but a direct, derivable consequence of the correlation structure among overlapping pairwise estimators themselves.

#### Chapter's Digest and Logical Bridge to Chapter 12

Chapter 11 showed that multi-asset stochastic volatility genuinely needs two new parameters beyond an asset/asset correlation matrix — $\chi_{S\sigma}$ and $\chi_{\sigma\sigma}$ — and that these are not academic refinements: real index skews are 25% steeper than component skews, the opposite sign from what the ubiquitous multi-asset local volatility model predicts, and the correlation swap's own fair strike is itself a genuinely model-dependent quantity that a pure local volatility view gets structurally wrong. In a pure local-volatility world, implied correlations backed out of different multi-asset payoffs used to cluster tightly around one number; with genuine multi-asset stochastic volatility, $\chi_{S\sigma}$ (measurable from the basket-to-component skew ratio) and $\chi_{\sigma\sigma}$ (measurable from the correlation-swap strike) become separate, independently-priced risks. This closes the book's excursion outward to multiple underlyings. Chapter 12, the book's final chapter, turns back to a single underlying but asks the remaining structural question left open since Chapter 2: can local volatility and stochastic volatility be combined into a single, mutually consistent model — a local-stochastic volatility model — that keeps the smile-calibration precision of local volatility while retaining the genuine forward-smile and vol-of-vol dynamics only a stochastic volatility model can supply?

---
### Chapter 12: Local-stochastic volatility models

The book's final chapter is explicitly a synthesis of its two main threads. The local volatility model (Chapter 2) is the simplest market model there is — one-dimensional, calibrated exactly to the vanilla smile — but offers zero control over the dynamics of implied volatilities it generates. Forward variance models (Chapter 7) sit at the opposite extreme: genuinely rich, tunable dynamics, but generally incapable of matching the full vanilla surface exactly. Local-stochastic ("mixed") volatility models are the natural next rung in the hierarchy: a market model with a Markovian representation in $(t,S)$ plus a handful of extra state variables, aiming to combine local volatility's exact calibration with stochastic volatility's genuine dynamical flexibility.

Section 12.2.1 defines the mixed-model ansatz directly: instantaneous volatility is the *product* of a stochastic-volatility factor and a local volatility function,

$$
\begin{equation*}
\sigma_{t}=\sqrt{\zeta_{t}^{t}} \sigma\left(t, S_{t}\right) \tag{12.1}
\end{equation*}
$$
[[Ch12#^eq12-1|Eq. (12.1)]]

with $\zeta_t^t$ driven by the two-factor model of Section 7.4 — deliberately relabeled $\zeta$ rather than $\xi$, since it is no longer itself a forward variance once multiplied by the local volatility component,

$$
\left\{\begin{align*}
d S_{t} & =(r-q) S_{t} d t+\sigma\left(t, S_{t}\right) \sqrt{\zeta_{t}^{t}} S_{t} d W_{t}^{S}  \tag{12.2}\\
d \zeta_{t}^{T} & =2 \nu \zeta_{t}^{T} \alpha_{\theta}\left((1-\theta) e^{-k_{1}(T-t)} d W_{t}^{1}+\theta e^{-k_{2}(T-t)} d W_{t}^{2}\right)
\end{align*}\right.
$$
[[Ch12#^eq12-2|Eq. (12.2)]]

The resulting pricing equation looks, on its face, like a straightforward generalization of Chapter 7's forward-variance pricing PDE, with the local volatility function simply inserted into the gamma and cross-gamma terms,

$$
\begin{align*}
& \frac{d P}{d t}+(r-q) S \frac{d P}{d S}+\frac{\zeta^{t} \sigma(t, S)^{2}}{2} S^{2} \frac{d^{2} P}{d S^{2}} \\
+ & \frac{1}{2} \int_{t}^{T} d u \int_{t}^{T} d u^{\prime} \nu\left(t, u, u^{\prime}, \zeta\right) \frac{d^{2} P}{\delta \zeta^{u} \delta \zeta^{u^{\prime}}}+S \sigma(t, S) \int_{t}^{T} d u \mu(t, u, \zeta) \frac{d^{2} P}{d S \delta \zeta^{u}}=r P \tag{12.3}
\end{align*}
$$
[[Ch12#^eq12-3|Eq. (12.3)]]

with $\mu,\nu$ the ordinary two-factor covariance functions carried over unchanged from Chapter 8,

$$
\begin{align*}
\mu(t, u, \xi)= & 2 \nu \xi^{u} \sqrt{\xi^{t}} \alpha_{\theta}\left[\rho_{S X^{1}}(1-\theta) e^{-k_{1}(u-t)}+\rho_{S X^{2}} \theta e^{-k_{2}(u-t)}\right]  \tag{12.4}\\
\nu\left(t, u, u^{\prime}, \xi\right)= & 4 \nu^{2} \xi^{u} \xi^{u^{\prime}} \alpha_{\theta}^{2}\left[(1-\theta)^{2} e^{-k_{1}\left(u+u^{\prime}-2 t\right)}+\theta^{2} e^{-k_{2}\left(u+u^{\prime}-2 t\right)}\right. \\
& \left.+\rho_{12} \theta(1-\theta)\left(e^{-k_{1}(u-t)} e^{-k_{2}\left(u^{\prime}-t\right)}+e^{-k_{2}(u-t)} e^{-k_{1}\left(u^{\prime}-t\right)}\right)\right]
\end{align*}
$$
[[Ch12#^eq12-4|Eq. (12.4)]]

Section 12.2.2 — posed as a genuinely open rhetorical question, "Is it a price?" — is the chapter's central methodological pivot, and it directly continues the book's opening theme from Chapter 1: a pricing equation only earns trust if it arises from a real replication argument with well-defined break-even levels. The forward-variance-model equation earned exactly this status in Chapter 7. But (12.3) is derived from it purely by the substitution

$$
\begin{equation*}
\sqrt{\zeta_{t}^{t}} \rightarrow \sqrt{\zeta_{t}^{t}} \sigma\left(t, S_{t}\right) \tag{12.5}
\end{equation*}
$$
[[Ch12#^eq12-5|Eq. (12.5)]]

— a substitution the chapter calls "seemingly innocuous," but one that is *not* itself derived from any replication argument. There is, therefore, no a priori guarantee that solving (12.3) actually produces a genuine price, in the sense that a delta-hedged position's P&L takes the ordinary gamma/theta form. This suspicion is left deliberately open here and confirmed later, with a stark verdict: most local-stochastic volatility models, it turns out, are not usable trading models at all.

Sections 12.2.3-12.2.5 set that open question aside temporarily to build the calibration machinery any mixed model needs regardless. Simulating the model reduces to simulating three processes — $S_t$ and the two OU factors,

$$
\left\{\begin{align*}
d S_{t} & =(r-q) S_{t} d t+\sqrt{\zeta_{t}^{t}} \sigma\left(t, S_{t}\right) S_{t} d W_{t}^{S}  \tag{12.6}\\
d X_{t}^{1} & =-k_{1} X_{t}^{1} d t+d W_{t}^{1} \\
d X_{t}^{2} & =-k_{2} X_{t}^{2} d t+d W_{t}^{2}
\end{align*}\right.
$$
[[Ch12#^eq12-6|Eq. (12.6)]]

with $\zeta_t^t$ read off the OU factors directly,

$$
\left\{\begin{array}{l}
\zeta_{t}^{t}=\zeta_{0}^{t} f\left(t, X_{t}^{1}, X_{t}^{2}\right)  \tag{12.7}\\
f\left(t, x_{1}, x_{2}\right)=e^{2 \nu \alpha_{\theta}\left[(1-\theta) x_{1}+\theta x_{2}\right]-\frac{(2 \nu)^{2}}{2}} \chi(t, t)
\end{array}\right.
$$
[[Ch12#^eq12-7|Eq. (12.7)]]

The remaining open task is calibrating the local volatility function $\sigma(t,S)$ itself. The general condition any diffusive model's instantaneous volatility must satisfy to match vanilla prices, already derived in Chapter 2,

$$
\begin{equation*}
E\left[\sigma_{t}^{2} \mid S_{t}=S\right]=\left.2 \frac{\frac{d C}{d T}+q C+(r-q) K \frac{d C}{d K}}{K^{2} \frac{d^{2} C}{d K^{2}}}\right|_{\substack{K=S \\ T=t}}=\sigma_{\mathrm{Mkt}}(t, S)^{2} \tag{12.8}
\end{equation*}
$$
[[Ch12#^eq12-8|Eq. (12.8)]]

translates, once $\sigma_t^2=\zeta_t^t\sigma(t,S)^2$ is substituted in, into a self-consistent equation for $\sigma(t,S)$ itself,

$$
\begin{equation*}
\sigma(t, S)^{2}=\frac{\sigma_{\mathrm{Mkt}}(t, S)^{2}}{E\left[\zeta_{t}^{t} \mid S_{t}=S\right]} \tag{12.9}
\end{equation*}
$$
[[Ch12#^eq12-9|Eq. (12.9)]]

with the unknown function appearing both explicitly on the left and implicitly, inside the conditional expectation, on the right. Solving it requires marching forward in time from a known joint density of $(S,\zeta)$: for a genuinely small number of factors this can be done by solving a forward Kolmogorov PDE directly, but this quickly becomes impractical beyond one factor, whereas the particle method of Henry-Labordère and Guyon — a Monte Carlo scheme based on simultaneously simulating interacting paths and re-estimating the conditional expectation from the empirical joint density at each time step,

$$
\begin{equation*}
\varphi_{e m}\left(t_{i+1}, S, \zeta\right)=\frac{1}{N} \Sigma_{k} \delta\left(S-S_{t_{i+1}}^{k}\right) \delta\left(\zeta-\zeta_{t_{i+1}}^{k}\right) \tag{12.23}
\end{equation*}
$$
[[Ch12#^eq12-23|Eq. (12.23)]]

— is entirely dimension-independent and immune to the curse of dimensionality, making it the practical method of choice whenever the underlying stochastic volatility model has more than one factor, exactly the two-factor case this chapter cares about.

Neither calibration technique, it is worth emphasizing, has anything to say about whether the resulting local volatility function is actually usable in the sense Section 12.2.2 raised. Calibration and usability are genuinely separate questions: a model can be calibrated flawlessly to the market smile via either the PDE or particle method, and still fail to produce a legitimate gamma/theta P&L once traded. The chapter deliberately keeps these two concerns apart before tackling the harder one next.

Section 12.3 finally delivers the payoff of Section 12.2.2's suspended question. The key distinction is between the pricing function $P^M(t,x)$ that takes the *local volatility function itself* as an explicit input — since it is not frozen, but recalibrated every time the underlying market observables move — and the pricing function $P(t,\widehat x)$ actually used in trading, parametrized directly by market observables (spot, implied volatilities) plus the underlying model's own state variables. Changing variables from the raw state $x$ to the observable-parametrized $\widehat x$ turns the mixed-model pricing PDE

$$
\begin{equation*}
\frac{d P^{\mathrm{M}}}{d t}+\left(\Sigma_{k} \mu_{k} \frac{d}{d x_{k}}+\frac{1}{2} \Sigma_{k l} a_{k l} \frac{d^{2}}{d x_{k} d x_{l}}\right) P^{\mathrm{M}}=0 \tag{12.24}
\end{equation*}
$$
[[Ch12#^eq12-24|Eq. (12.24)]]

into the equivalent equation directly in terms of observables,

$$
\begin{equation*}
\frac{d P}{d t}+\left(\Sigma_{i} \widehat{\mu}_{i} \frac{d}{d \widehat{x}_{i}}+\frac{1}{2} \Sigma_{i j} \widehat{a}_{i j} \frac{d^{2}}{d \widehat{x}_{i} d \widehat{x}_{j}}\right) P=0 \tag{12.25}
\end{equation*}
$$
[[Ch12#^eq12-25|Eq. (12.25)]]

with the transformed drift and covariance matrix given explicitly by the chain rule applied through the (implicit) recalibration map,

$$
\left\{\begin{array}{l}
\widehat{\mu}_{i}=\frac{d \widehat{x}_{i}}{d t}+\Sigma_{k} \mu_{k} \frac{d \widehat{x}_{i}}{d x_{k}}+\frac{1}{2} \Sigma_{k l} a_{k l} \frac{d^{2} \widehat{x}_{i}}{d x_{k} d x_{l}}  \tag{12.26}\\
\widehat{a}_{i j}=\Sigma_{k l} a_{k l} \frac{d \widehat{x}_{i}}{d x_{k}} \frac{d \widehat{x}_{j}}{d x_{l}}
\end{array}\right.
$$
[[Ch12#^eq12-26|Eq. (12.26)]]

Expanding the P&L of an *unhedged* short position over $\delta t$ in these observable variables, splitting them explicitly into market observables $O_i$ (spot, implied volatilities) and the underlying model's own state variables $\lambda_k$, gives the full decomposition

$$
\begin{align*}
P \& L= & -\Sigma_{i} \frac{d P}{d O_{i}}\left(\delta O_{i}-\widehat{\mu}_{i} \delta t\right)-\frac{1}{2} \Sigma_{i j} \frac{d^{2} P}{d O_{i} d O_{j}}\left(\delta O_{i} \delta O_{j}-\widehat{a}_{i j} \delta t\right) \\
& -\Sigma_{k} \frac{d P}{d \lambda_{k}}\left(\delta \lambda_{k}-\widehat{\mu}_{k} \delta t\right) \\
& -\frac{1}{2} \Sigma_{k l} \frac{d^{2} P}{d \lambda_{k} d \lambda_{l}}\left(\delta \lambda_{k} \delta \lambda_{l}-\widehat{a}_{k l} \delta t\right)-\Sigma_{i k} \frac{d^{2} P}{d O_{i} d \lambda_{k}}\left(\delta O_{i} \delta \lambda_{k}-\widehat{a}_{i k} \delta t\right) \tag{12.27}
\end{align*}
$$
[[Ch12#^eq12-27|Eq. (12.27)]]

Whether the hedge portfolio is parametrized by straight vanilla option prices,

$$
\begin{equation*}
P_{\mathrm{H}}=P-\left.\frac{d P}{d S}\right|_{\lambda, P_{K T}} S-\sum_{K T} \Delta_{K T} P_{K T} \tag{12.29}
\end{equation*}
$$
[[Ch12#^eq12-29|Eq. (12.29)]]

or equivalently by Black-Scholes implied volatilities,

$$
\begin{equation*}
P_{\mathrm{H}}=P-\left.\frac{d P}{d S}\right|_{\lambda, \widehat{\sigma}_{K T}} S-\sum_{K T} \Delta_{K T}\left(P_{K T}-\frac{d P_{K T}^{B S}}{d S} S\right) \tag{12.28}
\end{equation*}
$$
[[Ch12#^eq12-28|Eq. (12.28)]]

makes no difference to the resulting hedge ratios — canceling sensitivity to $S$ and vanilla prices is exactly equivalent to canceling sensitivity to $S$ and implied volatilities, exactly the market-model/sticky-strike delta equivalence already established for the local volatility model in Chapter 2. Delta- and vega-hedging removes the order-one $\delta O_i$ terms entirely, leaving the hedged position's P&L as

$$
\begin{align*}
P \& L_{\mathrm{H}}= & -\frac{1}{2} \Sigma_{i j} \frac{d^{2} P_{\mathrm{H}}}{d O_{i} d O_{j}}\left(\delta O_{i} \delta O_{j}-\widehat{a}_{i j} \delta t\right)  \tag{12.30a}\\
& -\Sigma_{k} \frac{d P_{\mathrm{H}}}{d \lambda_{k}}\left(\delta \lambda_{k}-\widehat{\mu}_{k} \delta t\right)  \tag{12.30b}\\
& -\frac{1}{2} \Sigma_{k l} \frac{d^{2} P_{\mathrm{H}}}{d \lambda_{k} d \lambda_{l}}\left(\delta \lambda_{k} \delta \lambda_{l}-\widehat{a}_{k l} \delta t\right)-\Sigma_{i k} \frac{d^{2} P_{\mathrm{H}}}{d O_{i} d \lambda_{k}}\left(\delta O_{i} \delta \lambda_{k}-\widehat{a}_{i k} \delta t\right) \tag{12.30c}
\end{align*}
$$
[[Ch12#^eq12-30a|Eq. (12.30a)]] [[Ch12#^eq12-30b|Eq. (12.30b)]] [[Ch12#^eq12-30c|Eq. (12.30c)]]

Because the mixed model is calibrated to hedge-instrument prices at every instant, the market observables never depend on the state variables ($dO_i/d\lambda_k=0$), which lets all $\lambda_k$-sensitivities of the hedged position be replaced with those of the unhedged position, giving the chapter's central technical result — the P&L splits cleanly into exactly two kinds of terms:

$$
\begin{align*}
P \& L_{\mathrm{H}}= & -\frac{1}{2} \Sigma_{i j} \frac{d^{2} P_{\mathrm{H}}}{d O_{i} d O_{j}}\left(\delta O_{i} \delta O_{j}-\widehat{a}_{i j} \delta t\right)  \tag{12.31a}\\
& -\Sigma_{k} \frac{d P}{d \lambda_{k}}\left(\delta \lambda_{k}-\widehat{\mu}_{k} \delta t\right)  \tag{12.31b}\\
& -\frac{1}{2} \Sigma_{k l} \frac{d^{2} P}{d \lambda_{k} d \lambda_{l}}\left(\delta \lambda_{k} \delta \lambda_{l}-\widehat{a}_{k l} \delta t\right)-\Sigma_{i k} \frac{d^{2} P}{d O_{i} d \lambda_{k}}\left(\delta O_{i} \delta \lambda_{k}-\widehat{a}_{i k} \delta t\right) \tag{12.31c}
\end{align*}
$$
[[Ch12#^eq12-31a|Eq. (12.31a)]] [[Ch12#^eq12-31b|Eq. (12.31b)]] [[Ch12#^eq12-31c|Eq. (12.31c)]]

Term (12.31a) is the genuine, expected gamma/theta P&L, with real financial break-even levels — this is what a legitimate market model is supposed to produce. Terms (12.31b) and (12.31c), by contrast, are unwanted P&L leakage driven purely by variations of the underlying stochastic-volatility model's own state variables — $V$ for a Heston-based mixed model, or the $\zeta^u$ curve for the two-factor case — variables that carry no financial significance whatsoever to a trader holding a delta- and vega-hedged position. These leakage terms were entirely absent from the local volatility model's own P&L decomposition in Chapter 2; they are a genuinely new pathology specific to mixed models.

Section 12.3.3 isolates exactly when this leakage vanishes. The condition is simply that the price must not depend on the state variables $\lambda_k$ at all,

$$
\begin{equation*}
\left.\frac{d P}{d \lambda_{k}}\right|_{S, \widehat{\sigma}_{K T}}=0, \forall k \tag{12.32}
\end{equation*}
$$
[[Ch12#^eq12-32|Eq. (12.32)]]

in which case the P&L collapses to the clean, purely observable-driven form

$$
\begin{equation*}
P \& L_{\mathrm{H}}=-\frac{1}{2} \Sigma_{i j} \frac{d^{2} P_{\mathrm{H}}}{d O_{i} d O_{j}}\left(\delta O_{i} \delta O_{j}-\widehat{a}_{i j} \delta t\right) \tag{12.33}
\end{equation*}
$$
[[Ch12#^eq12-33|Eq. (12.33)]]

Remarkably, condition (12.32) holds for the two-factor model: rescaling $\zeta_0^u\to\varphi^u\zeta_0^u$ and simultaneously rescaling $\sigma(u,S)\to\sigma(u,S)/\sqrt{\varphi^u}$ leaves the process for $S_t$ and its instantaneous volatility completely unchanged for arbitrary constants $\varphi^u$ — meaning price cannot depend on $\zeta^u$ at all, since it can be rescaled away without changing anything observable. The same holds whenever the underlying stochastic-volatility model's instantaneous volatility is itself lognormal, exactly the SABR-type case,

$$
\left\{\begin{align*}
d S_{t} & =(r-q) S_{t} d t+\sigma\left(t, S_{t}\right) \lambda_{t} S_{t} d W_{t}^{S}  \tag{12.35}\\
d \lambda_{t} & =\nu \lambda_{t} d W_{t}^{\lambda}
\end{align*}\right.
$$
[[Ch12#^eq12-35|Eq. (12.35)]]

But condition (12.32) fails outright for the Heston model, since price genuinely depends on $V$ directly, and it fails for the actively-used "Bloomberg model" too,

$$
\left\{\begin{align*}
d S_{t} & =(r-q) S_{t} d t+\sigma\left(t, S_{t}\right) \lambda_{t} S_{t} d W_{t}^{S}  \tag{12.34}\\
d \lambda_{t} & =-k\left(\lambda_{t}-\theta\right) d t+\xi(t) \lambda_{t} d W_{t}^{\lambda}
\end{align*}\right.
$$
[[Ch12#^eq12-34|Eq. (12.34)]]

Two models in active real-world trading use, by this analysis, are shown to be — strictly, by the book's own Chapter-1 standard of what counts as a genuine price — not usable trading models at all, generating exactly the leakage terms (12.31b)-(12.31c) rather than a clean gamma/theta breakdown.

This is a genuinely uncomfortable conclusion for a desk that has already built pricing and risk infrastructure around one of these models, and the chapter does not soften it: the leakage is not a numerical artifact that better calibration or finer time-stepping can fix, but a structural consequence of the model's own state variables carrying real price sensitivity that a delta/vega hedge cannot see or offset. A desk using such a model is, in effect, carrying an unhedged and unacknowledged risk every day, whether or not it shows up in the daily P&L explain.

Sections 12.4-12.6 turn to the admissible, two-factor-based mixed model's practical break-even levels. Since $dP/d\zeta^\tau=0$ for the admissible class, the OU factors can be normalized to $\zeta^\tau\equiv1$, giving a spot SDE

$$
\begin{equation*}
d S_{t}=(r-q) S_{t} d t+\sigma\left(t, S_{t}\right) \sqrt{f\left(t, X_{t}^{1}, X_{t}^{2}\right)} S_{t} d W_{t}^{S} \tag{12.36}
\end{equation*}
$$
[[Ch12#^eq12-36|Eq. (12.36)]]

Expanding the local volatility function around a time-dependent baseline $\bar\sigma(t)$ and the vol-of-vol factor to first order in $\nu$ simultaneously,

$$
\begin{equation*}
d S_{t}=(r-q) S_{t} d t+\left(\bar{\sigma}(t)+\delta \sigma\left(t, S_{t}\right)\right)\left(1+\frac{\nu}{2} g\left(t, X_{t}^{1}, X_{t}^{2}\right)\right) S_{t} d W_{t}^{S} \tag{12.37}
\end{equation*}
$$
[[Ch12#^eq12-37|Eq. (12.37)]]

decomposes the spot process itself into a baseline plus two independent order-one corrections — one from the local-volatility component, one from the stochastic-volatility component,

$$
\begin{equation*}
S_{t}=S_{t}^{0}+\delta S_{t}^{\mathrm{LV}}+\delta S_{t}^{\mathrm{SV}} \tag{12.38}
\end{equation*}
$$
[[Ch12#^eq12-38|Eq. (12.38)]]

which propagates directly into an identical split of option prices,

$$
\begin{equation*}
P=P_{0}+\delta P_{\mathrm{LV}}+\delta P_{\mathrm{SV}} \tag{12.39}
\end{equation*}
$$
[[Ch12#^eq12-39|Eq. (12.39)]]

and hence of implied volatility and the ATMF skew itself,

$$
\begin{align*}
\widehat{\sigma}_{K T} & =\widehat{\sigma}_{T}+\delta \widehat{\sigma}_{K T}^{\mathrm{LV}}+\delta \widehat{\sigma}_{K T}^{\mathrm{SV}}  \tag{12.40a}\\
\mathcal{S}_{T} & =\mathcal{S}_{T}^{\mathrm{LV}}+\mathcal{S}_{T}^{\mathrm{SV}} \tag{12.40b}
\end{align*}
$$
[[Ch12#^eq12-40a|Eq. (12.40a)]] [[Ch12#^eq12-40b|Eq. (12.40b)]]

Each piece was already computed, in isolation, in earlier chapters — the local-volatility contribution recycles Chapter 2's skew-averaging identity,

$$
\begin{equation*}
\mathcal{S}_{T}^{\mathrm{LV}}=\frac{1}{T} \int_{0}^{T} \frac{\widehat{\sigma}_{t}^{2} t}{\widehat{\sigma}_{T}^{2} T} \frac{\bar{\sigma}(t)}{\widehat{\sigma}_{T}} \alpha(t) d t \tag{12.41}
\end{equation*}
$$
[[Ch12#^eq12-41|Eq. (12.41)]]

which, expanding instead around a constant baseline volatility, simplifies to

$$
\begin{equation*}
\mathcal{S}_{T}^{\mathrm{LV}}=\frac{1}{T} \int_{0}^{T} \frac{t}{T} \alpha(t) d t \tag{12.42}
\end{equation*}
$$
[[Ch12#^eq12-42|Eq. (12.42)]]

while the stochastic-volatility contribution directly recycles Chapter 8's two-factor ATMF skew formula, applied to the initial variance curve $\xi_0^\tau=\bar\sigma^2(\tau)$,

$$
\begin{equation*}
\mathcal{S}_{T}^{\mathrm{SV}}=\frac{\nu \alpha_{\theta}}{\widehat{\sigma}_{T}^{3} T^{2}} \int_{0}^{T} d t \bar{\sigma}(t) \int_{t}^{T} d u \bar{\sigma}^{2}(u)\left[(1-\theta) \rho_{S X^{1}} e^{-k_{1}(u-t)}+\theta \rho_{S X^{2}} e^{-k_{2}(u-t)}\right] \tag{12.43}
\end{equation*}
$$
[[Ch12#^eq12-43|Eq. (12.43)]]

which for a constant baseline reduces to the familiar closed form,

$$
\begin{equation*}
\mathcal{S}_{T}^{\mathrm{SV}}=\nu \alpha_{\theta}\left[(1-\theta) \rho_{S X^{1}} \frac{k_{1} T-\left(1-e^{-k_{1} T}\right)}{\left(k_{1} T\right)^{2}}+\theta \rho_{S X^{2}} \frac{k_{2} T-\left(1-e^{-k_{2} T}\right)}{\left(k_{2} T\right)^{2}}\right] \tag{12.44}
\end{equation*}
$$
[[Ch12#^eq12-44|Eq. (12.44)]]

Summing both pieces gives the mixed model's combined ATMF skew, order one simultaneously in the local volatility component and in volatility of volatility,

$$
\begin{align*}
\mathcal{S}_{T}= & \frac{1}{T} \int_{0}^{T} \frac{\widehat{\sigma}_{t}^{2} t}{\widehat{\sigma}_{T}^{2} T} \frac{\bar{\sigma}(t)}{\widehat{\sigma}_{T}} \alpha(t) d t \\
& +\frac{\nu \alpha_{\theta}}{\widehat{\sigma}_{T}^{3} T^{2}} \int_{0}^{T} d t \bar{\sigma}(t) \int_{t}^{T} d u \bar{\sigma}^{2}(u)\left[(1-\theta) \rho_{S X^{1}} e^{-k_{1}(u-t)}+\theta \rho_{S X^{2}} e^{-k_{2}(u-t)}\right] \tag{12.45}
\end{align*}
$$
[[Ch12#^eq12-45|Eq. (12.45)]]

This total skew formula is, by itself, of limited direct use — the model is calibrated to market, so $\mathcal S_T$ can simply be read off the smile directly. Its real value is that it isolates $\mathcal S_T^{\text{SV}}$ specifically, which is what the rest of the section needs.

Section 12.4.2 turns to the ATMF volatility's own dynamics. The local-volatility component's contribution reuses Chapter 2's own SSR formula for that component,

$$
\begin{equation*}
\mathcal{R}_{T}^{\mathrm{LV}}=1+\frac{1}{T} \int_{0}^{T} \frac{\bar{\sigma}^{2}(t)}{\widehat{\sigma}_{t} \widehat{\sigma}_{T}} \frac{\mathcal{S}_{t}^{\mathrm{LV}}}{\mathcal{S}_{T}^{\mathrm{LV}}} d t \tag{12.46}
\end{equation*}
$$
[[Ch12#^eq12-46|Eq. (12.46)]]

while the stochastic-volatility component's contribution reuses Chapter 7's direct VS-volatility SDE, with coefficients

$$
\begin{equation*}
A_{i}=\frac{\int_{0}^{T} \xi_{0}^{\tau} e^{-k_{i} \tau} d \tau}{\int_{0}^{T} \xi_{0}^{\tau} d \tau}=\frac{\int_{0}^{T} \bar{\sigma}^{2}(\tau) e^{-k_{i} \tau} d \tau}{\int_{0}^{T} \bar{\sigma}^{2}(\tau) d \tau} \tag{12.47}
\end{equation*}
$$
[[Ch12#^eq12-47|Eq. (12.47)]]

Since the local-volatility component's own skew $\mathcal S_t^{\text{LV}}$ is not known in closed form, it is recovered directly from (12.40b) as $\mathcal S_t-\mathcal S_t^{\text{SV}}$, giving the combined ATMF-volatility dynamics entirely in terms of quantities already computed,

$$
\begin{align*}
d \widehat{\sigma}_{T}= & \left(\left(\mathcal{S}_{T}-\mathcal{S}_{T}^{\mathrm{SV}}\right)+\frac{1}{T} \int_{0}^{T} \frac{\bar{\sigma}^{2}(t)}{\widehat{\sigma}_{t} \widehat{\sigma}_{T}}\left(\mathcal{S}_{t}-\mathcal{S}_{t}^{\mathrm{SV}}\right) d t\right) d \ln S \\
& +\nu \alpha_{\theta} \widehat{\sigma}_{T}\left((1-\theta) A_{1} d W^{1}+\theta A_{2} d W^{2}\right) \tag{12.48}
\end{align*}
$$
[[Ch12#^eq12-48|Eq. (12.48)]]

which for a constant baseline simplifies to

$$
\begin{align*}
d \widehat{\sigma}_{T}= & \left(\left(\mathcal{S}_{T}-\mathcal{S}_{T}^{\mathrm{SV}}\right)+\frac{1}{T} \int_{0}^{T}\left(\mathcal{S}_{t}-\mathcal{S}_{t}^{\mathrm{SV}}\right) d t\right) d \ln S \\
& +\nu \alpha_{\theta} \widehat{\sigma}_{T}\left((1-\theta) \frac{1-e^{-k_{1} T}}{k_{1} T} d W^{1}+\theta \frac{1-e^{-k_{2} T}}{k_{2} T} d W^{2}\right) \tag{12.49}
\end{align*}
$$
[[Ch12#^eq12-49|Eq. (12.49)]]

Section 12.4.2.1 converts this directly into the mixed model's SSR. Starting from the ordinary SSR definition,

$$
\begin{equation*}
\mathcal{R}_{T}=\frac{1}{\mathcal{S}_{T}} \frac{\left\langle d \widehat{\sigma}_{T} d \ln S\right\rangle}{\left\langle(d \ln S)^{2}\right\rangle} \tag{12.50}
\end{equation*}
$$
[[Ch12#^eq12-50|Eq. (12.50)]]

writing $d\widehat\sigma_T$ as the sum of its LV and SV pieces gives a skew-weighted average of the two components' own SSRs,

$$
\begin{equation*}
\mathcal{R}_{T}=\frac{\mathcal{R}_{T}^{\mathrm{LV}} \mathcal{S}_{T}^{\mathrm{LV}}+\mathcal{R}_{T}^{\mathrm{SV}} \mathcal{S}_{T}^{\mathrm{SV}}}{\mathcal{S}_{T}^{\mathrm{LV}}+\mathcal{S}_{T}^{\mathrm{SV}}} \tag{12.51}
\end{equation*}
$$
[[Ch12#^eq12-51|Eq. (12.51)]]

which, substituting $\mathcal S_t^{\text{LV}}=\mathcal S_t-\mathcal S_t^{\text{SV}}$ throughout, resolves into the chapter's final working SSR formula for the mixed model,

$$
\begin{equation*}
\mathcal{R}_{T}=\mathcal{R}_{T}^{\mathrm{LV}}(\mathrm{Mkt})+\frac{\mathcal{S}_{T}^{\mathrm{SV}}}{\mathcal{S}_{T}}\left[\mathcal{R}_{T}^{\mathrm{SV}}-\mathcal{R}_{T}^{\mathrm{LV}}(\mathrm{SV})\right] \tag{12.52}
\end{equation*}
$$
[[Ch12#^eq12-52|Eq. (12.52)]]

where $\mathcal R_T^{\text{LV}}(\text{Mkt})$ and $\mathcal R_T^{\text{LV}}(\text{SV})$ are the local volatility model's own SSR formula evaluated with the *market* skew and with the stochastic-volatility component's skew, respectively,

$$
\begin{gather*}
\mathcal{R}_{T}^{\mathrm{LV}}(\mathrm{Mkt})=1+\frac{1}{T} \int_{0}^{T} \frac{\mathcal{S}_{t}}{\mathcal{S}_{T}} d t  \tag{12.53}\\
\mathcal{R}_{T}^{\mathrm{LV}}(\mathrm{SV})=1+\frac{1}{T} \int_{0}^{T} \frac{\mathcal{S}_{t}^{\mathrm{SV}}}{\mathcal{S}_{T}^{\mathrm{SV}}} d t  \tag{12.54}\\
\mathcal{S}_{T}^{\mathrm{SV}}=\nu \alpha_{\theta}\left[(1-\theta) \rho_{S X^{1}} \frac{k_{1} T-\left(1-e^{-k_{1} T}\right)}{\left(k_{1} T\right)^{2}}+\theta \rho_{S X^{2}} \frac{k_{2} T-\left(1-e^{-k_{2} T}\right)}{\left(k_{2} T\right)^{2}}\right]
\end{gather*}
$$
[[Ch12#^eq12-53|Eq. (12.53)]] [[Ch12#^eq12-54|Eq. (12.54)]]

and $\mathcal R_T^{\text{SV}}$ is exactly Chapter 9's own two-factor SSR formula,

$$
\begin{equation*}
\mathcal{R}_{T}^{\mathrm{SV}}=\frac{(1-\theta) \rho_{S X^{1}} \frac{1-e^{-k_{1} T}}{k_{1} T}+\theta \rho_{S X^{2}} \frac{1-e^{-k_{2} T}}{k_{2} T}}{(1-\theta) \rho_{S X^{1}} \frac{k_{1} T-\left(1-e^{-k_{1} T}\right)}{\left(k_{1} T\right)^{2}}+\theta \rho_{S X^{2}} \frac{k_{2} T-\left(1-e^{-k_{2} T}\right)}{\left(k_{2} T\right)^{2}}} \tag{12.55}
\end{equation*}
$$
[[Ch12#^eq12-55|Eq. (12.55)]]

Formula (12.52) is genuinely useful precisely because it exposes a tunable knob: since the market's overall skew $\mathcal S_T$ is fixed by calibration, the *split* of that skew between $\mathcal S_T^{\text{LV}}$ and $\mathcal S_T^{\text{SV}}$ — controlled by how the underlying two-factor model's parameters are chosen — directly moves the resulting SSR between the local-volatility limit and the stochastic-volatility limit, without changing the calibrated vanilla smile at all. The exact numerical value can also be obtained directly, without any of these approximations, by a first-order expansion of the Monte-Carlo-computed ATMF volatility in the state variables,

$$
\begin{equation*}
d \widehat{\sigma}_{T}=\frac{d \widehat{\sigma}_{T}}{d \ln S} d \ln S+\frac{d \widehat{\sigma}_{T}}{d X^{1}} d X^{1}+\frac{d \widehat{\sigma}_{T}}{d X^{2}} d X^{2} \tag{12.56}
\end{equation*}
$$
[[Ch12#^eq12-56|Eq. (12.56)]]

evaluated via one bump-and-reprice on the simulation grid, exactly the same numerical technique Chapter 9 used for the pure two-factor model, with no recalibration required since the break-even levels are those the model already generates with a fixed local volatility function.

Section 12.6 resolves the puzzle this decomposition exists to explain: numerically, mixed models calibrated to the *exact same* vanilla smile — but built on differently-parametrized underlying stochastic-volatility components — can generate markedly different SSRs and future skews, even though their static smile is identical. How can this be reconciled? The answer traces back to Chapter 8's general covariance-based skew identity,

$$
\begin{equation*}
\mathcal{S}_{T}=\frac{1}{2 \widehat{\sigma}_{T}^{3}(0) T} \int_{0}^{T} \frac{T-t}{T} \frac{\left\langle d \ln S_{t} d \widehat{\sigma}_{F_{T} T}^{2}(t)\right\rangle}{d t} d t \tag{12.57}
\end{equation*}
$$
[[Ch12#^eq12-57|Eq. (12.57)]]

which, at order one in both the local-volatility component and volatility of volatility, rewrites as

$$
\begin{equation*}
\mathcal{S}_{T}=\frac{1}{\widehat{\sigma}_{T}^{3}(0) T} \int_{0}^{T} \frac{T-t}{T} \widehat{\sigma}_{T}(t) \frac{\left\langle d \ln S_{t} d \widehat{\sigma}_{F_{T} T}(t)\right\rangle}{d t} d t \tag{12.58}
\end{equation*}
$$
[[Ch12#^eq12-58|Eq. (12.58)]]

The resolution is now direct: the *static* smile only pins down the *integrated*, weighted-average spot/volatility covariance appearing in (12.58) — it says nothing at all about how that total covariance is distributed, moment by moment, between the local-volatility mechanism and the stochastic-volatility mechanism. That distribution is exactly the extra degree of freedom mixed models were built to control, and it is precisely what a fixed vanilla smile cannot determine on its own.

This is worth restating plainly, since it resolves what could otherwise look like a contradiction: two mixed models, calibrated to literally the same option prices and therefore indistinguishable on any static risk report, can nonetheless generate opposite scenarios for how that smile evolves as the market moves. A risk manager comparing two desks' otherwise-identical vanilla books could easily miss this, since nothing in the static Greeks would flag it — only a genuinely dynamic stress test, or a direct SSR/future-skew comparison, would reveal the difference.

This is precisely why the whole apparatus of Chapters 8 and 9 — the covariance-based skew formula, the SSR, the static/dynamic classification — earns its place as more than academic machinery: it is the only lens fine enough to distinguish two models that a purely static vanilla comparison would otherwise judge identical.

The chapter's closing numerical illustration makes this control concrete and practically consequential. Pricing a strip of forward-start call spreads across the admissible mixed models shows that far-forward skew is generated almost entirely by the stochastic-volatility component — the local-volatility component's own contribution to forward skew dies off quickly with the forward-start date, exactly as one would expect from a mechanism ultimately tied to $(t,S)$-dependence alone. The practical upshot is unambiguous: it is the choice of underlying stochastic-volatility parametrization, not the calibrated local volatility function, that actually controls a mixed model's forward-smile risk — the local volatility component supplies the exact vanilla calibration, but the genuine dynamical content the whole chapter set out to add comes entirely from the stochastic-volatility half of the model.

#### Chapter's Digest and Logical Bridge

Chapter 12 closes the book by showing that combining local and stochastic volatility is not automatic: the "seemingly innocuous" substitution (12.5) does not by itself guarantee a genuine price, and most mixed models in active use — Heston-based and Bloomberg-model mixed constructions alike — fail the exact usability condition (12.32) and leak spurious P&L with no financial meaning. Only a specific, identifiable class — built on the two-factor model, or on any lognormal-instantaneous-volatility construction like SABR — passes this test, and for that admissible class the chapter delivers exactly the tool the rest of the book was building toward: a calibrated model that also lets a desk dial in the SSR and future-skew behavior independently of the static smile, resolving the apparent paradox that identical vanilla smiles can hide entirely different dynamics. The book's own conclusion is a deliberate scope boundary, not a triumphant final word: admissible mixed models are legitimate, genuinely usable trading models with real, if partial, control over break-even levels — but Chapter 7's forward variance models remain necessary whenever direct access to forward variances or VIX futures is required, since a mixed model's local-volatility component only ever depends on $(t,S)$, never on the variance curve itself. The two families are not competitors; they are tools for different jobs, and knowing which job is at hand is the practical judgment this entire book has been building toward.

---
