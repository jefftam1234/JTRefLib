---
tags: [quant/knowledge, quant/books, interest-rates]
---

## Callable Libor Exotics

Having discussed relatively simple vanilla securities in Chapters 16 and 17, we now move to the other extreme of the product spectrum and consider issues of valuation, calibration and risk management of callable Libor exotics (CLEs). CLEs were defined in Section 5.14 and constitute the most complicated class of interest rate derivatives traded in the market. The material in this chapter deals with CLEs in general; later chapters will take a more in-depth look at the idiosyncrasies of some of the most popular and/or challenging securities inside the CLE class.

### 18.1 Model Calibration for Callable Libor Exotics

Due to their inherent complexity, CLEs will have non-trivial dependencies on the dynamics of market rates and will require sophisticated term structure models for valuation and risk management. As discussed in earlier chapters, computation of non-vanilla prices in such models virtually always requires Monte Carlo simulation, which in turn introduces complications in determining the optimal exercise rules for CLEs. We shall spend a good portion of this chapter on a detailed discussion of algorithms for this particular problem. However, before a model can even be considered for CLE valuation, it needs to be calibrated, i.e. the volatility structure of the model needs to be parameterized to match available market information relevant for valuation.

At this point we should emphasize that CLEs and other exotic interest rate derivatives are different from many of vanilla securities considered previously in that their values are not observable in the market, a consequence of the fact that most exotic derivatives are sold to clients rather than traded between dealers. This fact leads to fundamental differences in the way models are used for valuation of vanilla vs. exotic derivatives. Whereas for many vanilla derivatives a model is primarily used to generate a hedging strategy and to perform inter- and extrapolation on market-observable prices, for
exotics the values are fundamentally derived from a model ${ }^{1}$. Given the absence of market prices of exotics themselves, models for these securities have to be calibrated indirectly, to other market (and perhaps also nonmarket) information that is deemed relevant for the class of exotics under consideration. Speaking very loosely, the purpose of a model for CLEs can therefore be characterized as performing a sophisticated extrapolation of information from a series of "spanning" vanilla markets to compute a meaningful exotics price, something we already alluded to in Section 14.5.5.

### 18.1.1 Risk Factors for CLEs

In order to make sure that a CLE model captures as much relevant market information as possible, it is important to understand what we can actually rely on when calibrating a model. There are essentially three sources of potentially relevant information. The first, and arguably the most important, source is the market prices of liquid vanilla interest rate derivatives. The second source is historical information about market quantities such as volatilities and correlations of various market rates that may not be available for observation directly from quoted market prices. The third, and somewhat more amorphous, source is the modeler's beliefs on what. constitutes reasonable behavior of the model parameters. The last category, for example, includes views on how time-stationary or how smooth model parameters should be.

To perform model calibration, we typically choose particular targets from all three sources above, and ultimately set model parameters in a way to match those targets in our model. In order to be able to choose relevant targets, it is important to identify risk factors that affect the valuation of the particular CLE in question. This part of the analysis is necessarily product specific, so let us look at a particular example to demonstrate some salient points. We choose a reasonably representative (albeit not very complicated) example, namely a callable inverse floater, or CIF, with exercise dates

$$
T_{1}<T_{2}<\ldots<T_{N-1}
$$

and structured coupon

$$
C_{n}=\min \left(\max \left(6 \%-L_{n}\left(T_{n}\right), 0 \%\right), 4 \%\right),
$$

where $L_{n}(t)$ is a forward Libor rate fixing at $T_{n}$. In the language of Section 5.14 , we see that the coupon strike is equal to $6 \%$, the cap is at $4 \%$ and the floor is at $0 \%$. It should be clear that the coupon can be decomposed into a portfolio of a long floorlet with strike $6 \%$ and a short floorlet with strike $2 \%$. Hence, the callable inverse floater can be thought of as a Bermudan-style option on a combination of floors and a Libor leg.

[^25]The underlying swap of this CIF consists of a collection of floorlets, which are vanilla derivatives with observable market prices. Market-implied volatilities ${ }^{2}$ as always serve as a convenient representation of these market values. For classification purposes that will be clear in a moment, we call such volatilities spot volatilities, as we can observe them on the valuation date, or "on the spot", from the market. Clearly, these volatilities affect the value of the exotic swap that underlies the CIF and hence have a direct impact on the value of the CIF. For that reason we should include them as targets for model calibration. Importantly, as caplet volatilities virtually always exhibit a volatility smile, volatilities for both $2 \%$ and $6 \%$ strikes are potentially relevant.

While the underlying exotic swap has no dependencies on spot volatilities of market rates other than those underlying the floorlets inherent in the CIF, the callability structure of the CIF nevertheless introduces additional dependence on other vanilla derivatives. To see this, suppose that we, say, ignored all exercise CIF dates but one, in which case the CIF would degenerate into a European-style option to enter the underlying swap. Even though the underlying swap is not vanilla (i.e. not a standard fixed-for-floating rate swap), it is clearly related to one and our single-call CIF is therefore related to the European swaption on that swap. More generally, the (multi-call) CIF will depend on the (spot) implied volatility of the swap rate that fixes on $T_{\imath}$, and runs for the period [ $T_{\imath}, T_{N}$ ], for each $i=1, \ldots, N-1$. It is less clear what strike we should use to define this implied volatility; while all strikes are in fact relevant, as shall be clear later, sometimes one may choose to simplify and just use the ATM swaption volatility.

To summarize our discussion so far, we have argued that the value of the CIF depends on market-implied spot Libor rate volatilities for all expiries until $T_{N-1}$ (of two different strikes), and spot swap rate volatilities for those swaptions for which expiry plus tenor is equal to $T_{N}$ - the so called core, diagonal, or coterminal swaptions. When calibrating a model to price this CIF, we would seek (as a minimum) to calibrate it to these market volatilities.

While there are no more obvious spot volatility targets that should be included in CIF calibration, there are other volatility-related quantities that will affect the valuation of the CIF and other similar CLEs. Let us imagine that we have arrived at exercise date $T_{n}, n<N-1$, and assume that the CIF has not been exercised in the past. At this point in time, the option holder needs to decide whether to exercise and receive the underlying value $U_{n}\left(T_{n}\right)$ (see (5.24)) or whether to continue to hold on to the CIF (with a view of perhaps exercising it later). The time $T_{n}$ value of the remaining part of the underlying swap depends on caplet volatilities as observed at time $T_{n}$. As these volatilities will be known only at time $T_{n}$, we call them forward volatilities of Libor rates. Likewise, the option to hold a Bermudan-style

[^26]option on the underlying $U_{n+1}\left(T_{n}\right)$, will depend on forward volatilities of swap rates, specifically core swaption volatilities observed at time $T_{n}$. As the exercise decision at time $T_{n}$ evidently will depend on the forward volatility structure at time $T_{n}$, the time 0 value of the CIF will depend on it as well.

At time 0 , forward volatilities are generally unknown random variables, but any model will impose certain dynamics on the volatility structure of interest rates, which will have a direct impact on the model value of the CIF (and other CLEs). Ideally, we should make sure that the model projections for forward volatilities are in line with market-implied information. The fact that such information is typically either not available or difficult to extract adds significantly to the challenge of model calibration for CLEs. Frequently, it will be necessary to lean on historical data and to impose exogenous assumptions that the model builder might feel are financially reasonable.

In Section 13.1.8.1 we identified inter-temporal correlations - i.e. correlations of core swap rates observed on different fixing dates - as value drivers for Bermudan swaptions. Inter-temporal correlations are closely linked to the forward volatilities discussed above, and the two concepts can often be used interchangeably. To demonstrate, note that the (inter-temporal) correlation between two core swap rates, $S_{n, N-n}\left(T_{n}\right)$ and $S_{m, N-m}\left(T_{m}\right)(n<m)$ is given by

$$
\begin{array}{r}
\operatorname{Corr}\left(S_{n, N-n}\left(T_{n}\right), S_{m, N-m}\left(T_{m}\right)\right)=\operatorname{Cov}\left(S_{n, N-n}\left(T_{n}\right), S_{m, N-m}\left(T_{m}\right)\right) \\
\times \operatorname{Var}\left(S_{n, N-n}\left(T_{n}\right)\right)^{-1 / 2} \operatorname{Var}\left(S_{m, N-m}\left(T_{m}\right)\right)^{-1 / 2}
\end{array}
$$

As $S_{m, N-m}\left(T_{m}\right)-S_{m, N-m}\left(T_{n}\right)$ is often only mildly dependent on $S_{m, N-m}\left(T_{n}\right)$ and $S_{n, N-n}\left(T_{n}\right)$, we can rewrite this as

$$
\begin{aligned}
\operatorname{Corr} & \left(S_{n, N-n}\left(T_{n}\right), S_{m, N-m}\left(T_{m}\right)\right) \\
\approx & \operatorname{Cov}\left(S_{n, N-n}\left(T_{n}\right), S_{m, N-m}\left(T_{m}\right)\right) \operatorname{Var}\left(S_{n, N-n}\left(T_{n}\right)\right)^{-1 / 2} \\
& \times\left(\operatorname{Var}\left(S_{m, N-m}\left(T_{n}\right)\right)+\operatorname{Var}\left(S_{m, N-m}\left(T_{m}\right)-S_{m, N-m}\left(T_{n}\right)\right)\right)^{-1 / 2} \\
\approx & \operatorname{Corr}\left(S_{n, N-n}\left(T_{n}\right), S_{m, N-m}\left(T_{n}\right)\right) \\
& \times\left(1+\frac{\operatorname{Var}\left(S_{m, N-m}\left(T_{m}\right)-S_{m, N-m}\left(T_{n}\right)\right)}{\operatorname{Var}\left(S_{m, N-m}\left(T_{n}\right)\right)}\right)^{-1 / 2}
\end{aligned}
$$

and we see that, with (spot) correlation $\operatorname{Corr}\left(S_{n, N-n}\left(T_{n}\right), S_{m, N-m}\left(T_{n}\right)\right)$ and (spot) variance $\operatorname{Var}\left(S_{m, N-m}\left(T_{n}\right)\right)$ fixed, the inter-temporal correlation is directly linked to the forward variance (square of forward volatility) $\operatorname{Var}\left(S_{m, N-m}\left(T_{m}\right)-S_{m, N-m}\left(T_{n}\right)\right)$. This equivalence will play a role when we discuss the local projection method in Section 18.4.

Another aspect of the model behavior worth emphasizing here is the volatility smile dynamics imposed by the model. In Sections 8.8 and 16.1.1 we considered the hedging impact of joint moves in rates and volatility smiles. For CLEs, volatility smile dynamics affect not only the hedging strategy, but
also the valuation of CLEs, due to the forward volatility effect we discussed earlier. We notice that for many CLEs the effective coupon strikes at the exercise boundary often end up being deeply in- or out-of-the-money, since interest rate levels at which option exercise is optimal are usually significantly different from the levels of rates at inception of the trade.

### 18.1.2 Model Choice and Calibration

As argued above, the value of a typical CLE depends strongly on volatility smile dynamics and on market (spot) implied volatilities for a wide selection of vanilla options, often across a range of strikes. As such the appropriate choice of a model for CLE valuation should typically involve the following criteria.

- Ability to calibrate to a large collection of vanilla options across expiries, tenors, and strikes.
- Reasonable and controllable dynamics of the volatility structure.
- Multi-factor interest rate dynamics, especially for CLEs on multi-rate underlyings.

A combination of these requirements often rules out simpler, lowdimensional models, especially for more complicated CLEs. We typically recommend using either Libor market, multi-factor quasi-Gaussian, or multifactor quadratic Gaussian models. While simpler low-dimensional models might succeed in fitting spot volatility information for selected Libor and swap rates, they often achieve this only by using values of model parameters that would imply unrealistic evolution of the volatility structure and, therefore, unreasonable pricing of CLEs (but see Section 18.4).

Another observation favors models that can calibrate to a large collection of European swaptions. A swap rate can be seen to be, approximately, a weighted average of Libor rates (see Section 14.4.2). Hence, an implied European swaption volatility contains some amount of information on marketconsistent correlations between Libor (and other swap) rates. Extracting this information is complicated, but in a sense is what a Libor market model calibrated to, say, the whole swaption grid (including caps) is designed to do. In the same spirit, a model calibrated to all swaptions and caps can be thought of as giving us the best available implied forward volatility structure, i.e. market-consistent information of what the most likely behavior of the volatility structure through time would be.

Another dimension to consider here is the volatility smile dynamics. Obviously, for the model to generate unspanned stochasticity in volatility movements, stochastic volatility (SV) version of the model would be required. Given the discussion in Section 8.8 and the importance for CLE pricing of controlling the dynamics of the volatility smile, SV models typically have a clear edge over their local volatility counterparts in CLE pricing applications.

With the choice of the model settled, calibration typically proceeds by identifying proper targets and fitting them with model parameters. For Libor market models we, in fact, have already covered relevant issues quite extensively, and the reader is advised to revisit Section 14.5.5. For other types of multi-factor models, the issues are similar and can be resolved in the same spirit. The mechanics of calibration for relevant models are covered in their respective chapters. For example, smile calibration of LM models is discussed in Section 15.2.

### 18.2 Valuation Theory

### 18.2.1 Preliminaries

While various measures could be used for CLE valuation, for concreteness we choose to use the spot Libor measure $\mathrm{Q}^{B}$ with the discrete money market numeraire $B(t)$ defined on a tenor structure

$$
0=T_{0}<T_{1}<\ldots<T_{N}, \quad \tau_{n}=T_{n+1}-T_{n} .
$$

For notational simplicity, we let $\mathrm{E}(\cdot)$ be the expected value operator in measure $\mathrm{Q}^{B}$.

We recall from Section 5.14 that a callable Libor exotic is a Bermudanstyle option on an exotic swap that specifies an exchange of structured coupons $C_{n}$ for Libor rates $L_{n}$, fixing on $T_{n}$ and paying on $T_{n+1}, n= 1, \ldots, N-1$. We denote by $X_{n}$ the net payment seen by the structured coupon receiver,

$$
\begin{equation*}
X_{n}=\tau_{n} \times\left(C_{n}-L_{n}\left(T_{n}\right)\right) \tag{18.1}
\end{equation*}
$$
^eq18-1

Also, we let $U_{n}(t)$ be the $n$-th exercise value, i.e. the value of all future payments if the callable Libor exotic is exercised at time $T_{n}$. Clearly

$$
\begin{equation*}
U_{n}(t)=B(t) \sum_{i=n}^{N-1} \mathrm{E}_{t}\left(B\left(T_{i+1}\right)^{-1} X_{i}\right) . \tag{18.2}
\end{equation*}
$$
^eq18-2

For completeness we set

$$
U_{N}(t) \equiv 0
$$

If a callable Libor exotic is exercised on $T_{n}$, the holder will receive $U_{n}\left(T_{n}\right)$, the present value of the remainder of the underlying exotic swap.

Finally, let $H_{n}(t)$ be the value at time $t$ of a callable Libor exotic where exercise opportunities have been restricted to the dates $\left\{T_{n+1}, \ldots, T_{N-1}\right\}$. In particular, $H_{0}(0)$ is then the time 0 value of the CLE. Each $H_{n}$ is called a hold value, since $H_{n}\left(T_{n}\right)$ is the value of the choice of not exercising on date $T_{n}$, i.e. continuing to "hold" the derivative. We must necessarily have

$$
\begin{equation*}
H_{0}(t) \geq H_{1}(t) \geq \ldots \geq H_{N-2}(t) \tag{18.3}
\end{equation*}
$$
^eq18-3

### 18.2.2 Recursion for Callable Libor Exotics

Let $\mathcal{T}_{n}$ be a set of all stopping time indices that exceed $n, n \in\{0, \ldots, N-1\}$, i.e. a set of random variables taking values in the set $\{n+1, \ldots, N\}$ such that for any $k$ and any $\xi \in \mathcal{T}_{n}$,

$$
\{\xi=k\} \in \mathcal{F}_{T_{h}}
$$

The sequence of random variables $H_{0}\left(T_{0}\right), \ldots, H_{N-1}\left(T_{N-1}\right)$ (as in Section 18.2.1) defines the Snell envelope for the sequence of (discounted) exercise values,

$$
\begin{equation*}
H_{n}\left(T_{n}\right)=B\left(T_{n}\right) \sup _{\xi \in \mathcal{T}_{n}} \mathrm{E}_{T_{n}}\left(B\left(T_{\xi}\right)^{-1} U_{\xi}\left(T_{\xi}\right)\right), \quad n=0, \ldots, N-1 \tag{18.4}
\end{equation*}
$$
^eq18-4

By the general theory of optimal stopping (see Chapter 1.10), the random time index $\eta_{n}$ that maximizes the right-hand side of (18.4) is given by

$$
\begin{equation*}
\eta_{n}(\omega)=\min \left\{k>n: U_{k}\left(T_{k}\right) \geq H_{k}\left(T_{k}\right)\right\} \wedge N, \tag{18.5}
\end{equation*}
$$
^eq18-5

and we set

$$
\eta \triangleq \eta_{0}=\min \left\{k \geq 1: U_{k}\left(T_{k}\right) \geq H_{k}\left(T_{k}\right)\right\} \wedge N .
$$

With this definition, the value of a callable contract can be re-written as

$$
\begin{equation*}
H_{0}(0)=\mathrm{E}\left(B\left(T_{\eta}\right)^{-1} U_{\eta}\left(T_{\eta}\right)\right)=\mathrm{E}\left(\sum_{n=\eta}^{N-1} B\left(T_{n+1}\right)^{-1} X_{n}\right) \tag{18.6}
\end{equation*}
$$
^eq18-6

The Hamilton-Jacobi-Bellman equation that corresponds to the optimal stopping problem (18.4) can be solved by backward induction. In particular, we have, for $n=N-1, \ldots, 1$,

$$
\begin{equation*}
H_{n-1}\left(T_{n-1}\right)=B\left(T_{n-1}\right) \mathrm{E}_{T_{n-1}}\left(B\left(T_{n}\right)^{-1} \max \left(H_{n}\left(T_{n}\right), U_{n}\left(T_{n}\right)\right)\right) \tag{18.7}
\end{equation*}
$$
^eq18-7

subject to the terminal condition $H_{N-1} \equiv 0$. The recursion starts at the final time $n=N-1$ and progresses backward in time until we obtain the desired time 0 security value $H_{0}(0)$.

The financial meaning of the recursion above is straightforward. If a callable contract $H_{0}$ has not been exercised up to and including time $T_{n}$, then it is worth the hold value $H_{n}\left(T_{n}\right)$. If the callable contract is exercised at time $T_{n}$ its value is equal to $U_{n}\left(T_{n}\right)$. Assuming optimal exercise, the value of the callable Libor exotic $H_{0}$ at time $T_{n}$ is the maximum of the two,

$$
\max \left(H_{n}\left(T_{n}\right), U_{n}\left(T_{n}\right)\right)
$$

The value of this payoff at time $T_{n-1}$ is then

$$
B\left(T_{n-1}\right) \mathrm{E}_{T_{n-1}}\left(B\left(T_{n}\right)^{-1} \max \left(H_{n}\left(T_{n}\right), U_{n}\left(T_{n}\right)\right)\right)
$$

But clearly, this is the value of the CLE that can only be exercised at times $T_{n}$ and beyond, i.e. of the CLE $H_{n-1}$, as specified in (18.7).

### 18.2.3 Marginal Exercise Value Decomposition

Before discussing techniques to numerically implement the valuation equations of the previous section, let us briefly review an important decomposition result for CLEs that follows from the recursion (18.7). After a slight rewrite, we obtain

$$
\begin{aligned}
& H_{n-1}\left(T_{n-1}\right)-B\left(T_{n-1}\right) \mathrm{E}_{T_{n-1}}\left(B\left(T_{n}\right)^{-1} H_{n}\left(T_{n}\right)\right) \\
& \quad=B\left(T_{n-1}\right) \mathrm{E}_{T_{n-1}}\left(B\left(T_{n}\right)^{-1}\left(U_{n}\left(T_{n}\right)-H_{n}\left(T_{n}\right)\right)^{+}\right)
\end{aligned}
$$

Note that

$$
B\left(T_{n-1}\right) \mathrm{E}_{T_{n-1}}\left(B\left(T_{n}\right)^{-1} H_{n}\left(T_{n}\right)\right)=H_{n}\left(T_{n-1}\right)
$$

so that

$$
H_{n-1}\left(T_{n-1}\right)-H_{n}\left(T_{n-1}\right)=B\left(T_{n-1}\right) \mathrm{E}_{T_{n-1}}\left(B\left(T_{n}\right)^{-1}\left(U_{n}\left(T_{n}\right)-H_{n}\left(T_{n}\right)\right)^{+}\right)
$$

Taking discounted expectations to time 0 we obtain

$$
H_{n-1}(0)-H_{n}(0)=\mathrm{E}\left(B\left(T_{n}\right)^{-1}\left(U_{n}\left(T_{n}\right)-H_{n}\left(T_{n}\right)\right)^{+}\right)
$$

and, summing up from $n=1$ to $N-1$,

$$
H_{0}(0)-H_{N-1}(0)=\sum_{n=1}^{N-1} \mathrm{E}\left(B\left(T_{n}\right)^{-1}\left(U_{n}\left(T_{n}\right)-H_{n}\left(T_{n}\right)\right)^{+}\right)
$$

Since $H_{N-1}(0)=0$ and $H_{0}(0)$ is the time 0 value of the CLE, we have established the following proposition, which is essentially inspired by an integral representation of an American option from Jamshidian [1992], as presented in Proposition 1.10.7.

Proposition 18.2.1. The time 0 value $H_{0}(0)$ of a callable Libor exotic is equal to the sum of European options on the difference between the exercise and hold values at all exercise dates,

$$
\begin{equation*}
H_{0}(0)=\sum_{n=1}^{N-1} \mathrm{E}\left(B\left(T_{n}\right)^{-1}\left(U_{n}\left(T_{n}\right)-H_{n}\left(T_{n}\right)\right)^{+}\right) \tag{18.8}
\end{equation*}
$$
^eq18-8

In Proposition 18.2.1 each of the terms $\mathrm{E}\left(B\left(T_{n}\right)^{-1}\left(U_{n}\left(T_{n}\right)-H_{n}\left(T_{n}\right)\right)^{+}\right.$) can be interpreted as a "marginal" exercise value, i.e. as the incremental value of having an exercise at time $T_{n}$. The total CLE value is then equal to the sum of the marginal exercise values.

### 18.3 Monte Carlo Valuation

If a model admits a low-dimensional Markovian representation, then PDE methods are available for valuation of CLEs, and the backward recursion (18.7) is easy to implement, see Section 2.7.4. The situation is more complicated with Monte Carlo based models; fortunately, a range of methods for approximate solutions of optimal exercise problems in Monte Carlo have been developed. The mechanics of the scheme have been broadly outlined in Section 3.5.4, and we now proceed to discuss implementation details for the CLE class.

### 18.3.1 Regression-Based Valuation of CLEs, Basic Scheme

We start with a basic regression-based method for estimating the value of a CLE. As we recall from the discussion of Section 3.5.4, the regression-based LS (for Least Squares) scheme builds on the idea that the expected value of a random variable conditioned on information at time $T$ can be calculated in a Monte Carlo simulation by regressing the random variable against simulated state variables of the model at time $T$.

To make matters precise, we introduce some notation. Let $\zeta(t)= \left(\zeta_{1}(t), \ldots, \zeta_{q}(t)\right)^{\top}$ be a $q$-dimensional vector process of regression variables, to be defined later. For a given Monte Carlo path $\omega$, let us denote the value of a random variable $X$ on that path by $X(\omega)$. Suppose $K$ paths $\omega_{1}, \ldots, \omega_{K}$ are generated. For a random variable $X$, we denote by

$$
\mathcal{R}_{T}(X)
$$

the results of regression of the $K$-dimensional vector $\left(X\left(\omega_{1}\right), \ldots, X\left(\omega_{K}\right)\right)$ on the $K \times q$ matrix of regression variable observations at time $T$, $\left(\zeta\left(T, \omega_{1}\right), \ldots, \zeta\left(T, \omega_{K}\right)\right)^{\top}$, i.e.

$$
\mathcal{R}_{T}(X)=\zeta(T)^{\top} \beta
$$

where the $q$-dimensional column vector $\beta$ is obtained by, for instance ${ }^{3}$, solving the minimization problem

$$
\begin{equation*}
\left\|\left(X\left(\omega_{1}\right)-\zeta\left(T, \omega_{1}\right)^{\top} \beta, \ldots, X\left(\omega_{K}\right)-\zeta\left(T, \omega_{K}\right)^{\top} \beta\right)\right\|^{2} \rightarrow \min \tag{18.9}
\end{equation*}
$$
^eq18-9

This least-squares problem can be solved in closed form as explained in Section 3.5.4; to link our discussion here to the results of that section we simply need to set $\zeta_{i}(t)=\psi_{\imath}(x(t)), i=1, \ldots, q$, where $\psi$ 's and $x$ 's are defined in Section 3.5.4. For future reference, we denote the solution vector $\beta$ to (18.9) by[^27]

$$
\mathcal{C}\left(\mathcal{R}_{T}(X)\right),
$$
with the notation meant to be read as "coefficients of the regression of $X$ ", so that

$$
\begin{equation*}
\mathcal{R}_{T}(X)=\mathcal{C}\left(\mathcal{R}_{T}(X)\right)^{\top} \zeta(T) . \tag{18.10}
\end{equation*}
$$
^eq18-10

As it turns out, there are several possible LS schemes for CLE valuation. The most basic one is based on the idea of simply replacing the conditional expected value operator $\mathrm{E}_{T}$ in (18.7) with the regression operator $\mathcal{R}_{T}$ introduced above, i.e. to write

$$
\begin{equation*}
\widetilde{H}_{n-1}\left(T_{n-1}\right)=\mathcal{R}_{T_{n-1}}\left(\frac{B\left(T_{n-1}\right)}{B\left(T_{n}\right)} \max \left(\widetilde{H}_{n}\left(T_{n}\right), U_{n}\left(T_{n}\right)\right)\right), \tag{18.11}
\end{equation*}
$$
^eq18-11

for $n=N-1, \ldots, 1$, where $\widetilde{H}_{n}$ is an approximation to the true hold value $H_{n}$. This approach was originally suggested in Carriëre [1996] and Tsitsiklis and Roy [2001]. While we shall later describe better LS schemes than (18.11), let us nevertheless take some care in documenting all the steps necessary to apply it. Some of the steps are shared with the standard (non-callable) Monte Carlo valuation algorithm, but we list them anyway for completeness.

1. Choose and calibrate a term structure model (such as the LM model).
2. Decide on what to use for the regression variables process $\zeta(t)$ (we will have more to say about this later).
3. Simulate $K$ paths $\omega_{1}, \ldots, \omega_{K}$. For the LM model in particular, each $\omega_{k}$ represents one simulated path of all core Libor rates.
4. For each path $\omega_{k}$ calculate simulated values of the numeraire $B\left(T_{n}, \omega_{k}\right)$, $n=1, \ldots, N-1$.
5. For each path $\omega_{k}$, calculate ${ }^{4}$ the value $U_{n}\left(T_{n}, \omega_{k}\right)$ of the underlying exotic swap on all exercise clates $n=1, \ldots, N-1$.
6. For each path $\omega_{k}$, calculate the values of the $q$-dimensional regression variables process $\zeta$ on the exercise dates, $\zeta\left(T_{n}, \omega_{k}\right), n=1, \ldots, N-1$.
7. Set $\widetilde{H}_{N-1} \equiv 0$.
8. For each $n=N-1, \ldots, 1$
a) Form a $K$-dimensional vector $V_{n}=\left(V_{n}\left(\omega_{1}\right), \ldots, V_{n}\left(\omega_{K}\right)\right)^{\top}$,

$$
\begin{equation*}
V_{n}\left(\omega_{k}\right)=\frac{B\left(T_{n-1}, \omega_{k}\right)}{B\left(T_{n}, \omega_{k}\right)} \max \left(\widetilde{H}_{n}\left(T_{n}, \omega_{k}\right), U_{n}\left(T_{n}, \omega_{k}\right)\right) \tag{18.12}
\end{equation*}
$$
^eq18-12

for $k=1, \ldots, K$.

[^28]b) Calculate

$$
\tilde{H}_{n-1}\left(T_{n-1}\right)=\mathcal{R}_{T_{n-1}}\left(V_{n}\right)
$$
by regressing $V_{n}\left(\omega_{k}\right)$ against the matrix of regression variables observed on date $T_{n-1}$, using (18.9).
9. At this point we have computed $\widetilde{H}_{0}\left(T_{0}\right)$ which is an estimate of the value of the CLE. Return it.

Note that the last iteration in Step $8(n=1)$ involves regression on the values of $\zeta\left(T_{0}\right), \widetilde{H}_{0}\left(T_{0}\right)=\mathcal{R}_{T_{0}}\left(V_{1}\right)$. As $T_{0}=0, \zeta\left(T_{0}\right)$ is not random, so the regression here degenerates into a simple average ${ }^{5}$ of $V_{1}\left(\omega_{1}\right), \ldots, V_{1}\left(\omega_{K}\right)$.

There are a number of shortcomings in the above scheme that make it poorly suited for industrial-strength pricing of CLEs. We list them below as they shall guide us in systematically building more refined versions of the algorithm.

1. In Step 5, there is an assumption that exercise values $U_{n}\left(T_{n}\right)$ can be computed in closed form from information available at time $T_{n}$. While this is possible for, say, simple Bermudan swaptions, more complicated exotic swaps will generally violate this assumption.
2. The use of "regressions upon regressions", i.e. the fact that we apply $\mathcal{R}_{T_{n-1}}$ in (18.11) to (a function of) a regressed value $\widetilde{H}_{n}\left(T_{n}\right)$, could lead to significant biases building up as the scheme marches backward in time.
3. In general we cannot state whether $\widetilde{H}_{n}$ 's are low- or high-biased estimates of $H_{n}$ 's.

### 18.3.2 Regression for Underlying

To improve our basic LS scheme, we start out by examining the assumption of the basic LS scheme that exercise value at time $t$ can be be evaluated in closed form from time $t$ state variables (typically the simulated yield curve and the state of stochastic volatility parameters). As mentioned earlier, vanilla fixedfor floating-rate swaps underlying Bermudan swaptions certainly satisfy this assumption, as their values can be obtained by discounting projected cash flows on a simulated yield curve at time $t$. In principle, a number of other exotic swaps could fit under the assumption as well. For example, for a callable inverse floater, the underlying swap is a collection of simple Libor rate options, the values of which could be calculated from simulated market data at times $T_{n}, n=1, \ldots, N-1$, by applying option pricing formulas developed earlier in this book. While such a scheme is indeed possible, it would come at a significant computational cost of having to invoke option valuation formulas multiple times for each simulated path, i.e. easily thousands of

[^29]times overall. Moreover, closed-form caplet/swaption/CMS option pricing is rarely exact in term structure interest rate models; embedded into a backwards recursion algorithm these errors may potentially build up and skew the pricing results for the CLE.

Fortunately, it turns out that the extension of the basic scheme to arbitrary underlying swaps is simple. Specifically, we can use formula (18.2) which states that the exercise value $U_{n}\left(T_{n}\right)$ is equal to the conditional expected value of all (net) coupons paid after $T_{n}$. As we have already introduced the regression operator as a numerical proxy for conditional expected value, all we need to do now is to approximate $U_{n}\left(T_{n}\right)$ with the regressed value of all (net) coupons paid after $T_{n}$. Of course a coupon paid at time $T_{\imath+1}$ is always ${ }^{6}$ measurable with respect to $\mathcal{F}_{T_{\imath+1}}$ or, equivalently, can be calculated from the knowledge of the simulated model state variables up to and including $T_{i+1}$. So, all coupon values are known once a given path is simulated, and we can extend the basic scheme to arbitrary underlyings with the following two modifications. First, we replace Step 5 with

5a. For each path $\omega_{k}$, calculate the values of all net coupons $X_{n}\left(T_{n}, \omega_{k}\right)$, $n=1, \ldots, N-1$.

Second, we must replace the formula (18.12) in Step 8 of the basic scheme with

$$
\begin{equation*}
V_{n}\left(\omega_{k}\right)=\frac{B\left(T_{n-1}, \omega_{k}\right)}{B\left(T_{n}, \omega_{k}\right)} \max \left(\tilde{H}_{n}\left(T_{n}, \omega_{k}\right), \tilde{U}_{n}\left(T_{n}, \omega_{k}\right)\right), \quad k=1, \ldots, K, \tag{18.13}
\end{equation*}
$$
^eq18-13

where

$$
\begin{equation*}
\widetilde{U}_{n}\left(T_{n}\right)=\mathcal{R}_{T_{n}}\left(\sum_{i=n}^{N-1} \frac{B\left(T_{n}\right)}{B\left(T_{i+1}\right)} X_{\imath}\right) . \tag{18.14}
\end{equation*}
$$
^eq18-14

We can write (18.14) in a backward-recursive format,

$$
\begin{equation*}
Y_{n}=\frac{B\left(T_{n}\right)}{B\left(T_{n+1}\right)}\left(X_{n}+Y_{n+1}\right), \quad \tilde{U}_{n}\left(T_{n}\right)=\mathcal{R}_{T_{n}}\left(Y_{n}\right), \quad n=N-1, \ldots, 1, \tag{18.15}
\end{equation*}
$$
^eq18-15

where we start from $Y_{N} \equiv 0$. In this form it fits nicely into the backward recursion of the basic LS scheme.

Interestingly, we can rewrite (18.2) (for $t=T_{n}$ ) as

$$
\begin{equation*}
U_{n}\left(T_{n}\right)=\mathrm{E}_{T_{n}}\left(B\left(T_{n}\right) B\left(T_{n+1}\right)^{-1}\left(X_{n}+U_{n+1}\left(T_{n+1}\right)\right)\right), \tag{18.16}
\end{equation*}
$$
^eq18-16

which gives raise to an alternative backward scheme for $\widetilde{U}_{n}\left(T_{n}\right)$,
$\tilde{U}_{n}\left(T_{n}\right)=\mathcal{R}_{T_{n}}\left(B\left(T_{n}\right) B\left(T_{n+1}\right)^{-1}\left(X_{n}+\tilde{U}_{n+1}\left(T_{n+1}\right)\right)\right), n=N-1, \ldots, 1$.

[^30]While (18.16) is trivially equivalent to (18.2) due to the additivity of the conditional expected values and the iterated conditional expectations property, (18.17) is not equivalent to (18.15). In (18.17), on step $n$ we add (discounted) coupon $X_{n}$ to an already-regressed value of future coupons $\widetilde{U}_{n+1}\left(T_{n+1}\right)$, whereas in (18.15) we sum up values of un-regressed coupons from $N-1$ to $n$ and then regress them to obtain $\widetilde{U}_{n}$. The difference between the two schemes originates with the fact that the regression operator does not satisfy an equivalent to the iterated conditional expectations property ${ }^{7}$ (see footnote 3 of Chapter 4). An inquisitive reader may ask whether one scheme is better (in the sense of producing smaller bias for the value of the CLE) than the other. The answer is not straightforward. While (18.15) avoids applying regression to the output of other regressions and thus could be expected to produce lower bias, empirical evidence suggests that no scheme is universally better than the other for all CLEs. We consequently recommend a flexible implementation that can use both schemes.

### 18.3.3 Valuing CLE as a Cancelable Note

Using regressions for the underlying is not the only approach that extends the basic scheme in Section 18.3.2 to arbitrary underlying swaps. An alternative scheme is based on the idea of representing a CLE as a cancelable note. To describe this in more detail, let us denote

$$
G_{n}(t)=H_{n}(t)-U_{n}(t)
$$

and obtain from (18.7) and (18.16) that

$$
\begin{align*}
& G_{n-1}\left(T_{n-1}\right) \\
& \quad=H_{n-1}\left(T_{n-1}\right)-U_{n-1}\left(T_{n-1}\right) \\
& \quad=\mathrm{E}_{T_{n-1}}\left(B\left(T_{n-1}\right) B\left(T_{n}\right)^{-1}\left(\max \left(H_{n}\left(T_{n}\right), U_{n}\left(T_{n}\right)\right)-\left(X_{n-1}+U_{n}\left(T_{n}\right)\right)\right)\right) \\
& \quad=\mathrm{E}_{T_{n-1}}\left(B\left(T_{n-1}\right) B\left(T_{n}\right)^{-1}\left(\left(H_{n}\left(T_{n}\right)-U_{n}\left(T_{n}\right)\right)^{+}-X_{n-1}\right)\right) \\
& \quad=\mathrm{E}_{T_{n-1}}\left(B\left(T_{n-1}\right) B\left(T_{n}\right)^{-1}\left(-X_{n-1}+G_{n}\left(T_{n}\right)^{+}\right)\right) \tag{18.18}
\end{align*}
$$
^eq18-18

$n=N, \ldots, 1$, where we for uniformity of notation have introduced a "fake" coupon at time zero,

$$
X_{0}=0 .
$$

We see that the value $G_{0}(0)$ is the value of the swap that pays (net) coupons $-X_{n}$ on dates $T_{n}, n=0, \ldots, N-1$, plus the right to cancel it, at zero cost, on any of the exercise dates $T_{1}, \ldots, T_{N-1}$. In fact, as explained in Section 5.14 , it is this structure (a callable structured note) that a bank usually sells to clients. While the representation of the cancelable note as a CLE plus the

[^31]underlying non-callable swap is often convenient for risk management, for the purposes of valuation it is, in fact, often useful to consider the original format.

The LS version of (18.18) is, naturally, given by

$$
V_{n}\left(\omega_{k}\right)=\frac{B\left(T_{n-1}, \omega_{k}\right)}{B\left(T_{n}, \omega_{k}\right)}\left(-X_{n-1}\left(\omega_{k}\right)+\tilde{G}_{n}\left(T_{n}, \omega_{k}\right)^{+}\right), \quad k=1, \ldots, K,
$$

followed by

$$
\begin{equation*}
\widetilde{G}_{n-1}\left(T_{n-1}\right)=\mathcal{R}_{T_{n-1}}\left(V_{n}\right), \tag{18.19}
\end{equation*}
$$
^eq18-19

for each $n=N-1, \ldots, 1$. The starting point is given by $\widetilde{G}_{N-1}\left(T_{N-1}\right) \equiv 0$. We trust the reader should have no problem amending the basic scheme to use the cancelable note representation.

Interestingly, by linearity (which does hold for the regression operator, unlike the tower property), we can rewrite (18.19) as

$$
\begin{aligned}
\widetilde{G}_{n-1}\left(T_{n-1}\right) & =-\mathcal{R}_{T_{"-1}}\left(B\left(T_{n-1}\right) B\left(T_{n}\right)^{-1} X_{n-1}\right) \\
+ & \mathcal{R}_{T_{"-1}}\left(B\left(T_{n-1}\right) B\left(T_{n}\right)^{-1} \widetilde{G}_{n}\left(T_{n}\right)^{+}\right), \quad n=N-1, \ldots, 1
\end{aligned}
$$

In other words, to get the value of the cancelable note at time $T_{n-1}$, start with its value $\widetilde{G}_{n}\left(T_{n}\right)$ at time $T_{n}$, apply the optimal cancelability condition $\widetilde{G}_{n}\left(T_{n}\right) \rightarrow \widetilde{G}_{n}\left(T_{n}\right)^{+}$, discount to $T_{n-1}$, regress on $\zeta\left(T_{n-1}\right)$, and then add the time $T_{n-1}$ value of the $(n-1)$-th coupon. If $X_{n-1}$ is actually $\mathcal{F}_{T_{n-1}}$ measurable, as is often the case (exceptions include range accrual coupons and coupons that depend on rates observed in-arrears, i.e. at the end of the observation period rather than the beginning), the scheme simplifies a bit more,

$$
\begin{align*}
& \widetilde{G}_{n-1}\left(T_{n-1}\right)=-P\left(T_{n-1}, T_{n}\right) X_{n-1} \\
& \quad+\mathcal{R}_{T_{n-1}}\left(B\left(T_{n-1}\right) B\left(T_{n}\right)^{-1} \widetilde{G}_{n}\left(T_{n}\right)^{+}\right), \quad n=N-1, \ldots, 1 \tag{18.20}
\end{align*}
$$
^eq18-20

Once the value $\widetilde{G}_{0}(0)$ has been calculated, the estimate of the CLE value $H_{0}(0)$ can be obtained via $H_{0}(0)=\widetilde{G}_{0}(0)+U_{0}(0)$, where $U_{0}$ may, if not available in closed form, be calculated via a (standard) Monte Carlo algorithm for (18.2). In estimating $U_{0}$, we would normally reuse the paths $\omega_{k}, k=1, \ldots, K$ that were used for the computation of $\widetilde{G}_{0}(0)$.

### 18.3.4 Using Regressed Variables for Decision Only

Our second criticism of the basic LS scheme of Section 18.3.1 focused on the issue that the values that are regressed at time $T_{n-1}$ themselves come as the result of a regression at time $T_{n}$. Such compounded regression could lead
to substantial biases. Interestingly, a small modification of the algorithm reduces this bias significantly.

Let us work with the cancelable note scheme (18.20), although the same idea can easily be applied to the basic CLE scheme. Going back to the expression (18.6) for the value of the CLE expressed as the sum of all coupons paid post optimal exercise, we note that a similar formula holds for cancelable notes,

$$
\begin{equation*}
G_{0}(0)=-\mathrm{E}\left(\sum_{n=0}^{\eta-1} B\left(T_{n+1}\right)^{-1} X_{n}\right) \tag{18.21}
\end{equation*}
$$
^eq18-21

(recall that $X_{0}=0$ ). The exercise index $\eta$ here is the same as in (18.6), since the optimal cancel time for the cancelable note $G$ is the same as the optimal exercise time for the CLE $H$. We rewrite the formula as

$$
G_{0}(0)=-E\left(\sum_{n=0}^{N-1} B\left(T_{n+1}\right)^{-1} X_{n} 1_{\{\eta>n\}}\right),
$$

and note that

$$
\begin{equation*}
1_{\{\eta>n\}}=\prod_{i=1}^{n} 1_{\left\{G_{2}\left(T_{1}\right)>0\right\}} . \tag{18.22}
\end{equation*}
$$
^eq18-22

Let us define $V_{n}$ 's recursively,
$V_{n}=B\left(T_{n}\right) B\left(T_{n+1}\right)^{-1}\left(-X_{n}+1_{\left\{G_{n+1}\left(T_{n+1}\right)>0\right\}} V_{n+1}\right), \quad n=N-1, \ldots, 0$,
with $V_{N}=0$. Then

$$
V_{0}=-\sum_{n=0}^{N-1} B\left(T_{n+1}\right)^{-1}\left(\prod_{\imath=1}^{n} 1_{\left\{G_{,}\left(T_{\imath}\right)>0\right\}}\right) X_{n}
$$

and, computing the expected value and using (18.22), we obtain that

$$
G_{0}(0)=\mathrm{E}\left(V_{0}\right)
$$

Moreover,

$$
\begin{equation*}
G_{n}\left(T_{n}\right)=\mathrm{E}_{T_{n}}\left(V_{n}\right) \tag{18.24}
\end{equation*}
$$
^eq18-24

We note that the recursion for $V_{n}$ involves the value of the cancelable note for exercise decisions only, through the indicators $1_{\left\{G_{2}\left(T_{i}\right)>0\right\}}$, whereas the coupon values that are added up, $X_{i}$ 's, are never regressed. Following Longstaff and Schwartz [2001], we can take advantage of this observation by defining a new approximation $\widehat{G}_{n}\left(T_{n}\right), n=1, \ldots, N-1$, to the true value of the cancelable note by

$$
\begin{align*}
& \widehat{G}_{n-1}\left(T_{n-1}\right)=B\left(T_{n-1}\right) B\left(T_{n}\right)^{-1}\left(-X_{n-1}+1_{\left\{\widetilde{G}_{\prime \prime}\left(T_{\prime \prime}\right)>0\right\}} \widehat{G}_{n}\left(T_{n}\right)\right),  \tag{18.25}\\
& \widetilde{G}_{n-1}\left(T_{n-1}\right)=\mathcal{R}_{T_{n-1}}\left(\widehat{G}_{n-1}\left(T_{n-1}\right)\right), \tag{18.26}
\end{align*}
$$
^eq18-25
^eq18-26

defined backwards for $n=N, \ldots, 1$. The first equation here comes from (18.23) and the second from (18.24). We emphasize that regression is only used to establish the exercise indicator functions in (18.25), i.e. whether and where to exercise on each path. In $(18.24), \widehat{G}_{0}(0)=\widehat{G}_{0}(0, \omega)$ is the (random) accumulation of discounted coupons up to exercise, and our estimate of the true value of the cancelable note $G_{0}$ is given by $\widetilde{G}_{0}(0)$, i.e the simple average of the realizations of $\widehat{G}_{0}\left(0, \omega_{1}\right), \ldots, \widehat{G}_{0}\left(0, \omega_{K}\right)$. Clearly, this is the numerical equivalent of the (unconditional) expected value in (18.21). We find that the scheme (18.25)-(18.26) typically has significantly less bias than the naive scheme (18.20).

While (18.25)-(18.26) is the scheme that we recommend for most applications, Egloff et al. [2007] introduce a "blend" of (18.25)-(18.26) and (18.20) with a tunable parameter that can be optimized over to select a scheme with the lowest bias. To briefly outline this idea, we note that the scheme (18.25)-(18.26) only uses un-regressed values of coupons while (18.20) always uses regressed values. A "blended" scheme uses the first few coupons that are unregressed, while the rest are regressed. For example, since we can write

$$
G_{0}(0)=-\mathrm{E}\left(\sum_{n=0}^{m-1} B\left(T_{n+1}\right)^{-1} X_{n} 1_{\{\eta>n\}}+B\left(T_{m}\right)^{-1} G_{m}\left(T_{m}\right) 1_{\left\{G_{m}\left(T_{m}\right)>0\right\}}\right)
$$

for any $m=1, \ldots, N-1$, we can replace (18.25) with

$$
\begin{aligned}
& \widehat{G}_{n-1}\left(T_{n-1}\right)=-\sum_{\imath=n-1}^{\mid n+m-1\lceil-1}\left(\frac{B\left(T_{n-1}\right)}{B\left(T_{i+1}\right)}\left(\prod_{j=n}^{i} 1_{\left\{\widetilde{G}_{i}\left(T_{1}\right)>0\right\}}\right) X_{i}\right) \\
+ & \mathcal{R}_{T_{n-1}}\left(\frac{B\left(T_{n-1}\right)}{B\left(T_{\prod_{n+m-1 \Gamma}}\right)}\left(\prod_{j=n}^{\mid n+m-1\lceil } 1_{\left\{\widetilde{G}_{\left.,\left(T_{,}\right)>0\right\}}\right.}\right) \widetilde{G}_{\left.\right|_{n+m-1 \Gamma}\left(T_{\eta_{n+m-1 \Gamma}}\right)}\right)
\end{aligned}
$$

coupled with (18.26), where we have denoted $l l \triangleq l \wedge N$. We would then run the regression algorithm for different $m$ 's (reusing the paths, of course) and choose the value of $m$ that would give us the highest value ${ }^{8}$. We recover (18.20) for $m=1$ and (18.25)-(18.26) for $m=+\infty$.

### 18.3.5 Regression Valuation with Boundary Optimization

The regression-based method can sometimes be improved with methods similar to those of the parametric boundary optimization discussed in Sec-

[^32]tion 3.5.2, an idea popularized by Bender et al. [2006]. Starting with, for example, scheme (18.25)-(18.26), we fix a collection of trigger thresholds $\mathbf{h}=\left(h_{1}, \ldots, h_{N-1}\right)$ and rewrite (18.25)-(18.26) as (note dependence on the triggers in the indicator functions)
$\widehat{G}_{n-1}\left(T_{n-1}, \mathbf{h}\right)=B\left(T_{n-1}\right) B\left(T_{n}\right)^{-1}\left(-X_{n-1}+1_{\left\{\tilde{G}_{n}\left(T_{n}, \mathbf{h}\right)>h_{n}\right\}} \widehat{G}_{n}\left(T_{n}, \mathbf{h}\right)\right)$, $\widetilde{G}_{n-1}\left(T_{n-1}, \mathbf{h}\right)=\mathcal{R}_{T_{n-1}}\left(\widehat{G}_{n-1}\left(T_{n-1}, \mathbf{h}\right)\right)$,
$n=N-1, \ldots, 1$. For a given set $\mathbf{h}$, the scheme produces an estimate $\widetilde{G}_{0}(0, h)$ of the value of the cancelable note. Then, iterating over $h$, we find the optimal (highest) value of $\widetilde{G}_{0}$ and return this as our improved estimate for the cancelable note. Of course, if our original regression was fundamentally sound, we should see the optimal value of $h$ being very close to $(0,0, \ldots, 0)$, in which case the trigger iteration adds no value. To the extent, however, that the original LS method produces significantly sub-optimal exercise decisions (e.g. due to a poor choice of regression variables), the trigger iteration may lead to pick-up of substantial additional value.

The search for the optimal value of the triggers can be efficiently organized as a sequence of $N-1$ one-dimensional optimizations, along the same lines as the algorithm in Section 3.5.3. In particular, we see that $\widehat{G}_{n-1}\left(T_{n-1}, \mathrm{~h}\right)$ depends on $h_{n}, \ldots, h_{N-1}$ only. Moreover, if for a given $n, h_{n}$ maximizes the value of $\widetilde{G}_{0}(0, h)$ then it also maximizes the value of the cancelable note with first $n-1$ exercise dates removed, i.e. the value $\mathcal{R}_{0}\left(\widehat{G}_{n-1}\left(T_{n-1}, \mathrm{~h}\right)\right)$ (recall that $\mathcal{R}_{0}(X)=K^{-1} \sum X\left(\omega_{k}\right)$, i.e. just the average of path values of $X)$. Hence, we find the optimal value of the $n$-th trigger $h_{n}^{*}$ via

$$
h_{n}^{*}=\underset{h_{n}}{\operatorname{argmin}} \mathcal{R}_{0}\left(\widehat{G}_{n-1}\left(T_{n-1},\left(h_{n}, h_{n+1}^{*}, \ldots, h_{N-1}^{*}\right)\right)\right), \quad n=N-1, \ldots, 1,
$$

where, slightly abusing notation, we use ( $h_{n}, h_{n+1}^{*}, \ldots, h_{N-1}^{*}$ ) to denote a vector h with the last $N-n-1$ elements fixed to the optimal values found on previous steps (and first $n-1$ elements irrelevant). The optimization problem above is easy to solve, but we remind the reader of our comments from Section 3.5.3, where we noted that for a finite-path simulation, the objective functions in each of optimization problems will not be smooth, so one should avoid the use of a derivatives-based numerical optimizer.

### 18.3.6 Lower Bound via Regression Scheme

All the variations of the regression scheme developed so far produce an estimate of the value of the CLE $H_{0}$ (or, equivalently, the value of the corresponding cancelable note $G_{0}$ ) but the bias of the estimate is generally unknown. On one hand, the exercise decisions used in the schemes are necessarily suboptimal, as we use estimates, rather than actual values, of
hold/underlying variables to define them. This, in isolation, suggests that our estimate should be biased low. But on the other hand, our schemes use the same set of sample paths to estimate the exercise decision as to calculate the values of the security if it is exercised. This could lead to an upward bias in the estimate as some amount of future information can affect our decision to exercise, leading to a "perfect foresight" high-bias, see (3.110).

For risk assessment and for price quotation, it is typically desirable to know the sign of biases in computed security prices. To control the sign of the bias in CLE valuation by regression methods, we can follow the advice of Section 3.5 and use the LS regression scheme to estimate the exercise decision rule only, while using an independent simulation to calculate the value of the CLE given that exercise rule. A typical implementation algorithm is outlined below.

1. Run the basic regression-based scheme in Section 18.3.1 or any of the alternatives in Sections 18.3.2-18.3.5.
2. The output of the regression are the regression coefficients for the hold and exercise values at all exercise times, $\mathcal{C}\left(\widetilde{H}_{n}\left(T_{n}\right)\right), \mathcal{C}\left(\widetilde{U}_{n}\left(T_{n}\right)\right)$, $n=1, \ldots, N-1$. See (18.10).
3. Simulate additional $K^{\prime}$ paths, $\omega_{1}^{\prime}, \ldots, \omega_{K^{\prime}}^{\prime}$, that are independent from the paths used in the regression scheme.
4. For each path $\omega_{k}^{\prime}$, calculate the values of the $q$-dimensional regression variables process $\zeta$ on the exercise dates, $\zeta\left(T_{n}, \omega_{k}^{\prime}\right), n=1, \ldots, N-1$.
5. For each path $\omega_{k}^{\prime}$, calculate an estimate of the exercise index $\widetilde{\eta}$ by

$$
\begin{align*}
& \tilde{\eta}\left(\omega_{k}^{\prime}\right)=\min \left\{n \geq 1: \mathcal{C}\left(\tilde{U}_{n}\left(T_{n}\right)\right)^{\top} \zeta\left(T_{n}, \omega_{k}^{\prime}\right)\right. \\
& \left.\quad \geq \mathcal{C}\left(\tilde{H}_{n}\left(T_{n}\right)\right)^{\top} \zeta\left(T_{n}, \omega_{k}^{\prime}\right)\right\} \wedge N . \tag{18.27}
\end{align*}
$$
^eq18-27

In simple terms the exercise is based on the regression estimates of the exercise and hold values obtained in the basic scheme, applied to the new values of the regression variables $\zeta$.
6. Calculate the CLE value as the Monte Carlo value of a knock-in discrete barrier option based on the exercise rule $\widetilde{\eta}$,

$$
\begin{equation*}
H_{0}(0) \approx \frac{1}{K^{\prime}} \sum_{k=1}^{K^{\prime}}\left(\sum_{n=\widetilde{\eta}\left(\omega_{k}^{\prime}\right)}^{N-1} B\left(T_{n+1}, \omega_{k}^{\prime}\right)^{-1} X_{n}\left(\omega_{k}^{\prime}\right)\right) \tag{18.28}
\end{equation*}
$$
^eq18-28

The guaranteed low bias of this two-stage scheme comes at an additional cost of simulating and evaluating payoffs on extra $K^{\prime}$ paths. For the record, we often choose $K^{\prime} \approx 10,000$ to 100,000 and $K \approx K^{\prime} / 4$ to $K^{\prime} / 2$. So the cost is not inconsiderable, on average slowing down a valuation by a factor of up to 2 or so. Importantly, this additional cost is not incurred in the evaluation of most risk sensitivities, since we can reuse the exercise
boundary in calculations of first-order sensitivities, as we explain in Chapter 24. However, if performance is an issue, it is worth pointing out that we often find the values obtained by the scheme (18.25)-(18.26) (with only a single batch of paths simulated) to be pretty close to those from the more costly two-stage scheme (18.28). In particular, it appears that in many practical situations (18.25)-(18.26) is biased low, even though it is not guaranteed to be so.

It is easy to see the strong connection between the scheme (18.25)-(18.26) and (18.28). If, instead of the independent paths in the second stage of (18.28) we used the same paths as in the regression scheme, i.e. $K^{\prime}=K$ and $\omega_{k}^{\prime}=\omega_{k}$ for $k=1, \ldots, K$, then the values produced by the two scheme would be exactly the same. We leave the verification of this simple fact to the reader.

### 18.3.7 Iterative Improvement of Lower Bound

Consider the cancelable note that pays (net) coupons $-X_{n}$ at $T_{n+1}, n= 0, \ldots, N-1$. Suppose we are given some exercise policy, generally suboptimal. An exercise policy is, in essence, a stopping time index $\alpha$ that specifies when the note is canceled. For technical reasons we want this exercise strategy to be specified not just for the original cancelable note but also for cancelable notes with the first $k, k=1, \ldots, N-2$ exercise dates removed. This is most conveniently expressed as a collection of stopping times $\boldsymbol{\alpha}=\left(\alpha_{0}, \ldots, \alpha_{N-1}\right)$, with $\alpha_{0}=\alpha$, that satisfy the following exercise policy consistency conditions,

$$
n+1 \leq \alpha_{n} \leq N, \quad \alpha_{N-1} \equiv N,
$$

and

$$
\alpha_{n}>n+1 \Rightarrow \alpha_{n}=\alpha_{n+1}, \quad n=0, \ldots, N-2 .
$$

Let us consider the part of the note that includes coupons $-X_{n}, \ldots,-X_{N-1}$ only. Then the value of this note, exercised per stopping time $\alpha_{k}, k \geq n$, is given at time $T_{n}$ by

$$
\begin{aligned}
G_{n}^{\alpha_{k}}\left(T_{n}\right) & =-B\left(T_{n}\right) \mathrm{E}_{T_{n}}\left(\sum_{i=n}^{\alpha_{k}-1} B\left(T_{i+1}\right)^{-1} X_{i}\right) \\
& =-B\left(T_{n}\right) \mathrm{E}_{T_{n}}\left(\sum_{i=n}^{N-1} B\left(T_{i+1}\right)^{-1} X_{i} 1_{\left\{\alpha_{k}>i\right\}}\right)
\end{aligned}
$$

We note that $\alpha_{k} \geq k+1$, so the coupons $-X_{n}, \ldots,-X_{k}$ are always paid. The optimal exercise policy $\boldsymbol{\eta} \triangleq\left(\eta_{0}, \ldots, \eta_{N-1}\right)$, as defined by (18.5), satisfies the consistency conditions, and we of course have

$$
G_{n}^{\eta_{n}}\left(T_{n}\right)=G_{n}\left(T_{n}\right), \quad n=0, \ldots, N-1,
$$

where on the right-hand side we have the actual values of the remaining parts of the cancelable note. The approximations to the optimal exercise rule
we developed in previous sections, such as (18.27), satisfy the consistency conditions as well.

As pointed out previously, for any exercise policy $\boldsymbol{\alpha}$ we have

$$
G_{n}^{\alpha_{n}}\left(T_{n}\right) \leq G_{n}\left(T_{n}\right), \quad n=0, \ldots, N-1,
$$

so an exercise policy gives us a way to obtain a lower bound, namely $G_{0}^{\alpha_{11}}(0)$, for the time 0 cancelable note value $G_{0}(0)$. It turns out that the general theory of optimal stopping tells us how to improve a given exercise policy, i.e. how to find another exercise policy that would be better, in the sense of producing a higher lower bound. The improvements could be iterated, eventually (after $N-1$ iterations) converging to the optimal exercise policy irrespective of the starting point. This "policy iteration" method was first applied to the pricing of callable derivatives by Bender et al. [2006] and Kolodko and Schoenmakers [2006]. To demonstrate how policy iteration works, for a given policy $\boldsymbol{\alpha}$, let us define its improvement $\widehat{\boldsymbol{\alpha}}=\left(\widehat{\alpha}_{0}, \ldots, \widehat{\alpha}_{N-1}\right)$ by

$$
\begin{equation*}
\widehat{\alpha}_{n}=\min \left\{k>n: \max _{j \geq k} G_{k}^{\alpha_{\prime}}\left(T_{k}\right) \leq 0\right\} \wedge N, \quad n=0, \ldots, N-1 . \tag{18.29}
\end{equation*}
$$
^eq18-29

To understand this definition, first note that $G_{k}^{\alpha i}\left(T_{k}\right)$ is the time- $T_{k}$ value of the coupon stream from $T_{k}$ onwards, exercised per stopping time $\alpha_{j}$, i.e. we basically add up the (discounted) coupons $-X_{k}, \ldots,-X_{j}$ and thereafter the remaining coupons subject to the original exercise rule. Hence, the improved exercise rules states that we should exercise (cancel the note) at the first date $T_{k}$ for which holding on to the note under the original exercise policy does not make sense, specifically if the maximum of the remaining value of the note over all original exercise rules that go strictly past $T_{k}$ is non-positive.

The proof that $\widehat{\boldsymbol{\alpha}}$ is an improvement over $\boldsymbol{\alpha}$, i.e. that

$$
G_{0}^{\widehat{\alpha}_{10}}(0) \geq G_{0}^{\alpha_{10}}(0)
$$

can be found in Kolodko and Schoenmakers [2006]. The improvement could be applied to the policy $\widehat{\boldsymbol{\alpha}}$ as well, and this can be iterated multiple times. After $N-1$ iterations the exercise policy converges to the optimal one, as also proven in Kolodko and Schoenmakers [2006].

While in theory we can find the optimal exercise strategy by iteration starting from any initial policy - such as a trivial policy $\alpha_{n}=n+1$ for any $n=0, \ldots, N-1$ - performing more than a few iterations is impractical, as shall be discussed shortly. In practice, a sensible strategy would apply just one round of improvements to an already decent exercise policy, such as (18.27) obtained by the basic scheme of Section 18.3.1 or its variations.

The high numerical cost of the policy iteration stems from the presence of terms like $G_{k}^{\alpha \prime}\left(T_{k}\right)$ in (18.29). These are the time $T_{k}$ values of the coupon stream exercised per some rule, and are rarely, if ever, available in closed
form. One might guess that we could estimate the required conditional expected values via regressions, as we have done for other quantities. Alas, using regressed values in (18.29) in place of true conditional expected values generally leads to no policy improvements. Consider, for example, the scheme (18.25)-(18.26), with the corresponding exercise policy defined by

$$
\alpha_{n}=\min \left\{k>n: \widetilde{G}_{k}\left(T_{k}\right) \leq 0\right\} \wedge N, \quad n=0, \ldots, N-1 .
$$

Then the iterated exercise policy would base the exercise decision on $\max _{j \geq k} \widetilde{G}_{k}^{\alpha_{\prime}}\left(T_{k}\right)$. However, as $\alpha_{n}$ 's are constructed to be the optimal stopping policy for $\widetilde{G}_{n}\left(T_{n}\right)$ 's, we always have

$$
\widetilde{G}_{k}^{\alpha_{i}}\left(T_{k}\right) \leq \widetilde{G}_{k}^{\alpha_{k}}\left(T_{k}\right)=\widetilde{G}_{k}\left(T_{k}\right),
$$

and

$$
\max _{j \geq k} \widetilde{G}_{k}^{\alpha \prime}\left(T_{k}\right)=\widetilde{G}_{k}\left(T_{k}\right)
$$

for any $k$. Thus the "improved" policy $\widehat{\boldsymbol{\alpha}}$ computed from (18.29) would coincide with $\alpha$, the original policy.

As regression methods cannot be used, an alternative is to resort to "bruteforce" nested simulation to compute unbiased estimates for the conditional expectations $G_{k}^{\alpha_{\prime}}\left(T_{k}\right)$. The basic idea is simple: if we need to estimate $\mathrm{E}_{T}(X)$ for some random variable $X$ then, for each Monte Carlo path $\omega$, we simulate a number of additional sub-paths with the simulated model state at time $T$ on path $\omega$ as a starting point, and then estimate the conditional expected value on path $\omega$ by averaging the values of $X$ realized on all sub-paths ${ }^{9}$. For policy iteration, such sub-paths must be launched for each (outer) simulated path $\omega_{k}, k=1, \ldots, K$, for each exercise time $T_{n}, n=1, \ldots, N-1$, so the computational expense is quite considerable even with a modest number of sub-paths. Bender et al. [2006] recommend using control variates to speed up valuation, and Beveridge and Joshi [2009] list a number of additional suggestions to improve computational performance. Nevertheless, there is no doubt that to keep Monte Carlo simulation error of the conditional expectation low enough for the policy iteration scheme to be effective, a substantial numerical effort is required. As such, we do not recommend routine application of policy iteration in the pricing of CLEs. This, of course, will require that the regression estimates of the exercise policy ${ }^{10}$ are sufficiently accurate as is, something that to a large extent depends on how careful we are in choosing the regression variable vector $\zeta(t)$. We turn to this topic in Section 18.3.9, but first we need to address the fundamental problem of how to test whether a given exercise policy is close to optimal in the first place. One useful approach to this problem is to use the estimated

[^33]exercise policy to construct an upper bound for the option price, which, if close to the lower bound, will give us confidence that our exercise policy is close to optimal. Section 18.3.8 below discusses this technique in detail.

### 18.3.8 Upper Bound

Computing a lower bound for a CLE price is straightforward: pick some exercise policy and price the CLE by Monte Carlo methods. Assuming our computation of the exercise policy did not "cheat" by using information from the Monte Carlo trials subsequently used for valuation, the resulting price estimates will always have a non-positive bias, as the chosen exercise strategy will almost certainly be suboptimal. The lower bound algorithm of Section 18.3.6 was based on precisely such a strategy. To complement a computed lower bound for the CLE, we are now interested in using the lower bound exercise policy to construct an upper bound for the CLE price. Taken together with the lower bound, the upper bound can be used to construct a valid confidence interval for the CLE price; this, in turn, will allow us to assess the quality of the exercise policy.

### 18.3.8.1 Basic Ideas

Our strategy to construct an upper bound for CLEs will draw directly on the duality results in Section 1.10.2 and the generic upper bound simulation ideas in Section 3.5.5. To formulate these results in the CLE setting, let us start with the description (18.4)-(18.6):

$$
\begin{equation*}
H_{0}(0)=\sup _{\xi \in \mathcal{T}_{0}} \mathrm{E}\left(B\left(T_{\xi}\right)^{-1} U_{\xi}\left(T_{\xi}\right)\right)=\mathrm{E}\left(B\left(T_{\eta}\right)^{-1} U_{\eta}\left(T_{\eta}\right)\right) . \tag{18.30}
\end{equation*}
$$
^eq18-30

Following Section 1.10.2, let $\mathcal{K}$ denote the space of adapted martingales $M$ for which $\sup _{\xi \in \mathcal{T}_{0}} \mathrm{E}\left|M\left(T_{\xi}\right)\right|<\infty$. For any martingale $M \in \mathcal{K}$, (1.71) demonstrates that

$$
\begin{equation*}
H_{0}(0) \leq M(0)+\mathrm{E}\left(\max _{n=1, \ldots, N-1}\left(\frac{U_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}-M\left(T_{n}\right)\right)\right) . \tag{18.31}
\end{equation*}
$$
^eq18-31

Also, we know from the duality result (1.72) that this upper bound will become an equality provided that $M$ is chosen to be the martingale component of the (supermartingale) deflated value process of the CLE. To emphasize this result, set

$$
\begin{equation*}
V_{\mathrm{CLE}}(t)=B(t) \mathrm{E}_{t}\left(B\left(T_{n}\right)^{-1} \max \left(U_{n}\left(T_{n}\right), H_{n}\left(T_{n}\right)\right)\right), \quad t \in\left(T_{n-1}, T_{n}\right], \tag{18.32}
\end{equation*}
$$
^eq18-32

and use the Doob-Meyer decomposition of Section 1.10.2 to write $V_{\text {CLE }}(t) / B(t)=m_{\text {CLE }}(t)-A(t)$, where $m_{\text {CLE }}(t)$ is a martingale and $A(t)$ an increasing predictable process with $A(0)=0$. Then setting $M(t)=m_{\text {CLE }}(t)$ yields

$$
\begin{equation*}
m_{\mathrm{CLE}}(0)+\mathrm{E}\left(\max _{n=1, \ldots, N-1}\left(\frac{U_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}-m_{\mathrm{CLE}}\left(T_{n}\right)\right)\right)=H_{0}(0) . \tag{18.33}
\end{equation*}
$$
^eq18-33

If the underlying model is driven by a vector-valued $\mathrm{Q}^{B}$-Brownian motion $W(t)$, the martingale representation theorem (Theorem 1.1.4) shows that any martingale $M$ in (18.31) must be of the form

$$
\begin{equation*}
M(t)=\int_{0}^{t} \sigma(s)^{\top} d W(s) \tag{18.34}
\end{equation*}
$$
^eq18-34

for some adapted vector-process $\sigma(t)$ satisfying the usual conditions required for the stochastic integral to be proper martingale. Clearly, however, if $\sigma(t)$ is chosen arbitrarily, the resulting upper bound computed from (18.31) is likely to be very loose, and probably not particularly useful. While (18.33) is of little immediate practical use (since we do not know the process $V_{\text {CLE }}(t) / B(t)$ ), it does suggest that for a chosen martingale $M(t)$ in (18.31) to produce a tight upper bound, it needs to be "close" to $m_{\text {CLE }}(t)$.

Several strategies have been proposed for constructing a good martingale $M(t)$. When working in a simple model setup on simple payouts, sometimes one can make inspired guesses for what $M(t)$ should be. For instance, in a one-dimensional Black-Scholes model, Rogers [2001] shows that using the numeraire-deflated European put option price (which is analytically known) as a guess for $M(t)$ generates good bounds for a Bermudan put option price. This approach, however, does not easily generalize to the CLE setting with its more complicated model and exercise payouts.

### 18.3.8.2 Nested Simulation (NS) Algorithm

Andersen and Broadie [2004] propose a general strategy for generating upper bounds, starting from any approximation $\tilde{\eta}$ to the optimal exercise strategy $\eta$. Typically, this approximation would originate from an LS regression, e.g. as in Section 18.3.6, or from an optimization of a parametric formulation of the exercise rule, as in Section 3.5.2 or Section 19.6.2. In a nutshell, the algorithm in Andersen and Broadie [2004] uses nested simulation - also known as "simulation within a simulation" and already mentioned in Section 18.3.7to construct an estimate of the low-bound value process $\widetilde{V}_{\text {CLE }}(t)$ generated from $\widetilde{\eta}$. The martingale component of the numeraire-deflated value of this process is then used as $M(t)$ in (18.31).

To outline the basic nested simulation (NS) algorithm in further detail, let us work on the exercise time line $\left\{T_{1}, T_{2}, \ldots, T_{N-1}\right\}$ and define

$$
\begin{align*}
M\left(T_{n+1}\right)-M\left(T_{n}\right) & =\frac{\tilde{V}_{\mathrm{CLE}}\left(T_{n+1}\right)}{B\left(T_{n+1}\right)}-\mathrm{E}_{T_{n}}\left(\frac{\tilde{V}_{\mathrm{CLE}}\left(T_{n+1}\right)}{B\left(T_{n+1}\right)}\right) \\
& =\mathrm{E}_{T_{n+1}}\left(\frac{U_{\tilde{\eta}_{n}}\left(T_{\tilde{\eta}_{n}}\right)}{B\left(T_{\tilde{\eta}_{n}}\right)}\right)-\mathrm{E}_{T_{n}}\left(\frac{U_{\tilde{\eta}_{n}}\left(T_{\tilde{\eta}_{n}}\right)}{B\left(T_{\tilde{\eta}_{n}}\right)}\right) \tag{18.35}
\end{align*}
$$
^eq18-35

with $M(0)=\widetilde{V}_{\text {CLE }}(0)$ and $M\left(T_{1}\right)=\widetilde{V}_{\text {CLE }}\left(T_{1}\right) / B\left(T_{1}\right)$. In the second equality of (18.35), we use $\widetilde{\eta}_{n}$ to denote the restriction of our approximate exercise policy exercise to the index set $\{n+1, \ldots, N\}$, such that $\widetilde{\eta}_{0}=\widetilde{\eta}$. For instance, if we use the algorithm in Section 18.3.6, we have (compare to (18.27))

$$
\begin{equation*}
\tilde{\eta}_{n}=\min \left\{k \geq n+1: U_{k}\left(T_{k}\right) \geq \mathcal{C}\left(\tilde{H}_{k}\left(T_{k}\right)\right)^{\top} \zeta\left(T_{k}\right)\right\} \wedge N . \tag{18.36}
\end{equation*}
$$
^eq18-36

For convenience, let us define an $\mathcal{F}_{T_{n+1}}$-measurable exercise indicator

$$
\iota\left(T_{n+1}\right)=1_{\left\{\tilde{\eta}_{n}=n+1\right\}}
$$

which will be one at time $T_{n+1}$ if our exercise policy indicates that the CLE should be exercised, and zero otherwise. We can also define hold values

$$
\frac{\tilde{H}_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}=\mathrm{E}_{T_{n}}\left(\frac{U\left(T_{\tilde{\eta}_{n}}\right)}{B\left(T_{\tilde{\eta}_{n}}\right)}\right)
$$

in which case (18.35) can be rewritten as

$$
\begin{align*}
M\left(T_{n+1}\right)-M\left(T_{n}\right)= & \iota\left(T_{n+1}\right) \frac{U_{n+1}\left(T_{n+1}\right)}{B\left(T_{n+1}\right)} \\
& +\left(1-\iota\left(T_{n+1}\right)\right) \frac{\widetilde{H}_{n+1}\left(T_{n+1}\right)}{B\left(T_{n+1}\right)}-\frac{\widetilde{H}_{n}\left(T_{n}\right)}{B\left(T_{n}\right)} \tag{18.37}
\end{align*}
$$
^eq18-37

Notice that

$$
\begin{aligned}
\frac{\widetilde{V}_{\mathrm{CLE}}\left(T_{n+1}\right)}{B\left(T_{n+1}\right)}-\frac{\widetilde{V}_{\mathrm{CLE}}\left(T_{n}\right)}{B\left(T_{n}\right)}= & \frac{\widetilde{V}_{\mathrm{CLE}}\left(T_{n+1}\right)}{B\left(T_{n+1}\right)}-\mathrm{E}_{T_{n}}\left(\frac{\widetilde{V}_{\mathrm{CLE}}\left(T_{n+1}\right)}{B\left(T_{n+1}\right)}\right) \\
& -\left(\frac{\widetilde{V}_{\mathrm{CLE}}\left(T_{n}\right)}{B\left(T_{n}\right)}-\mathrm{E}_{T_{n}}\left(\frac{\widetilde{V}_{\mathrm{CLE}}\left(T_{n+1}\right)}{B\left(T_{n+1}\right)}\right)\right) \\
= & M\left(T_{n+1}\right)-M\left(T_{n}\right)-\left(A\left(T_{n+1}\right)-A\left(T_{n}\right)\right)
\end{aligned}
$$

where we have denoted

$$
\begin{aligned}
A\left(T_{n+1}\right)-A\left(T_{n}\right) & =\frac{\tilde{V}_{\mathrm{CLE}}\left(T_{n}\right)}{B\left(T_{n}\right)}-\mathrm{E}_{T_{n}}\left(\frac{\tilde{V}_{\mathrm{CLE}}\left(T_{n+1}\right)}{B\left(T_{n+1}\right)}\right) \\
& =\iota\left(T_{n}\right)\left\{\frac{U_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}-\frac{\widetilde{H}_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}\right\}
\end{aligned}
$$

with $A(0)=0$. The second equality follows from the fact that $\widetilde{H}_{n}\left(T_{n}\right) / B\left(T_{n}\right)=\mathrm{E}_{T_{n}}\left(\widetilde{V}_{\mathrm{CLE}}\left(T_{n+1}\right) / B\left(T_{n+1}\right)\right)$. Therefore, we have the following decomposition

$$
\frac{\widetilde{V}_{\mathrm{CLE}}\left(T_{n}\right)}{B\left(T_{n}\right)}=M\left(T_{n}\right)-A\left(T_{n}\right)
$$

Notice that the process $A$ is not an increasing process (and $\tilde{V}_{\text {CLE }}$ therefore not a supermartingale), since we cannot guarantee that $\widetilde{H}_{n}\left(T_{n}\right) \geq U\left(T_{n}\right)$ : whenever an incorrect exercise decision is made, $A$ decreases.

For the purpose of computing an upper bound, the hold values $\widetilde{H}_{n}\left(T_{n}\right)$ and $\widetilde{H}_{n+1}\left(T_{n+1}\right)$ in (18.37) cannot be estimated by regression; doing so will introduce unknown biases which will destroy the martingale property of $M$ and, in turn, invalidate the inequality in (18.31). Instead, following Andersen and Broadie [2004] we can launch at times $T_{n}$ and $T_{n+1}$ Monte Carlo simulations to estimate the two expectations in (18.35). Notice that these "inner" Monte Carlo simulations will be nested inside a main "outer" simulation trial that generate sample paths of $U, B$, and $M$, as needed to estimate the expectation on the right-hand side of (18.31).

Now, we insert the martingale defined by (18.37) into the right-hand side of (18.31), which gives rise to a high-biased estimate $H_{0}^{\mathrm{hi}}(0)$ for the CLE value,

$$
\begin{equation*}
H_{0}^{\mathrm{hi}}(0)=\widetilde{H}_{0}(0)+\Delta \geq H_{0}(0) \tag{18.38}
\end{equation*}
$$
^eq18-38

where the duality gap $\Delta$ is defined as

$$
\begin{equation*}
\Delta=\mathrm{E}(D), \quad D \triangleq \max _{n=1, \ldots, N-1}\left(\frac{U_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}-M\left(T_{n}\right)\right) \tag{18.39}
\end{equation*}
$$
^eq18-39

The hold value $\widetilde{H}_{0}(0)=\widetilde{V}_{\text {CLE }}(0)$ can be estimated bias-free from the given exercise strategy $\widetilde{\eta}$ by standard Monte Carlo methods (see Section 18.3.6), so we focus on providing an estimate of the duality gap $\Delta$. The following NS algorithm can be used for this purpose.

1. Simulate $K_{U}$ paths $\omega_{1}, \ldots, \omega_{K_{I I}}$.
2. For each path $\omega_{k}$ calculate simulated values of the numeraire $B\left(T_{n}, \omega_{k}\right)$, $n=1, \ldots, N-1$.
3. For each path $\omega_{k}$, calculate the value $U_{n}\left(T_{n}, \omega_{k}\right)$ of the underlying exotic swap on all exercise dates $n=1, \ldots, N-1$.
4. For each path $\omega_{k}$ and each $T_{n}, n=1, \ldots, N-2$, launch $K_{\text {nest }}$ independent sub-paths $\left\{\omega_{k, 1}^{n}, \ldots, \omega_{k, K_{11} \cdots \text { st }}^{n}\right\}$ to time $T_{N}$ and estimate hold values $\widetilde{H}_{n}\left(T_{n}, \omega_{k}\right)$ as

$$
\widetilde{H}_{n}\left(T_{n}, \omega_{k}\right) \approx \widehat{H}_{n}\left(T_{n}, \omega_{k}\right) \triangleq \frac{1}{K_{\text {nest }}} \sum_{j=1}^{K_{\text {neat }}} \frac{U_{\widetilde{\gamma}(j, k, n)}\left(T_{\tilde{\gamma}(j, k, n)}\right)}{B\left(T_{\widetilde{\gamma}(j, k, n)}\right)}
$$

where, in slightly labored notation, $\widetilde{\gamma}(j, k, n)$ is the first exercise date for sub-path $j$, date $T_{n}$, and "outer" path $k$ :

$$
\widetilde{\gamma}(j, k, n)=\min \left\{l>n: \iota\left(T_{l}, \omega_{k, j}^{n}\right)=1\right\} \wedge N .
$$

5. For each path $\omega_{k}$, use $\widehat{H}_{n}\left(T_{n}, \omega_{k}\right), n=1, \ldots, N-1$, to form martingale estimates $\widehat{M}\left(T_{n}, \omega_{k}\right)$ by substituting $\widehat{H}_{n}\left(T_{n}, \omega_{k}\right)$ for $\widetilde{H}_{n}\left(T_{n}, \omega_{k}\right)$ in (18.37).
6. For each path $\omega_{k}$, compute pathwise duality gaps as

$$
\widehat{D}\left(\omega_{k}\right)=\max _{n=1, \ldots, N-1}\left(\frac{U_{n}\left(T_{n}, \omega_{k}\right)}{B\left(T_{n}, \omega_{k}\right)}-\widehat{M}\left(T_{n}, \omega_{k}\right)\right) .
$$

7. Estimate $\Delta$ as $\widehat{\Delta}$, where

$$
\widehat{\Delta}=\frac{1}{K_{U}} \sum_{k=1}^{K_{\prime \prime}} \widehat{D}\left(\omega_{k}\right) .
$$

### 18.3.8.3 Bias and Computational Cost of NS Algorithm

Having outlined the basic NS algorithm above, our first order of business is to establish formally that the estimator for $H_{0}^{\mathrm{hi}}(0)$ resulting from the NS algorithm is, in fact, biased high. Our primary concern here is the effect of the usage in Step 4 of nested Monte Carlo estimators $\widehat{H}_{n}\left(T_{n}, \omega_{k}\right)$ in place of the true hold values $H_{n}\left(T_{n}, \omega_{k}\right)$. The key result is the following.

Proposition 18.3.1. In the NS algorithm, the estimator for the duality gap is biased high, i.e.

$$
\mathrm{E}\left(\frac{1}{K_{U}} \sum_{k=1}^{K_{11}} \widehat{D}\left(\omega_{k}\right)\right) \geq \mathrm{E}(D) .
$$

Proof. We drop $\omega_{k}$ throughout the proof, such that $\widehat{H}_{n}\left(T_{n}\right)=\widehat{H}_{n}\left(T_{n}, \omega_{k}\right)$ and so forth. By construction, $\widehat{H}_{n}\left(T_{n}\right)$ is an unbiased estimator for $\widetilde{H}_{n}\left(T_{n}\right)$, i.e.

$$
\widehat{H}_{n}\left(T_{n}\right)=\widetilde{H}_{n}\left(T_{n}\right)+e_{n},
$$

where $e_{n}$ is a pure-noise error term with mean 0 and standard deviation proportional to $1 / \sqrt{K_{\text {nest }}}$. It follows from (18.37) that Step 5 in the NS algorithm will compute

$$
\begin{aligned}
& \widehat{M}\left(T_{n+1}\right)-\widehat{M}\left(T_{n}\right) \\
& \quad=\iota\left(T_{n+1}\right) \frac{U_{n+1}\left(T_{n+1}\right)}{B\left(T_{n+1}\right)} \\
& \quad+\left(1-\iota\left(T_{n+1}\right)\right) \frac{\widetilde{H}_{n+1}\left(T_{n+1}\right)+e_{n+1}}{B\left(T_{n+1}\right)}-\frac{\widetilde{H}_{n}\left(T_{n}\right)+e_{n}}{B\left(T_{n}\right)} \\
& \quad=M\left(T_{n+1}\right)-M\left(T_{n}\right)+\left(1-\iota\left(T_{n+1}\right)\right) \frac{e_{n+1}}{B\left(T_{n+1}\right)}-\frac{e_{n}}{B\left(T_{n}\right)} .
\end{aligned}
$$

By induction, it follows that, for $n=1, \ldots, N-1$,

$$
\widehat{M}\left(T_{n}\right)=M\left(T_{n}\right)+q_{n}
$$

where $q_{n}=q_{n}\left(\omega_{k}\right)$ is a random variable with zero mean (the explicit expression for $q_{n}$ is irrelevant) for all $n$. Hence,

$$
\widehat{D}=\max _{n=1, \ldots, N-1}\left(\frac{U_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}-M\left(T_{n}\right)-q_{n}\right)
$$

On the path $\omega_{k}$, let $\bar{n}\left(\omega_{k}\right)$ be the date index at which $U_{n}\left(T_{n}, \omega_{k}\right) / B\left(T_{n}, \omega_{k}\right)- M\left(T_{n}, \omega_{k}\right)$ attains it maximum. We can therefore write (again dropping $\omega_{k}$ 's)

$$
\begin{aligned}
& \mathrm{E}\left(\max _{n=1, \ldots, N-1}\left(\frac{U_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}-M\left(T_{n}\right)-q_{n}\right)\right) \\
& \quad \geq \mathrm{E}\left(\frac{U_{\bar{n}}\left(T_{\bar{n}}\right)}{B\left(T_{\bar{n}}\right)}-M\left(T_{\bar{n}}\right)-q_{\bar{n}}\right) \\
& \quad=\mathrm{E}\left(\frac{U_{\bar{n}}\left(T_{\bar{n}}\right)}{B\left(T_{\bar{n}}\right)}-M\left(T_{\bar{n}}\right)\right) \\
& \quad=\mathrm{E}\left(\max _{n=1, \ldots, N-1}\left(\frac{U_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}-M\left(T_{n}\right)\right)\right)
\end{aligned}
$$

where the first equality follows from the zero mean of $q_{\bar{n}}$.
Proposition 18.3.1 demonstrates that our estimate of the duality gap is biased high for finite values of the sub-path sample size $K_{\text {nest }}$. As such, our finite sample estimate for $H_{0}^{\mathrm{hi}}(0)$ will itself be biased high and have a mean that is above the true CLE value, as desired. By increasing $K_{\text {nest }}$ we can reduce the bias originating from the finite sample size of sub-paths, and thereby tighten the upper bound. Of course, even in the limit $K_{\text {nest }} \rightarrow \infty$ we will still produce an estimator that is biased high; the size of this bias will reflect the quality of our exercise strategy choice, and will only vanish in the (unlikely) event that we manage to use precisely the optimal exercise strategy, i.e. when $\widetilde{\eta}=\eta$.

While the basic NS algorithm is quite straightforward, the need for nested simulations makes it numerically expensive: the worst-case workload will be proportional to

$$
\begin{equation*}
K_{\text {nest }} \cdot K_{U} \cdot N^{2} \tag{18.40}
\end{equation*}
$$
^eq18-40

For comparison, with $K^{\prime}$ simulation trials the lower bound simulation in Section 18.3.6 has a workload proportional to $K^{\prime} \cdot N$, plus the work required to estimate the exercise rule in a pre-simulation. In many cases the inner simulations of the NS algorithm can be stopped quickly (due to exercise of the CLE), so in practice the dependence on $N$ in (18.40) is often less than quadratic and sometimes close to linear. Finally, $K_{\text {nest }}$ can often be set to a number much smaller than $K_{U}$ without significantly affecting the quality of the upper bound, and even very small values of $K_{\text {nest }}$ (e.g., 50-100 or less) may yield informative results. If one additionally takes advantage of the
algorithm refinements discussed in Section 18.3.8.6, it is typically possible to execute an upper bound computation for a long-dated (say, 30 years) CLE in a few minutes. This is still relatively time-consuming, so upper bound computations are often most useful in practice as a way to test the quality of a postulated exercise strategy. We expand on this topic in the next section.

### 18.3.8.4 Confidence Intervals and Practical Usage

For concreteness, assume that the algorithm in Section 18.3.6 has been used to compute a $K^{\prime}$-path lower bound estimate of $V_{\text {CLE }}(0)=H_{0}(0)$. Let this estimate be denoted $\widehat{V}_{\text {CLE }}$, and let its recorded sample standard deviation be $\widehat{s}_{L}$. Also assume that an independent simulation of the NS algorithm with $K_{U}$ outer simulation trials has produced an estimate $\widehat{\Delta}$ for the duality gap, with sample standard deviation $\widehat{s}_{U}$. Asymptotically, a $100(1-\gamma) \%$ confidence interval for the true price $V_{\text {CLE }}(0)$ must be tighter than

$$
\begin{equation*}
\left[\widehat{V}_{\mathrm{CLE}}-u_{\gamma / 2} \frac{\widehat{s}_{L}}{\sqrt{K^{\prime}}} ; \widehat{V}_{\mathrm{CLE}}+\widehat{\Delta}+u_{\gamma / 2} \sqrt{\frac{\widehat{s}_{U}^{2}}{K_{U}}+\frac{\widehat{s}_{L}^{2}}{K^{\prime}}}\right], \tag{18.41}
\end{equation*}
$$
^eq18-41

where $\Phi\left(u_{\gamma / 2}\right)=1-\gamma / 2$. As already mentioned in Section 3.5.6, the confidence interval is conservative ${ }^{11}$ because of the low bias in the sample estimator $\widehat{V}_{\text {CLE }}$ (i.e., $\mathrm{E}\left(\widehat{V}_{\text {CLE }}\right) \leq H_{0}(0)$ ) and the high bias in $\widehat{V}_{\text {CLE }}+\widehat{\Delta}$, which originates in part from the nature of the upper bound, and in part from the earlier mentioned additional high bias introduced by the finite sample size of the inner simulations (see Proposition 18.3.1).

As noted in Section 3.5.6, it is not uncommon that the upper and lower bounds for the option price often are roughly symmetric around the true value, so in the event that we have computed both bounds, the obvious point estimate

$$
\begin{equation*}
\widehat{V}_{\mathrm{CLE}}+\frac{1}{2} \widehat{\Delta} \tag{18.42}
\end{equation*}
$$
^eq18-42

will often give better price estimates than either the upper or lower bound alone.

Upper bound simulation algorithms can typically be expected to be both more involved and/or more expensive than lower bound simulation methods. In many cases, the best use of the upper bound simulation algorithm will therefore be to test whether postulated lower bound exercise strategies are tight or not. Specifically, starting from some guess for the exercise strategy, we can produce confidence intervals using (18.41) to test whether the lower

[^34]bound estimate is of good quality, in which case the confidence interval can be made tight by using large values of $K_{U}$ and $K^{\prime}$ (as well as the number of inner simulation trials, $K_{\text {nest }}$ ). In case the lower bound estimator is deemed unsatisfactory, we can iteratively refine it, by altering the choice of basis functions, say, until the confidence interval is tight. Importantly, such tests can often be done at a high level, covering entire classes of payouts and/or models. Once an exercise strategy has been validated for a particular product or model, day-to-day pricing of callable securities can be done by the lower bound method, with only occasional runs of the upper bound method needed (e.g., if market conditions change markedly). If upper bound methods are predominantly used in this fashion, the fact that they may sometimes be computationally intensive ${ }^{12}$ becomes less punitive.

### 18.3.8.5 Non-Analytic Exercise Values

The observant reader might have noticed that the NS algorithm outlined in Section 18.3.8.2 assumes (in Steps 3, 4, and 6) that exercise values are directly computable at each maturity and each state of the world. A similar assumption was made in the basic LS regression algorithm of Section 18.3.1, but relaxed later, in Section 18.3.2. To handle CLEs that involve exotic swap underlyings for which the values are not easily computable, we can modify the upper bound simulation algorithm in at least two ways.

Our first approach is straightforward, and based on the representations

$$
\frac{U_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}=\mathrm{E}_{T_{n}}\left(\sum_{i=n}^{N-1} \frac{X_{\imath}}{B\left(T_{i+1}\right)}\right)
$$

and

$$
\begin{align*}
\frac{\widetilde{H}_{n}\left(T_{n}\right)}{B\left(T_{n}\right)} & =\mathrm{E}_{T_{n}}\left(\frac{U_{\tilde{\eta}_{n}}\left(T_{\tilde{\eta}_{n}}\right)}{B\left(T_{\tilde{\eta}_{n}}\right)}\right)=\mathrm{E}_{T_{n}}\left(\frac{\mathrm{E}_{T_{\tilde{\eta}_{n}}}\left(\sum_{i=\tilde{\eta}_{n}}^{N-1} \frac{X_{i}}{B\left(T_{i+1}\right)}\right)}{B\left(T_{\tilde{\eta}_{n}}\right)}\right) \\
& =\mathrm{E}_{T_{n}}\left(\frac{\sum_{i=\tilde{\eta}_{n}}^{N-1} \frac{X_{i}}{B\left(T_{i+1}\right)}}{B\left(T_{\tilde{\eta}_{n}}\right)}\right) \tag{18.43}
\end{align*}
$$
^eq18-43

where the last equality follows from the optional sampling theorem. The relevant expectations can be computed bias-free by launching a nested simulation at time $T_{n}$, generating sample paths from time $T_{n}$ to $T_{N-1}$. For instance, for the path $\omega_{k}$ and date $T_{n}$, we would write[^35]

$$
U_{n}\left(T_{n}, \omega_{k}\right) \approx \widehat{U}_{n}\left(T_{n}, \omega_{k}\right)=\frac{1}{K_{\text {nest }}} \sum_{j=1}^{K_{\text {nest }}} \sum_{i=n}^{N-1} \frac{X_{i}\left(\omega_{k, j}^{n}\right)}{B\left(T_{i+1}, \omega_{k, j}^{n}\right)},
$$
where $\left\{\omega_{k, 1}^{n}, \ldots, \omega_{k, K_{\text {nest }}}^{n}\right\}$ is the set of "inner" sub-paths spawned at time $T_{n}$ for the "outer" path $\omega_{k}$. An estimator $\widehat{H}_{n}\left(T_{n}, \omega_{k}\right)$ for $\widetilde{H}_{n}\left(T_{n}, \omega_{k}\right)$ may be computed from (18.43) the same way, except that only net coupons $X_{i}\left(\omega_{k, j}^{n}\right)$ after exercise will be counted in the accumulation of deflated coupons. With these estimators used in Steps 4 and 6, the NS algorithm in Section 18.3.8.2 may proceed as before. As the nested simulation will produce bias-free (but noisy) estimates for the true exercise and hold values, it follows from Proposition 18.3.1 that the resulting upper bound estimator will still be guaranteed to be biased high.

While using nested simulation in the manner described above does not add to the order of the computational complexity of the upper bound algorithm (it remains as in (18.40)), the need to construct exercise values by simulating coupon streams will obviously add additional noise to the basic algorithm, which in turn will increase the finite sample bias and also widen the confidence interval (18.41) for a given computational budget.

In our second approach to dealing with non-analytic exercise values, we follow the alternative route also taken in Section 18.3.3, and focus on pricing the cancelable note $G_{0}(0)$, where

$$
G_{0}(0)=H_{0}(0)-U_{0}(0)
$$

As a starting point, consider the expression (18.21), which we may rewrite as

$$
\begin{align*}
G_{0}(0) & =\sup _{\xi \in \mathcal{T}_{0}} \mathrm{E}\left(-\sum_{n=0}^{\xi-1} B\left(T_{n+1}\right)^{-1} X_{n}\right) \\
& =\sup _{\xi \in \mathcal{T}_{0}} \mathrm{E}\left(B\left(T_{\xi}\right)^{-1} J\left(T_{\xi}\right)\right), \quad J\left(T_{\xi}\right) \triangleq-\sum_{n=0}^{\xi-1} \frac{B\left(T_{\xi}\right)}{B\left(T_{n+1}\right)} X_{n} \tag{18.44}
\end{align*}
$$
^eq18-44

Financially, the quantity $J\left(T_{\xi}\right)$ can be interpreted as the payout at the exercise date $T_{\xi}$ from re-investing all pre-exercise coupons into the numeraire $B$. Effectively, this formulation removes all pre-exercise cash flows from the expectation for $G_{0}(0)$, making (18.44) structurally identical to (18.30), but with the exercise value $U_{\xi}\left(T_{\xi}\right)$ replaced by $J\left(T_{\xi}\right)$. Of course, while $U_{\xi}\left(T_{\xi}\right)$ may be difficult to compute at time $T_{\xi}$ in a Monte Carlo simulation, $J\left(T_{\xi}\right)$ is not.

To construct an upper bound for $G_{0}(0)$, we follow the same principles that lead to (18.35) and construct a martingale $M$ as

$$
\begin{equation*}
M\left(T_{n+1}\right)-M\left(T_{n}\right)=\mathrm{E}_{T_{n+1}}\left(\frac{J\left(T_{\tilde{\eta}_{n}}\right)}{B\left(T_{\tilde{\eta}_{n}}\right)}\right)-\mathrm{E}_{T_{n}}\left(\frac{J\left(T_{\tilde{\eta}_{n}}\right)}{B\left(T_{\tilde{\eta}_{n}}\right)}\right) \tag{18.45}
\end{equation*}
$$
^eq18-45

where $\widetilde{\eta}$ is a given exercise strategy. With $\widetilde{G}_{0}(0)$ being the lower bound value for $G_{0}(0)$ computed by using $\widetilde{\eta}$ as the exercise strategy, an upper bound for $G_{0}(0)$ is (compare to (18.38))

$$
G_{0}^{\mathrm{hi}}(0)=\widetilde{G}_{0}(0)+\Delta_{G} \geq G_{0}(0)
$$

where the duality gap $\Delta_{G}$ is defined as

$$
\Delta_{G}=\mathrm{E}\left(D_{G}\right), \quad D_{G} \triangleq \max _{n=1, \ldots, N-1}\left(\frac{J\left(T_{n}\right)}{B\left(T_{n}\right)}-M\left(T_{n}\right)\right) .
$$

When using Monte Carlo simulation to estimate $\Delta_{G}$, we can use nested simulation to establish a bias-free estimator for the martingale in (18.45), in the same way as was done for the NS algorithm. By the arguments in Proposition 18.3.1, the resulting estimator for $\Delta_{G}$ will be biased high, i.e. our upper bound is valid. Confidence intervals for $G_{0}(0)$ (and for $H_{0}(0)= \left.G_{0}(0)+U_{0}(0)\right)$ can be constructed using the principles of Section 18.3.8.4.

### 18.3.8.6 Improvements to NS Algorithm

In Broadie and Cao [2008], the authors outline a number of improvements to both upper and lower bound simulations. As some of the proposed techniques are fairly involved, we cannot give this paper full justice, but contend ourselves with listing one relatively straightforward trick from Broadie and Cao [2008]. The reader interested in additional techniques should consult Broadie and Cao [2008] directly.

Let us return to the setting of Section 18.3.8.2, and assume that an easily computable lower limit $\underline{H}\left(T_{n}\right)$ exists for the hold value at time $T_{n}$,

$$
\begin{equation*}
H_{n}\left(T_{n}\right) \geq \underline{H}\left(T_{n}\right), \quad n=1, \ldots, N-1 . \tag{18.46}
\end{equation*}
$$
^eq18-46

We comment on how to choose $\underline{H}(\cdot)$ later. Let us also assume that the inequality (18.46) is honored by the given exercise policy $\widetilde{\eta}$, i.e. we assume that exercise will never take place when $U_{n}\left(T_{n}\right)<\underline{H}\left(T_{n}\right)$. When using a regression approach, we can ensure that this assumption is true by simply writing (compare to (18.36))

$$
\begin{equation*}
\tilde{\eta}_{n}=\min \left\{k \geq n+1: U_{k}\left(T_{k}\right) \geq \max \left(\underline{H}\left(T_{k}\right), \mathcal{C}\left(\tilde{H}_{k}\left(T_{k}\right)\right)^{\top} \zeta\left(T_{k}\right)\right)\right\} \wedge N . \tag{18.47}
\end{equation*}
$$
^eq18-47

Modifying an exercise policy to accommodate bounds such as (18.46) is sometimes known as policy fixing, see Broadie and Glasserman [2004].

Before stating our next result, let us make the additional assumption that $\underline{H}(t) / B(t)$ is a submartingale in measure $\mathrm{Q}^{B}$ for $t \leq T_{N}$; this is, for instance, the case when $\underline{H}(t)$ is chosen to be either zero or the price of an asset that pays no cash flows before $T_{N}$. With this assumption, we have

$$
\begin{equation*}
\frac{\tilde{H}_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}=\mathrm{E}_{T_{n}}\left(\frac{U_{\tilde{\eta}_{n}}\left(T_{\tilde{\eta}_{n}}\right)}{B\left(T_{\tilde{\eta}_{n}}\right)}\right) \geq \mathrm{E}_{T_{n}}\left(\frac{\underline{H}\left(T_{\tilde{\eta}_{n}}\right)}{B\left(T_{\tilde{\eta}_{n}}\right)}\right) \geq \frac{\underline{H}\left(T_{n}\right)}{B\left(T_{n}\right)} . \tag{18.48}
\end{equation*}
$$
^eq18-48

We use this result to show the following proposition, adapted from Broadie and Cao [2008].

Proposition 18.3.2. Let the martingale $M(t)$ be defined as in (18.35) and (18.37), and assume that $\underline{H}(t) / B(t)$ is a submartingale. Assume also that $U_{k}\left(T_{k}\right)<\underline{H}\left(T_{k}\right)$ for all $k=l, \ldots, n$, with $1 \leq l \leq n$. Then we have,

$$
M\left(T_{n}\right)=M\left(T_{l-1}\right)-\tilde{H}_{l-1}\left(T_{l-1}\right) / B\left(T_{l-1}\right)+\tilde{H}_{n}\left(T_{n}\right) / B\left(T_{n}\right)
$$

and

$$
\begin{equation*}
U_{n}\left(T_{n}\right) / B\left(T_{n}\right)-M\left(T_{n}\right)<\tilde{H}_{l-1}\left(T_{l-1}\right) / B\left(T_{l-1}\right)-M\left(T_{l-1}\right) . \tag{18.49}
\end{equation*}
$$
^eq18-49

In particular, if $l=1$ then

$$
\begin{equation*}
M\left(T_{n}\right)=\tilde{H}_{n}\left(T_{n}\right) / B\left(T_{n}\right), \quad U_{n}\left(T_{n}\right) / B\left(T_{n}\right)<M\left(T_{n}\right) . \tag{18.50}
\end{equation*}
$$
^eq18-50

Proof. First, we notice that exercise will never take place on the interval $\left[T_{l}, T_{n}\right]$ such that, from (18.37),

$$
\begin{aligned}
M\left(T_{n}\right) & =M\left(T_{n-1}\right)+\frac{\widetilde{H}_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}-\frac{\widetilde{H}_{n-1}\left(T_{n-1}\right)}{B\left(T_{n-1}\right)} \\
& =M\left(T_{l-1}\right)+\sum_{j=l}^{n}\left(\frac{\widetilde{H}_{j}\left(T_{j}\right)}{B\left(T_{j}\right)}-\frac{\widetilde{H}_{j-1}\left(T_{j-1}\right)}{B\left(T_{j-1}\right)}\right) \\
& =M\left(T_{l-1}\right)-\frac{\widetilde{H}_{l-1}\left(T_{l-1}\right)}{B\left(T_{l-1}\right)}+\frac{\widetilde{H}_{n}\left(T_{n}\right)}{B\left(T_{n}\right)} .
\end{aligned}
$$

Using (18.48), we also have

$$
\begin{equation*}
\frac{\widetilde{H}_{n}\left(T_{n}\right)}{B\left(T_{n}\right)} \geq \frac{\underline{H}\left(T_{n}\right)}{B\left(T_{n}\right)}>\frac{U_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}, \tag{18.51}
\end{equation*}
$$
^eq18-51

which proves (18.49). The results for the special case where $l=1$ follows from the fact that $M(0)=\widetilde{H}_{0}\left(T_{0}\right)$.

Recalling (18.39), it is clear from (18.50) that there is no contribution to the upper bound increment $D$ before the underlying process enters the region where $U_{n}\left(T_{n}\right) \geq \underline{H}\left(T_{n}\right)$. Moreover, whenever the option is inside a region where the exercise value is less than the lower limit $\underline{H}, M I$ does not depend on the actual path taken inside the region. This allows for the following straightforward modification to the NS algorithm: whenever at time $T_{n}$ we observe that $U_{n}\left(T_{n}\right)<\underline{H}\left(T_{n}\right)$, we simply skip the nested simulations and proceed to the next date $T_{n+1}$; otherwise we launch $K_{\text {nest }}$ sub-simulations
and update $M\left(T_{n}\right)$ according to Proposition 18.3.2. When computing $\widehat{D}$ in Step 6, we can ignore all dates where $U_{n}\left(T_{n}\right)<\underline{H}\left(T_{n}\right)$. For options that are deeply out-of-the money, a substantial number of nested simulations can be avoided in this fashion, leading to significant improvements in computational performance.

The refinement in the Broadie and Cao [2008] algorithm hinges upon our ability to establish a meaningful lower limit submartingale $\underline{H}$. While this must generally be done on a case-by-case basis, one obvious possibility in the standard CLE setting is to set the lower limit to zero. This is a relatively coarse bound, and it may be tempting to sharpen it by observing that, for any $t \geq 0$ and $n, m \geq 0$,

$$
H_{n}(t) \geq H_{n+m}(t) \geq U_{n+m+1}(t)
$$

which, assuming we can calculate $U_{n+m+1}(t)$ 's analytically (admittedly a strong assumption), shows that a lower limit can be computed as

$$
\underline{H}\left(T_{n}\right)=\max _{k=n+1, \ldots, N-1}\left(U_{k}\left(T_{n}\right)\right) .
$$

However, while this choice of the lower limit ${ }^{13} \underline{H}(\cdot)$ can certainly be used in policy fixing to improve the lower bound, it is generally not a submartingale and therefore cannot be used in upper bound computations.

### 18.3.8.7 Other Upper Bound Algorithms

To avoid the need for nested simulation, it is tempting to return to the representation in (18.34) and contemplate whether one can estimate the optimal choice of $\sigma(t)$ by extracting it from the empirical volatility of (the martingale component of) $\widetilde{V}_{\text {CLE }}(t) / B(t)$. One advantage to this approach is that any errors in the estimation of $\sigma(t)$ will not affect the martingale property of $M(t)$ in (18.34). Starting again from a postulated exercise strategy $\widetilde{\eta}$, Belomestny et al. [2007] use this observation to construct a regression on a set of basis functions to produce an estimate for the function $\sigma(t)$. By applying regression techniques this way, the authors are able to construct a true martingale process $M(t)$, which can be turned into a valid upper bound through (18.31). While the resulting algorithm involves no nested simulation, it requires considerable care in its implementation, in part because the optimal integrand $\sigma(t)$ can be expected to be considerably less regular than the optimal martingale $M(t)$ itself. This, in turn, requires additional thought in the selection of appropriate basis functions for the regression. One possibility advocated in Belomestny et al. [2007] is to include, whenever available, exact or approximate expressions for the diffusion term

[^36]in dynamics of several still-alive European options underlying the Bermudan option. This strategy is akin to that of Rogers [2001], and its feasibility depends on the pricing problem at hand. In cases where it does apply, the authors of Belomestny et al. [2007] demonstrate that their method gives good results, with the upper bound often being nearly as tight as that of the nested algorithm in Andersen and Broadie [2004]. They also show how to use their technique to develop a variance-reduced version of the algorithm in Andersen and Broadie [2004].

Additional techniques for computing upper bounds can be found in Broadie and Glasserman [1997] and Glasserman and Yu [2002]. The latter is based purely on regression, but requires strong conditions on regression basis functions, that may be hard to check in practice.

### 18.3.9 Regression Variable Choice

The single most critical determinant of the performance of the regressionbased valuation methods is the choice of the regression variables ${ }^{14} \zeta(t)= \left(\zeta_{1}(t), \ldots, \zeta_{q}(t)\right)^{\top}$. Recall that the values of these variables at time $T$ partly ${ }^{15}$ serve to approximate the information contained in the sigma-algebra $\mathcal{F}_{T}$. Some of this information is relevant to the valuation of a given security and some is not. The closer $\zeta(T)$ approximates the information in $\mathcal{F}_{T}$ that is relevant for the security, the better the regression method performs, i.e. the smaller is the bias of lower bound estimates of the security value.

### 18.3.9.1 State Variables Approach

For Markovian models, all information observed at (but not before) time $T$ is encoded in the Markovian state variables, say some ( $d$-dimensional) vector $x(T)$. As discussed earlier in Section 3.5.4, for such models it is often natural take the $\zeta(T)$ 's to be deterministic functions of the state variables; these functions should approximately span the set of all functions of the state variables. A good choice is to use monomials of $x(T)$, i.e. functions of the type $\prod_{i=1}^{d} x_{i}(T)^{p_{i}}$, which corresponds to using a polynomial basis in the LS regression. One caveat applies: the filtration $\sigma(x(T))$ generated by the state variables at $x(T)$ is clearly smaller than $\mathcal{F}_{T}$, since the former consists of model information observed at time $T$ only, while the latter contains information for all times from 0 to $T$. For some derivatives the reduction from $\mathcal{F}_{T}$ to $\sigma(x(T))$ is irrelevant - this includes all callable securities whose coupons are not path-dependent. On the other hand, for securities such

[^37]as callable snowballs (Section 5.14.4) whose coupons depend on the whole history of interest rate evolution, the regression variables must be augmented with some that carry information from the past. To do this in a productindependent way, we can, for instance, include state variables from previous times in $\zeta(T)$.

Markovian interest rate models that would allow for such essentially product-independent choice of regression variables include quasi-Gaussian models and quadratic Gaussian models, at least if the number of Markov state variables is not too high. While Libor market (LM) models are Markovian in the set of all Libor forward rates (and any stochastic volatility variables), the dimensionality of the state vector is typically so large that the regression will suffer from numerical problems, an issue we shall discuss further later. One general way to reduce the dimension of the state variable vector is to perform a principal components analysis. To demonstrate this idea in a LM model setting, suppose that we wish to synthesize $d$ "state variables" by forming the first $d$ principal components of the vector of Libor rates. At a given time $T_{n}$ we have a vector of (centered) still-alive forward Libor rates

$$
A=\left(L_{n}\left(T_{n}\right)-\mathrm{E}\left(L_{n}\left(T_{n}\right)\right), \ldots, L_{N-1}\left(T_{n}\right)-\mathrm{E}\left(L_{N-1}\left(T_{n}\right)\right)\right)^{\top} .
$$

We assume $N-n \geq d$. The term covariance matrix of the rates,

$$
c=\mathrm{E}\left(A A^{\top}\right)
$$

can be estimated using methods from Chapter 14. Then, by principal components analysis (see Section 3.1.3) we can find an $(N-n) \times d$ matrix $D$ such that $D D^{\top}$ is the closest (in the Frobenius norm) rank- $d$ approximation to $c$. Then

$$
A \approx D x
$$

for some $d$-dimensional vector $x=x\left(T_{n}\right)$, which we recover from the leastsquares (regression) problem

$$
\|A-D x\|^{2} \rightarrow \min ,
$$

with the solution

$$
x=\left(D^{\top} D\right)^{-1} D^{\top} A
$$

Doing this for each $T_{n}$ gives us a set of $d$ approximate state variables for regression on each exercise date.

### 18.3.9.2 Explanatory Variables

While the state variable approach (with or without principal components dimension reduction) is appealing in its independence of security specifics, it has several shortcomings. First, as mentioned earlier, situations may arise where the information carried in the state variables is inadequate for the
security in question (e.g. snowballs or other path-dependent CLEs). Second, there may be cases where there is too much information in the state variables, leading to an overabundance of regression variables. For instance, a standard Bermudan swaption turns out (see Section 19.6.1) to primarily depend on the overall level of interest rates, so if, say, the number $d$ of PCA state variables used in an LM model is much larger than 1 , many of the regression variables will have little or no explanatory power.

Having too many regression variables not only adds more work to the numerical scheme, it ultimately tends to reduce the quality of the CLE price estimate. Indeed, for a finite budget of simulated paths used in regression, having too many regression variables ${ }^{16}$ will induce errors in the regression coefficient estimates and, therefore, in the exercise rule and, ultimately, in the value of the callable security. Moreover, if some of the regression variables add no explanatory power to the regression, their inclusion in the regression may lead to spurious noise and further issues with the estimation quality of the exercise rule ${ }^{17}$.

In light of the comments above, in practice it is hard to avoid a careful analysis of each security type, to ensure that neither too little or too much information is contained within the set of regression variables. In such a "product-specific" approach to choosing regression variables, we would aim to choose regression variables in a way that maximizes their explanatory power for the particular security in question.

Similar to the state variables approach, we find it convenient to specify regression variables $\zeta(t)$ as simple functions of so-called explanatory variables, which can be thought of as product-specific analogs to state variables of the model. Let us use $x(t)=\left(x_{1}(t), \ldots, x_{d}(t)\right)^{\top}$ to denote these explanatory variables; the fact that we recycle the notation for $x$ and $d$ from Section 18.3.9.1 should not lead to any confusion. When constructing regression variables from explanatory variables, we typically choose monomials of explanatory variables

$$
\begin{equation*}
\prod_{i=1}^{d} x_{i}(t)^{p_{i}}, \tag{18.52}
\end{equation*}
$$
^eq18-52

up to a low order $r$, say 3 to 4 , so that $\sum_{i=1}^{d} p_{i} \leq r$. This, of course, amounts to fitting exercise and hold values with polynomials of degree $r$ in explanatory variables. It is worth noting that while higher-order polynomials could give a closer fit, they could also lead to unexpected behavior outside of the range of the points used in fitting. As exercise boundaries may lie far away

[^38]from typical simulation scenarios, this could lead to poor estimation of the exercise rule.

When selecting explanatory variables, we should look for variables that have high explanatory power in the regression of hold and exercise values of a given CLE. Inevitably, this process is trade-type specific and combines elements of both science and art. While trial and error is always needed, we can make a few recommendations. Generally, we always prefer using financially meaningful explanatory variables such as various market rates or values of (simple) market instruments, and find it convenient to distinguish three classes of potential explanatory variables. The first class contains the market variables that directly drive the values of coupons. For example, if a coupon is linked to a CMS spread between two rates, we should include the spread as a potential explanatory variable. The second class contains those variables that describe the past; these variables are only relevant for CLEs with path-dependent coupons. For instance, recall that in a snowball the value of a current coupon depends on the values of the previous coupons, so we would recommend including the current coupon in the set of potential explanatory variables. Finally, the last class contains those market variables that are thought to drive the exercise decision. Inspired by the case of a simple Bermudan swaption, we would often include a variable responsible for the overall level of the yield curve. In particular, for the regression on date $T_{n}$ we would typically include a swap rate that fixes at $T_{n}$ and spans all periods to the last exercise (a so-called "core" swap rate). We might also include a variable that reflects the slope of the yield curve on the exercise date; a front Libor rate, i.e. a Libor rate that fixes on the exercise date, is a good choice for this.

After collecting potential explanatory variables from all three classes, we would typically proceed to analyze the list and try to reduce it to a manageable number of variable, such as 2 to 4 , that we would then use in our regression method. While nothing replaces careful analysis, the selection process can be automated somewhat and done during the actual regression - we discuss this in more details later on, in Section 18.3.10.1.

Whatever variables we choose, we should be careful to always choose variables for regression on $T_{n}$ that are $\mathcal{F}_{T_{n}}$-measurable. In plain speak, variables should not be "future-looking", but should be computable by using only the state of the model as observed up to and at time $T_{n}$. This requirement, while seemingly technical, is critical to the success of the regression algorithm, as one must not be allowed to "see into the future" when making decisions about exercise.

In selecting explanatory variables we should be thinking about qualitative impact of various variables on exercise/hold values but, fortunately, we do not need to be quantitatively exact in capturing the effects. As long as the general influence is accounted for correctly, the fitting of parametric functions will generally take care of choosing the best scalings and/or linear
combinations of explanatory variables ${ }^{18}$. For example, one may decide (as we often recommend) to include in the set of explanatory variables the level of the yield curve, as measured by a core swap rate, as well as the slope of the yield curve, as measured by the difference between the swap rate and some short-tenor Libor rate. For the latter variable, one does not need to explicitly include the difference of the rates as an explanatory variable, since the shorttenor Libor rate will work just as well - the properly weighted difference of the two rates will implicitly be used in the fitting of the polynomials.

### 18.3.9.3 Explanatory Variables with Convexity

Sometimes the quality of exercise boundary estimation is improved by such seemingly minor changes as using core swap values in the set of explanatory variables, rather than core swap rates. In practical applications, we always recommend trying both and seeing which one is better. The reasons for this effect are not always entirely clear, but originate with subtle convexity differences between the two choices. It turns out that convexity, i.e. nonlinear dependence between exercise values and the explanatory variables, can be quite important. In particular, instead of using simple rates such as a core swap rate or a front Libor rate, we sometimes find it useful to use functions of them as explanatory variables. While this may seem superfluous - after all we will ultimately be applying polynomial functions to explanatory variables to get our regression variables - it turns out that using functional mappings more finely tuned to the features of the trade is often beneficial. The effect is due to the fact that for many CLEs the underlying swap consists of coupons that are options on rates. By matching, roughly, the resulting convexity in value of coupons relative to underlying rates, we can often improve the regression fit. Sometimes this can be accomplished by simply using $\left(S\left(T_{n}\right)-K\right)^{+}$, instead of some rate $S\left(T_{n}\right)$, as an explanatory variable for a CLE whose coupons are strike- $K$ options on the rate $S\left(T_{n}\right)$ (or some other, relatively similar, rate). For a more refined approach we could try to roughly estimate the option value of the remaining coupons as of the exercise date $T_{n}$. Even if coupons do not have any optionality, such as for Bermudan swaptions, we may wish to use European options on the core swap rates as explanatory variables to better fit hold values, which are always convex due to callability.

If one uses a model that allows for closed-form valuation of European options on rates at any future time, then, for suitably simple CLEs, it may be possible to use the exact exercise value as an explanatory variable. This, however, only works for a subset of models, and the approach can significantly impair the speed of valuation, since one needs to compute European option prices a large number of times (one per path per exercise time per underlying

[^39]option). Fortunately, it is not necessary to be particularly precise in matching the curvature of the exercise values with the explanatory variables, and a rough estimate will typically do just fine. For example one can use the Black formula with an approximate volatility to value the options in the underlying, even in non-Black models. Also, if the underlying is a strip of options, one does not need to value all options in the underlying, but just one, e.g. the one in the middle of the strip.

To give an example of an implementation using option-like explanatory variables, consider a callable capped floater. The structured coupon at time $T_{n}$ is here given by

$$
C_{n}=\min \left(L_{n}\left(T_{n}\right)+s, c\right)
$$

received against a Libor rate payment,

$$
X_{n}=\tau_{n} \times\left(C_{n}-L_{n}\left(T_{n}\right)\right)
$$

We note that the net coupon $X_{n}$ can be written

$$
\begin{aligned}
X_{n} & =\tau_{n} \times\left(\min \left(L_{n}\left(T_{n}\right)+s, c\right)-L_{n}\left(T_{n}\right)\right) \\
& =\tau_{n} \times s-\tau_{n} \times\left(L_{n}\left(T_{n}\right)-(c-s)\right)^{+}
\end{aligned}
$$

so the payment is a combination of a fixed rate payment with rate $s$ and a call option on the Libor rate with strike $c-s$.

For the first explanatory variable we may use an approximate value of the exercise value on each exercise date

$$
x_{1}\left(T_{n}\right)=\sum_{i=n}^{N-1} \tau_{i} \times P\left(T_{n}, T_{\imath+1}\right) \times\left(s-c_{\mathrm{B}}\left(T_{n}, L_{i}\left(T_{n}\right) ; T_{i}, c-s ; \sigma_{n, i}\right)\right)
$$

where the Black formula $c_{\mathrm{B}}(\ldots)$ is defined in Remark 7.2.8. The volatility $\sigma_{n, i}$ is the log-normal volatility of the Libor rate $L_{i}(t)$ over the interval $\left[T_{n}, T_{i}\right]$ as given by the model. To reiterate, very crude approximations can be used here. For the displaced log-normal Libor market model (see (14.12)) one may use $\sigma_{n, i}=\bar{\lambda}$, where $\bar{\lambda}$ is some average of relative forward Libor model volatilities over time and Libor rates, and for the stochastic volatility version (see (14.16)) one may use

$$
\begin{equation*}
\sigma_{n, i}=\bar{\lambda} \sqrt{z\left(T_{n}\right)} \tag{18.53}
\end{equation*}
$$
^eq18-53

For the second explanatory variable, we can use a front Libor rate fixing on an exercise date, as we recommended before. Alternatively, we can use just the approximate value of the first coupon payment,

$$
x_{2}\left(T_{n}\right)=\tau_{n} \times P\left(T_{n}, T_{n+1}\right) \times\left(s-c_{\mathrm{B}}\left(T_{n}, L_{n}\left(T_{n}\right), c-s ; T_{n} ; \sigma_{n, n}\right)\right)
$$

or the last one,

$$
x_{2}\left(T_{n}\right)=\tau_{n} \times P\left(T_{n}, T_{N}\right) \times\left(s-c_{\mathrm{B}}\left(T_{n}, L_{N-1}\left(T_{n}\right) ; T_{N-1}, c-s ; \sigma_{n, N-1}\right)\right) .
$$

In conjunction with the first explanatory variable, either choice of $x_{2}\left(T_{n}\right)$ will capture the effect of changes in the interest rate curve slope.

In the example above, we proposed using an explanatory variable dependent on the stochastic volatility process (see (18.53)). It turns out that this is quite important to do in general, even for CLEs that do not have optionality in the underlying, such as Bermudan swaptions. The stochastic volatility process $z(t)$ can either be incorporated into an explanatory variable in the way of (18.53), or used as a separate explanatory variable.

### 18.3.10 Regression Implementation

While the selection of the regression variables is the primary key to the success of a regression-based method for CLE valuation, the details of implementation of the numerical algorithm for performing regressions are important as well. The basic regression algorithm does not involve much more than matrix inversion, see Section 3.5.4, but there are a number of ways in which the algorithm could be made more robust. We discuss them in this section.

For future reference, let us quickly fix notations. Throughout, we shall consider a particular date $T_{n}$ and assume that paths $\omega_{1}, \ldots, \omega_{K}$ have been simulated, allowing us to collect a $K \times q$ matrix $Z$ of simulated values of the regression variables $\zeta\left(T_{n}\right)$, with

$$
Z_{k, j}=\zeta_{j}\left(T_{n}, \omega_{k}\right), \quad k=1, \ldots, K, \quad j=1, \ldots, q .
$$

Recall that each $\zeta_{j}\left(T_{n}\right)$ is typically obtained as a particular monomial (18.52) applied to the vector of explanatory variables (or state variables, per Section 18.3.9.1) observed at time $T_{n}$. We also assume that we have available a $K$-dimensional vector of simulated values $Y=\left(Y_{1}, \ldots, Y_{K}\right)^{\top}$ that we would like to regress, for example the vector $Y_{k}=V_{n}\left(\omega_{k}\right), k=1, \ldots, K$, from (18.12) in the basic regression scheme of Section 18.3.1. The goal is to find a $q$-dimensional vector $\beta$ such that $Z \beta$ approximates $Y$ in some sense, e.g. in the least-squares sense

$$
\begin{equation*}
\|Y-Z \beta\|^{2} \rightarrow \min . \tag{18.54}
\end{equation*}
$$
^eq18-54

A (naive) solution to this problem is, of course, given by

$$
\begin{equation*}
\beta=\left(Z^{\top} Z\right)^{-1} Z^{\top} Y . \tag{18.55}
\end{equation*}
$$
^eq18-55

### 18.3.10.1 Automated Explanatory Variable Selection

As we explained in Section 18.3.9.2, we often find ourselves in a situation where we have many potential candidates for the explanatory variables but want to prune the set to keep only the most relevant variables for regression.

One approach here is to analyze potential candidates for explanatory variables and, based on experience, choose the ones that subjectively appear to be the most relevant. Another approach is to try to extract a subset of variables based on numerical criteria of regression fit. It turns out that the latter approach is quite common in econometrics circles, allowing us to draw on known techniques. In the problem of automatic econometric model selection, one considers a given time series of data - says investment returns of a hedge fund - and tries to choose macro variables, such as equity returns, oil prices, etc., that contribute the most to explaining the time series data. Conceptually, the solution to the problem is simple: one tries regressing the time series on subsets of potential regression variables and observes which subset of variables provides the best fit. The details of how these trials are conducted are, however, quite important, as a brute-force test of all variable combinations would often be impractical: even with just 10 variables to choose from, the number of potential subsets of variables to check is $2^{10}-1$. The literature on the subject of automatic econometric model selection is quite extensive ${ }^{19}$ so we do not go into much detail here, but merely scratch the surface by demonstrating a simple algorithm.

In an algorithm for automatic pruning of explanatory variables, the first question that needs to be answered is how to measure the quality of fit of a given set of regression variables. In the context of linear regression (18.54), an often-used measure ${ }^{20}$ is $R^{2}$ ("R-squared") which measures the variance of the residual $Y-Z \beta$ relative to the variance of $Y$

$$
R^{2}=1-\frac{\|Y-Z \beta\|^{2}}{\|Y\|^{2}}
$$

With the solution (18.55), $R^{2}$ is equal to

$$
R^{2}=\frac{Y^{\top} Z \beta}{Y^{\top} Y}
$$

The higher $R^{2}$, the better the fit is deemed to be. So, in our algorithm, we might choose different subsets of explanatory variables and examine the resulting values of the regression $R^{2}$. The results could, say, be used to choose a subset of variables that would give us the highest $R^{2}$, subject to a constraint on the maximum size $d$ of the variable vector.

Many search paths through the collection of subsets of explanatory variables are possible. For demonstration, let us consider a special case of the General-to-Specific (GETS) approach ${ }^{21}$. Here we start with all potential

[^40]explanatory variables included in the regression and calculate the baseline $R^{2}$. Then we remove each explanatory variable in turn (always returning previously-removed variables to the set after each regression) and calculate the $R^{2}$. Having tried removing all variables, we then select the one that gave us the smallest reduction in $R^{2}$ and throw it away. With the potential set of explanatory variables reduced by one, we repeat the procedure and again exclude a variable that gives us the smallest effect on the $R^{2}$. We continue until we have either reached a pre-determined number of explanatory variables, or until we have reached some maximum allowed reduction in $R^{2}$. With each regression typically being pretty quick to execute, this approach does not affect the overall valuation time much, yet generally gives good results. As one can imagine, however, substantially more sophisticated approaches for variable pruning are possible; we refer the reader to Campos et al. [2005] as a starting point.

While we motivated our discussion of explanatory variable selection with a problem in econometrics, we note that our regression problem is not exactly the same as that typically faced in econometric analysis. While econometricians often try to explain virtually all changes in their time series through changes in explanatory variables - i.e. they seek values of $R^{2}$ that are close to one - the problems we are interested in here would typically be characterized by values of $R^{2}$ that are much lower than one. This is a consequence of the fact that even the full information set available at time $T$ cannot explain all changes in hold/exercise values after time $T$. In fact, if in the regression we obtain high values of $R^{2}$, this would indicate that there is something wrong with our choice of explanatory variables and they most likely are "future-looking".

### 18.3.10.2 Suboptimal Point Exclusion

The main point of performing a regression on exercise/hold values is to determine an exercise rule. If for a given path, at a given point in time, we can prove that the exercise can never be optimal then, arguably, we should exclude this point from the regression that defines the exercise rule.

Interestingly, it turns out that in some general cases we can indeed establish situations where exercise is never optimal. Consider a cancelable note on a coupon stream $-X_{1}, \ldots,-X_{N-1}$, see Section 18.3.3. From (18.18),

$$
\begin{align*}
G_{n-1}\left(T_{n-1}\right) & =\mathrm{E}_{T_{n-1}}\left(B\left(T_{n-1}\right) B\left(T_{n}\right)^{-1}\left(-X_{n-1}+\left(G_{n}\left(T_{n}\right)\right)^{+}\right)\right) \\
& \geq-\mathrm{E}_{T_{n-1}}\left(B\left(T_{n-1}\right) B\left(T_{n}\right)^{-1} X_{n-1}\right) \tag{18.56}
\end{align*}
$$
^eq18-56

so if, for a given path $\omega,-\left(\mathrm{E}_{T_{n-1}}\left(B\left(T_{n-1}\right) B\left(T_{n}\right)^{-1} X_{n-1}\right)\right)(\omega)$ is positive, then $G_{n-1}\left(T_{n-1}, \omega\right)$ is positive and it cannot be optimal to cancel the note at $T_{n-1}$ on that path ${ }^{22}$. Fortunately, in many cases, $X_{n-1}$ is $\mathcal{F}_{T_{n-1}}$-measurable and[^41]

$$
\mathrm{E}_{T_{n-1}}\left(B\left(T_{n-1}\right) B\left(T_{n}\right)^{-1} X_{n-1}\right)=P\left(T_{n-1}, T_{n}\right) X_{n-1},
$$
a quantity easily computable at time $T_{n-1}$. Since we know the scenarios where it is never optimal to exercise, we can do two things. First, as a simple application of the policy fixing rule (18.47), we can outright forbid exercise in the suboptimal scenarios, even if our regression-based rule instructs us to exercise. Second, as suggested by Beveridge and Joshi [2008], we can exclude those paths $\omega_{k}$ for which $X_{n-1}\left(\omega_{k}\right)$ is negative from participating in the regression fit at time $T_{n-1}$, potentially improving the fit in the region that matters for the exercise rule.

The idea of excluding "uninteresting" paths could be taken further, albeit based on more practical than theoretical considerations. For example, we may decide to exclude (or at least de-emphasize) paths that are very deeply in or out of the money, since we want most precision around the exercise frontier itself. We formalize this idea in the next section.

### 18.3.10.3 Two Step Regression

As mentioned in our discussion on convexity in Section 18.3.9.3, simple polynomial functions may not capture to sufficient precision the functional dependence of regressed hold and exercise variables on regression variables. The functional mapping of Section 18.3.9.3 gives us one way of addressing this shortcoming. Another idea is to fit polynomials separately in different regions of values of explanatory variables. For example, with a single explanatory variable, we can split the range of explanatory variable values into a few intervals, and fit polynomials separately in each interval. With many explanatory variables, simple subdivision of the space into intervals quickly becomes impractical. In this case, however, Beveridge and Joshi [2009] suggest a somewhat similar idea of space subdivision based on the moneyness of the derivative in question. They point out that global regression fit of polynomials often works well for values of the regressed variables "in the wings", while for points near the "at-the-money" of the CLE (whatever that might mean), the richer functional structure of the hold and exercise values is often not well-approximated by globally fit polynomials. Beveridge and Joshi [2009] propose a two-step scheme, that we describe in the context of cancelable notes of Section 18.3.3 and, in particular, for the scheme (18.25)(18.26). First, for a given $T_{n}$, we regress the variables $Y_{k}=\widehat{G}_{n}\left(T_{n}, \omega_{k}\right)$ to obtain values $\zeta\left(T_{n}, \omega_{k}\right)^{\top} \beta$, an approximation to $\widetilde{G}_{n}\left(T_{n}, \omega_{k}\right), k=1, \ldots, K$. Subsequently we can choose those $\omega_{k}$ for which $\widetilde{G}_{n}\left(T_{n}, \omega_{k}\right)$ is close to zero, which we consider the definition of "at-the-money" (ATM) for a cancelable note. Specifically, for some threshold value $\epsilon$ we set

$$
\mathcal{W}^{\epsilon}=\left\{\omega_{k}:\left|\zeta\left(T_{n}, \omega_{k}\right)^{\top} \beta\right| \leq \epsilon\right\}
$$

Excluding paths outside of $\mathcal{W}^{\epsilon}$, we now perform a separate fit of $\widehat{G}_{n}\left(T_{n}, \omega\right)$, $\omega \in \mathcal{W}^{\epsilon}$, with the same polynomials; since the set of regression points is
different, a new regression coefficient vector $\beta^{\epsilon}$ is obtained. Finally, we set

$$
\tilde{G}_{n}\left(T_{n}, \omega_{k}\right)=\left\{\begin{array}{l}
\zeta\left(T_{n}, \omega_{k}\right)^{\top} \beta, \quad \omega_{k} \notin \mathcal{W}^{\epsilon} \\
\zeta\left(T_{n}, \omega_{k}\right)^{\top} \beta^{\epsilon}, \quad \omega_{k} \in \mathcal{W}^{\epsilon}
\end{array}\right.
$$

i.e. we use the values of the original regression for non-ATM points $\omega_{k} \notin \mathcal{W}^{\epsilon}$, and a new regression for ATM points $\omega_{k} \in \mathcal{W}^{\epsilon}$ only. The value of $\epsilon$ could be set to be some (small) fraction of the notional of the derivative, say $5 \%$; ultimately some experimentation here may be required.

Rather than a binary division of regression points into "near" and "far" groups, one could also imagine using some kind of smooth kernel that weighs points in the regression differently, depending on how close one is perceived to be to the exercise point, e.g. by how far away from zero $\widetilde{G}_{n}\left(T_{n}, \omega_{k}\right)$ is. One can also refine it through iterations, where the procedure is repeated multiple times.

### 18.3.10.4 Robust Implementation of Regression Algorithm

In Section 18.3.1 we assumed that the regression operator $\mathcal{R}_{T}$ at time $T=0$ reduces to averaging of the values of the random variable it is applied to, because the values of the regression variables $\zeta(0)$ are the same for all paths. This is indeed the case for the true solution of the problem (18.54). If we assume that all the rows of $Z$ are the same and equal to a vector $\varsigma^{\top}$, then the objective function in (18.54) is equal to

$$
\sum_{k=1}^{K}\left(Y_{k}-\varsigma^{\top} \beta\right)^{2}
$$

Differentiating with respect to $\beta_{j}$ and setting the derivative to zero we obtain

$$
\sum_{k=1}^{K}\left(Y_{k}-\varsigma^{\top} \beta\right) \varsigma_{j}=0
$$

giving us a regressed solution

$$
\varsigma^{\top} \beta=\frac{1}{K} \sum_{k=1}^{K} Y_{k}
$$

i.e. the average of $Y_{k}$ 's as advertised. The solution, however, is not given by (18.55) as the inverse of $Z^{\top} Z=K \varsigma \varsigma^{\top}$ here will not exist ${ }^{23}$.

This simple example highlights the danger of using the textbook regression solution (18.55) to the problem (18.54). Beyond the degenerate case

[^42]described above, similar issues will arise whenever the matrix $Z^{\top} Z$ is illconditioned, i.e. close to singular. This can happen either due to outright user error (e.g. an inexperienced user accidentally entering the same explanatory variable twice) or due to subtle near-linear dependencies between explanatory variables. In such cases, the regression problem becomes ill-posed and the numerical solution of the regression problem will be unstable. To counteract this, the user will normally have to add additional structure to the regression problem in order for a robust solution with desirable properties to exist.

To stabilize an ill-posed regression problem, we should first contemplate what would constitute desirable properties for the vector $\beta$ provided that the regression data matrix $Z$ imposes insufficient constraints on its behavior. A standard approach is to give preference to solutions for $\beta$ with smaller norms, a choice we can motivate by the observation that if some of the regression coefficients are not constrained by the data, they should be set to zero to ensure that the corresponding regression variables (monomials) do not contribute to the fit. With this in mind, we choose a scalar regularization weight $w_{\text {reg }}>0$ and replace (18.54) with

$$
\begin{equation*}
\|Y-Z \beta\|^{2}+w_{\text {reg }}\|\beta\|^{2} \rightarrow \min . \tag{18.57}
\end{equation*}
$$
^eq18-57

Intuitively, we should choose the weight $w_{\text {reg }}$ small enough so that the extra term in the objective function does not interfere with the regression objective, yet sufficiently large that the regularization term performs its function. Numerous methods for data-driven selection of $w_{\text {reg }}$ have been published in the literature, including the $L$-curve method in Hansen [1992], the discrepancy principle in Morozov [1966], and generalized cross-validation in Craven and Wahba [1979]. Andersen [2005] contains a review of many of these methods in the setting of yield curve construction. Whatever method is used, the size of $w_{\text {reg }}$ should obviously reflect the relative scale of the numbers used in the regression. To examine the scaling issue a bit further, notice that the quadratic (in $\beta$ ) term in the objective function is given by

$$
\beta^{\top}\left(Z^{\top} Z+w_{\mathrm{reg}} I\right) \beta
$$

where $I$ is the identity $q \times q$ matrix. The sum of squares of all elements in $I$ and $Z^{\top} Z$ equals $q$ and

$$
\operatorname{tr}\left(\left(Z^{\top} Z\right)^{\top}\left(Z^{\top} Z\right)\right)
$$

respectively. Consequently it is natural to write

$$
\begin{equation*}
w_{\mathrm{reg}}=\epsilon\left(\frac{1}{q} \operatorname{tr}\left(Z^{\top} Z Z^{\top} Z\right)\right)^{1 / 2} \tag{18.58}
\end{equation*}
$$
^eq18-58

where $\epsilon$ is a new scale-free constant to be determined. While ideally we should rely on one of the data-driven approaches above, in a pinch we can always try to set $\epsilon$ equal to a small number, such as $10^{-4}$.

The formal solution to (18.57) is given by

$$
\begin{equation*}
\beta=\left(Z^{\top} Z+w_{\mathrm{reg}} I\right)^{-1} Z^{\top} Y \tag{18.59}
\end{equation*}
$$
^eq18-59

Note that the matrix $Z^{\top} Z+w_{\text {reg }} I$ is of full rank even if $Z^{\top} Z$ is not, so the matrix inversion in (18.59) is always well-defined, even when $Z$ is illconditioned. The resulting method (which we have used before, in Section 6.4.3) is often called Tikhonov regularization or ridge regression.

Tikhonov regularization is attractive because it retains a fair amount of intuition as to what happens to the regression coefficients as a result of regularization. We should, however, mention other regularization alternatives, in particular the pseudo-inverse or truncated singular value decomposition (TSVD) method. To briefly outline this approach, let us rewrite (18.55) as a system of linear equations on $\beta$,

$$
\begin{equation*}
M \beta=Z^{\top} Y, \tag{18.60}
\end{equation*}
$$
^eq18-60

where $M=Z^{\top} Z$. The SVD method, see Press et al. [1992], allows us to decompose the $q \times q$ matrix $M$ into a product of three matrices,

$$
M=U \Sigma V^{\top}
$$

where $U$ (not to be confused with the exercise value notation from earlier in the chapter) and $V$ are $q \times q$ orthogonal matrices (i.e. $U^{\top}=U^{-1}, V^{\top}=V^{-1}$ ) and $\Sigma$ is a diagonal $q \times q$ matrix. The diagonal elements of $\Sigma$ are called singular values and are ordered by their absolute value (highest first). The decomposition applies even to singular matrices $M$. In particular, if $M$ is of rank $r, r \leq q$, then only the first $r$ diagonal elements of $\Sigma$ will be non-zero.

The pseudo-inverse of the matrix $M$ is defined by

$$
M^{+}=V \Sigma^{+} U^{\top}
$$

where $\Sigma^{+}$is a diagonal matrix with elements

$$
\Sigma_{i, i}^{+}= \begin{cases}1 / \Sigma_{i, i}, & \Sigma_{i, i} \neq 0 \\ 0, & \Sigma_{i, i}=0\end{cases}
$$

Then we have

$$
M M^{+}=M^{+} M=\operatorname{diag}\left((1, \ldots, 1,0, \ldots, 0)^{\top}\right)
$$

where there are as many 1's on the diagonal as the rank of he matrix $M$.
The pseudo-inverse allows us to define a solution to (18.60) that always exists. For numerical stability, it is common to modify the solution slightly by choosing a truncation cut-off value $\epsilon>0$ and defining a diagonal matrix $\Sigma^{\epsilon}$ by

$$
\Sigma_{i, i}^{\epsilon}= \begin{cases}1 / \Sigma_{i, i}, & \left|\Sigma_{i, i}\right| \geq \epsilon\left|\Sigma_{1,1}\right| \\ 0, & \left|\Sigma_{i, i}\right|<\epsilon\left|\Sigma_{1,1}\right|\end{cases}
$$

Then the solution to (18.60) and, ultimately, the regression problem (18.54), is given by

$$
\beta=V \Sigma^{\epsilon} U^{\top} Z^{\top} Y .
$$

A possible choice for $\epsilon$ is $\epsilon=10^{-6}$.
To understand better the intuition behind TSVD, let us highlight an interesting connection between Tikhonov regularization and TSVD. Specifically, it can be shown that the Tikhonov solution (18.59) to the regression problem can be written as

$$
\beta=V \Sigma^{\text {Tikhonov }} U^{\top} Z^{\top} Y .
$$

where $\Sigma^{\text {Tikhonov }}$ is a diagonal $q \times q$ matrix computed from the singular value matrix $\Sigma$ as

$$
\Sigma_{i, i}^{\text {Tikhonov }}=\frac{\Sigma_{i, i}}{\Sigma_{i, i}^{2}+w_{\text {reg }}}
$$

We recognize this as a smoothed version of the cut-off matrix $\Sigma^{\epsilon}$ above, with the Tikhonov factor $w_{\text {reg }}$ determining how much small singular values get dampened out.

Singular values of widely different magnitudes in the matrix $Z^{\top} Z$ generally cause numerical problems in inversion of the matrix, something that the Tikhonov and TSVD method can help rectify. However, widely different scales of singular values do not necessarily arise only when explanatory variables are poorly chosen (e.g. highly dependent among themselves), but can also emerge if explanatory variables used are themselves of diffcrent scales. For example, if one variable is a swap rate measured in the units of a few percent and another is the value of the current coupon measured in the unit of millions of dollars, such scale discrepancy could lead to numerical problems in the regression. The problem is exacerbated by our choice of polynomials as basis functions, as one million to the power of, say, 4 is obviously quite different from one percentage point to the 4th power.

Fortunately, such scaling issues are easy to rectify, as we only need to rescale all variables to the same base before applying the regression. So, instead of the matrix $Z$ we would use $\widetilde{Z}$, whose elements are given by

$$
\begin{equation*}
\widetilde{Z}_{i, j}=\frac{Z_{i, j}-\bar{Z}_{\cdot, j}}{\left(K^{-1} \sum_{k=1}^{K}\left(Z_{k, j}-\bar{Z}_{\cdot, j}\right)^{2}\right)^{1 / 2}}, \quad i=1, \ldots, K, \quad j=1, \ldots, q \tag{18.61}
\end{equation*}
$$
^eq18-61

where $\bar{Z}_{\cdot, j}=K^{-1} \sum_{k=1}^{K} Z_{k, j}$. This transformation sets all columns in the matrix $Z$ to have zero mean and unit (empirical) standard deviation. The only caveat with (18.61) is that, if applied to a column of $Z$ with constant numbers - a column that is always present as we typically include a constant function in our regression - a division by 0 would occur. One obvious workaround is to simply avoid scaling the constant columns.

Once the standard deviation, as required by (18.61), of each column of $Z$ is calculated, then we can apply another simple method to increase the robustness: we can watch out for points that are well outside a "reasonable" range - say outside of 10 standard deviations for a given column - and exclude them from our regression.

### 18.4 Valuation with Low-Dimensional Models

Libor market models are often our preferred choice for valuation and risk management of callable Libor exotics. For some CLEs, however, we can use simpler and faster models without sacrificing the benefits of proper calibration and good model dynamics. The trick here is to calibrate a simpler model in a special way, an approach we call the local projection method.

### 18.4.1 Single-Rate Callable Libor Exotics

The local projection method builds on the calibration discussion of Section 18.1 by calibrating a "local", low-dimensional model to the volatility information that we identified as important to the CLE valuation. Information obtained from the market is used directly, and the rest is extracted from a "global", fully-calibrated model such as the Libor market model. The success of the method depends on our ability to identify the relevant volatility information, and how well the local model can calibrate to this information. As a low-dimensional model has only a limited number of parameters, it can only be successfully calibrated for a CLE that depends on a relatively small subset of all available market information.

Callable Libor exotics most amendable to the local projection method are those that have, for each $n, n=1, \ldots, N-1$, coupon $C_{n}$ that is a function of at most a single market rate. We denote such structures single-rate CLEs. Examples include Bermudan swaptions, callable inverse floaters, callable CMS capped floaters and fixed-rate callable range accruals. Excluded are CLEs whose coupons depend on spreads between CMS rates, floating-rate callable range accruals, and similar.

The main attraction of using a low-dimensional model is the ability to apply PDE methods for valuation. We have already briefly discussed a relevant pricing scheme, see Section 2.7.4, and the mechanics of the valuation algorithm typically present no special difficulties unless the underlying CLE has path-dependent features. Section 18.4.5 discusses certain PDE techniques that can be used if the path-dependency is sufficiently weak.

### 18.4.2 Calibration Targets for the Local Projection Method

To start, let us focus on the heart of the local projection method, namely the choice of calibration targets for the local model. Let $\left\{S_{n}^{1}(t)\right\}_{n=1}^{N-1}$ be the
strip of swap rates that define coupons of a CLE, so that each $C_{n}$ depends on $S_{n}^{1}\left(T_{n}\right), n=1, \ldots, N-1$. We assume that the swap rate $S_{n}^{1}$ has $\mu(n)$ periods, so that $S_{n}^{1}(t)=S_{n, \mu(n)}(t)$, where $S_{n, m}(t)$ is the standard notation for a swap rate fixing on $T_{n}$ and covering $m$ periods, see e.g. (4.10). For example, for a callable inverse floater we have $S_{n}(t)=L_{n}(t)$ (a Libor rate) and $\mu(n)=1$; and for a callable CMS capped floater $\mu(n) \equiv k$, where $k$ is the number of periods for the underlying CMS rate. In addition, we define a second strip of swap rates $\left\{S_{n}^{2}(t)\right\}_{n=1}^{N-1}$ to be the core, or coterminal, swap rate strip, i.e. $S_{n}^{2}(t)=S_{n, N-n}(t), n=1, \ldots, N-1$.

The underlying of the CLE, an exotic swap, can be expressed as a strip of options where the $n$-th option is written on the rate $S_{n}^{1}\left(T_{n}\right)$. Thus, for a model to reprice the underlying correctly, it should be calibrated to the market (spot) volatilities of the first swap rate strip, i.e. the implied volatilities of the European swaptions defined by $\left\{S_{n}^{1}\left(T_{n}\right)\right\}_{n=1}^{N-1}$. The ability to match the underlying exotic swap (i.e. the non-callable CLE) is certainly a prerequisite for any reasonable model for a callable CLE, but, as we have already seen in Section 18.1.1, we also need to consider that the volatilities and inter-temporal correlations of core swap rates $\left\{S_{n}^{2}(\cdot)\right\}_{n=1}^{N-1}$ will affect the value of the callability feature of the CLE. In light of this, as a starting point (to be refined later, see Section 18.4.4) we suggest that any low-dimensional model be calibrated to the following targets:

- The underlying volatilities, or swap rate volatilities for $\left\{S_{n}^{1}\left(T_{n}\right)\right\}_{n=1}^{N-1}$ that correspond to strikes relevant for the coupons $C_{n}$ or, in a pinch, to at-the-money strikes.
- The core volatilities, or swap rate volatilities for $\left\{S_{n}^{2}\left(T_{n}\right)\right\}_{n=1}^{N-1}$. The choice of swaption strikes used to define core volatilities is often not straightforward, but at-the-money strikes is a common choice. In some cases more advanced methods for strike selection are available, see e.g. Section 19.3.
- The core correlations, or inter-temporal correlations for $\left\{S_{n}^{2}\left(T_{n}\right)\right\}_{n=1}^{N-1}$.

While volatilities of swap rates are directly observable from the market, the inter-temporal correlations are not. This is where we can draw on the LM (or similar) global model; once it has been calibrated to the market as a whole, we can calculate the required correlations from the global model. In a nutshell, the role of the global model is to serve as our "correlation extractor". The important point here is that by including dynamic information such as inter-temporal correlations as calibration targets, the local model not only captures the static information about interest rate volatilities at valuation time, but also the transition densities and dynamics of the volatility structure, as seen by a global, fully calibrated and, presumably, realistic model.

### 18.4.3 Review of Suitable Local Models

The one-dimensional quasi-Gaussian (qG) model developed in Chapter 13 is a natural candidate to consider for the role of the local model in the local projection method for single-rate CLEs. A simple, yet useful, special case setup is based on the version of the qG model with linear local volatility, see Section 13.1.6. The volatility structure of such a qG model is controlled by several time-dependent functions, including the volatility function, the skew function and the mean reversion function. If convexity in the volatility smile is deemed important, the model could be upgraded to the stochastic volatility version in Section 13.2.

Let us first look at the volatility structure specification; we will consider skew and stochastic volatility parameter selection later on. With the volatility function and the mean reversion function discretized over the tenor structure $\left\{T_{n}\right\}_{n=0}^{N}$, the qG model has $2(N-1)$ independent parameters for volatility calibration. As discussed in Section 13.1.7, the volatility parameters can be used to calibrate the model to term volatilities for one of the swap rate strips. The mean reversion can be used to either match the term volatilities for the second swap rate strip (Section 13.1.8.2), or the inter-temporal correlations (Section 13.1.8.3). As such, the one-factor qG model is not large enough to match all three sets of calibration targets identified in Section 18.4.2 above. In some situations, however, this might be acceptable as some securities may turn out to depend only weakly on one of the three calibration targets. For example, for shorter-dated CLEs, inter-temporal correlations may not affect the CLEs value all that much. Likewise, if the underlying has options on the rates $\left\{S_{n}^{1}\left(T_{n}\right)\right\}_{n=1}^{N-1}$ that are deep in or out of the money, this set of calibration targets can potentially be dropped. Of course, all such decisions must be supported by extensive testing, which fortunately is easy to do as we always have the global (LM or similar) model to benchmark against.

For derivatives where all three sets of calibration targets are important, a one-factor qG model will not suffice ${ }^{24}$, and we ideally need to move to models with more stochastic factors. One particularly simple choice is here the two-factor Gaussian model, see Section 12.1.4, which has enough degrees of freedom to match all volatility targets. In this model, some of the time-dependent parameters can be chosen to be constant to make the dynamics of the volatility structure implied by the model more realistic, see Lemma 12.1.11.

The disadvantage of the two-factor Gaussian model is, of course, its lack of control over the volatility smile, so calibration to the volatility targets will require us to identify a single strike per swaption. Improved smile fits can be accomplished by using two-factor versions of either the quasi-Gaussian model (Chapter 13), the affine model, or the quadratic Gaussian model (Chapter

[^43]12). While all these models are different in some regards, the underlying philosophy and calibration methods will be quite similar.

For models that are sufficiently rich to incorporate volatility skew/smile effects (such as local volatility or stochastic volatility qG models), we also need to select the market information to which we wish to calibrate skew and smile parameters. Normally we would extend one of the swaption strips, $\left\{S_{n}^{1}\right\}$ or $\left\{S_{n}^{2}\right\}$, to multiple strikes for this purpose (for mechanics of calibration see e.g. Section 13.2.3). The choice of the strip is typically driven by an analysis of relative importance of the two sets of smiles to the value of the CLE. It is difficult to state any firm general guidelines here, but we can observe that it is often fairly easy to match the underlying exotic swap value by a judicious choice of a single strike per maturity in the swaption strip $\left\{S_{n}^{1}\right\}$. On the other hand, it is often difficult to establish which strikes are the most relevant for the "callability" value. Given this, it is often reasonable to use whatever skew/smile parameters we have at our disposal to improve the broad fit of implied core swaption volatilities (the strip $\left\{S_{n}^{2}\right\}$ ) at multiple strikes per maturity. If we only have skew, but not smile, parameters, we can use these to match two volatilities at each maturity, or to match the slope of the volatility smile at a given strike. The latter could be important if the underlying structured coupons are not simple European options but, for example, of digital or range-accrual type, in which case it is the slope of the volatility smile, and not the overall level, that drives the underlying value. In this case we might, in fact, want to use skew to calibrate the underlying, rather than the callability, value.

### 18.4.4 Defining a Suitable Analog for Core Swap Rates

When we in Section 18.4.2 looked for the elements of the volatility structure that are relevant for a callable Libor exotic security, we argued that the callability value is driven by the volatilities of core swap rates $S_{n}^{2}(t)=S_{n, N-n}(t)$, since a CLE is related to a standard Bermudan swaption. This argument, clearly, has limitations of its applicability. For instance, in Section 19.4 we study Bermudan swaptions on amortizing swaps and show that the most relevant European swaptions in this case are not the standard core European swaptions, but swaptions with tenors based on the durations of the underlying amortizing swap.

In light of this, let us try to refine the selection of the volatility targets relevant for the callability option of a CLE. As a starting point we can use the idea that the local model should match the values of European options on exercise values $U_{n}\left(T_{n}\right), n=1, \ldots, N-1$. While market values of options on $U_{n}\left(T_{n}\right)$ could be hard to come by, we can linearize the underlying $U_{n}\left(T_{n}\right)$ of the CLE and use the resulting rate as a replacement for the core swap rate that should be used in volatility and correlation calibration.

Using a LM model as a backdrop for our analysis, the exercise values $U_{n}\left(T_{n}\right)$, for each $n=1, \ldots, N-1$, are functions of the vector of primary

Libor rates

$$
\mathbf{L}\left(T_{n}\right)=\left(L_{n}\left(T_{n}\right), \ldots, L_{N-1}\left(T_{n}\right)\right)^{\top}
$$

observed on the date $T_{n}$, i.e.,

$$
U_{n}\left(T_{n}\right)=f_{n}\left(\mathbf{L}\left(T_{n}\right)\right), \quad n=1, \ldots, N-1 .
$$

Linearizing this expression, we obtain

$$
U_{n}\left(T_{n}\right) \approx f_{n}(\mathbf{L}(0))+\nabla f_{n}(\mathbf{L}(0))\left(\mathbf{L}\left(T_{n}\right)-\mathbf{L}(0)\right),
$$

where $\nabla f_{n}(x)$ is the (row vector) gradient of $f_{n}(x)$. Hence, the value of the European option on the underlying,

$$
\mathrm{E}\left(B\left(T_{n}\right)^{-1}\left(U_{n}\left(T_{n}\right)\right)^{+}\right)
$$

can be approximated with

$$
\begin{equation*}
\mathrm{E}\left(B\left(T_{n}\right)^{-1}\left(\nabla f_{n}(\mathbf{L}(0)) \mathbf{L}\left(T_{n}\right)-\left(\nabla f_{n}(\mathbf{L}(0)) \mathbf{L}(0)-f_{n}(\mathbf{L}(0))\right)\right)^{+}\right) . \tag{18.62}
\end{equation*}
$$
^eq18-62

We can therefore argue that the most relevant "interest rate" is

$$
R_{n}\left(T_{n}\right)=\nabla f_{n}(\mathbf{L}(0)) \mathbf{L}\left(T_{n}\right)=\sum_{j} w_{n, j} L_{j}\left(T_{n}\right), \quad w_{n, j} \triangleq \frac{\partial f_{n}}{\partial L_{j}}(\mathbf{L}(0)) .
$$

Being a linear combination of Libor rates, $R_{n}\left(T_{n}\right)$ is not, strictly speaking, a market swap rate. However, the volatility of the rate $R_{n}\left(T_{n}\right)$ can be approximated in a Libor market model (along the same lines as in Section 14.4.2), as well as in local models we may wish to use. Therefore, we can easily use the term volatilities of $R_{n}\left(T_{n}\right), n=1, \ldots, N-1$, as volatility targets in place of core swap rates.

The underlying $U_{n}(t)$ typically consists of options on market rates. The derivatives $\partial f_{n} / \partial L_{j}$ can then be computed quite easily with, say, Black-type approximations to option values. Volatilities that should be used in these calculations are the forward (as observed at time $T_{n}$ ) volatilities. Needless to say, high degree of precision is not necessary in these calculations.

To see how consistent the method defined above with our recommendations in Section 18.4.2, let us apply it to a standard Bermudan swaption. For a payer swap with coupon $K$ we have

$$
\begin{aligned}
U_{n}\left(T_{n}\right) & =\sum_{i=n}^{N-1} \tau_{i}\left(L_{i}\left(T_{n}\right)-K\right) P\left(T_{n}, T_{i+1}\right) \\
& =\sum_{i=n}^{N-1} \tau_{i}\left(L_{i}\left(T_{n}\right)-K\right)\left(\prod_{k=n}^{i} \frac{1}{1+\tau_{k} L_{k}\left(T_{n}\right)}\right) .
\end{aligned}
$$

Hence,

$$
f_{n}(x)=\sum_{i=n}^{N-1} \tau_{\imath}\left(x_{i}-K\right)\left(\prod_{k=n}^{i} \frac{1}{1+\tau_{k} x_{k}}\right),
$$

and

$$
\frac{\partial f_{n}}{\partial L_{j}}(x)=\tau_{j} \prod_{k=n}^{j} \frac{1}{1+\tau_{k} x_{k}}-\frac{\tau_{j}}{1+\tau_{j} x_{j}} \sum_{i=j}^{N-1} \tau_{\imath}\left(x_{i}-K\right)\left(\prod_{k=n}^{\imath} \frac{1}{1+\tau_{k} x_{k}}\right) .
$$

Thus

$$
w_{n, j}=\tau_{j} P\left(0, T_{n}, T_{j+1}\right)\left(1-U_{j}(0) / P\left(0, T_{j}\right)\right) .
$$

If we compare these weights with those obtained by decomposition of the swap rate into a sum of Libor rates via the "freezing" techniques in Section 14.4.2, we see that they are roughly the same, up to a constant scaling. Thus

$$
R_{n}\left(T_{n}\right)=\sum_{j=n}^{N-1} w_{n, j} L_{j}\left(T_{n}\right)
$$

is quite close to the (scaled) core swap rate $S_{n, N-n}\left(T_{n}\right)$. Therefore, for standard Bermudan swaptions, using $R_{n}\left(T_{n}\right)$ to define volatility calibration targets should be approximately consistent with the standard method of using core swap rates. When applied to non-standard (e.g., amortizing) Bermudan swaptions, this method produces results that are similar to what we propose later in Section 19.4.

It should be clear that the choice of calibration targets has carries significant impact on the value of a CLE in a local model. Equally important, it also defines the "basis" for vegas (volatility sensitivities), i.e. the set of swaption volatilities to which the CLE is sensitive; hedging of volatility exposure in local model would therefore, as a practical matter, only be done with the swaptions included in the calibration strips. Using a particular swaption for calibration implies the dependence of the CLE value to the volatility of that swaption; conversely, omitting a swaption from the calibration set makes the CLE value (in the local model) insensitive to its volatility. This, of course, is not wholly realistic as even simple CLEs (including Bermudan swaptions) would, when priced in a global model, typically show sensitivity to volatilities of all swaptions whose total maturity is no greater than the CLE maturity ${ }^{25}$. In some sense, vegas from a local model (to a subset of swaptions) can be thought of as an aggregation of vegas from the global model. Some traders may in fact prefer such an aggregated view as it (seemingly) simplifies the job of vega hedging.

[^44]
### 18.4.5 PDE Methods for Path-Dependent CLEs

As mentioned earlier, one attraction of the local projection method is the fact that the resulting model state can often be represented by a low-dimensional state vector $x(t)$. If the dimension of $x(t)$ is less than 3 or 4 , this will often allow us to state the CLE value as the solution to a PDE, a problem that can be attacked by the finite difference methods in Chapter 2. While callability is easy to handle (see Section 2.7.4), most path-dependent CLEs are outside the scope of finite difference methods. Exceptions do exist, however, if the path-dependency is sufficiently mild. We show some examples of this below. Of course, even in those cases where a PDE solution is technically possible, one should contemplate whether a local projection model is fundamentally suitable for the path-dependent derivative in question. In particular, the basic single-rate CLE calibration strategies may need adjustment to better capture the path-dependent feature of the payout.

### 18.4.5.1 CLEs Accreting at Coupon Rate

One particular class of path-dependent CLEs that is amenable to PDE methods has its path-dependency confined to the CLE notional only, see Piterbarg [2002]. A prime example of such a CLE is a callable Libor exotic accreting at a coupon rate, see Section 5.14.5, which is the example we consider here. Recall that a CLE is defined by its structured coupon $C_{n}$ that is fixed at time $T_{n}$ and paid at time $T_{n+1}, n=1, \ldots, N-1$. In the standard CLE, the notional of the coupon is constant, or at least deterministic, and has been factored out from the definition in Section 18.2.1. For a couponaccreting CLE, the notional to which the coupon rate and the Libor rate are applied at time $T_{n+1}$ is equal to the notional at time $T_{n}$ times an accretion factor that depends on $C_{n}$.

Formally, we replace (18.1) with

$$
X_{n}=D_{n} \tau_{n}\left(C_{n}-L_{n}\left(T_{n}\right)\right)
$$

where $D_{1}=1$ and

$$
\begin{equation*}
D_{n}=D_{n-1} \times\left(1+\tau_{n-1} C_{n-1}\right), \quad n=2, \ldots, N-1 . \tag{18.63}
\end{equation*}
$$
^eq18-63

A coupon-accreting CLE is defined as a Bermudan-style option to enter, on date $T_{n}$, the remaining part of the underlying, i.e. an exotic swap with the value (at time $t \leq T_{n}$ ),

$$
U_{n}(t)=B(t) \sum_{\imath=n}^{N-1} \mathrm{E}_{t}\left(B\left(T_{i+1}\right)^{-1} X_{i}\right)
$$

The backward-induction scheme (18.7), (18.16) is trivially extended,

$$
\begin{align*}
U_{n}\left(T_{n}\right)= & \tau_{n} D_{n}\left(C_{n}-L_{n}\left(T_{n}\right)\right) P\left(T_{n}, T_{n+1}\right) \\
& +B\left(T_{n}\right) \mathrm{E}_{T_{n}}\left(B\left(T_{n+1}\right)^{-1} U_{n+1}\left(T_{n+1}\right)\right)  \tag{18.64}\\
H_{n}\left(T_{n}\right)= & B\left(T_{n}\right) \mathrm{E}_{T_{n}}\left(B\left(T_{n+1}\right)^{-1} \max \left(U_{n+1}\left(T_{n+1}\right), H_{n+1}\left(T_{n+1}\right)\right)\right) \tag{18.65}
\end{align*}
$$
^eq18-64
^eq18-65

As written, the scheme cannot be directly implemented in a PDE solver because of path-dependency under the expected value operator in (18.64). However, by employing the method of similarity reduction, the scheme can be rewritten in a way amendable to a PDE representation.

Dividing both sides of (18.64)-(18.65) by $D_{n}$, and using the fact that $D_{n}$ is $\mathcal{F}_{T}$,-measurable, we get

$$
\begin{aligned}
U_{n}^{\prime}\left(T_{n}\right)= & \tau_{n}\left(C_{n}-L_{n}\left(T_{n}\right)\right) P\left(T_{n}, T_{n+1}\right) \\
& +B\left(T_{n}\right) \mathrm{E}_{T_{n}}\left(B\left(T_{n+1}\right)^{-1} \frac{D_{n+1}}{D_{n}} U_{n+1}^{\prime}\left(T_{n+1}\right)\right) \\
H_{n}^{\prime}\left(T_{n}\right)= & B\left(T_{n}\right) \mathrm{E}_{T_{n}}\left(B\left(T_{n+1}\right)^{-1} \frac{D_{n+1}}{D_{n}} \max \left(U_{n+1}^{\prime}\left(T_{n+1}\right), H_{n+1}^{\prime}\left(T_{n+1}\right)\right)\right)
\end{aligned}
$$

where we have defined

$$
U_{n}^{\prime}\left(T_{n}\right)=\frac{U_{n}\left(T_{n}\right)}{D_{n}}, \quad H_{n}^{\prime}\left(T_{n}\right)=\frac{H_{n}\left(T_{n}\right)}{D_{n}}, \quad n=1, \ldots, N-1
$$

From (18.63),

$$
\frac{D_{n+1}}{D_{n}}=1+\tau_{n} C_{n}
$$

where $1+\tau_{n} C_{n}$ is $\mathcal{F}_{T_{n}}$-measurable. Therefore, the factor $D_{n+1} / D_{n}$ can be pulled out from inside the expected value operator, to give us

$$
\begin{align*}
U_{n}^{\prime}\left(T_{n}\right)= & \tau_{n}\left(C_{n}-L_{n}\left(T_{n}\right)\right) P\left(T_{n}, T_{n+1}\right) \\
& +\left(1+\tau_{n} C_{n}\right) B\left(T_{n}\right) \mathrm{E}_{T_{n}}\left(B\left(T_{n+1}\right)^{-1} U_{n+1}^{\prime}\left(T_{n+1}\right)\right)  \tag{18.66}\\
H_{n}^{\prime}\left(T_{n}\right)= & \left(1+\tau_{n} C_{n}\right) B\left(T_{n}\right) \\
& \times \mathrm{E}_{T_{n}}\left(B\left(T_{n+1}\right)^{-1} \max \left(U_{n+1}^{\prime}\left(T_{n+1}\right), H_{n+1}^{\prime}\left(T_{n+1}\right)\right)\right) \tag{18.67}
\end{align*}
$$
^eq18-66
^eq18-67

These equations are used for $n=N-1, \ldots, 0$, with $U_{N}^{\prime}\left(T_{N}\right)=H_{N}^{\prime}\left(T_{N}\right)=0$. We have the following result.

Proposition 18.4.1. For a coupon-accreting $C L E, U_{n}^{\prime}\left(T_{n}\right)$ and $H_{n}^{\prime}\left(T_{n}\right)$ can, for each $n=0, \ldots, N$, be written as deterministic functions of the model state variables model at time $T_{n}$.

Proof. The proof is by induction. The statement is trivially true for $n=N$. To prove the induction step, we assume it is true for $n+1$. We note that

$$
B\left(T_{n}\right) \mathrm{E}_{T_{n}}\left(B\left(T_{n+1}\right)^{-1} U_{n+1}^{\prime}\left(T_{n+1}\right)\right)
$$

and

$$
B\left(T_{n}\right) \mathrm{E}_{T_{n}}\left(B\left(T_{n+1}\right)^{-1} \max \left(U_{n+1}^{\prime}\left(T_{n+1}\right), H_{n+1}^{\prime}\left(T_{n+1}\right)\right)\right)
$$

can be written as functions of the state variables at time $T_{n}$ by applying the PDE rollback scheme to

$$
U_{n+1}^{\prime}\left(T_{n+1}\right), \quad \max \left(U_{n+1}^{\prime}\left(T_{n+1}\right), H_{n+1}^{\prime}\left(T_{n+1}\right)\right)
$$

(which are functions of the state variables at time $T_{n+1}$ by the induction hypothesis). The accreting factor ( $1+\tau_{n} C_{n}$ ) and the marginal coupon $\tau_{n}\left(C_{n}-L_{n}\left(T_{n}\right)\right) P\left(T_{n}, T_{n+1}\right)$ are functions of the state variables at time $T_{n}$ as well. The proposition is proved.

The new scheme (18.66)-(18.67) in fact looks just like (18.64)-(18.65) for a unit-notional CLE, with one modification: on each backward induction step, the values of $U^{\prime}$ and $H^{\prime}$ are rescaled by the "marginal" accreting notional $\left(1+\tau_{n} C_{n}\right)$. The key fact here is that this factor is known at time $T_{n}$.

### 18.4.5.2 Snowballs

Certain other path-dependent CLEs can be valued by PDE methods by introducing extra state variables, along the lines of Section 2.7.5. Among the more popular CLEs for which this method is applicable are the snowball swaps and callables, see Chapter 5. In a snowball, the structured coupon at time $T_{n}$ is a function of the structured coupon at time $T_{n-1}$ (and rates at time $T_{n}$ ). As an example, recall the basic structure from Section 5.13.5 with the coupon at time $T_{n}$ given by

$$
C_{n}=\left(C_{n-1}+s_{n}-g_{n} \times L_{n}\left(T_{n}\right)\right)^{+}, \quad n=1, \ldots, N-1,
$$

with $C_{0}$ being a fixed initial coupon.
We can value snowball swaps and callables by PDE method after the introduction of an extra state variable $I(t)$ defined to be the current coupon,

$$
I(t)=\sum_{n=0}^{N-1} 1_{\left\{t \in\left\{T_{n}, T_{n+1}\right\}\right\}} C_{n}
$$

The backward recursion for the exercise value at time $T_{n}$ then reads,

$$
\begin{equation*}
U_{n}\left(T_{n}\right)=\tau_{n}\left(I\left(T_{n}\right)-L_{n}\left(T_{n}\right)\right) P\left(T_{n}, T_{n+1}\right)+\mathrm{E}_{T_{n}}\left(\frac{B\left(T_{n}\right)}{B\left(T_{n+1}\right)} U_{n+1}\left(T_{n+1}\right)\right) \tag{18.68}
\end{equation*}
$$
^eq18-68

Hence we obtain the following continuity condition for a given value of $I\left(T_{n}\right)=I$ (where $T_{n}$ - is a time immediately prior to $T_{n}$ ),

$$
\begin{aligned}
U_{n}\left(T_{n}-, I\right)=\tau_{n}\left(I-L_{n}\left(T_{n}\right)\right) P & \left(T_{n}, T_{n+1}\right) \\
& +U_{n}\left(T_{n},\left(I+s_{n}-g_{n} \times L_{n}\left(T_{n}\right)\right)^{+}\right)
\end{aligned}
$$

combined with the following one-period rollback scheme,

$$
\begin{equation*}
U_{n}\left(T_{n}, I\right)=\mathrm{E}_{T_{n}}\left(\frac{B\left(T_{n}\right)}{B\left(T_{n+1}\right)} U_{n+1}\left(T_{n+1}-, I\right)\right) \tag{18.69}
\end{equation*}
$$
^eq18-69

$n=N-1, \ldots, 1$. For the hold value, the continuity condition is

$$
\begin{equation*}
H_{n}\left(T_{n}-, I\right)=H_{n}\left(T_{n},\left(I+s_{n}-g_{n} \times L_{n}\left(T_{n}\right)\right)^{+}\right) \tag{18.70}
\end{equation*}
$$
^eq18-70

where

$$
\begin{equation*}
H_{n}\left(T_{n}, I\right)=\mathrm{E}_{T_{"}}\left(\frac{B\left(T_{n}\right)}{B\left(T_{n+1}\right)} \max \left(U_{n+1}\left(T_{n+1}-, I\right), H_{n+1}\left(T_{n+1}-, I\right)\right)\right) \tag{18.71}
\end{equation*}
$$
^eq18-71

The PDE scheme may be implemented by discretizing $I$ over an appropriate range, solving PDEs (18.69)-(18.71) on each $I$-plane, and interfacing the solutions between $I$-planes at times $T_{n}, n=1, \ldots, N-1$ using (18.68)(18.70). The details of implementation follow the general plan of Section 2.7.5, and we do not repeat them here.

