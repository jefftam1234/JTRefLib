---
tags: [quant/knowledge, quant/books, interest-rates]
---

## Bermudan Swaptions

After our general discussion of callable Libor exotics in the previous chapter, we now turn our attention to an important subset of the generic CLE class, the Bermudan swaptions. Bermudan swaptions are among the most liquid exotic interest rate derivatives, and the demands they place on accuracy, fidelity and performance of term structure models have driven many advances in interest rate modeling. While the ideas and methods from the previous chapter all apply to Bermudan swaptions, the simpler structure of Bermudan swaptions compared to general CLEs allows us to considerably deepen our analysis of valuation and risk management methods.

### 19.1 Definitions

As defined in Section 5.12, a Bermudan swaption is a callable Libor exotic with the coupon paying a fixed rate, $C_{n}=k, n=1, \ldots, N-1$. Alternatively, we can consider it a Bermudan-style option to enter a simple fixed-for-floating swap. The fixed rate $k$ is often referred to as the strike of the Bermudan swaption. Exercise dates of a Bermudan swaption are typically ${ }^{1}$ a subset of a tenor structure $\left\{T_{n}\right\}_{n=0}^{N}$ that defines the underlying swap. In a standard structure, exercise is restricted to the dates $\left\{T_{n}\right\}_{n=s}^{N-1}$ where $s \geq 1$; as we explained in Section 5.12, the period up to $T_{s}$ is known as the lockout or no-call period. Recall that a Bermudan swaption on, say, a 10 year swap with a 2 year lockout period (at inception) is known as a " 10 no-call 2 ", or "10nc2", Bermudan swaption. For convenience (and without loss of generality) we assume in most of this chapter that all $\left\{T_{n}\right\}_{n=1}^{N-1}$ are, in fact, exercise dates. If the Bermudan swaption is exercised at time $T_{n}$, the exercise value, for a payer swap, is given by

$$
\begin{equation*}
U_{n}(t)=\sum_{i=n}^{N-1} \tau_{i} P\left(t, T_{i+1}\right)\left(L_{\imath}(t)-k\right) \tag{19.1}
\end{equation*}
$$
^eq19-1

[^45]where $k$ is the fixed rate. We note that $U_{n}(t)$ can be written as

$$
U_{n}\left(T_{n}\right)=A_{n}\left(T_{n}\right)\left(S_{n}\left(T_{n}\right)-k\right),
$$
where $A_{n}(t) \triangleq A_{n, N-n}(t)$ is the annuity, and $S_{n}(t) \triangleq S_{n, N-n}(t)$ is the swap rate for the swap into which one can exercise at time $T_{n}$ (see notations (4.8), (4.10)). The definition of hold values carries over unchanged from Chapter 18.

### 19.2 Local Projection Method

As Bermudan swaptions are liquid and their volume is relatively high, the performance advantages of PDE methods over Monte Carlo simulation lead many market participants to value Bermudan swaptions in low-factor Markovian models, using either finite difference grids or trees. A sound framework for the usage of low-dimensional Markovian models is provided by the local projection method for single-rate CLEs that we discussed in Section 18.4. The method takes a particularly simple form for Bermudan swaptions, as the underlying swaps have no optionality and only the volatility parameters of core swap rates $\left\{S_{n}(\cdot)\right\}_{n=1}^{N-1}$ are relevant. As we discussed before (in Section 13.1.8.1), we can view a Bermudan swaption as the option to choose the "best" among a collection of swap values observed on different dates. This implies that a Bermudan swaption value is driven by core volatilities, or volatilities of the core swap rates $\left\{S_{n}\left(T_{n}\right)\right\}_{n=1}^{N-1}$, and core correlations, or inter-temporal correlations of the core swap rates $\left\{S_{n}\left(T_{n}\right)\right\}_{n=1}^{N-1}$. Alternatively, of course, we can think of forward volatilities in place of inter-temporal correlations as the source of "exotic" risk in Bermudan swaptions, see Section 18.1.1.

The relative simplicity of the dependence of Bermudan swaptions on the volatility structure allows us to use models as simple as the one-factor Gaussian model (Section 10.1.2) for valuation and risk management. The time-dependent volatility is typically calibrated to core swaption volatilities, while the mean reversion is calibrated to inter-temporal correlations of core swap rates (see Section 13.1.8.3); these correlations could, for instance, be extracted from an LM model. In practice, it is not unusual to skip the last step - since Bermudan swaptions have been traded well before LM models (or other practical multi-factor models) were invented, a market practice has developed whereby the mean reversion of a Gaussian (or similar) model is used essentially as a free parameter, rather than implied from a global model. This practice continues today, with mean reversions often set to match the "market" prices of Bermudan swaptions that are sometimes observable, or quasi-market prices such as independently-produced averages of dealer-submitted prices of a few typical structures ${ }^{2}$. Another fairly popular

[^46]choice would set mean reversions to match caplet volatilities, although using caplets for mean reversion calibration is rather arbitrary and can sometimes lead to odd mean reversion curves (see for instance the discussion in Section 10.1.2.3). The practice can, however, perhaps be justified if caplets are used as hedging instruments for inter-temporal correlation or forward volatility; technical details are available in Section 13.1.8.2.

Turning to the issue of volatility smile, we recall that the Gaussian model basically has no control over it, and the model can only be calibrated to one ${ }^{3}$ volatility per expiry $T_{n}, n=1, \ldots, N-1$. A one-factor quasi-Gaussian (qG) model with local volatility (Section 13.1) would constitute an improvement, as it will also allow to capture the slopes of volatility smiles of core swaptions, in addition to volatilities at specific strikes. Finally, the stochastic volatility version of the one-dimensional qG model (Section 13.2) would essentially allow for (best-fit) calibration to all core swaption volatility smiles across all strikes. On balance, the local volatility qG model is probably sufficient for effective risk management of Bermudan swaptions, although we would of course choose the SV version if available computing power permits. Finally, a two-factor quadratic Gaussian model of Section 12.3 is also a viable choice for Bermudan swaption pricing.

While we are on the subject of the model choice, let us briefly comment on the discussion around what number of factors is appropriate for a Bermudan swaption model. While the usage of single-factor, or essentially single-factor models such as the qG model, for Bermudan swaption valuation is widespread, some argue that single-factor models significantly underprice Bermudan swaptions. The basic claim is that higher de-correlation in rates has a positive impact on Bermudan swaption prices (as a Bermudan swaption is a "best-of" option on swap rates) and two- and multi-factor models intrinsically are able to de-correlate rates more than a single-factor model (where the instantaneous correlations between moves in all forward rates is always one). There are a number of flaws in this argument, starting with the fact that the correlations relevant for Bermudan swaptions are the inter-temporal correlations, which can be easily manipulated in a one-factor model through the choice of mean reversion. In addition, when comparing one- and multifactor rates models, it is obviously important that calibration to European swaptions is unaffected by changes in the number of factors. Careful analysis in Andersen and Andreasen [2001] of Bermudan swaption pricing in a oneand a two-factor Gaussian models shows that, if the models are calibrated in consistent fashion to core European swaptions, the two-factor model price is in fact slightly lower than the one-factor price. Experiments with LM models with different numbers of Brownian motions, but all calibrated to the full swaption grid, confirm this analysis.

[^47]The reason for the slight decrease in Bermudan swaption price as a function of the number of factors may seem puzzling at first and is, indeed, a rather subtle effect. While there are numerous factors in play (see Andersen and Andreasen [2001] for a full analysis), one important observation is that forward volatility in a low-factor model generally is higher than in a multifactor model, as long as both models are in calibration with the European swaption market. A technical explanation for this phenomenon can be found in Appendix 19.A of this chapter; loosely speaking the effect stems from the fact that a one-factor model will imply a lower time 0 instantaneous forward volatility term structure than will a multi-factor model, a relationship that must be reversed as time progresses to preserve overall variance.

### 19.3 Smile Calibration

Let us now discuss the issues of smile calibration in more detail. For concreteness, we consider the linear local volatility version of the quasi-Gaussian model from Section 13.1. As we have seen before, the model has enough flexibility to match the level and the slope of the volatility smile for each of the core swap rates. The market volatility smile is, of course, not close to linear, so we often seek to match the volatility of a particular strike exactly, or as closely as possible, while roughly matching the overall slope of the volatility smile.

The simplest approach to choosing the strikes that define core volatilities for calibration involves using at-the-money (ATM) strikes for each core swap rate. This rather crude approach is still in use for (we assume) historical reasons, as Bermudan swaptions started trading well before pronounced market smiles developed in interest rate markets, and probably even before the non-ATM points of the swaption volatility cube became liquid enough to keep track of them. Proponents of this approach sometimes rely on hedging arguments, as volatility exposure of a Bermudan swaption is often vega hedged using the most liquid European swaptions - which happen to be at-the-money. Yet another possible justification for the ATM strike choice notes that using the same volatilities for Bermudan swaptions of different strikes (seemingly) ensures consistency of valuation across Bermudan swaptions of different strikes. In reality, however, the ATM strike choice leads to inconsistent valuation between European and Bermudan swaptions: if one uses a Bermudan swaption model calibrated to ATM swaptions, and applies it to a Bermudan swaption with a non-ATM strike and just a single exercise date, the value is going to be different from the value of the same derivative priced as a European swaption. Clearly, this is a strongly undesirable feature of the ATM strike calibration idea.

To ensure consistency between the Bermudan swaption and its underlying core European swaptions, it suffices to set the calibration strike equal to that of the Bermudan swaption itself. That is, if the fixed rate of the

Bermudan swaption is $k$, then for each expiry $T_{n}$ one uses the volatility of the appropriate core European swaption that corresponds to the strike $k$. This method automatically ensures that a European swaption valued as a single-exercise Bermudan swaption has exactly the same value in the model as in the market. The fact that all swaptions we can exercise into are priced exactly is intuitively appealing, and also ensures that certain rational bounds for the Bermudan swaption price will not be violated. Indeed, if $V_{\text {swaption, } n}(0 ; k)$ is the price of a $k$-strike, $T_{n}$-expiry swaption on a swap that matures at time $T_{N}$, then clearly ${ }^{4}$ the $k$-strike Bermudan swaption price $V_{\text {Berm }}(0 ; k)$ at time 0 must satisfy (compare to (18.3))

$$
\begin{equation*}
V_{\text {Berm }}(0 ; k) \geq \max _{n=1, \ldots, N-1} V_{\text {swaption }, n}(0 ; k), \tag{19.2}
\end{equation*}
$$
^eq19-2

where we as mentioned earlier have assumed that exercise can take place at all $T_{n}, n \geq 1$. If our model fundamentally matches all swaption prices inside the max-operator on the right-hand side of (19.2), then pricing the Bermudan swaption in, say, a finite difference grid will always return a Bermudan swaption that satisfies (19.2). We notice as an aside that the (non-negative) difference between the left- and right-hand sides of (19.2) is sometimes known in trader jargon as the Bermudanality of the Bermudan swaption.

While enforcing consistency with European swaptions is useful, the idea of at-the-strike calibration is not a panacea, as Bermudan and European swaptions can behave quite differently. For instance, Bermudan swaptions have other "interesting", e.g. high-convexity, points in the swap rate dimension than just the underlying strike $k$, the most important of which is the exercise boundary, i.e. the swap rate level (for each time $T_{n}$ ) at which the decision to exercise the swaption switches to the decision to hold. The importance of this point is clearly seen from the marginal exercise value decomposition (18.8), as it corresponds to "strikes" of European options in the representation of the value of a Bermudan swaption as a sum of European options. Hence, a third calibration option available is to use the swap rate volatilities that correspond to the exercise boundary on each of the exercise dates. As was the case for the Bermudan strike method, this method makes valuation of single-exercise Bermudan swaptions consistent with the valuation of (equivalent) European swaptions, since the exercise boundary for a European swaption coincides with its strike. For the same reason, any weighted average of the strike and the exercise boundary would also produce a consistent scheme. We do point out, however, that using calibration strikes other than that of the Bermudan swaption may lead to violations of (19.2).

[^48]To provide a bit more detail on the idea of using the exercise boundary to select calibration strikes, let us first observe that for any given model, one can determine the exercise boundary as a function of the state variables. To be able to calibrate to European swaption volatilities with the strike at the exercise boundary, one has to be able to translate this "model" exercise boundary into a value of the corresponding core swap rate. Strictly speaking, this can be done unambiguously only in single-factor models, such as the Gaussian model. For the one-factor qG model with its two state variables the boundary is, in fact, represented by a line in a two-dimensional plane of possible values of the $x$ and $y$ state variables; each of the points on this line corresponds, potentially, to a different value of the core swap rate. However, the dependence of the exercise boundary on $y$ is rather mild, a fact that should come as no surprise if one recalls the "auxiliary" nature of the $y$ state variable, see Chapter 13. Hence, for the qG model, we can use the expected value of $y\left(T_{n}\right)$ when converting the exercise boundary for the state variable $x\left(T_{n}\right)$ into the strike for the swap rate $S_{n}\left(T_{n}\right)$.

The choice of exercise boundary for calibration is, unfortunately, rather inconvenient from the implementation point of view because the exercise boundary information is not available until after the valuation algorithm has been run. One can try a recursive scheme where one uses some (e.g., strikecalibrated) volatilities for an initial calibration, values a Bermudan swaption, calculates the exercise boundary, looks up the core volatilities for calibration at this boundary, calibrates the model again, values the Bermudan swaption, and so on. Such a procedure can in fact diverge; thus, one is forced to limit the number of iterations artificially, potentially resulting in unstable risk sensitivities and other problems. Moreover, this scheme in general consumes more computational resources due to the multiple valuations required. For all these reasons the at-the-boundary calibration is of limited use and the at-the-strike volatility calibration method is probably the most reasonable in practice, combining ease of implementation and consistency with European swaptions.

### 19.4 Amortizing, Accreting, and Other Non-Standard Bermudan Swaptions

A standard (or vanilla or bullet) Bermudan swaption is characterized by the fact that the notionals on which coupons are paid are all identical, as in (19.1) (where the notionals of all coupons are 1). A relatively popular extension involves making the notional of a Bermudan swaption time-dependent and deterministic, with the exercise value given by

$$
\begin{equation*}
\tilde{U}_{n}(t)=\sum_{i=n}^{N-1} R_{i} \tau_{i} P\left(t, T_{i+1}\right)\left(L_{i}(t)-k\right) \tag{19.3}
\end{equation*}
$$
^eq19-3

(compare to (19.1)). Here $R_{i}$ is the notional for the $i$-th coupon, $i= 1, \ldots, N-1$.

If the notional increases with the coupon index, the Bermudan swaption is said to be accreting; if it decreases, it is said to be amortizing. Other profiles are possible but are much less common. Amortizing Bermudan swaptions are often used as hedges for pools of mortgages, with the amortization feature mimicking prepayments on the pool. Accreting Bermudan swaptions, on the other hand, often appear as a result of issuing "zero coupon" structured notes, i.e. notes with the repayment notional growing over time but paying no coupons during the life of the note, as explained in Section 19.4.6.

Since the notionals in Bermudan swaptions of type (19.3) are still deterministic (even if time-varying), their valuation in a properly calibrated model does not present any particular technical difficulties ${ }^{5}$ - but what constitutes "properly calibrated", however, is not always obvious. Of course, in models with global calibration (e.g. LM models), calibration for non-standard Bermudan swaptions is no different from calibration for standard Bermudan swaptions, as model calibration is product-independent by definition. On the other hand, for models requiring local calibration, such as a one-factor Gaussian or a quasi-Gaussian model, calibration for non-standard Bermudan swaptions will require additional analysis. We consider this problem in the next few sections, but it is worth pointing out that, in the opinion of some, making notionals time-dependent pushes Bermudan swaptions across the boundary that separates those securities for which local models are acceptable to use from those for which globally-calibrated multi-factor models are required.

Before commencing on a more detailed analysis, let us first briefly try to understand the basic complications involved in local model calibration for non-vanilla Bermudan swaptions. As established previously, a locallycalibrated model should, as a minimum, be calibrated to the volatilities of core swap rates. For an amortizing Bermudan swaption, say, a core swap rate would correspond to an amortizing swap. Volatilities of amortizing swaps can be extracted from amortizing European swaptions, but the liquidity of such swaptions is significantly poorer than for vanilla European swaptions in fact, amortizing European swaptions are about as illiquid as Bermudan swaptions themselves. In practice, if one wishes to calibrate a local model to core (amortizing) swaptions, one may need to use a "pre-processing" step to extract amortizing European swaption prices from a model calibrated to liquid vanilla European swaptions, as we do later in Section 19.4.4. Alternatively, one needs to choose a different set of calibration targets in the first place, see Section 19.4.3.

[^49]
### 19.4.1 Relationship Between Non-Standard and Standard Swap Rates

Regardless of the calibration method ultimately used, it is useful to understand the relationship between non-standard and standard swaps (and, hence, swap rates). To be consistent with (19.3), consider a swap that starts at $T_{n}$, ends at $T_{N}$, and has a notional schedule $\left\{R_{i}\right\}$. The time $t$ value of such a swap is given by $\widetilde{U}_{n}(t)$ in (19.3). To make some of the formulas below simpler, let us extend the notional schedule by one period and set $R_{N}=0$. We denote the annuity and the swap rate that correspond to this non-standard swap by $\widetilde{S}_{n}(\cdot)$ and $\widetilde{A}_{n}(\cdot)$, so that

$$
\widetilde{A}_{n}(t)=\sum_{i=n}^{N-1} R_{i} \tau_{i} P\left(t, T_{i+1}\right), \quad \widetilde{S}_{n}(t)=\widetilde{A}_{n}(t)^{-1} \sum_{i=n}^{N-1} R_{i} \tau_{i} P\left(t, T_{i+1}\right) L_{i}(t)
$$

We would like to decompose the non-standard swap into a linear combination of standard swaps. Such a decomposition is, however, not unique and could be done in a multitude of ways, potentially using any of the standard swaps with starting date on or after $T_{n}$, and final payment date on or before $T_{N}$. To narrow down the problem, we note that we here are ultimately interested in establishing the volatility of the non-standard rate $\widetilde{S}_{n}(\cdot)$ over the period $\left[0, T_{n}\right]$. Since the values of standard European swaptions provide us with the information on the volatilities of swap rates only over the period from time 0 to their start dates, we should focus only on standard swaps that start on $T_{n}$; this choice makes the decomposition unique.

Let us denote the value of a standard swap starting at $T_{n}$ and covering $m$ periods by $V_{n, m}(t)$, and the corresponding annuity and swap rate by $S_{n, m}(t)$ and $A_{n, m}(t)$ (see (4.8), (4.10)). In light of the discussion above, we want to find weights $\left\{v_{n, m}\right\}, m=1, \ldots, N-n$, such that

$$
\tilde{U}_{n}\left(T_{n}\right)=\sum_{m=1}^{N-n} v_{n, m} V_{n, m}\left(T_{n}\right)
$$

Note that only swaps starting at time $T_{n}$ are used in the right-hand side of this expression. Matching terms to (19.3) we obtain

$$
\begin{equation*}
v_{n, m}=R_{n+m-1}-R_{n+m}, \quad m=1, \ldots, N-n, \tag{19.4}
\end{equation*}
$$
^eq19-4

so that

$$
\widetilde{U}_{n}\left(T_{n}\right)=\sum_{m=1}^{N-n}\left(R_{n+m-1}-R_{n+m}\right) V_{n, m}\left(T_{n}\right) .
$$

After some algebraic manipulations, we obtain the following relationship for the swap rates,

$$
\begin{equation*}
\widetilde{S}_{n}\left(T_{n}\right)=\sum_{m=1}^{N-n} w_{n, m}\left(T_{n}\right) S_{n, m}\left(T_{n}\right) \tag{19.5}
\end{equation*}
$$
^eq19-5

where (recall that we set $R_{N}=0$ )

$$
w_{n, m}\left(T_{n}\right)=\left(R_{n+m-1}-R_{n+m}\right) \frac{A_{n, m}\left(T_{n}\right)}{\widetilde{A}_{n}\left(T_{n}\right)} .
$$

While the swap weights $v_{n, m}$ are deterministic, the swap rate weights $w_{n, m}\left(T_{n}\right)$ are not. For a qualitative discussion, however, we note that the weights $w_{n, m}\left(T_{n}\right)$ can be approximated reasonably well by their values at time 0 ,

$$
\begin{equation*}
\widetilde{S}_{n}\left(T_{n}\right) \approx \sum_{m=1}^{N-n} w_{n, m}(0) S_{n, m}\left(T_{n}\right), \quad w_{n, m}(0)=\left(R_{n+m-1}-R_{n+m}\right) \frac{A_{n, m}(0)}{\widetilde{A}_{n}(0)} . \tag{19.6}
\end{equation*}
$$
^eq19-6

From the expression (19.6) it follows that the volatility of the non-standard swap rate is a function of volatilities of all standard swap rates with a given expiry ( $T_{n}$ in our case), and of their correlations. Putting correlations aside for a moment, observe that to price a non-standard Bermudan swaption, in principle one needs to calibrate the model to volatilities of standard rates with all expiries ( $T_{1}, \ldots, T_{N-1}$ ) and all maturities, something a lowdimensional local model will virtually never be able to do. Below, we discuss two possible strategies for going forward.

### 19.4.2 Same-Tenor Approach

Perhaps the simplest calibration approach for non-standard Bermudan swaptions is to simply pretend that they are standard Bermudan swaptions and set up the model calibration accordingly. So, for expiry $T_{n}$, one would choose a European swaption on the $(N-n)$-period swap as the calibration instrument. While easy to implement, the merits of this approach are obviously somewhat wanting, and we do not recommend it. Nevertheless, it is instructive to investigate the issues that would come up if we adopted this scheme. As an example, consider an amortizing Bermudan swaption and a one-factor model calibrated to standard swaptions of the same tenor as the core amortizing swaptions. As a thought experiment, suppose that we increase mean reversion in the model, while keeping it calibrated to our calibration swaption set. In this case, the core amortizing swaption prices would increase, a simple consequence of our decomposition (19.6) and the fact that shorter-tenor (standard) swap rate volatilities increase as a function of mean reversion when volatility of a longer-tenor swap rate is kept fixed, see the discussion in Section 13.1.8.1. As a consequence, mean reversion would affect not only the inter-temporal correlations that are important for Bermudan swaptions, but would also affect the volatilities of core swap rates. In the context of a local projection method, we would then face a dilemma as to which targets to calibrate the mean reversion to: the inter-temporal correlations or the prices of amortizing European swaptions (the latter, just
like the former, would be available from a global model). Of course it would be highly unlikely that both calibration targets would imply the same mean reversion.

Volatility smile calibration presents another challenge for the same-tenor approach. For instance, if one chooses a particular strike of the non-standard swaption to calibrate to (e.g. the fixed rate of the non-standard swap), which strike for the standard swaption would that correspond to?

### 19.4.3 Representative Swaption Approach

The idea of the representative swaption approach is to choose a standard swap that approximates the non-standard swap in some reasonable sense, and then to calibrate the Bermudan swaption model to the market-implied volatilities of swaptions on these standard swaps, one per exercise date.

One can define a "representative" swap in many ways. For example, a fairly simple $P V B P$ matching method chooses the standard swap whose PVBP (Present Value of a Basis Point, see Section 5.5) matches the PVBP of the non-standard swap most closely. In this case, for expiry $T_{n}$, we would choose the tenor $\mu_{n}$ of the standard swap by

$$
\begin{equation*}
\mu_{n}=\underset{m}{\operatorname{argmin}}\left\{\left|R_{n} A_{n, m}(0)-\widetilde{A}_{n}(0)\right|\right\} \tag{19.7}
\end{equation*}
$$
^eq19-7

(note that we use the notional of the first period of the swap $\tilde{U}_{n}$ to scale the PVBP of the standard swap). While somewhat simplistic, the method actually turns out to be reasonably robust for some non-standard Bermudan swaptions. We proceed to improve it and make more rigorous, which will also help us identify not just the right tenor, but also the right strike for the standard calibration swaptions.

We work in the context of a one-factor Gaussian model to demonstrate the main idea, although the method is not tied to a particular model. Let us fix a start date $T_{n}$ and note that the value of (any) swap at time $T_{n}$ is a function of the Gaussian short rate state $x=x\left(T_{n}\right)$ on that date. Let $\tilde{U}_{n}(x)$ be the value of the non-standard swap $\tilde{U}_{n}\left(T_{n}\right)$, as a function of the short rate state. Note that $\tilde{U}_{n}(x)$ depends on the mean reversion, but not the volatility parameter of the model, as follows from bond reconstruction formulas (Proposition 10.1.7). Define $V(x ; R, q, m)$ to be the value of a standard swap starting on $T_{n}$ as a function of $x$, with constant notional $R$, fixed rate $q$, and covering $m$ periods. In departure from our normal conventions, we here allow $m$ to be any real number and not necessarily an integer; we interpret a value of, say, $m=5.3$ as 5 full periods plus three tenths of the sixth period. The rationale for allowing for fractional periods will become clear shortly.

Now, we have three parameters that define the standard swap, $R, q$, and $m$. In the payoff matching method, we choose the three parameters by
matching the level, slope and curvature of the swap payoffs as functions of the state variable,

$$
\begin{align*}
V\left(x_{0} ; R, q, m\right) & =\widetilde{U}_{n}\left(x_{0}\right)  \tag{19.8}\\
\frac{\partial}{\partial x} V\left(x_{0} ; R, q, m\right) & =\frac{\partial}{\partial x} \widetilde{U}_{n}\left(x_{0}\right)  \tag{19.9}\\
\frac{\partial^{2}}{\partial x^{2}} V\left(x_{0} ; R, q, m\right) & =\frac{\partial^{2}}{\partial x^{2}} \widetilde{U}_{n}\left(x_{0}\right) \tag{19.10}
\end{align*}
$$
^eq19-8
^eq19-9
^eq19-10

where the expansion point $x_{0}$ is the expected value of $x\left(T_{n}\right)$ (or close to it). In the parameterization of Section 10.1.2, it suffices to set $x_{0}=0$. The system of equation (19.8)-(19.10) is easy to solve (numerically) in the one-factor Gaussian model; let us denote the solution by $R_{n}^{*}, q_{n}^{*}, m_{n}^{*}$.

Even though we fix the parameters of the standard swap by local conditions around $x_{0}$, numerical experiments show that the swap that solves (19.8)-(19.10) tends to match that of the non-standard swap across a large range of state values $x\left(T_{n}\right)$, suggesting considerable robustness. In addition, even though the functions $V$ and $\widetilde{U}_{n}$ depend on mean reversion, numerical experiments show that the best-fit parameters $R_{n}^{*}, q_{n}^{*}, m_{n}^{*}$ are only mildly sensitive to mean reversion.

Incorporating the payoff matching method into a volatility calibration routine is quite easy, since the choice of the best-fitting standard swaps is independent of volatility, which allows us to identify the calibration targets before we commence on the volatility calibration. Moreover, the strike of the calibration swaptions are produced automatically as part of the payoff matching routine, facilitating calibration to a particular point of the observed volatility smile.

Before discussing application of the representative swaption idea to accreters and amortizers, let us briefly motivate our usage of fractional swap tenors. If tenors are restricted to an integer number of periods, then a perturbation of the market data, e.g. when shifting a yield curve to calculate an interest rate delta, could potentially alter the solved-for number of periods $m$ by plus or minus one period. Hence, restricting $m$ to be an integer would potentially introduce a discontinuity in the calibrated model parameter and therefore in the value of the Bermudan swaption - as a function of market data. Such a discontinuity would be purely artificial and, as explained at length in Chapter 23, highly undesirable for stability of risk sensitivities. By allowing fractional tenors, we eliminate these problems. Of course, to complete the volatility calibration, we need to know implied volatilities of swaptions with fractional tenors, but these could be obtained by (smoothly!) interpolating implied volatilities of swaptions with integer-valued tenors.

Now, let us see how the representative swaption method works for an amortizing Bermudan swaption, i.e. a Bermudan swaption with decreasing notionals $R_{1} \geq R_{2} \geq \ldots \geq R_{N-1}$. We note that, according to (19.4), all weights $v_{n, m}$ in the decomposition of the amortizing swap into standard
swaps are positive, $v_{n, m} \geq 0, m=1, \ldots, N-n$. It follows that both the PVBP matching and payoff matching methods produce a standard swap whose tenor is some average of tenors of the standard swaps in the basket; in particular, the resulting standard swap have a final maturity that is shorter than the amortizing swap. This is an intuitive result, and leads to an amortizing Bermudan swaption being sensitive to interest rate volatilities of standard swaps $V_{n, m_{n}^{*}}, n=1, \ldots, N-1$, with $n+m_{n}^{*} \leq N$ for any $n=1, \ldots, N-1$.

The situation is different for accreting Bermudan swaptions, i.e. when $R_{1}<R_{2}<\ldots<R_{N-1}$. According to (19.4),

$$
\begin{array}{ll}
v_{n, m}<0, & m=1, \ldots, N-n-1, \\
v_{n, m}>0, & m=N-n .
\end{array}
$$

So, an accreting swap decomposes into a standard swap of the maximum tenor $N-n$ and maximum notional $R_{N-1}$, minus a basket of swaps of smaller tenors. The PVBP $\widetilde{A}_{n}(0)$ of an accreting swap is larger than the PVBP of a standard swap of the matching tenor (times starting notional $R_{n}$ ) so the PVBP matching method would calculate an optimal tenor $m^{*}$ that is longer than the tenor of the amortizing swap, $m^{*}>N-n$. The same would be true of the payoff matching method as well. This is, of course, rather problematic, as our calibration method would suggest that an accreting Bermudan swaption is sensitive to volatilities of swaptions with total length (the sum of expiry and tenor) exceeding the final maturity of the Bermudan swaption. This would be in direct contradiction to what, say, a globally-calibrated LM model would suggest, as in the latter the price of any derivative is fully determined by the volatility structure of Libor rates that fix before the final maturity of the derivative, and this volatility structure, in turn, is fully determined by the volatilities of swaptions with total length less than the final maturity.

The reader may ask why we are getting reasonable results for amortizing swaptions and unreasonable ones for accreting swaptions. A bit of reflection reveals that the discrepancy originates with the single-factor assumption that we made implicitly in the PVBP matching method, and explicitly in the payoff matching method. For the amortizing swap, the decomposition resulted in a basket of standard swaps with positive weights, a basket that can be reasonably well-hedged with a single swap of average tenor. In the accreting case, our decomposition resulted in a spread position in standard swaps: long a long-tenor swap and short a basket of short-tenor swaps. Our one-factor methods suggest hedging this spread position with a single (very) long-dated swap - perfectly reasonable in a one-factor world, but not in actual reality.

While sensibly hedging the spread position in an accreting swap with a single standard swap is not possible, things improve markedly if we allow usage of two standard swaps in the hedge. In particular, we may then take
as one of our hedges a long position in the $T_{N}$-maturity swap with maximum ( $R_{N-1}$ ) notional, and construct the second swap hedge by PVBP (or payoff) matching the remaining short basket of swaps $\sum_{m=1}^{N-n-1} v_{n, m} V_{n, m}$. As all weights in the short basket are of the same (negative) sign, a single standard swap would often provide a good hedge. Thus, to get a reasonable calibration scheme for an accreting Bermudan swaption, we would need to calibrate to two standard European swaptions per expiry, both of which would have their final payment date on or before the final payment date of the Bermudan swaption. Of course, we would find it difficult to accurately calibrate a one-factor Markovian model to two swaptions per expiry, and would likely need to move on to more elaborate models with additional factors.

In conclusion, we note that the two-swaps approach works universally for Bermudan swaptions with arbitrary notional schedules. To apply it, for each $n$ we would combine swaps with positive weights $v_{n, m}$ into one basket, and swaps with negative weights $v_{n, m}$ into another basket. Then we would represent each basket by one standard swap by the procedures discussed above, yielding the calibration swaption targets for that expiry.

### 19.4.4 Basket Approach

The discussion above suggests that the pricing of at least some Bermudan swaptions with non-standard notional schedules is best done in multi-factor models. Still, one-factor Markovian models are highly popular due to their performance advantages, and the desire to use them even in situations where they might be overstretched is often considerable. Consequently, rather creative ways of using one-factor models for non-standard Bermudan swaptions have been developed, resulting in a family of approaches that we here all categorize as basket methods.

The basket methods generally split the valuation of a non-standard Bermudan swaption into two stages. During the first stage, some model is used to calculate values of core non-standard European swaptions. During the second stage, a one-factor model is calibrated to these values of nonstandard core European swaptions, and subsequently used to compute the value of the non-standard Bermudan swaption. Various method differ in how the values of non-standard European swaptions are calculated. One perfectly sound method uses a globally calibrated model such as the LM model for the task, resulting in a local projection method for non-standard Bermudan swaptions. We have discussed the local projection method in various flavors often enough, so we trust the reader with filling in remaining details. Instead, we review some alternatives for how to execute the first stage of the basket method.

For concreteness, let us focus on the first non-standard European swaption underlying the Bermudan swaption, i.e. the option expiring at $T_{1}$ on a swap that covers $N-1$ periods. As follows from (19.4), this non-vanilla European swaption can be interpreted as an option on a basket of standard swaps,
all starting on $T_{1}$ but with different maturities. Hence, to compute the price of the non-standard swaption, one can use a model calibrated to the volatilities of options on such swaps, as well as relevant swap rate correlations. Notice that the standard swaptions involved here all form a "row" of the swaption grid (see Section 5.10), as they all share the same expiry but have different tenors. A one-factor mean-reverting Gaussian (or quasi-Gaussian or quadratic Gaussian) model can be calibrated to this set of European swaptions, although the calibration will be different from our standard procedure. In particular, the prices of all swaption targets depend on the model volatility function over the same interval $\left[0, T_{1}\right]$ and, thus, the short rate volatility function (e.g. $\sigma_{r}(\cdot)$ in the notation of Proposition 10.1.7) cannot be used if we want to match each swaption volatility exactly. Upon reflection, it should be clear that we instead can use the time-dependent mean reversion function ( $\varkappa(\cdot)$ in the notation of Proposition 10.1.7) as our main calibration "knob", since the pricing of a $T_{1}$-expiry swaption on a swap that covers $m$ periods will depend on the mean reversion function over the period $\left[T_{1}, T_{1+m}\right]$. Hence, a sequential mean reversion calibration is possible: after calibrating to the first $m$ standard swaptions, the ( $m+1$ )-th is matched by changing the mean reversion function ${ }^{6}$ over the time interval $t \in\left[T_{m+1}, T_{m+2}\right]$, for $m=0, \ldots, N-2$. The (constant) level of volatility over the first period could be set arbitrarily, as its scaling effect would be compensated by the mean reversion calibration. However, it is advisable to keep it at a "typical" value of, say, $1 \%$ so the calibrated mean reversions would also remain in a "typical" region, as the numerical implementation of the model might not cope well with extreme values of mean reversion.

To summarize, the basket method for a mean-reverting one-factor short rate model works like this. First, for each row of the swaption grid that corresponds to an exercise date $T_{n}, n=1, \ldots, N-1$, of the non-vanilla Bermudan swaption, we fit separate instances of the one-factor model by sequentially calibrating the mean reversion function to all relevant swaptions in the row. For each $T_{n}$, the relevant instance of the model is then used to compute the price of the $T_{n}$-expiry non-vanilla European swaption that the Bermudan swaption can be exercised into. Finally, we calibrate the model once more by setting its short rate volatility function to match the prices of the non-standard European swaptions established in the previous step. In the final calibration, we would typically keep the mean reversion fixed, either at a user-specified level or (ideally) at a level that makes inter-temporal correlations of core swap rates match those coming from a global model; see Section 19.2 for additional discussion.

So far we have side-stepped the issue of what strikes to choose for various calibrations in the basket scheme above. It is fair to say that this choice is a non-trivial problem. Various ad-hoc schemes could be imagined, such as using

[^50]the standard swaption of the same relative moneyness as the non-standard one, but they are rarely entirely satisfactory. Using a quasi-Gaussian or quadratic Gaussian model (or, even better, a stochastic volatility extension of these models) is obviously preferable to, say, using a Gaussian model, as the former models will allow us to calibrate to the volatility smile at more than a single strike, thereby alleviating somewhat the strike selection problem.

Another issue that we should touch on concerns correlations. By using a one-factor model for establishing non-standard European swaption values, we are implicitly assuming high ${ }^{7}$ correlations between standard swaps in the portfolio that replicates the exercise value of the swaption. This is not necessarily as constraining as it might appear to be, as in reality these swaps are indeed highly correlated. Still, we may want to contemplate methods to somehow incorporate into our procedure observations about correlation extracted, say, from historical analysis. One possible route for this would be to apply approaches inspired by basket valuation methods from equities modeling. For example, we could (rather crudely) value a non-standard European swaption by the Black formula on the (non-standard) swap rate whose volatility is obtained by moment matching ${ }^{8}$. This method would need to approximate the weights in the decomposition (19.5) as being deterministic (although they are not). We can also use more advanced approaches, such as the copula methods, or even SV methods, that we developed for multirate derivatives in Chapter 17, allowing us to incorporate volatility smile information into the basket valuation.

A few final comments on the method developed in this section are in order. First, it is worth pointing out that most of the basket methods are consistent with the way standard Bermudan swaptions are valued. Specifically, if we apply these methods to standard Bermudan swaptions, we obtain the same price as if we had used the "standard" valuation method of Section 19.2. Second, notice that the method produces volatility sensitivities for nonstandard Bermudan (and European, for that matter) swaptions that tend to be intuitive and reasonable. In particular, each underlying European non-standard swaption will show sensitivities only to the correct row of the swaption grid and to ordinary swaptions on swaps with maturities that do not exceed the maturity of the swap in the corresponding non-standard swaption. For an accreting European swaption in particular with, say, an expiry $T_{n}$ and swap maturity $T_{N}$, the sensitivities will be negative for all standard swaptions with expiry $T_{n}$ and swap maturities $T_{i}, i=n+1, \ldots, N-1$, and will be positive for a swaption with expiry $T_{n}$ and swap maturity $T_{N}$, thus faithfully representing the accreting swap rate as a "spread".

[^51]
### 19.4.5 Super-Replication for Non-Standard Bermudan Swaptions

The replication method of Proposition 8.4.13 links the value of a European option with a non-standard payoff to that of a portfolio of standard European options. Not only does that give us a way to value a non-standard derivative, it also allows us to fully hedge it in a model-independent way. Such static replication results are extremely convenient, and much research in derivatives pricing theory have been directed towards the search for static hedges for exotic derivatives, see e.g. Andersen et al. [2002]. Unfortunately, the availability of truly model-free static replication methods for non-European options is an exceedingly rare phenomenon and no such results are known to exist for Bermudan swaptions. Interestingly, however, Bermudan swaptions with non-standard notional schedules can be super-replicated in a modelindependent way, in the sense that for any given non-standard Bermudan swaption we can find a portfolio of standard Bermudan swaptions that would dominate the value of the non-standard Bermudan swaption in all states of the world. Morcover, in some cases the difference in value between the nonstandard Bermudan swaption and its super-replicating portfolio can be quite small. While not as convenient as a replicating portfolio, the super-replicating portfolio has several practical uses. First, the value of the portfolio provides a hard no-arbitrage bound for the value of the non-standard Bermudan swaption, and any modeling procedure (including calibration rules, choice of strikes, etc.) can be checked against this bound. Second, if the upper bound provided by the portfolio is known to be relatively tight (as is the case for, e.g., amortizing Berniudan swaptions, as we shall see shortly), then the super-replicating portfolio can be used directly for valuation purposes. perhaps amended with a small ad-hoc adjustment. More importantly, the super-replicating portfolio can be used as a robust hedge that requires little rebalancing over time.

The easiest way to demonstrate the construction of a super-replicating portfolio is by example. Consider first an amortizing Bermudan swaption. For concreteness, assume it is a 10 year ( 10 y ) Bermudan swaption with exercises every year, starting in year 1. Suppose the initial notional is 10 and it decreases by 1 every year. The super-replicating portfolio then consists of nine standard Bermudan swaptions, all with unit notional: a 10 no-call 1 , a 9 no-call $1, \ldots$, a 3 no-call 1 , and a 2 no-call 1 . To see that this is indeed a super-replicating portfolio, suppose the amortizing Bermudan swaption is exercised at year 5 . Then the option holder receives an amortizing 5 year swap, with a starting notional of 5 that decreases by 1 every year; the value of this swap is equal to the sum of standard swaps of tenors $5 \mathrm{y}, 4 \mathrm{y}, \ldots, 1 \mathrm{y}$, each with notional 1. But clearly the value of each of these standard swaps is dominated by the value of a corresponding standard Bermudan swaption in the super-replicating basket.

Another way to test that the super-replicating strategy dominates the amortizing Bermudan swaption is to impose a particular exercise strategy on
the portfolio of standard Bermudan swaptions. Specifically, we simply require that all still-alive (i.e. non-expired) standard Bermudan swaptions shall be exercised at the same time as when the amortizing Bermudan swaption is exercised. A little thought shows that the exercise value obtained from the super-replicating basket is then exactly the same as from the amortizing Bermudan swaption. Hence, the value of the portfolio of standard Bermudan swaptions with this specific exercise rule enforced must precisely equal the value of the amortizing Bermudan swaption. However, as the chosen exercise strategy will generally be sub-optimal for each of the standard Bermudan swaptions in the portfolio, the true value (i.e. the value obtained with optimal exercise) of the basket of standard Bermudan swaptions will be higher than the amortizing Bermudan swaption, and will dominate its value in all states of the world.

To give another example, consider a 10 year accreting Bermudan swaption with an initial notional of 1 that increases by 1 every year. Assuming that exercise can take place annually starting in year 1, the super-replicating portfolio will now consist of a collection of nine 10 year standard unit notional Bermudan swaptions exercisable annually, with lockout periods of $1,2, \ldots, 9$ years, respectively. To check that this hedge works, let us again assume the accreter is exercised in year 5. Then the holder would receive a swap that can be decomposed into a spot-starting 5 year standard swap, a 4 year standard swap starting in 1 y , a 3 year standard swap starting in 2 years and so on, all with notional 1. Again, the value of each of these swaps is dominated by a corresponding standard Bermudan swaption in the super-replicating portfolio.

Super-replicating portfolios for Bermudan swaptions with arbitrary notional schedule always exist, but are rather tedious to write down explicitly (for example, see a related algorithm in Evers and Jamshidian [2005]). The basic idea, however, is quite simple: for any exercise opportunity one needs to ensure that each of the standard swaps into which the exercise value can be decomposed (see (19.4)) is matched by a standard Bermudan swaption in the super-replicating portfolio.

Let us show some numerical results as a way to examine the tightness of the upper value bound produced by the super-replicating portfolio. For the numerical experiments, we throughout use a one-factor quasi-Gaussian model with some reasonable, representative parameters and a yield curve flat at $6 \%$. In Table 19.1 we show values for Bermudan (European) amortizing swaptions of different maturities (with 1y lockout) against the values of corresponding super-replicating portfolios of standard Bermudan (European) swaptions; all contracts are receivers (options on receive-fixed swaps) with $6 \%$ strike. The notionals of amortizing swaps decrease linearly from the initial notional indicated in the table by 1 every year. We notice that the upper bounds produced by the super-replicating portfolio are here quite tight, something that appears to hold generally for amortizing Bermudan swaptions.

| Maturity | 10 y | 10 y | 30 y | 30 y |
| :--- | :---: | :---: | :---: | :---: |
| Initial notional | 10 | 10 | 30 | 30 |
| Bermudan/European | E | B | E | B |
| Amortizer value | 0.606 | 0.630 | 2.180 | 2.700 |
| Portfolio value | 0.614 | 0.650 | 2.230 | 2.830 |

Table 19.1. The value of an amortizing Bermudan or European swaption vs. the value of a super-replicating portfolio of standard Bermudan/European swaptions for different maturities and contract types.

As a second test we look at a particular amortizing Bermudan swaption across a range of strikes, see Table 19.2. We consider a 30 year amortizing Bermudan swaption with the initial notional of 30 , and compare it to the super-replicating portfolio. For reference, the vega (change in Bermudan swaption value to $1 \%$ change of Black volatilities of European swaptions) for the $6 \%$ Bermudan swaption is about 0.1 , with European swaption implied volatilities around $12 \%$. Again, the results from the super-replicating portfolio are quite close to the real option values.

| Strike | $3 \%$ | $4 \%$ | $5 \%$ | $6 \%$ | $7 \%$ | $8 \%$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| Amortizer value | 0.267 | 0.632 | 1.372 | 2.700 | 4.640 | 7.017 |
| Portfolio value | 0.355 | 0.766 | 1.534 | 2.830 | 4.710 | 7.052 |

Table 19.2. The value of a 30 year amortizing Bermudan swaption vs. the value of a super-replicating portfolio of standard Bermudan swaptions for different strikes.

Next, we look at the results for accreting swaptions. In Table 19.3 we have the results for accreting receivers with $6 \%$ strike with notional accreting at $6 \%$ relative rate, across different contract types and maturities. Clearly, the super-replicating portfolio here is substantially more expensive than the accreting Bermudan swaption.

| Maturity | 10 y | 10 y | 30 y | 30 y |
| :--- | :---: | :---: | :---: | :---: |
| Initial notional | 1 | 1 | 1 | 1 |
| Bermudan/European | E | B | E | B |
| Accreter value | 0.096 | 0.113 | 0.121 | 0.255 |
| Portfolio value | 0.123 | 0.143 | 0.280 | 0.380 |

Table 19.3. The value of an accreting Bermudan or European swaption vs. the value of a super-replicating portfolio of standard Bermudan/European swaptions for different maturities and contract types.

Finally, we look at a particular 30 year accreting Bermudan swaption with initial notional of 1 across different strikes, see Table 19.4. For each contract, the notional compounds at the rate given by the fixed rate.

| Strike | $3 \%$ | $4 \%$ | $5 \%$ | $6 \%$ | $7 \%$ | $8 \%$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| Accreter value | 0.021 | 0.053 | 0.121 | 0.255 | 0.495 | 0.880 |
| Portfolio value | 0.021 | 0.063 | 0.164 | 0.380 | 0.773 | 1.400 |

Table 19.4. The value of a 30 year accreting Bermudan swaption vs. the value of a super-replicating portfolio of standard Bermudan swaptions for different strikes.

Judging by Tables 19.1-19.4, it appears that the super-replicating portfolio tends to produce a tighter bound for amortizers than for accreters. This observation, as it turns out, is not tied to the particular structures and market data used in the tables, but is true in general. To understand why this is the case, recall that the super-replicating portfolio for a non-standard Bermudan swaption will have the exact same value as the non-standard Bermudan swaption if each standard Bermudan swaption in the portfolio is exercised at the same time as the non-standard Bermudan swaption. In other words, the tightness of the upper bound for amortizers therefore suggests that it is optimal to exercise all standard Bermudan swaptions in the superreplicating portfolio at about the same time (for an amortizing Bermudan swaption, we recall that this portfolio consists of Bermudan swaptions with identical lockout periods, but different maturities). However, according to Proposition 19.7.1 proven later in this chapter, arbitrage arguments can be used to show that, at any exercise date, if one does not exercise the standard Bermudan swaption with the shortest tenor among remaining in the basket (i.e. with a remaining 1 year swap in the example above), then one should never exercise any of the remaining standard Bermudan swaptions (with tenors $2,3, \ldots$ years). In light of this result, the tightness of the super-replication bound is therefore not surprising. The same argument does not hold for the accreters, because for these structures the super-replicating portfolio consists of standard Bermudan swaptions with different lockout periods, rather than different underlying swap tenors. As such, it is obviously not reasonable to assume that the standard Bermudans will be optimally exercised at the same time.

Finally, let us briefly comment on the lower bound for the value of a nonstandard Bermudan swaption. While we know of no general results, the carry argument of Section 19.7.2 allows us to show that an amortizing Bermudan swaption with a final notional (i.e. the notional for the final exercise date) of 1 is bounded from below by the standard Bermudan swaption of notional 1, as long as we keep the exercise schedule, strike, etc. unchanged. This holds
in a model-independent way. Unfortunately, this lower bound is typically quite loose, except for Bermudan swaptions that amortize slowly.

### 19.4.6 Zero-Coupon Bermudan Swaptions

Let us momentarily turn to the question of where accreting Bermudan swaptions come from in the first place. Consider an accreting (receiver) Bermudan swaption with the notional defined by

$$
\begin{equation*}
R_{i}=\prod_{j=0}^{i-1}\left(1+\tau_{j} k\right), \quad i=1, \ldots, N-1 \tag{19.11}
\end{equation*}
$$
^eq19-11

According to (19.3), its $n$-th exercise value at time $T_{n}$ is given by

$$
\begin{align*}
\tilde{U}_{n}\left(T_{n}\right) & =\sum_{i=n}^{N-1} R_{i} \tau_{i} P\left(T_{n}, T_{i+1}\right)\left(k-L_{i}\left(T_{n}\right)\right) \\
& =\sum_{i=n}^{N-1}\left(\prod_{j=0}^{i-1}\left(1+\tau_{j} k\right)\right)\left(P\left(T_{n}, T_{i+1}\right)\left(1+\tau_{i} k\right)-P\left(T_{n}, T_{i}\right)\right) \\
& =\sum_{i=n}^{N-1}\left(R_{i+1} P\left(T_{n}, T_{i+1}\right)-R_{i} P\left(T_{n}, T_{i}\right)\right) \\
& =R_{N} P\left(T_{n}, T_{N}\right)-R_{n} \tag{19.12}
\end{align*}
$$
^eq19-12

where we have used the defining relation $L_{i}\left(T_{n}\right)=\left(P\left(T_{n}, T_{i}\right)-\right. P\left(T_{n}, T_{i+1}\right) /\left(\tau_{i} P\left(T_{n}, T_{i+1}\right)\right)$.

Now consider a contract in which an investor gives ${ }^{9}$ the dealer 1 at time $T_{0}$, and the dealer promises to pay the investor the amount of $R_{N}= \prod_{j=0}^{N-1}\left(1+\tau_{j} k\right)$ at $T_{N}$. The payment of $R_{N}$ at $T_{N}$ can be seen as the value of the original investment compounded at the fixed rate $k$ over the time period $\left[T_{0}, T_{N}\right]$. The contract is essentially a zero-coupon bond with a discretely compounding rate of $k$. Suppose now that the dealer is granted a Bermudanstyle option to cancel the zero-coupon note at any time $T_{n}, n=1, \ldots, N-1$, in return for paying the investor the accumulated amount to that date, i.e. $R_{n}$. For reasons that should be obvious, the embedded Bermudan option is called a zero-coupon Bermudan swaption. The payoff to the dealer upon exercise at time $T_{n}$ of this option is evidently equal to i) an immediate outflow of $R_{n}$; and ii) release from the obligation to pay $R_{N}$ at time $T_{N}$, the value of which is $P\left(T_{n}, T_{N}\right) R_{N}$. In other words, the exercise value of a zero-coupon Bermudan swaption, $U_{\mathrm{ZC}, n}\left(T_{n}\right)$, is equal to

$$
U_{\mathrm{ZC}, n}\left(T_{n}\right)=-R_{n}+R_{N} P\left(T_{n}, T_{N}\right)
$$

[^52]which allows us to deduce from (19.12) that

$$
\widetilde{U}_{n}\left(T_{n}\right)=U_{\mathrm{ZC}, n}\left(T_{n}\right), \quad n=1, \ldots, N-1 .
$$

Therefore, the value of a zero-coupon Bermudan swaption with rate $k$ is equal to the value of an accreting Bermudan swaption with the notional accreting at rate $k$, as in (19.11). Our earlier discussion of accreting Bermudan swaptions therefore holds unchanged for zero-coupon Bermudan swaptions.

### 19.4.7 American Swaptions

Having a time-varying notional schedule is not the only non-standard feature that can be attached to Bermudan swaptions. One relatively common deviation from the standard contract permits the option holder to exercise on any business date, after a lockout period. Not surprisingly, such swaptions are called American swaptions. These are fairly popular in the US as hedges for mortgage bonds, presumably because the American exercise feature might be considered a better hedge for the prepayment behavior of mortgage borrowers.

The coupon-paying nature of the underlying swap makes the definition of an American swaption somewhat complicated. If exercise takes place during a period $\left[T_{n}, T_{n+1}\right]$, then the option holder receives the swap starting at $T_{n+1}$, i.e. $U_{n+1}(\cdot)$, as well as an "exercise fee" equal to the difference between the Libor rate effective for the period $\left[T_{n}, T_{n+1}\right]$ (i.e. $L_{n}\left(T_{n}\right)$ ) and the fixed rate, times the notional and times the remaining time to $T_{n+1}$ (in the appropriate day count convention). In mathematical notations, the exercise value per unit notional at $t, t \in\left(T_{n}, T_{n+1}\right]$, is given by

$$
U_{n}^{A}(t)=\left(L_{n}\left(T_{n}\right)-k\right)\left(T_{n+1}-t\right)+U_{n+1}(t) .
$$

We emphasize that, as a rule, the time $t$ fee is set to the "accrued current coupon" $\left(L_{n}\left(T_{n}\right)-k\right)\left(T_{n+1}-t\right)$, not its discounted value from $t$ to the payment date $T_{n+1}$. This choice, as is true of many others related to contract specification, is made by those who write term sheets (documents outlining details of derivatives contracts) rather than by those responsible for valuation algorithms, and implies that the exercise value will be discontinuous in time,

$$
\begin{equation*}
U_{n}^{A}\left(T_{n}+\right) \neq U_{n}^{A}\left(T_{n}\right)=U_{n}\left(T_{n}\right) . \tag{19.13}
\end{equation*}
$$
^eq19-13

Odd as it is, this discontinuity is not the main issue with American swaptions. From a valuation standpoint the biggest problem with American swaptions is the fact that the exercise value at $t$ (if not equal to one of the coupon dates) depends on the value of the Libor rate at $T_{n}<t$ and is thus path-dependent.

For Monte Carlo based valuation methods, the path dependence is not a problem, as the value of the Libor rate is known on each path when estimating the exercise boundary or using it in valuation. In a PDE setting, matters
are more complicated. Before describing possible methods, let us comment on the somewhat prevalent view that one can approximate an American swaption with a Bermudan swaption with high frequency of exercise dates.

### 19.4.7.1 American Swaptions vs. High-Frequency Bermudan Swaptions

Let us choose a particular period [ $T_{n}, T_{n+1}$ ] and consider it subdivided into $M$ periods

$$
T_{n}=t_{0}<t_{1}<\ldots<t_{\Lambda I}=T_{n+1}
$$

Then, consider a Bermudan swaption with exercises at $\left\{t_{i}\right\}_{i=1}^{M I}$, versus an American swaption that can be exercised on the same dates. Also, for simplicity assume that the "exercise fee" for an American is in fact, properly discounted, such that the discontinuity in (19.13) is removed. Concentrating only on the exercise value contributions $u_{n}^{A}(t), u_{n}^{B}(t)$ paid in the period $\left[T_{n}, T_{n+1}\right]$, so that

$$
U_{n}^{A}\left(t_{m}\right)=u_{n}^{A}\left(t_{m}\right)+U_{n+1}\left(t_{m}\right), \quad U_{n}^{B}\left(t_{m}\right)=u_{n}^{B}\left(t_{m}\right)+U_{n+1}\left(t_{m}\right)
$$

a standard Bermudan swaption exercised at $t_{m}$ gives the holder an exercise value

$$
u_{n}^{B}\left(t_{m}\right)=\sum_{i=m}^{M-1}\left(t_{i+1}-t_{i}\right) P\left(t_{m}, t_{i+1}\right)\left(\frac{P\left(t_{m}, t_{i}\right)-P\left(t_{m}, t_{i+1}\right)}{\left(t_{i+1}-t_{i}\right) P\left(t_{m}, t_{i+1}\right)}-k\right)
$$

Here

$$
\frac{P\left(t_{m}, t_{i}\right)-P\left(t_{m}, t_{i+1}\right)}{\left(t_{i+1}-t_{i}\right) P\left(t_{m}, t_{i+1}\right)}
$$

is just the time- $t_{m}$ forward Libor rate for the period $\left[t_{l}, t_{i+1}\right]$. Using the standard rearrangement of terms, we obtain

$$
\begin{align*}
u_{n}^{B}\left(t_{m}\right)=L\left(t_{m}, t_{m}, T_{n+1}\right) & \left(T_{n+1}-t_{m}\right) P\left(t_{m}, T_{n+1}\right) \\
& \quad-k \sum_{i=m}^{M-1}\left(t_{i+1}-t_{i}\right) P\left(t_{m}, t_{i+1}\right) \tag{19.14}
\end{align*}
$$
^eq19-14

with

$$
L\left(t_{m}, t_{m}, T_{n+1}\right)=\frac{1-P\left(t_{m}, T_{n+1}\right)}{\left(T_{n+1}-t_{m}\right) P\left(t_{m}, T_{n+1}\right)}
$$

where we have used the full notation (4.2) for forward Libor rates. The first term represents payment of a Libor rate covering the period $\left[t_{n}, T_{n+1}\right]$, and the second is a collection of fixed-rate payments on a schedule $\left\{t_{m}, \ldots, t_{M}=\right. \left.T_{n+1}\right\}$.

Let us contrast this with the exercise value of a true American swaption exercised at $t_{m}$,

$$
\begin{align*}
u_{n}^{A}\left(t_{m}\right)=L\left(T_{n}, T_{n}, T_{n+1}\right)\left(T_{n+1}\right. & \left.-t_{m}\right) P\left(t_{m}, T_{n+1}\right) \\
& -k\left(T_{n+1}-t_{m}\right) P\left(t_{m}, T_{n+1}\right) \tag{19.15}
\end{align*}
$$
^eq19-15

with

$$
L\left(T_{n}, T_{n}, T_{n+1}\right)=\frac{1-P\left(T_{n}, T_{n+1}\right)}{\left(T_{n+1}-T_{n}\right) P\left(T_{n}, T_{n+1}\right)}
$$

Clearly, there are two differences between $u_{n}^{B}$ and $u_{n}^{A}$. The first one is the difference in the fixed leg, the payment of the annuity $k \sum_{i=m}^{M-1}\left(t_{i+1}-\right. \left.t_{i}\right) P\left(t_{m}, t_{i+1}\right)$ versus a single bullet payment $k\left(T_{n+1}-t_{i}\right) P\left(t_{m}, T_{n+1}\right)$. The difference is in the timing of discounting and is normally small, as $T_{n+1}-T_{n}$ would often be equal to 3 months in the US. Even if one deems this to be an issue, any perceptible difference could be eliminated almost fully by imposing appropriately chosen deterministic exercise fees. The second difference, on the other hand, is much greater, and concerns the difference of which Libor rate is applied to the period $\left[t_{m}, T_{n+1}\right]$. For the Bermudan swaption, it is $L\left(t_{m}, t_{m}, T_{n+1}\right)$ and for the American, it is $L\left(T_{n}, T_{n}, T_{n+1}\right)$. Not only will the two Libor rates have different forward values, their different fixing dates ( $t_{m}$ vs. $T_{n}$ ) will affect the amount of volatility each Libor rate experiences over its lifetime. These effects can yield quite significant valuation differences, especially for steeper yield curves and shorter maturities.

### 19.4.7.2 The Proxy Libor Rate Method

Besides highlighting the fallacy of using a high-frequency Bermudan swaption as a proxy for an American, the analysis above also hints at proper remediation. Indeed, it should be clear that the approximation of an American swaption with a Bermudan swaption suffers the most not from the mismatched exercise frequency, but from the difference in the exercise values in the two contracts, as the Bermudan approximation uses a Libor rate that has the wrong forward value and the wrong volatility. The idea behind the proxy Libor rate method involves correcting the forward value/volatility as appropriate, while removing the path-dependence of exercise that hinders the application of backward-induction methods.

Continuing with the notations of the previous section, we define the proxy Libor rate $\widetilde{L}\left(t_{m}, t_{m}, T_{n+1}\right)$ by

$$
\begin{align*}
& \widetilde{L}\left(t_{m}, t_{m}, T_{n+1}\right)=L\left(0, T_{n}, T_{n+1}\right) \\
& +\frac{\operatorname{Stdev}\left(L\left(T_{n}, T_{n}, T_{n+1}\right)\right)}{\operatorname{Stdev}\left(L\left(t_{m}, t_{m}, T_{n+1}\right)\right)}\left(L\left(t_{m}, t_{m}, T_{n+1}\right)-L\left(0, t_{m}, T_{n+1}\right)\right) \tag{19.16}
\end{align*}
$$
^eq19-16

Here $\operatorname{Stdev}(X)$ is defined as the Normal (or basis-point, see Remark 7.2.9 term volatility of the rate $X$, which we may compute from any particular model used for American swaption valuation. The proxy Libor rate enjoys the following properties.

- Its expected value under the $T_{n+1}$-forward measure is equal to the forward $L\left(0, T_{n}, T_{n+1}\right)$, i.e. the forward of the rate used in the "real" American swaption.
- Its (term) volatility is equal to that of $L\left(T_{n}, T_{n}, T_{n+1}\right)$.
- It is a function of the yield curve at time $t_{m}$ and, unlike for $L\left(T_{n}, T_{n}, T_{n+1}\right)$, its value is available in a backward-induction scheme at $t_{m}$.

Having defined the proxy Libor rate, we define a Bermudan swaption whose exercise value approximates the exercise value of an American swaption (compare to (19.14)),

$$
\begin{aligned}
u_{n}^{B}\left(t_{m}\right) \approx \widetilde{L}\left(t_{m}, t_{m}, T_{n+1}\right)\left(T_{n+1}-t_{m}\right) & P\left(t_{m}, T_{n+1}\right) \\
& -k \sum_{\imath=m}^{M-1}\left(t_{i+1}-t_{i}\right) P\left(t_{m}, t_{i+1}\right)
\end{aligned}
$$

Assuming that the underlying model is Markovian and low-dimensional, a Bermudan swaption with these exercise values can easily be evaluated in a finite difference lattice, and is a close proxy for the "real" American swaption.

### 19.4.7.3 The Libor-as-Extra-State Method

While having an accurate approximation is a good step forward, it may still be of use to have an exact lattice-based valuation method for American swaptions, especially when assessing the accuracy of various approximations. One approach for this is the "extra state variable" method, where pathdependence is dealt with by back-propagating values of the security in all possible states of the path-dependent state variable, and then applying update conditions between different "slices", see Sections 2.7.5 and 18.4.5. For our purposes here, we define the extra state variable by

$$
I(t)=\sum_{n=0}^{N-1} 1_{\left\{t \in\left[T_{n}, T_{n+1}\right)\right\}} L\left(T_{n}, T_{n}, T_{n+1}\right)
$$

Clearly, the payoff of an American swaption at time $t$ could be expressed as a function of the yield curve observed at $t$, and the state variable $I(t)$. To elaborate, let $H_{n}(t, I)$ be the value of an American swaption with exercise dates after $T_{n}$ at time $t$ given $I(t)=I$. Focusing on the period $\left[T_{n}, T_{n+1}\right)$ and assuming the exercise dates fall on the grid $\left\{t_{m}\right\}$ as defined in Section 19.4.7.1, we have the following recursion,

$$
\begin{aligned}
& H_{n}\left(t_{m}, I\right) \\
& \quad=\mathrm{E}_{t_{m}}\left(\frac{B\left(t_{m}\right)}{B\left(t_{m+1}\right)} \max \left(u_{n}^{A}\left(t_{m+1}, I\right)+U_{n+1}\left(t_{m+1}\right), H_{n}\left(t_{m+1}, I\right)\right)\right)
\end{aligned}
$$

for $m=M-1, \ldots, 0$, where we have defined (see (19.15))

$$
u_{n}^{A}\left(t_{m}, I\right)=(I-k)\left(T_{n+1}-t_{m}\right) P\left(t_{m}, T_{n+1}\right)
$$

As $I\left(T_{n}\right)=L\left(T_{n}, T_{n}, T_{n+1}\right)$ is a function of the state variables of the model at $T_{n}$, we have that

$$
H_{n}\left(T_{n}\right)=H_{n}\left(t_{0}, L\left(T_{n}, T_{n}, T_{n+1}\right)\right)
$$

is independent of $I$ and only a function of the state variables of the model at $T_{n}$. We then use $H_{n}\left(T_{n}\right)$ to start a recursion for the $(n-1)$ period. The full American value is given by $H_{0}(0)$.

### 19.4.8 Mid-Coupon Exercise

Half-way between standard Bermudan swaptions and American swaptions lie Bermudan swaptions that allow exercises on the standard tenor dates $\left\{T_{n}\right\}$ plus a select few - often just one - extra dates per coupon period. At this point, we should note that while we so far for convenience have assumed that fixed and floating payments take place on the same schedule, in reality swaps in many currencies (including USD and EUR.) pay floating coupons more frequently than fixed coupons. Taking the US as an example, standard conventions specify that floating rate payments occur every three months (and are linked to three-month Libor rates), while fixed rate payments are made every six months. Exercise dates are often chosen to coincide with floating rate fixing dates, i.e. are spaced three months apart. Having an exercise take place in the middle of a fixed-rate coupon period is not a problem, however, as the value of the remaining part of the fixed-rate coupon is trivial to estimate on the exercise date ${ }^{10}$. Less common are exercise dates in the middle of a floating-rate coupon period. For such contracts we face the same issue as with American swaptions: the exercise value on the exercise date is path-dependent, as it is linked to a fixing of the Libor rate that occurs prior to the exercise date.

PDE pricing of structures with exercise taking place inside a floating rate period involves the same issues as those discussed in Section 19.4.7 for American swaptions, and remediation follows the same path. As an approximation we can use an expression like (19.16) to replace the Libor rate with a proxy Libor rate setting on the exercise date, with the proxy rate constructed to have the same forward value and volatility as the real rate. For exact valuation, the state-variable approach of Section 19.4.7 can be used.

[^53]
### 19.5 Flexi-Swaps

A Bermudan swaption can be interpreted as a fixed-floating swap with zero notional and a (single) option to increase the notional to a given level on any of the exercise dates. Likewise, a cancelable swap can be seen as a swap of full notional with an option to decrease the notional to zero on any of the exercise dates; once the option is exercised, the right goes away. More flexibility in choosing swap notionals is afforded in a so-called flexi-swap (also known as a chooser swap or a band swap), a swap with multiple options to change the notional on a given set of exercise dates, subject to certain constraints. Flexi-swaps are related to the flexi-caps discussed in Section 2.7.6 and are most often used as hedges for so-called balance-guarantee swaps, i.e. swaps with a notional linked to a pool of mortgages. The ability to gradually decrease the notional in a flexi-swap on each exercise date allows the option holder to mimic (random) prepayments in the mortgage pool. We briefly consider flexi-swap valuation in this section.

With a tenor structure $\left\{T_{n}\right\}_{n=0}^{N}$ and a collection of net coupons $X_{n}$ with unit notional, fixing at $T_{n}$ and paying at $T_{n+1}, n=0, \ldots, N-1$, we define a flexi-swap to be a contract that pays a net coupon $X_{n} R_{n}$ at time $T_{n+1}$, where the starting notional $R_{0}$ is fixed up-front, and time- $T_{n}$ notional $R_{n}$ is chosen by the holder of the option at time $T_{n}$ (so that $R_{n}$ is $\mathcal{F}_{T_{n}}$-measurable) for each $n=1, \ldots, N-1$, subject to some constraints. The constraint set for the decision at time $T_{n}$ may include

- Global deterministic bounds, e.g. $R_{n} \in\left[g_{n}^{\text {lo }}, g_{n}^{\text {hi }}\right]$.
- Local bounds that are functions of the current notional, e.g. $R_{n} \in \left[l_{n}^{\mathrm{lo}}\left(R_{n-1}\right), l_{n}^{\mathrm{hi}}\left(R_{n-1}\right)\right]$.
- Bounds that are function of market data $x_{n}$ (such as Libor and swap rates) at time $T_{n}$, e.g. $R_{n} \in\left[m_{n}^{\text {lo }}\left(x_{n}\right), m_{n}^{\text {hi }}\left(x_{n}\right)\right]$.

In a general flexi-swap the constraint set for time $T_{n}$ is the intersection of the global, local and market constraint sets; let us denote this set $\mathcal{C}_{n}\left(R_{n-1}, x_{n}\right)$, so that $R_{n} \in \mathcal{C}_{n}\left(R_{n-1}, x_{n}\right) \subset \mathbb{R}$. It is common to require that the notional may only decrease, so that $\mathcal{C}_{n}\left(R_{n-1}, x_{n}\right) \subset\left[0, R_{n-1}\right]$. Of course, the larger the constraint set is for each date (i.e. the fewer the constraints enforced), the more expensive the flexi-swap will be.

The valuation of a flexi-swap may proceed by backward induction, while keeping track of "current notional". To demonstrate, let $V_{n}(t, R)$ be the time $t$ of the part of the flexi-swap paying strictly after $T_{n}$, given that $R_{n}=R$. At time $T_{n-1}$, the flexi-swap value must be equal to the discounted expected value of the maximum value at time $T_{n}$, with the maximum taken over all possible choices of the notional. This observation allows us to write down the backward recursion equation,

$$
\begin{align*}
& V_{n-1}\left(T_{n-1}, R\right)=P\left(T_{n-1}, T_{n}\right) X_{n-1} R \\
& \quad+B\left(T_{n-1}\right) \mathrm{E}_{T_{n-1}}\left(B\left(T_{n}\right)^{-1} \max _{R^{\prime} \in \mathcal{C}_{n}\left(R, x_{n}\right)}\left\{V_{n}\left(T_{n}, R^{\prime}\right)\right\}\right) \tag{19.17}
\end{align*}
$$
^eq19-17

for $n=N, \ldots, 1$, with the terminal condition $V_{N}\left(T_{N}, R\right) \equiv 0$. The time 0 actual value of the flexi-swap is given by $V_{0}\left(T_{0}, R_{0}\right)$. The recursion (19.17) can be implemented in a PDE model by introducing an extra state variable to keep track of the current notional $R$, along the lines of Sections 2.7.5 and 19.4.7.3.

### 19.5.1 Purely Global Bounds

Using ideas similar to those from Section 19.4.5, Evers and Jamshidian [2005] demonstrate that a flexi-swap with purely global deterministic bounds can be decomposed exactly into a portfolio of Bermudan swaptions. While theoretically interesting, the replication is sometimes awkward in practice as a typical flexi-swap will decompose into hundreds of Bermudan swaptions; valuing them all one by one is rarely more efficient than just applying the recursion (19.17). On the other hand, when using a local model such as a onefactor qG model, valuing the Bermudan swaptions one by one would allow us to tailor calibration to each individual Bermudan swaption, leading to increased precision in the value of the flexi-swap. The choice of the valuation method will be dictated by the trade-off between calibration accuracy and performance.

### 19.5.2 Purely Local Bounds

Flexi-swaps that involve non-global constraints generally cannot be replicated with portfolios of Bermudan swaptions. However, in the practically relevant special case of purely local bounds of scaling type we may obtain a more efficient valuation formula involving no state variables beyond those driving the yield curve. Specifically, let us assume that only local constraints are enforced and that these are given by the current notional multiplied by lower and upper multipliers, i.e.

$$
l_{n}^{\mathrm{lo}}\left(R_{n-1}\right)=\lambda_{n}^{\mathrm{lo}} R_{n-1}, \quad l_{n}^{\mathrm{hi}}\left(R_{n-1}\right)=\lambda_{n}^{\mathrm{hi}} R_{n-1}, \quad 0 \leq \lambda_{n}^{\mathrm{lo}} \leq \lambda_{n}^{\mathrm{hi}}
$$

for $n=1, \ldots, N-1$. To simplify the valuation method (19.17), we make the critical observation that the value of the flexi-swap scales linearly in notional, i.e.

$$
\begin{equation*}
V_{n}\left(T_{n}, R\right)=R V_{n}\left(T_{n}, 1\right) \tag{19.18}
\end{equation*}
$$
^eq19-18

for any $T_{n}, R$. This follows from the fact that all coupons scale linearly with notional $R$, as do all constraints. In particular, as there are no global constraints, our exercise decision at any time $T_{n}$ is independent of the absolute size of the notional.

An important corollary to (19.18) is that on any step $n$, the optimal notional choice is of the "all or nothing" type (known in control theory as "bang-bang"). This follows from the fact that

$$
\begin{aligned}
\max _{R^{\prime} \in\left[\lambda^{\prime l},: R, \lambda^{l \prime \prime} R\right]}\left\{V_{n}\left(T_{n}, R^{\prime}\right)\right\} & =\max _{x \in\left[\lambda_{\prime \prime}^{\prime l}, \lambda_{\prime \prime \prime}^{\prime l \prime}\right]}\left\{V_{n}\left(T_{n}, x R\right)\right\} \\
& =\max _{x \in\left[\lambda_{\prime \prime \prime}^{\prime l}, \lambda_{\prime \prime \prime}^{\prime \prime \prime}\right]}\left\{V_{n}\left(T_{n}, R\right) x\right\}
\end{aligned}
$$

whereby the function being maximized, $V_{n}\left(T_{n}, R\right) x$, is linear in the maximization variable $x$. Hence the maximum is attained at the boundary of the interval,

$$
\begin{equation*}
\max _{R^{\prime} \in\left[\lambda_{\prime \prime \prime}^{\prime \mathrm{lo}} R, \lambda_{\prime \prime \prime}^{\mathrm{ln}} R\right]}\left\{V_{n}\left(T_{n}, R^{\prime}\right)\right\}=\max \left(V_{n}\left(T_{n}, \lambda_{n}^{\mathrm{lo}} R\right), V_{n}\left(T_{n}, \lambda_{n}^{\mathrm{lin}} R\right)\right) . \tag{19.19}
\end{equation*}
$$
^eq19-19

We can use the two observations above to simplify the valuation algorithm. Rewriting (19.17) with the help of (19.19) we get

$$
\begin{aligned}
& V_{n-1}\left(T_{n-1}, R\right)=P\left(T_{n-1}, T_{n}\right) X_{n-1} R \\
& \quad+B\left(T_{n-1}\right) \mathrm{E}_{T_{n-1}}\left(B\left(T_{n}\right)^{-1} \max \left(V_{n}\left(T_{n}, \lambda_{n}^{\mathrm{lo}} R\right), V_{n}\left(T_{n}, \lambda_{n}^{\mathrm{hi}} R\right)\right)\right)
\end{aligned}
$$

Dividing through by $R$, using (19.18), and introducing the abbreviated notation $V_{n}\left(T_{n}, 1\right)=V_{n}\left(T_{n}\right)$, we obtain the valuation equation

$$
\begin{aligned}
V_{n-1}\left(T_{n-1}\right)= & P\left(T_{n-1}, T_{n}\right) X_{n-1} \\
& +B\left(T_{n-1}\right) \mathrm{E}_{T_{n-1}}\left(B\left(T_{n}\right)^{-1} \max \left(V_{n}\left(T_{n}\right) \lambda_{n}^{\mathrm{lo}}, V_{n}\left(T_{n}\right) \lambda_{n}^{\mathrm{hi}}\right)\right) .
\end{aligned}
$$

Clearly, if $V_{n}\left(T_{n}\right)$ is positive, then

$$
\max \left(V_{n}\left(T_{n}\right) \lambda_{n}^{\mathrm{lo}}, V_{n}\left(T_{n}\right) \lambda_{n}^{\mathrm{hi}}\right)=V_{n}\left(T_{n}\right) \max \left(\lambda_{n}^{\mathrm{lo}} \cdot \lambda_{n}^{\mathrm{hi}}\right)=V_{n}\left(T_{n}\right) \lambda_{n}^{\mathrm{hi}},
$$

and if $V_{n}\left(T_{n}\right)$ is negative, then

$$
\max \left(V_{n}\left(T_{n}\right) \lambda_{n}^{\mathrm{lo}}, V_{n}\left(T_{n}\right) \lambda_{n}^{\mathrm{hi}}\right)=V_{n}\left(T_{n}\right) \min \left(\lambda_{n}^{\mathrm{lo}}, \lambda_{n}^{\mathrm{hi}}\right)=V_{n}\left(T_{n}\right) \lambda_{n}^{\mathrm{lo}} .
$$

In total, we therefore have the ultimate valuation recursion

$$
\begin{align*}
& V_{n-1}\left(T_{n-1}\right)=P\left(T_{n-1}, T_{n}\right) X_{n-1} \\
& +B\left(T_{n-1}\right) \mathrm{E}_{T_{\prime \prime-1}}\left(B\left(T_{n}\right)^{-1} V_{n}\left(T_{n}\right)\left(\lambda_{n}^{\mathrm{lo}} 1_{\left\{V_{\prime \prime}\left(T_{\prime \prime}\right)<0\right\}}+\lambda_{n}^{\mathrm{hi}} 1_{\left\{V_{\prime \prime}\left(T_{\prime \prime}\right)>0\right\}}\right)\right) . \tag{19.20}
\end{align*}
$$
^eq19-20

The value at time 0 is given by $R_{0} V_{0}\left(T_{0}\right)$; only one PDE plane is required to calculate it, unlike for (19.17) which requires multiple $R$-planes. We note that the standard cancelable swap valuation recursion is recovered with $\lambda_{n}^{\mathrm{lo}}=0, \lambda_{n}^{\mathrm{hi}}=1$.

### 19.5.3 Marginal Exercise Value Decomposition

In it instructive to see what the marginal exercise value decomposition of Section 18.2.3 looks like for a flexi-swap with local bounds. We rewrite (19.20) in a slightly different way,

$$
\begin{aligned}
V_{n-1}\left(T_{n-1}\right)= & P\left(T_{n-1}, T_{n}\right) X_{n-1} \\
& +B\left(T_{n-1}\right) \mathrm{E}_{T_{n-1}}\left(B\left(T_{n}\right)^{-1}\left[\lambda_{n}^{\mathrm{lo}} V_{n}\left(T_{n}\right)^{-}+\lambda_{n}^{\mathrm{hi}} V_{n}\left(T_{n}\right)^{+}\right]\right) .
\end{aligned}
$$

Then, since

$$
V_{n}\left(T_{n}\right)^{+}=V_{n}\left(T_{n}\right)-V_{n}\left(T_{n}\right)^{-}
$$

we have

$$
\begin{aligned}
& V_{n-1}\left(T_{n-1}\right)-\lambda_{n}^{\mathrm{hi}} V_{n}\left(T_{n-1}\right) \\
& \quad=P\left(T_{n-1}, T_{n}\right) X_{n-1}+B\left(T_{n-1}\right) \mathrm{E}_{T_{n-1}}\left(B\left(T_{n}\right)^{-1}\left(\lambda_{n}^{\mathrm{hi}}-\lambda_{n}^{\mathrm{lo}}\right)\left(-V_{n}\left(T_{n}\right)^{-}\right)\right) \\
& \quad=P\left(T_{n-1}, T_{n}\right) X_{n-1}+B\left(T_{n-1}\right) \mathrm{E}_{T_{n-1}}\left(B\left(T_{n}\right)^{-1}\left(\lambda_{n}^{\mathrm{hi}}-\lambda_{n}^{\mathrm{lo}}\right)\left(-V_{n}\left(T_{n}\right)\right)^{+}\right)
\end{aligned}
$$

and, taking discounted expected values to time 0 ,

$$
\begin{aligned}
& V_{n-1}(0)-\lambda_{n}^{\mathrm{hi}} V_{n}(0)=\mathrm{E}\left(B\left(T_{n}\right)^{-1} X_{n-1}\right) \\
& \quad+\mathrm{E}\left(B\left(T_{n}\right)^{-1}\left(\lambda_{n}^{\mathrm{hi}}-\lambda_{n}^{\mathrm{lo}}\right)\left(-V_{n}\left(T_{n}\right)\right)^{+}\right)
\end{aligned}
$$

This holds for $n=1, \ldots, N$. Weighting the $n$-th equality with

$$
\alpha_{n}^{\mathrm{hi}} \triangleq \prod_{i=1}^{n-1} \lambda_{i}^{\mathrm{hi}}
$$

(with $\alpha_{1}^{\mathrm{hi}}=1$ ), summing all terms, and observing that

$$
\sum_{n=1}^{N} \alpha_{n}^{\mathrm{hi}}\left(V_{n-1}(0)-\lambda_{n}^{\mathrm{hi}} V_{n}(0)\right)=V_{0}(0)-\alpha_{N}^{\mathrm{hi}} \lambda_{N}^{\mathrm{hi}} V_{N}(0)=V_{0}(0)
$$

we obtain

$$
\begin{align*}
V_{0}(0)=\sum_{n=1}^{N} \alpha_{n}^{\mathrm{hi}} \mathrm{E} & \left(B\left(T_{n}\right)^{-1} X_{n-1}\right) \\
& +\sum_{n=1}^{N} \alpha_{n}^{\mathrm{hi}}\left(\lambda_{n}^{\mathrm{hi}}-\lambda_{n}^{\mathrm{lo}}\right) \mathrm{E}\left(B\left(T_{n}\right)^{-1}\left(-V_{n}\left(T_{n}\right)\right)^{+}\right) \tag{19.21}
\end{align*}
$$
^eq19-21

More generally

$$
\begin{align*}
V_{n}\left(T_{n}\right)= & \frac{1}{\alpha_{n+1}^{\mathrm{hi}}} \sum_{i=n+1}^{N} \alpha_{i}^{\mathrm{hi}} \mathrm{E}\left(B\left(T_{i}\right)^{-1} X_{i-1}\right) \\
& +\frac{1}{\alpha_{n+1}^{\mathrm{hi}}} \sum_{i=n+1}^{N} \alpha_{i}^{\mathrm{hi}}\left(\lambda_{i}^{\mathrm{hi}}-\lambda_{i}^{\mathrm{lo}}\right) \mathrm{E}\left(B\left(T_{i}\right)^{-1}\left(-V_{i}\left(T_{i}\right)\right)^{+}\right) . \tag{19.22}
\end{align*}
$$
^eq19-22

### 19.5.4 Narrow Band Limit

The decomposition (19.21) turns out to be useful to study the flexi-swap when the notional range $\left|\lambda_{n}^{\mathrm{hi}}-\lambda_{n}^{\mathrm{lo}}\right|$ is small, which is often the case as clients look for cheaper means to hedge their balance-guarantee swaps (recall that narrow range implies less optionality and lower cost). Let $\epsilon=\left(\lambda_{n}^{\mathrm{hi}}-\lambda_{n}^{\mathrm{lo}}\right) / \lambda_{n}^{\mathrm{hi}}$ be small, and denote by $U_{n}$ the value of all coupons fixing on or after $T_{n}$ weighted by $\alpha_{i}^{\text {hi }}$ so that

$$
U_{n}(t)=B(t) \sum_{i=n}^{N-1} \alpha_{i+1}^{\mathrm{hi}} \mathrm{E}_{t}\left(B\left(T_{i+1}\right)^{-1} X_{i}\right)
$$

This is the value of the portion of the (amortizing) swap after $T_{n}$ assuming that on each exercise date the option holder always chooses the multiplier $\lambda_{i}^{\mathrm{hi}}$. Then, it follows from (19.22) that to first order in $\epsilon$,

$$
\begin{equation*}
\alpha_{n+1}^{\mathrm{hi}} V_{n}\left(T_{n}\right)-U_{n}\left(T_{n}\right)=O(\epsilon) \tag{19.23}
\end{equation*}
$$
^eq19-23

and we obtain from (19.21) that

$$
\begin{aligned}
V_{0}(0)= & U_{0}(0)+\epsilon \sum_{n=1}^{N} \alpha_{n+1}^{\mathrm{hi}} \mathrm{E}\left(B\left(T_{n}\right)^{-1}\left(-V_{n}\left(T_{n}\right)\right)^{+}\right) \\
= & U_{0}(0)+\epsilon \sum_{n=1}^{N} \mathrm{E}\left(B\left(T_{n}\right)^{-1}\left(-U_{n}\left(T_{n}\right)\right)^{+}\right) \\
& +\epsilon \sum_{n=1}^{N} \mathrm{E}\left(B\left(T_{n}\right)^{-1}\left[\alpha_{n+1}^{\mathrm{hi}}\left(-V_{n}\left(T_{n}\right)\right)^{+}-\left(-U_{n}\left(T_{n}\right)\right)^{+}\right]\right) .
\end{aligned}
$$

The last line is of the second order in $\epsilon$ per (19.23), so that

$$
\begin{equation*}
V_{0}(0)=U_{0}(0)+\epsilon \sum_{n=1}^{N} \mathrm{E}\left(B\left(T_{n}\right)^{-1}\left(-U_{n}\left(T_{n}\right)\right)^{+}\right)+O\left(\epsilon^{2}\right) . \tag{19.24}
\end{equation*}
$$
^eq19-24

We recognize the terms in the sum above as European swaptions on the amortizing swap $-U_{n}\left(T_{n}\right)$, so the value of a narrow-band flexi-swap with local constraints is approximately equal to the underlying amortizing swap plus a strip of European swaptions on the remaining parts of the reverse of the underlying amortizing swap.

### 19.6 Monte Carlo Valuation

With our discussion so far, we have demonstrated that low-dimensional models, if appropriately calibrated, can be used effectively for Bermudan swaption valuation with PDE or tree-based methods. Still, it is sometimes useful to be able to price Bermudan swaptions with Monte Carlo methods, e.g. to compare prices computed in a low-dimensional model to those of a larger globally-calibrated model (such as the LM model). The mechanics of Monte Carlo valuation follow our discussion in Section 18.3 closely, so here we merely point out the simplifications made possible by the simpler structure of Bermudan swaptions, as compared to the general class of callable Libor exotics (CLEs).

### 19.6.1 Regression Methods

Bermudan swaptions can be valued by Monte Carlo simulation in straightforward fashion, using the general regression-based methods of Section 18.3. There are, however, a number of shortcuts that are worth pointing out. First, observe that the exercise values of a Bermudan swaption can be calculated directly off the yield curve at the time of exercise, whereby there is no need to use regression methods to estimate exercise values. It follows that we can use the simple algorithm of Section 18.3.1 directly.

Another advantage that Bermudan swaptions enjoy over more complex CLEs is the relative ease with which good explanatory variables can be selected. It is clear that for regressing the hold value of a Bermudan swaption on a given exercise date, the value of the underlying swap is important, suggesting that the overall level of interest rates on each exercise date - as represented by either the swap rate or the value of the swap starting on the exercise date and maturing on the final date of the Bermudan swaption should always be included in the set of explanatory variables. The slope of the yield curve on each of the exercise dates turns out to be relevant as well, as it is actually the forward-starting swap - that is, the swap that underlies the European swaption expiring on the next exercise date - that impacts the hold value, and the difference in value between a spot-starting and a forward-starting swaps clearly originates with the slope of the yield curve. We can either include the forward starting swap as an additional explanatory variable on each exercise date or, better yet, include the spot Libor rate for the next period on each exercise date. The latter suggestion achieves nearly the same result as the spot-starting swap could be decomposed into an FRA for the next period and a forward-starting swap. Note that empirical evidence shows that it is not advisable to use the forward-starting swap alone as the sole explanatory variable per exercise date as it appears both the level and the slope of the yield curve should be represented in the set, especially in the setting of a multi-factor model. We investigate and compare various concrete exercise strategies in Section 19.6.2 below.

Another observation that can be fruitfully explored in the LS regression algorithm is the fact that prices of European options on the underlying swaps, i.e. European swaptions, could be calculated (or approximated in a computationally efficient manner) in most models of interest. As mentioned in Section 18.3.9.3, usage of (proxies of) European swaptions allows us to better incorporate convexity in the hold values in the regressions, improving the final value estimate. Additionally, we can draw on all tricks and enhancements from Section 18.3.10, including, in particular, policy improvement based on the carry argument of Section 18.3.10.2 (we extend the carry results for Bermudan swaptions in Section 19.7.2).

### 19.6.2 Parametric Boundary Methods

One hallmark of the LS regression approach is its "semi-automatic" nature: once we have identified some potentially meaningful variables and assumed a particular form of the regression basis functions, we let the regression algorithm work its magic to sort out a reasonable exercise strategy. In contrast, to be effective, the boundary optimization technique introduced in Section 3.5.2 requires more careful thought about the functional form of the exercise boundary. As there is considerable intuition to be gained from the results of the boundary optimization technique, let us spend some time on the application of this method to Bermudan swaption pricing. The material in this section draws on results in Andersen [2000a], where the reader can look up many additional details and numerical results that we do not list in our brief treatment here.

### 19.6.2.1 Sample Exercise Strategies for Bermudan Swaptions

Perhaps the simplest exercise strategy for a Bermudan swaption (and for many other options with early exercise rights) is to "exercise when the option is sufficiently deep in the money". Mathematically speaking, if $\iota(\cdot)$ is the exercise indicator function ${ }^{11}$ then our first proposed strategy is

$$
\begin{equation*}
\text { Strategy I: } \iota\left(T_{n}\right)=1_{\left\{V_{n, N-n}\left(T_{n}\right) \geq h_{I}\left(T_{n}\right)\right\}} \text {, } \tag{19.25}
\end{equation*}
$$
^eq19-25

where $V_{n, N-n}\left(T_{n}\right)$ is the underlying swap value (see Section 19.4.1) and $h_{I}(\cdot)$ is some unknown deterministic function. Assuming that the Bermudan swaption is of the payer type, the swap value $V_{n, N-n}\left(T_{n}\right)$ in (19.25) can be computed directly from the yield curve as

$$
V_{n, N-n}\left(T_{n}\right)=A_{n, N-n}\left(T_{n}\right)\left(S_{n, N-n}\left(T_{n}\right)-k\right) .
$$

It is clear that the function $h_{I}(\cdot)$ must be strictly non-negative, a constraint that should be checked and enforced in the search for $h_{I}$. As discussed in

[^54]Section 3.5.2, the search for the $N-1$ values $h_{I}\left(T_{N-1}\right), h_{I}\left(T_{N-2}\right), \ldots, h_{I}\left(T_{1}\right)$ can be conducted in backwards fashion from a set of Monte Carlo pre-trials, starting from the known condition $h_{I}\left(T_{N-1}\right)=0$ and using the fact that a Bermudan swaption with first exercise date $T_{1}$ will have the same optimal exercise indicator function at time $T_{n}$ as will a Bermudan swaption with first exercise date $T_{n}$ (as long as both are written on swaps with identical coupons and terminal maturity, of course). For each value of $n$, establishing $h_{I}\left(T_{n}\right)$ involves a one-dimensional optimization only, to be done either by outright sorting or by a derivatives-free one-dimensional optimizer. We emphasize that all pre-trials should be cached for numerical efficiency; when using Strategy I, for each path it suffices to store at every date $T_{n}, n=1, \ldots, N-1$, the intrinsic swap value as well as the numeraire (e.g. the spot numeraire $B\left(T_{n}\right)$ if working in the spot measure), for a total of $2(N-1)$ double-precision numbers per path.

As the function $h_{I}\left(T_{n}\right)$ tends to be decreasing roughly linearly as a function of $T_{n}$, it often will suffice to assume that $h_{I}(\cdot)$ is piecewise linear on the interval [ $T_{1}, T_{N-1}$ ], with a low-dimensional number $b$ of break-points $t_{1}<t_{2}<\ldots<t_{b}$, satisfying $t_{1}=T_{1}$ and $t_{b}=T_{N-1}$. The $b-1$ values $h_{I}\left(t_{1}\right), \ldots, h_{I}\left(t_{b-1}\right)$ can be found by a series of one-dimensional optimizations as described earlier, with the values of $h_{I}(\cdot)$ at coupon dates $T_{1}, T_{2}, \ldots, T_{N-1}$ easily computed by linear interpolation. The piecewise linear representation of the exercise rule not only improves numerical efficiency by reducing the number of optimizations to be performed from $N-1$ to $b-1$, but also makes the overall algorithm more robust by assigning more explanatory value to each quantity that is optimized over. Indeed, the fewer parameters that have to be estimated by optimization, the less Monte Carlo pre-trials are necessary to get a smooth, noise free estimation of the exercise boundary. Andersen [2000a] demonstrates that even very low values of $b$ (e.g. 2 or 3) will often suffice, a consequence of the well-known fact that prices of options with early exercise rights tend to be quite insensitive to the precise location of the exercise barrier.

We shall show some numerical results for Strategy I in (19.25) shortly, but let us first introduce some more advanced strategies. Recalling that Bermudan swaptions cannot be worth less than the most expensive core European swaption (see the bound (19.2)), it is reasonable to contemplate the application of a policy improvement step in (19.25) to enforce this constraint. Let us therefore define

$$
V_{\text {swaption, } M(n)}^{\max }\left(T_{n}\right)=\max _{i=n+1, \ldots, M(n)} V_{\text {swaption, } 2}\left(T_{n} ; k\right)
$$

where $M(n)$ is some $n$-dependent upper bound for the number of European swaptions to include in the max-operation (and $k$ is the strike, see (19.2)). Then, a second exercise strategy is

Strategy II: $\iota\left(T_{n}\right)=\left\{\begin{array}{l}1, V_{n . N-n}\left(T_{n}\right)>\max \left(h_{I I}\left(T_{n}\right), V_{\text {swaption }, M(n)}^{\text {max }}\left(T_{n}\right)\right), \\ 0, \text { otherwise } .\end{array}\right.$

Setting $M(n)=N-1$ for all $n$ would ensure that our strategy never breaks the hard value bound (19.2), but could also make computation of the strategy computationally expensive, particularly in models where European swaption pricing requires non-trivial work. As typically only the first few European swaptions are candidates for the maximum in (19.26), to cut down on numerical work ${ }^{12}$ it may make sense to write

$$
\begin{equation*}
M(n)=\min (N-1, n+1+m) \tag{19.27}
\end{equation*}
$$
^eq19-27

where $m$ is some relatively small integer, e.g. 1 or 2 .
A strategy related to (19.26) is

$$
\text { Strategy III: } \quad \iota\left(T_{n}\right)=\left\{\begin{array}{l}
1, V_{n, N-n}\left(T_{n}\right)>h_{I I I}\left(T_{n}\right)+V_{\text {swaption }, M(n)}^{\max }\left(T_{n}\right),  \tag{19.28}\\
0, \text { otherwise } .
\end{array}\right.
$$
^eq19-28

Strategy III replaces the absolute trigger condition of Strategy I with a relative one, where exercise takes place when the intrinsic value is sufficiently high relative to the most expensive core European swaption. To some extent a Bermudan swaption can be viewed as a multi-factor best-of option (that is, an option to choose the most expensive of several assets, see Section 19.2), and Strategy III allows one to impose the well-known condition that exercise never takes place when the underlying assets are too close to each other, irrespective of their magnitudes ${ }^{13}$. Notice that if we enforce that $h_{I I I}$ be strictly non-negative, Strategy III would automatically enforce the policy improvement condition in (19.26). By considering multiple component swaption, Strategies II and III effectively embed more information about the detailed state of the yield curve into the exercise decision than Strategy I. Strategies II and III can thus be expected to be most useful in a multi-factor model ${ }^{14}$. Note that both Strategies II and III can be modified the same way as Strategy I to allow for a piecewise linear representation of the trigger functions $h_{I I}$ and $h_{I I I}$ on some low-dimensional grid $\left\{t_{i}\right\}_{i=1}^{b}$. Also note that the storage requirements for pre-simulations of Strategies II and III will involve $3(N-1)$ number per pre-simulation, as we must store on each exercise date $T_{n}$ i) the numeraire value; ii) the intrinsic swap value; and iii) the maximum core European swaption value.

[^55]Strategies I-III all involve only sequences of one-dimensional optimizations to uncover the scalar functions $h_{I}, h_{I I}, h_{I I I}$; all optimizations start with the boundary condition $h_{I}\left(T_{N-1}\right)=h_{I I}\left(T_{N-1}\right)=h_{I I I}\left(T_{N-1}\right)=0$. Higher-dimensional strategies are possible, too, although they rarely seem worth the extra effort. In Andersen [2000a], the strategies (19.25) and (19.28) are combined into

$$
\text { Strategy IV: } \iota\left(T_{n}\right)=\left\{\begin{array}{l}
1, V_{n, N-n}\left(T_{n}\right)>\max \left(h_{I V}^{1}\left(T_{n}\right), h_{I V}^{2}\left(T_{n}\right)\right. \\
\left.\quad+V_{\text {swaption.M(n) }}^{\max }\left(T_{n}\right)\right) \\
0, \text { otherwise },
\end{array}\right.
$$

where now two functions, $h_{I V}^{1}$ and $h_{I V}^{2}$, have to be determined by optimization. In Andersen [2000a] it is found that this strategy results in no statistically significant pick-up in Bermudan value compared to the simpler strategies above.

### 19.6.2.2 Some Numerical Tests

To test the exercise strategies outlined above, we shall use simple one- and two-factor log-normal LM models. Specifically, we consider a setting where 3 month Libor rates satisfy

$$
d L_{k}(t) / L_{k}(t)=O(d t)+\lambda_{k}(t)^{\top} d W(t)
$$

where $W(t)$ is a vector Brownian motion and where the drift term depends on the probability measure, see Chapter 14 . We consider two settings of $\lambda_{k}(t)$ :

$$
\begin{array}{ll}
\text { Scenario } \mathrm{A}: & \lambda_{k}(t)=20 \% \\
\text { Scenario } \mathrm{B}: & \lambda_{k}(t)=\left(15 \%, 15 \%-\sqrt{0.009\left(T_{k}-t\right)}\right)^{\top}
\end{array}
$$

In Tables 19.5 and 19.6 are numerical results for various Bermudan swaptions, using several of the strategies outlined earlier. We used an initial forward curve that was flat at $10 \%$ (quarterly compounded). In computing the lower bounds in the tables, we used 5,000 pre-trials to establish the trigger functions $h_{I}, h_{I I}$, and $h_{I I I}$ on a time line with $b=4$ break-points; in Strategies II and III, we used $M(n)=N-1$ for all $n .50,000$ independent pricing paths were subsequently drawn to compute the lower bounds for each strategy. The tables also include upper bound duality results, computed from Strategy I using the nested simulation algorithm in Section 18.3.8.2, with $K_{U}=750$ outer paths and $K_{\text {nest }}=300$ inner paths. The $95 \%$ confidence interval (CI) listed in the tables were computed as outlined in Section 18.3.8.4.

To comment on the tables, we first notice from Table 19.5 that the duality gap computed from Strategy I is never more than $1-2$ basis points, leading

| Type | Strike | Strategy I | Strategy II | Strategy III | $\widehat{\Delta}$ |  | 95\% CI |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| $15 \mathrm{M} / 3 \mathrm{M}$ | 8\% | 184.6 (0.1) | 184.6 (0.1) | 184.6 (0.1) | 0.02 |  | 184.5-184.8 |
| $15 \mathrm{M} / 3 \mathrm{M}$ | 10\% | 49.1 (0.1) | 49.1 (0.1) | 48.9 (0.1) | 0.02 |  | 48.7-49.2 |
| $15 \mathrm{M} / 3 \mathrm{M}$ | 12\% | 8.9 (0.1) | 8.9 (0.1) | 8.7 (0.1) | 0.004 |  | 8.5-8.9 |
| 3Y/1Y | 8\% | 355.6 (0.4) | 355.6 (0.4) | 355.1 (0.4) | 0.07 |  | 354.3-355.9 |
| 3Y/1Y | 10\% | 157.8 (0.5) | 157.8 (0.5) | 156.8 (0.5) | 0.2 |  | 156.0-158.0 |
| 3Y/1Y | 12\% | 61.8 (0.4) | 61.8 (0.4) | 61.0 (0.3) | 0.04 |  | 60.2-61.7 |
| 6Y/1Y | 8\% | 807.2 (0.9) | 807.2 (0.9) | 808.0 (0.9) | 0.23 |  | 805.9-809.8 |
| 6Y/1Y | 10\% | 417.8 (0.9) | 417.8 (0.9) | 416.9 (0.9) | 0.63 |  | 13.7-418.0 |
| 6Y/1Y | 12\% | 212.7 (0.9) | 212.7 (0.9) | 212.6 (0.9) | 0.33 |  | 11.4-215.2 |
| 11Y/1Y | 8\% | 1381.6 (1.6) | 1381.6 (1.6) | 1380.2 (1.6) | 1.33 |  | 1378.5-1386.3 |
| 11Y/1Y | 10\% | 812.9 (1.4) | 812.9 (1.4) | 813.2 (1.4) | 1.26 |  | 810.1-817.1 |
| 11Y/1Y | 12\% | 495.8 (1.5) | 495.8 (1.5) | 496.7 (1.4) | 0.71 |  | 495.3-502.1 |
| 6Y/3Y | 8\% | 493.2 (0.8) | 493.7 (0.8) | 493.3 (0.8) | 0.08 |  | 492.3-495.7 |
| 6Y/3Y | 10\% | 293.6 (0.9) | 294.6 (0.9) | 293.0 (0.9) | 0.65 |  | 292.4-296.7 |
| 6Y/3Y | 12\% | 170.3 (0.8) | 170.3 (0.8) | 169.9 (0.8) | 0.53 |  | 168.9-172.8 |

Table 19.5. Upper and lower bound results for the one-factor model in Scenario A. The initial forward curve is flat at $10 \%$, quarterly compounded. All values are computed using Euler-style discretization and are reported in upfront basis points; numbers in parentheses are sample Monte Carlo errors. "Type" refers to the maturity/lockout period of the Bermudan swaption. " $\widehat{\Delta}$ " is the upper-lower duality gap estimate and " $95 \% \mathrm{CI}$ " is the $95 \%$ confidence interval for the Bermudan swaption price. The computational setup is described in more detail in the text.

| Type | Strike | Strategy I | Strategy II | Strategy III | $\widehat{\Delta}$ | 95\% CI |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| $15 \mathrm{M} / 3 \mathrm{M}$ | 8\% | 184.0 (0.0) | 184.0 (0.0) | 184.0 (0) | 0.05 | 183.9-184.1 |
| $15 \mathrm{M} / 3 \mathrm{M}$ | 10\% | 43.3 (0.1) | 43.4 (0.1) | 43.2 (0.1) | 0.06 | 43.1-43.6 |
| $15 \mathrm{M} / 3 \mathrm{M}$ | 12\% | 5.6 (0.1) | 5.6 (0.1) | 5.6 (0.1) | 0.01 | 5.5-5.7 |
| 3Y/1Y | 8\% | 339.7 (0.2) | 339.8 (0.2) | 339.4 (0.2) | 0.4 | 339.2-340.6 |
| 3Y/1Y | 10\% | 125.8 (0.3) | 125.9 (0.3) | 125.7 (0.3) | 0.7 | 125.1-127.2 |
| 3Y/1Y | 12\% | 36.9 (0.2) | 36.8 (0.2) | 36.6 (0.2) | 0.2 | 36.4-37.6 |
| 6Y/1Y | 8\% | 750.2 (0.6) | 749.6 (0.6) | 751.6 (0.6) | 3.7 | 749.0-755.2 |
| 6Y/1Y | 10\% | 317.0 (0.7) | 315.9 (0.7) | 319.4 (0.7) | 5.0 | 315.6-323.5 |
| 6Y/1Y | 12\% | 127.7 (0.6) | 128.0 (0.6) | 129.2 (0.6) | 2.6 | 126.5-131.6 |
| 11Y/1Y | 8\% | 1247.3 (1.2) | 1250.9 (1.2) | 1253.7 (1.3) | 18.1 | 1245.1-1269.0 |
| 11Y/1Y | 10\% | 620.8 (1.1) | 627.1 (1.1) | 633.2 (1.3) | 20.8 | 618.4-645.0 |
| 11Y/1Y | 12\% | 327.1 (1.2) | 331.8 (1.1) | 337.0 (1.2) | 14.8 | 324.7-345.0 |
| 6Y/3Y | 8\% | 444.7 (0.6) | 444.4 (0.6) | 445.2 (0.6) | 0.8 | 443.6-446.6 |
| 6Y/3Y | 10\% | 226.9 (0.7) | 227.2 (0.7) | 227.5 (0.7) | 1.2 | 225.5-229.5 |
| 6Y/3Y | 12\% | 107.1 (0.6) | 107.1 (0.6) | 107.6 (0.6) | 0.8 | 105.9-109.0 |

Table 19.6. Upper and lower bound results for two-factor model in Scenario B. All values are in upfront basis points; numbers in parentheses are sample Monte Carlo errors. Labels are identical to those of Table 19.5.

us to conclude that Strategy 1 very accurately captures the correct exercise decision for the model setup in Table 19.5. Supporting this conclusion is the fact that Strategies II and III lead to no statistically significant increase in the Bermudan swaption value. In the two-factor scenario in Table 19.6, the duality gaps are, not surprisingly, wider than for the one-factor case, although still relatively small for most of the contracts examined. Reasonably significant spreads, in the order of 15 to 20 basis points, can be observed for the 11 year contract with 1 year lockout. Intuitively, for the correlation effects introduced by the two-factor model to matter, the exercise period must be quite long; otherwise, even a two-factor model would imply near-perfect correlation of the different swaps the option holder can exercise into. The suboptimality of exercise based on Strategy I for the 11 no-call 1 Bermudan swaption is also reflected in the fact that the more complicated Strategies II and, especially, III here pick up significant additional value relative to Strategy I. In fact, Strategy III produces prices that lie close to the average of the upper and lower bound, suggesting that this strategy is likely quite close to optimal. Using Strategy III (rather than Strategy I) to form an upper bound confirms this: the duality gap for the 11 year contract with 1 year lockout is reduced to $7.3,6.3$, and 3.5 basis points for coupons of $8 \%$, $10 \%$, and $12 \%$, respectively.

While one should not read too much into the limited set of test data presented above, our results do suggest that for models without stochastic volatility, Strategy I is sufficient for short-dated Bermudan swaptions and for models with high forward rate correlation. For longer-dated structures and for multi-factor models, Strategy III is a safer bet. In a LS regression setting, this reinforces our observations of Section 19.6.1 on the importance of including variables that represent both the level and the slope of the yield curve on exercise dates.

### 19.6.2.3 Additional Comments

A number of papers in the literature elaborate on the analysis in Andersen [2000a]. For instance, in an LM model setting Jensen and Svenstrup [2003] conclude that for Strategy III just setting $m=1$ in (19.27) typically yields Bermudan swaption values that are indistinguishable from those computed using $M(n)=N-1$. Jensen and Svenstrup [2003] also compare the parametric boundary optimization technique against an LS regression where the basis functions include the first two powers of the intrinsic swap value and the spot numeraire, as well as their cross product. For an LM model without stochastic volatility, the parametric boundary technique with Strategy III is found to slightly outperform this particular setup of the LS regression. A similar conclusion is reached in Pedersen [1999], where more details on the LS regression for Bermudan swaptions can also be found.

The analysis in Andersen [2000] (and our discussion in the previous section) concerns itself only with models that contain no stochastic volatility
component. Jensen and Svenstrup [2003] examine an LM model with stochastic volatility and conclude that in this case a rather small, but economically significant, duality gap opens up for Strategy III, especially when the volatility of variance is large and/or the mean reversion speed of volatility is low. Not surprisingly, an LS regression where the variance level itself is included in the set of regressors manages to lower this duality gap. In general, for models with stochastic volatility, explicitly specifying the functional form of the exercise boundary seems difficult, and the best approach is typically to use a regression approach to uncover it, as we described in Section 18.3.9.3.

### 19.7 Other Topics

### 19.7.1 Robust Bermudan Swaption Hedging with European Swaptions

As we explain in more detail in Chapter 22, risk management of exotic derivatives such as Bermudan swaptions generally involves both delta hedging (offsetting sensitivity to the yield curve by dynamically trading in swaps) and vega hedging (offsetting sensitivity of Bermudan swaptions to changes in volatility). Vega hedging of Bermudan swaptions and other CLEs is typically done by trading in European swaptions, but as transaction costs for options can be relatively high, dealers would prefer hedges that do not require frequent rebalancing. A good example of such a hedge would be the static hedging of CMS-linked derivatives with European swaptions at multiple strikes, as specified by the replication method in Section 16.6.1. The resulting hedge not only needs no rebalancing, but is model-independent (up to the annuity mapping function selection).

For Bermudan swaptions we are not aware of any known modelindependent static hedge position in European swaptions, but some insights can be gained from the marginal exercise value decomposition of Section 18.2.3. Even though each of the European options in the decomposition, $\mathrm{E}\left(B\left(T_{n}\right)^{-1}\left(U_{n}\left(T_{n}\right)-H_{n}\left(T_{n}\right)\right)^{+}\right.$), is not a standard European swaption, the replication method of Section 16.6.1 tells us that it can easily be represented as a static position of European swaptions over a continuum of strikes. This position, however, is not model-independent, as each payoff $\left(U_{n}\left(T_{n}\right)-H_{n}\left(T_{n}\right)\right)^{+}$is sensitive (through $H_{n}\left(T_{n}\right)$ ) to the model volatilities and, in particular, the forward volatilities produced by the model. Also, as volatilities inevitably change over time, the effective payoffs are liable to change, as is therefore the composition of the European swaption portfolio that the Bermudan swaption is decomposed into. This implies that the hedge portfolio will need rebalancing over time, which of course rules it out as a truly static hedging portfolio.

Another approach that could be pursued is the semi-static decomposition of barrier options into European swaptions developed by Andersen et al.
[2002]. While developed specifically for barrier options, the technique of this paper also applies to Bermudan swaptions, as these can be interpreted as barrier options with a knock-in barrier set to the optimal exercise boundary (in fact, we already used this representation in developing valuation algorithms). Unfortunately, this line of attack also fails to produce a static hedge, for the same reasons as for the marginal exercise value decomposition: the hedge portfolio depends strongly on the model-specific volatility structure, and is also likely to need rebalancing over time as volatility moves around randomly.

While no theoretically airtight static hedge for Bermudan swaptions is known (as far as we are aware), various pragmatic strategies - often collectively known as the portfolio replication approach - have been more successful. Let us briefly summarize the main idea here, while acknowledging the fact that it can be implemented in many different ways. As an aside, we note that the approach can be applied to many exotic derivatives, although its performance would ultimately depend on the specific risk characteristics of the specific derivative under consideration.

We start by identifying a universe of potential hedging instruments. For a Bermudan swaption of given maturity $T_{N}$, we typically would select all European swaptions with expiry + tenor less than or equal $T_{N}$, with strikes chosen to span a reasonably wide range. Having identified the hedging instruments, we formulate market data scenarios that we would want the hedge portfolio to cover. These would typically be scenarios of joint moves of the yield curve and the volatility surface. For Bermudan swaptions, it is probably sufficient to choose parallel moves in the yield curve, moves by a pre-specified amount within a given range, although one can add more complicated ones, e.g. the yield curve twists and "bends" suggested by principal components analysis (see Section 14.3.1). Similar types of volatility scenarios could be used, such as parallel and non-parallel shifts across all swaption expiries and maturities.

Suppose we have defined $M$ scenarios and chosen $K$ hedging instruments. Let $\Delta V_{\text {Berm }}$ denote the $M$-dimensional vector of value changes of the Bermudan swaption in all $M$ scenarios, and let the vector of value changes of the $k$-th hedging instrument be $\Delta V_{k}, k=1, \ldots, K$. Then, on the last step of the portfolio replication method, we look for a vector of weights $\chi=\left(\chi_{1}, \ldots, \chi_{K}\right)^{\top}$ such that the portfolio of hedges defined by these weights immunizes the changes in values of the Bermudan swaption in all scenarios. This is usually formalized as a least-squares optimization problem,

$$
\left\|\Delta V_{\text {Berm }}-\sum_{k=1}^{K} \chi_{k} \Delta V_{k}\right\|^{2} \rightarrow \min .
$$

Variations are possible, including weighting different scenarios differently or adding additional terms to the objective function to express user preferences,
such as minimizing the total notional of all swaptions, penalizing excessive use of deep out-of-the-money swaptions, and so forth.

We generally find the portfolio replication method to be an effective risk management tool for Bermudan swaptions. Anecdotal evidence provided by traders suggests that the method often outperforms the standard delta/vega hedging approach relying on "local" sensitivities. It also provides a relatively straightforward way of dealing with the well-known gamma-theta mismatch that plagued many dealers' Bermudan swaption portfolios in the 1980s and 1990s. As it turns out, if one uses the volatility sensitivity information from certain simple models to hedge the vega (volatility sensitivity) of a Bermudan swaption by selling European swaptions, the resulting position will sometimes be short gamma (second order yield curve sensitivity) and short time decay (theta). This, however, runs counter to the "standard" Black-Scholes theory in which the gamma and theta balance each other: a long gamma position is always short time decay and makes money in volatile markets and loses money in calm markets; and a short-gamma, long-theta position does the opposite. The unenviable position of being short gamma and short theta will tend to loose money in all markets, volatile and calm, and historically resulted in a number of Bermudan swaption book disasters over the years. The portfolio replication method can help resolve the gamma-theta problem by effectively hedging the exposure to forward volatility or to inter-temporal correlation with European swaptions across multiple expiries and tenors ${ }^{15}$ - something a globally calibrated LM model would do, for example. On the other hand, vega hedging positions computed in short rate models calibrated according to the views held at the time (where mean reversion was often considered superfluous and either excluded from consideration or linked, directly or indirectly, to volatility as in, e.g., the BDT model from Section 11.1.1) would generally suggest incorrect European swaption hedges for the unobservable volatility positions. For example if one does not explicitly link mean reversion to market values of (off-diagonal) swaptions, the forward volatility exposure could either be not hedged at all, or might (wrongly) be linked to the diagonal European swaptions.

### 19.7.2 Carry and Exercise

In Section 18.3.10.2 we showed that if on a given exercise date the next (net) coupon of a cancelable note is positive, then it is not optimal to cancel the note on that date. The net coupon of a derivative security is often referred to as its carry, so we can state that carry-positive cancelable notes should never be exercised.

[^56]Cancelable notes and callable Libor exotics are, of course, intricately linked (see for example Section 18.3.3), so it should come as no surprise that carry-based restrictions on exercise decisions exist for Bermudan swaptions. In fact a result more general than (18.56) holds for cancelable notes and CLEs, but we present it in this chapter, since only for Bermudan swaptions is this more general result actually easy to apply.

We start by recalling that the $n$-th exercise value of a Bermudan payer swaption can be written as

$$
U_{n}(t)=A_{n}(t)\left(S_{n}(t)-k\right) .
$$

A simple relation follows,

$$
\begin{aligned}
U_{n}(t) & =\left[P\left(t, T_{n}\right)-P\left(t, T_{n+1}\right)-k \tau_{n} P\left(t, T_{n+1}\right)\right]+U_{n+1}(t) \\
& =\tau_{n} P\left(t, T_{n+1}\right)\left[\frac{P\left(t, T_{n}\right)-P\left(t, T_{n+1}\right)}{\tau_{n} P\left(t, T_{n+1}\right)}-k\right]+U_{n+1}(t) \\
& =A_{n, 1}(t)\left[S_{n, 1}(t)-k\right]+U_{n+1}(t)
\end{aligned}
$$

and, more generally,

$$
\begin{equation*}
U_{n}(t)=A_{n, m}(t)\left[S_{n, m}(t)-k\right]+U_{n+m}(t), \quad m \geq 1, \tag{19.29}
\end{equation*}
$$
^eq19-29

where we used the notation $A_{n, m}(t)$ and $S_{n, m}(t)$ for the annuity and the swap rate, respectively, for a swap starting at $T_{n}$ and covering $m$ periods. Clearly, for the hold values we have ( $m \geq 1$ ),

$$
H_{n}(t) \geq H_{n+m-1}(t) \geq U_{n+m}(t)
$$

hence it follows from (19.29) that

$$
U_{n}(t) \leq A_{n, m}(t)\left[S_{n, m}(t)-k\right]+H_{n}(t), \quad m \geq 1 .
$$

Taking a minimum over all $m$ we obtain the following result.
Proposition 19.7.1. For a given Bermudan payer swaption and a given exercise date $T_{n}$, we have

$$
\begin{equation*}
U_{n}\left(T_{n}\right)-H_{n}\left(T_{n}\right) \leq \min _{m \geq 1}\left\{A_{n, m}\left(T_{n}\right)\left[S_{n, m}\left(T_{n}\right)-k\right]\right\}, \tag{19.30}
\end{equation*}
$$
^eq19-30

and so if any of the swaps that start at $T_{n}$ and have maturities up to the final maturity of the Bermudan swaption have negative value at $T_{n}$, it is never optimal to exercise at time $T_{n}$.

Proof. If there exists $m$ such that $S_{n, m}\left(T_{n}\right)-k<0$ then by (19.30), the exercise value is strictly less than the hold value, and the exercise is not optimal.

As annuity factors are always positive, it follows from Proposition 19.7.1 above that a Bermudan payer swaption should never be exercised if any swap rate of any still-alive swap is less than the fixed coupon $k$. A similar result holds for Bermudan receiver swaptions - we trust the reader can derive it himself.

### 19.7.3 Fast Pricing via Exercise Premia Representation

There are situations when the speed of valuation of Bermudan swaptions is key yet the accuracy could be sacrificed. One example is robust hedging of Bermudan swaptions in Section 19.7.1 where high-precision pricing is not particularly important. Another is the calculation of a credit value adjustment (CVA), an adjustment to the value of a derivative that takes into account the possibility that the counterparty could default on its payments. While CVA calculations are outside the scope of this book (see Gregory [2009] for a good overview), in essence the evaluation of CVA requires prices of a portfolio of Bermudan swaptions at many future dates under many simulated market conditions. Again, speed of valuation here is very important. One can speed up valuation by using a simple model, such as the one-factor Gaussian model of Section 10, but even higher performance is often desired. In this section we consider a useful approximation based on the representation of a Bermudan swaption as a stream of coupons paid in the exercise region, an adaptation of the representation we developed in Section 1.10.3 for American options.

Recall the definition of $\iota(\cdot)$, the exercise indicator, from Section 19.6.2.1 and Section 18.3.8.2. Let us denote by $V(t)$ the value of the Bermudan swaption at time $t$ as in (18.32).

Proposition 19.7.2. The following holds for any $n=1, \ldots, N-1$,

$$
\begin{align*}
& \mathrm{E}_{T_{n}}\left(\frac{V\left(T_{n}\right)}{B\left(T_{n}\right)}-\frac{V\left(T_{n+1}\right)}{B\left(T_{n+1}\right)}\right)=\mathrm{E}_{T_{n}}\left(\iota\left(T_{n}\right) B\left(T_{n+1}\right)^{-1} X_{n}\right) \\
& \quad+\mathrm{E}_{T_{n}}\left(\iota\left(T_{n}\right)\left(1-\iota\left(T_{n+1}\right)\right) B\left(T_{n+1}\right)^{-1}\left(U_{n+1}\left(T_{n+1}\right)-H_{n+1}\left(T_{n+1}\right)\right)\right) \tag{19.31}
\end{align*}
$$
^eq19-31

where we have used the convention that $\iota\left(T_{N}\right) \equiv 1$ and $V\left(T_{N}\right) \equiv 0$. In particular,

$$
\begin{align*}
& V(0)=\mathrm{E}\left(\sum_{n=1}^{N-1} \iota\left(T_{n}\right) B\left(T_{n+1}\right)^{-1} X_{n}\right) \\
+ & \mathrm{E}\left(\sum_{n=1}^{N-1} \iota\left(T_{n}\right)\left(1-\iota\left(T_{n+1}\right)\right) B\left(T_{n+1}\right)^{-1}\left(U_{n+1}\left(T_{n+1}\right)-H_{n+1}\left(T_{n+1}\right)\right)\right) . \tag{19.32}
\end{align*}
$$
^eq19-32

Proof. We have that

$$
V\left(T_{n}\right)=\iota\left(T_{n}\right) U_{n}\left(T_{n}\right)+\left(1-\iota\left(T_{n}\right)\right) H_{n}\left(T_{n}\right)
$$

and

$$
\begin{aligned}
& \frac{V\left(T_{n}\right)}{B\left(T_{n}\right)}-\frac{V\left(T_{n+1}\right)}{B\left(T_{n+1}\right)} \\
& =\iota\left(T_{n}\right)\left(\frac{U_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}-\frac{V\left(T_{n+1}\right)}{B\left(T_{n+1}\right)}\right)+\left(1-\iota\left(T_{n}\right)\right)\left(\frac{H_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}-\frac{V\left(T_{n+1}\right)}{B\left(T_{n+1}\right)}\right) .
\end{aligned}
$$

Taking expected value conditioned on $\mathcal{F}_{T}$ and using the fact that $\iota\left(T_{n}\right)$ is $\mathcal{F}_{T_{\text {n }}}$-measurable and

$$
H_{n}\left(T_{n}\right)=B\left(T_{n}\right) \mathrm{E}_{T_{n}}\left(\frac{V\left(T_{n+1}\right)}{B\left(T_{n+1}\right)}\right)
$$

by definition, we get

$$
\mathrm{E}_{T_{n}}\left(\frac{V\left(T_{n}\right)}{B\left(T_{n}\right)}-\frac{V\left(T_{n+1}\right)}{B\left(T_{n+1}\right)}\right)=\iota\left(T_{n}\right) \mathrm{E}_{T_{n}}\left(\frac{U_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}-\frac{V\left(T_{n+1}\right)}{B\left(T_{n+1}\right)}\right) .
$$

Moreover,

$$
\begin{aligned}
\iota\left(T_{n}\right) & \left(\frac{U_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}-\frac{V\left(T_{n+1}\right)}{B\left(T_{n+1}\right)}\right) \\
= & \iota\left(T_{n}\right)\left(\frac{U_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}\right. \\
& \left.-\frac{1}{B\left(T_{n+1}\right)}\left(\iota\left(T_{n+1}\right) U_{n+1}\left(T_{n+1}\right)+\left(1-\iota\left(T_{n+1}\right)\right) H_{n+1}\left(T_{n+1}\right)\right)\right) \\
= & \iota\left(T_{n}\right)\left(\frac{U_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}-\frac{U_{n+1}\left(T_{n+1}\right)}{B\left(T_{n+1}\right)}\right) \\
& +\iota\left(T_{n}\right)\left(1-\iota\left(T_{n+1}\right)\right) \frac{1}{B\left(T_{n+1}\right)}\left(U_{n+1}\left(T_{n+1}\right)-H_{n+1}\left(T_{n+1}\right)\right)
\end{aligned}
$$

Then, since

$$
\mathrm{E}_{T_{n}}\left(\frac{U_{n}\left(T_{n}\right)}{B\left(T_{n}\right)}-\frac{U_{n+1}\left(T_{n+1}\right)}{B\left(T_{n+1}\right)}\right)=\mathrm{E}_{T_{"}}\left(\frac{X_{n}}{B\left(T_{n+1}\right)}\right)
$$

(see (18.2)), the result (19.31) follows. Finally, the result (19.32) follows by summing up equalities (19.31) for $n=1, \ldots, N-1$ and taking the (unconditional) expected value.

Let us consider the second term on the right-hand side of (19.31). It represents the contribution of those paths that are in the exercise region at time $T_{n}\left(\iota\left(T_{n}\right)=1\right)$ and are in the hold region at time $T_{n+1}\left(\iota\left(T_{n+1}\right)=0\right)$. One can argue that there are "not too many" of such paths, especially if $T_{n}$ and $T_{n+1}$ are relatively close. Moreover, quantities that are actually evaluated for those paths, the differences between exercise and hold values $U_{n+1}\left(T_{n+1}\right)-H_{n+1}\left(T_{n+1}\right)$, will be small because the exercise value is close to the hold value on the border between exercise and hold regions (by definition of the exercise boundary). Indeed, in the continuous-exercise limit these
terms simply disappear, as should be clear from comparing (19.32) to (1.77). These considerations lead us to suggest an approximation to the value of a Bermudan swaption in which we simply disregard the second sum on the right-hand side of (19.32):

Corollary 19.7.3. The value of a Bermudan swaption (or, indeed, any callable Libor exotic) is approximately equal to the sum of (net) coupons that are paid only in the exercise region, i.e.

$$
\begin{align*}
V(0) & \approx \mathrm{E}\left(\sum_{n=1}^{N-1} \iota\left(T_{n}\right) B\left(T_{n+1}\right)^{-1} X_{n}\right) \\
& =\mathrm{E}\left(\sum_{n=1}^{N-1} 1_{\left\{U_{n}\left(T_{n}\right) \geq H_{n}\left(T_{n}\right)\right\}} B\left(T_{n+1}\right)^{-1} X_{n}\right) . \tag{19.33}
\end{align*}
$$
^eq19-33

The error of approximation is given by the second term in (19.32); the error will decrease as the frequency of exercise of the Bermudan swaption is lowered.

At this point the reader may recall that we have already derived a similar-looking representation for CLEs, namely the marginal exercise value decomposition of Proposition 18.2.1, where we can rewrite (18.8) as

$$
\begin{equation*}
V(0)=\mathrm{E}\left(\sum_{n=1}^{N-1} 1_{\left\{U_{n}\left(T_{n}\right) \geq H_{n}\left(T_{n}\right)\right\}} B\left(T_{n}\right)^{-1}\left(U_{n}\left(T_{n}\right)-H_{n}\left(T_{n}\right)\right)\right) . \tag{19.34}
\end{equation*}
$$
^eq19-34

Not surprisingly, (19.33) could also be derived from (19.34) if we observe that, when the Bermudan is "deep in the money" at time $T_{n}$ then

$$
\begin{aligned}
H_{n}\left(T_{n}\right) & =B\left(T_{n}\right) \mathrm{E}_{T_{n}}\left(B\left(T_{n+1}\right)^{-1} \max \left(U_{n+1}\left(T_{n+1}\right), H_{n+1}\left(T_{n+1}\right)\right)\right) \\
& \approx B\left(T_{n}\right) \mathrm{E}_{T_{n}}\left(B\left(T_{n+1}\right)^{-1} U_{n+1}\left(T_{n+1}\right)\right) \\
& =U_{n+1}\left(T_{n}\right)=U_{n}\left(T_{n}\right)-P\left(T_{n}, T_{n+1}\right) X_{n},
\end{aligned}
$$

and, therefore, we see that the marginal exercise value decomposition implies that

$$
\begin{aligned}
V(0) & \approx \mathrm{E}\left(\sum_{n=1}^{N-1} \iota\left(T_{n}\right) B\left(T_{n}\right)^{-1} P\left(T_{n}, T_{n+1}\right) X_{n}\right) \\
& =\mathrm{E}\left(\sum_{n=1}^{N-1} \iota\left(T_{n}\right) B\left(T_{n+1}\right)^{-1} X_{n}\right)
\end{aligned}
$$

which is (19.33) in Corollary 19.7.3.
Everything we have discussed so far is valid for general CLEs. Let us now specialize our setup to Bermudan swaptions, with the net coupon $X_{n}$
given by $\tau_{n}\left(L_{n}\left(T_{n}\right)-k\right)$. The coupon is a function of the Libor rate $L_{n}\left(T_{n}\right)$; critically, in pretty much all one-factor models (and certainly in the onefactor Gaussian model) the exercise boundary at time $T_{n}$ could also be parameterized by the same Libor rate and expressed in the form

$$
\iota\left(T_{n}\right)=1_{\left\{L_{n}\left(T_{n}\right) \geq h\left(T_{n}\right)\right\}}
$$

for some deterministic function $h(\cdot)$. Then we can rewrite (19.33) as

$$
\begin{equation*}
V(0) \approx \sum_{n=1}^{N-1} P\left(0, T_{n+1}\right) \tau_{n} \mathrm{E}^{T_{n+1}}\left(\left(L_{n}\left(T_{n}\right)-k\right) 1_{\left\{L_{n}\left(T_{n}\right) \geq h\left(T_{n}\right)\right\}}\right) . \tag{19.35}
\end{equation*}
$$
^eq19-35

Each term on the right-hand side can be expressed as a combination of caplets and digital caplets on the Libor rate $L_{n}$. Notice that the expected value is taken under the $T_{n+1}$-forward measure, a measure under which the Libor rate is a martingale. In the one-factor Gaussian model the distribution of $L_{n}\left(T_{n}\right)$ under the $T_{n+1}$-forward measure is well approximated by the Gaussian distribution, and each term could be evaluated rapidly with just a few applications of the Bachelier formula (7.16). Similar approximations could be derived in many other models. For the Libor market model, in particular, the distributions of Libor rates is often known exactly.

The exercise boundary function $h(\cdot)$ in (19.35) is not known a-priori and needs to be found as part of valuation. This can be done efficiently in a backward induction algorithm that utilizes the representation (19.35) for a Bermudan swaption at future times. In particular, we can find the value $h\left(T_{k}\right)$ by solving

$$
\left.H_{k}\left(T_{k}\right)\right|_{L_{k}\left(T_{k}\right)=h\left(T_{k}\right)}=\left.U_{k}\left(T_{k}\right)\right|_{L_{k}\left(T_{k}\right)=h\left(T_{k}\right)},
$$

where $H_{k}\left(T_{k}\right)$ is calculated by an analog to (19.35) with $h\left(T_{k+1}\right), \ldots, h\left(T_{N-1}\right)$ already determined from previous steps. The recursion can be accelerated further by search for the exercise boundary only on the subset of exercise dates, with the missing points filled by interpolation (see Ju [1998]). The final algorithm turns out to be quick and robust, and is well-suited for situations where valuation speed is the primary consideration.

## 19.A Appendix: Forward Volatility and Correlation

When European swaption prices are kept fixed, increasing correlation between forward rates (e.g. by moving from a two- to a one-factor model) will tend to increase forward volatility. There are a number of ways to explain this effect. One approach relies on swaption "triangles" (see Section 20.3), another on the BGM/HJM model formalism. In this appendix we explore the latter.

First, a bit of notation. Let Libor rates $L_{k}(t)=L_{k}\left(t, T_{k}, T_{k+1}\right)$ satisfy

$$
d L_{k}(t)=O(d t)+\sigma_{k}(t) d W_{k}(t), \quad k=1, \ldots, N-1
$$

where $\sigma_{k}$ 's are scalar and the $W_{k}$ 's are correlated Brownian motions. We assume that our calibration is global (see Section 14.5) and therefore the model calibrates properly to caplets (or, equivalently, swaptions on shorttenor swaps). If follows that the quantities

$$
\begin{equation*}
\int_{0}^{T_{k}} \sigma_{k}(u)^{2} d u \tag{19.36}
\end{equation*}
$$
^eq19-36

must be invariants, independent of the correlation between the $W_{k}$ 's.
Consider now the market for short-expiry swaptions - which must also be calibrated in our setup (which is global) - and let us study two different settings of the average correlation between the $W_{k}$ 's, "high" and "low", indicated by appropriate subscripts in what follows. To match the short-expiry, short-maturity swaptions (i.e., caplets), we fundamentally need

$$
\sigma_{1}^{\mathrm{hi}}(0)=\sigma_{1}^{\mathrm{lo}}(0)
$$

On the other hand, to match long-tenor (but still short-expiry) swaptions, we need

$$
\begin{equation*}
\sigma_{k}^{\mathrm{hi}}(0)<\sigma_{k}^{\mathrm{lo}}(0), \quad k=2, \ldots, N-1, \tag{19.37}
\end{equation*}
$$
^eq19-37

a relationship that follows from the fact that swaption volatilities are effectively volatilities of sums of Libor rates. Specifically, as the volatility of a sum increases in correlation, we need to lower the volatility of the "components of the sum" (that is, the Libor rates) to preserve swaption volatilities.

Let us pick some $k>1$. If we look at (19.36) and (19.37), it is obvious that to satisfy both conditions simultaneously $\sigma_{k}^{\mathrm{hi}}(t)$ must ${ }^{16}$ ultimately "overtake" $\sigma_{k}^{\mathrm{lo}}(t)$ as $t$ is increased from 0 to $T_{k}$. As this holds for all $k$, it is clear that forward volatilities of both caps and swaptions will, as promised above, be higher in the high-correlation model than in the low-correlation model.

## 19.B Appendix: A Primer on Moment Matching

## 19.B. 1 Basics

Let there be given $d \log$-normal random variables $X_{1}, \ldots, X_{d}$, with known distribution parameters $m_{i}$ and $s_{i}$ :

$$
\ln \left(X_{i}\right) \sim \mathcal{N}\left(m_{i}, s_{i}^{2}\right), \quad i=1, \ldots, d
$$

[^57]We assume that the $d \times d$ correlation matrix $\rho$ of logarithms is known,

$$
\rho_{i, j} \triangleq \operatorname{Corr}\left(\ln \left(X_{i}\right), \ln \left(X_{j}\right)\right), \quad i, j=1, \ldots, d
$$

From standard results for log-normal variables, the first two moments of the $X_{\imath}$ can be computed as

$$
\begin{align*}
\mathrm{E}\left(X_{\imath}\right) & =\exp \left(m_{i}+\frac{s_{\imath}^{2}}{2}\right)  \tag{19.38}\\
\operatorname{Var}\left(X_{i}\right) & =\mathrm{E}\left(X_{\imath}\right)^{2}\left(\exp \left(s_{i}^{2}\right)-1\right) \tag{19.39}
\end{align*}
$$
^eq19-38
^eq19-39

Also,

$$
\begin{aligned}
\mathrm{E}\left(X_{1} X_{j}\right) & =\mathrm{E}\left(\exp \left(\ln \left(X_{i}\right)+\ln \left(X_{j}\right)\right)\right) \\
& =\mathrm{E}\left(\exp \left(m_{i}+s_{i} Z+m_{j}+s_{j}\left(\rho_{i, j} Z+\sqrt{1-\rho_{i, j}^{2}} Y\right)\right)\right)
\end{aligned}
$$

where $Z$ and $Y$ are independent standard Gaussian variables, and where we have used the Cholesky decomposition. Therefore, using the result (19.38),

$$
\begin{equation*}
\mathrm{E}\left(X_{i} X_{j}\right)=\mathrm{E}\left(X_{i}\right) \mathrm{E}\left(X_{j}\right) \exp \left(\rho_{i, j} s_{i} s_{j}\right) . \tag{19.40}
\end{equation*}
$$
^eq19-40

We note in passing that therefore (see e.g. (17.66))

$$
\begin{align*}
\operatorname{Cov}\left(X_{i}, X_{j}\right) & =\mathrm{E}\left(X_{i} X_{j}\right)-\mathrm{E}\left(X_{i}\right) \mathrm{E}\left(X_{j}\right) \\
& =\mathrm{E}\left(X_{i}\right) \mathrm{E}\left(X_{j}\right)\left(\exp \left(\rho_{i, j} s_{i} s_{j}\right)-1\right) . \tag{19.41}
\end{align*}
$$
^eq19-41

Suppose now that we are interested in approximating the moments of the weighted sum

$$
\begin{equation*}
\widehat{X}=\sum_{i=1}^{d} w_{i} X_{i} \tag{19.42}
\end{equation*}
$$
^eq19-42

where the $w_{2}$ 's are given positive constants. Clearly

$$
\begin{equation*}
\mathrm{E}(\widehat{X})=\sum_{i=1}^{d} w_{i} \mathrm{E}\left(X_{\imath}\right) \tag{19.43}
\end{equation*}
$$
^eq19-43

with $\mathrm{E}\left(X_{1}\right)$ given in (19.38). Also,

$$
\begin{align*}
\mathrm{E}\left(\widehat{X}^{2}\right) & =\sum_{i=1}^{d} \sum_{j=1}^{d} w_{i} w_{j} \mathrm{E}\left(X_{i} X_{j}\right) \\
& =\sum_{i=1}^{d} \sum_{j=1}^{d} w_{i} w_{j} \mathrm{E}\left(X_{i}\right) \mathrm{E}\left(X_{j}\right) \exp \left(\rho_{\imath, j} s_{i} s_{j}\right) \\
& =\sum_{i=1}^{d} w_{i}^{2} \mathrm{E}\left(X_{i}\right)^{2} e^{s_{\imath}^{2}}+2 \sum_{\imath=1}^{d} \sum_{j=i+1}^{d} w_{\imath} w_{j} \mathrm{E}\left(X_{i}\right) \mathrm{E}\left(X_{j}\right) e^{\rho_{i, j} s_{\imath} s_{j}} . \tag{19.44}
\end{align*}
$$
^eq19-44

In many applications, we are interested in representing $\widehat{X}$ as being approximately log-normal, i.e. we would like to write

$$
\ln (\widehat{X}) \sim \mathcal{N}\left(m_{\widehat{X}}, s_{\widehat{X}}^{2}\right)
$$

Using a moment-matching principle, we would determine $m_{\widehat{X}}$ and $s_{\widehat{X}}$ from the equations

$$
\begin{aligned}
\exp \left(m_{\hat{X}}+\frac{s_{\widehat{X}}^{2}}{2}\right) & =\mathrm{E}(\widehat{X}) \\
\exp \left(m_{\widehat{X}}+\frac{s_{\widehat{X}}^{2}}{2}\right)^{2} \exp \left(s_{\widehat{X}}^{2}\right) & =\mathrm{E}\left(\widehat{X}^{2}\right)
\end{aligned}
$$

which can be solved to yield

$$
\begin{equation*}
s_{\widehat{X}}=\sqrt{\ln \left(\mathrm{E}\left(\widehat{X}^{2}\right)\right)-\ln \left(\mathrm{E}(\widehat{X})^{2}\right)}, \quad m_{\widehat{X}}=\ln (\mathrm{E}(\widehat{X}))-\frac{s_{\widehat{X}}^{2}}{2} . \tag{19.45}
\end{equation*}
$$
^eq19-45

In these formulas $\mathrm{E}(\widehat{X})$ and $\mathrm{E}\left(\widehat{X}^{2}\right)$ should be computed from formulas (19.43) and (19.44), respectively.

Note that if we are willing to relax the requirement that $\widehat{X}$ be approximately log-normal, we can obtain more accurate approximations. A popular choice here is to assume that $\widehat{X}$ is approximately displaced log-normal. This introduces one more degree of freedom in the matching distribution (the displacement parameter) which, together with the mean and variance, could be used to match three, rather than two, moments. We leave the details of this for the reader to work out, and in the examples below we stick with simple log-normal moment matching.

## 19.B. 2 Example 1: Asian Option in BSM Model

Let $I(t)$ be some asset following the simple process

$$
\begin{equation*}
d I(t) / I(t)=-b(t) d t+\sigma(t) d W(t) \tag{19.46}
\end{equation*}
$$
^eq19-46

where $W(t)$ is a scalar Brownian motion in the risk-neutral measure. For certain weights $w_{i}$, we form the weighted average

$$
M\left(T_{d}\right)=\sum_{i=1}^{d} w_{i} I\left(T_{i}\right),
$$

on some schedule $0<T_{1}<T_{2}<\ldots<T_{d}$. An Asian option pays

$$
V_{\text {Asian }}\left(T_{p a y}\right)=\left(M\left(T_{d}\right)-K\right)^{+}, \quad T_{p a y} \geq T_{d}
$$

where typically the weights are $w_{i}=1 / d$ for all $i$. Standing at time 0 , we wish to use moment-matching to model the $T_{d}$-observed average $M\left(T_{d}\right)$ as a log-normal variable. From (19.46), it is clear that $I\left(T_{i}\right)$ is a log-normal random variable, since

$$
I\left(T_{i}\right)=I(0) l\left(T_{i}\right) \exp \left(-\frac{1}{2} v\left(T_{i}\right)^{2} T_{i}+\int_{0}^{T_{i}} \sigma(u) d W(u)\right)
$$

where we have defined

$$
l\left(T_{i}\right) \triangleq \exp \left(-\int_{0}^{T_{2}} b(u) d u\right), \quad v\left(T_{i}\right)^{2} \triangleq T_{i}^{-1} \int_{0}^{T_{2}} \sigma(u)^{2} d u
$$

If we define $X_{\imath}=I\left(T_{i}\right)$, it follows that, in the notation of Section 19.B.1,

$$
\begin{aligned}
m_{i} & =\ln \left(l\left(T_{i}\right)\right)-\frac{1}{2} v\left(T_{i}\right)^{2} T_{i} \\
s_{i} & =v\left(T_{i}\right) \sqrt{T_{i}}
\end{aligned}
$$

and that $M\left(T_{d}\right)=\widehat{X}$, where $\widehat{X}$ is defined in (19.42). To use the results (19.45) it only remains to find the correlation matrix $\rho$. But clearly

$$
\begin{aligned}
\rho_{i, j} & =\operatorname{Corr}\left(\int_{0}^{T_{2}} \sigma(u) d W(u), \int_{0}^{T_{1}} \sigma(u) d W(u)\right) \\
& =\frac{\int_{0}^{\min \left(T_{1}, T_{1}\right)} \sigma(u)^{2} d u}{\sqrt{\int_{0}^{T_{1}} \sigma(u)^{2} d u} \sqrt{\int_{0}^{T_{1}} \sigma(u)^{2} d u}} \\
& =\frac{\min \left(v\left(T_{i}\right)^{2} T_{i}, v\left(T_{j}\right)^{2} T_{j}\right)}{v\left(T_{i}\right) v\left(T_{j}\right) \sqrt{T_{2} T_{j}}} \\
& =\frac{\min \left(v\left(T_{2}\right) \sqrt{T_{i}}, v\left(T_{j}\right) \sqrt{T_{j}}\right)}{\max \left(v\left(T_{i}\right) \sqrt{T_{i}}, v\left(T_{j}\right) \sqrt{T_{j}}\right)}
\end{aligned}
$$

Applying (19.43), (19.44), and finally (19.45) then allows us to write, approximately,

$$
\ln \left(M\left(T_{d}\right)\right) \sim \mathcal{N}\left(m_{\widehat{X}}, s_{\widehat{X}}^{2}\right)
$$

for computed constants $m_{\widehat{X}}$ and $s_{\widehat{X}}$. Assuming deterministic interest rates, standard Black-Scholes arguments (see Section 1.9) allow us to finally approximate the time 0 option price as

$$
\begin{align*}
V_{\text {Asian }}(0) & \approx P\left(0, T_{\text {pay }}\right)\left(e^{\left.m_{\widehat{X}}+\frac{1}{2} s_{\widehat{X}}^{2} \Phi\left(d_{+}\right)-K \Phi\left(d_{-}\right)\right)}\right.  \tag{19.47}\\
d_{ \pm} & =\frac{\ln \left(e^{m_{\widehat{X}}+\frac{1}{2} s_{\widehat{X}}^{2}} / K\right) \pm \frac{1}{2} s_{\widehat{X}}^{2}}{s_{\widehat{X}}}=\frac{m_{\widehat{X}}+\frac{1}{2} s_{\widehat{X}}^{2}-\ln (K) \pm \frac{1}{2} s_{\widehat{X}}^{2}}{s_{\widehat{X}}}
\end{align*}
$$
^eq19-47

where $\Phi$ is the Gaussian CDF and $P\left(0, T_{p a y}\right)$ is a risk-free discount factor to time $T_{p a y}$. We note that we may, of course, rewrite this expression in the perhaps slightly more convenient form

$$
\begin{align*}
V_{\text {Asian }}(0) & \approx P\left(0, T_{\text {pay }}\right)\left(\mathrm{E}\left(M\left(T_{d}\right)\right) \Phi\left(d_{+}\right)-K \Phi\left(d_{-}\right)\right)  \tag{19.48}\\
d_{ \pm} & =\frac{\ln \left(\mathrm{E}\left(M\left(T_{d}\right)\right) / K\right) \pm \frac{1}{2} s_{\hat{X}}^{2}}{s_{\hat{X}}}
\end{align*}
$$
^eq19-48

where

$$
\mathrm{E}\left(M\left(T_{d}\right)\right)=\sum_{i=1}^{d} w_{i} \mathrm{E}\left(I\left(T_{i}\right)\right)=\sum_{i=1}^{d} w_{i} I(0) l\left(T_{i}\right)
$$

Note that this form does not require us to compute $m_{\widehat{X}}$, as only $s_{\widehat{X}}$ is needed.

## 19.B. 3 Example 2: Basket Option in BSM Model

Consider $d$ risk-neutral processes

$$
d I_{i}(t) / I_{i}(t)=-b_{i}(t) d t+\sigma_{i}(t) d W_{i}(t), \quad i=1, \ldots, d
$$

where we assume that $\left\langle d W_{i}(t), d W_{j}(t)\right\rangle=\rho_{i, j} d t$. Also consider the payout of a basket option

$$
V_{\text {basket }}\left(T_{p a y}\right)=(\widehat{I}(T)-K)^{+}, \quad T_{p a y} \geq T
$$

where

$$
\widehat{I}(T)=\sum_{i=1}^{d} w_{i} I_{i}(T)
$$

with the understanding that all basket weights $w_{i}$ are positive. In the framework of Section 19.B.1, we now set $X_{i}=I_{i}(T)$, such that, at time 0, $X_{i}$ is log-normal with parameters

$$
m_{i}=\ln \left(l_{i}(T)\right)-\frac{1}{2} v_{i}(T)^{2} T, \quad s_{i}=v_{i} T \sqrt{T}
$$

where we have defined, for $T>0$,

$$
l_{i}(T) \triangleq \exp \left(-\int_{0}^{T} b_{i}(u) d u\right), \quad v_{i}(T)^{2} \triangleq T^{-1} \int_{0}^{T} \sigma_{i}(u)^{2} d u
$$

In the notation of Section 19.B.1, clearly $\widehat{I}(T)=\widehat{X}$ and we may proceed as in Example 1 above to find $m_{\widehat{X}}$ and $s_{\widehat{X}}$, at which point the formula (19.47) (or (19.48)) will price the basket option at time 0 .

