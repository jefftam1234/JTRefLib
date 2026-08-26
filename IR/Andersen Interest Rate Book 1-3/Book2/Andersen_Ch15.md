---
tags: [quant/knowledge, quant/books, interest-rates]
---

## The Libor Market Model II

For the sake of cohesion, our discussion of LM models in Chapter 14 silently skipped over a number of practical issues. Chief amongst these is the problem of how to construct an entire (continuous) discount curve from knowledge of only a finite set of simply compounded Libor rates. This surprisingly subtle issue shall be discussed in this chapter, along with a select set of other advanced pricing and calibration topics in LM modeling. For instance, we provide a number of extensions to the stochastic volatility setup of Chapter 14 and also show how to construct swaption pricing formulas more accurate and more general than those in Chapter 14. We cover the generic problem of evolving separate discount and forward curves, and also include brief discussions of the so-called swap market models and of LM models with "near-Markov" structure. The latter topic shall be taken up again in Chapter 25.

### 15.1 Interpolation

The simulation schemes that we developed so far (see Section 14.6) allow us to obtain at any time $t$ a vector of forward Libor rates on a pre-specified tenor structure. As should be obvious, and as pointed out previously in Sections 14.1.2 and 14.2.3, this information is not sufficient to recover the full interest rate yield curve at time $t$. At the very least, to be able to compute $P\left(t, T_{n}\right)$ for all $n$ (see (14.3)), we need to additionally establish the front stub discount factor $P\left(t, T_{q(t)}\right)$. In addition, as many actual security payoffs dictate that we calculate $P(t, T)$ for an arbitrary $T$, the back stub (forward) discount factor $P\left(t, T, T_{q(T)}\right)=P\left(t, T_{q(T)}\right) / P(t, T)$ will also be required, since

$$
P(t, T)=P\left(t, T_{q(t)}\right) \times\left(\prod_{i=q(t)}^{q(T)-1}\left(1+\tau_{i} L_{i}(t)\right)^{-1}\right) / P\left(t, T, T_{q(T)}\right)
$$

Both the front and back stubs cannot, in general, be expressed as a function of Libor rates on a fixed tenor structure.

There are a number of approaches that could be employed to obtain the front and back stub in a simulation. We start with the back stub as it is somewhat easier to handle.

### 15.1.1 Back Stub, Simple Interpolation

Let us fix the discount factor maturity time $T$, and set $m=q(T)$ such that

$$
T_{m-1} \leq T<T_{m}
$$

Observe that as $T \rightarrow T_{m-1}$ the back stub $P\left(t, T, T_{m}\right)=P\left(t, T_{m}\right) / P(t, T)$ converges to $P\left(t, T_{m-1}, T_{m}\right)$, a discount factor that can be calculated from Libor rates as $P\left(t, T_{m-1}, T_{m}\right)=\left(1+L_{m-1}(t) \tau_{m-1}\right)^{-1}$. At the other limit. when $T \rightarrow T_{m}$, the back stub converges to 1 . Hence, it seems reasonable to approximate $P\left(t, T, T_{m}\right)$ by interpolating between these two known extremes. This idea gives rise to a number of plausible schemes that we now proceed to describe.

A particularly simple idea is to apply linear interpolation directly to bond prices, resulting in the scheme

$$
\begin{equation*}
P\left(t, T, T_{m}\right)=\frac{T-T_{m-1}}{T_{m}-T_{m-1}}+\frac{T_{m}-T}{T_{m}-T_{m-1}} P\left(t, T_{m-1}, T_{m}\right) \tag{15.1}
\end{equation*}
$$
^eq15-1

Using $P\left(t, T_{m-1}, T\right)=P\left(t, T_{m-1}, T_{m}\right) / P\left(t, T, T_{m}\right)$ as the interpolation variable instead, another linear interpolation scheme arises:

$$
P\left(t, T_{m-1}, T\right)=\frac{T_{m}-T}{T_{m}-T_{m-1}}+\frac{T-T_{m-1}}{T_{m}-T_{m-1}} P\left(t, T_{m-1}, T_{m}\right)
$$

or

$$
\begin{equation*}
P\left(t, T, T_{m}\right)=\frac{P\left(t, T_{m-1}, T_{m}\right)}{\frac{T_{m-}-T}{T_{m-}-T_{m-1}}+\frac{T-T_{m-1}}{T_{m-}-T_{m-1}} P\left(t, T_{m-1}, T_{m}\right)} \tag{15.2}
\end{equation*}
$$
^eq15-2

Alternatively, we can apply piecewise constant interpolation to continuously compounded instantaneous forward rates $f(t, u)$ for $u \in\left[T_{m-1}, T_{m}\right]$, yielding

$$
-\frac{1}{T_{m}-T} \ln P\left(t, T, T_{m}\right)=-\frac{1}{T_{m}-T_{m-1}} \ln P\left(t, T_{m-1}, T_{m}\right)
$$

or, in terms of forward bond prices,

$$
\begin{equation*}
P\left(t, T, T_{m}\right)=P\left(t, T_{m-1}, T_{m}\right)^{\frac{T_{m}-T}{T_{m-}-T_{m-1}}} \tag{15.3}
\end{equation*}
$$
^eq15-3

Yet another interpolation scheme is obtained by constant interpolation of simply compounded rates,

$$
\frac{1}{T_{m}-T}\left(\frac{1}{P\left(t, T, T_{m}\right)}-1\right)=\frac{1}{T_{m}-T_{m-1}}\left(\frac{1}{P\left(t, T_{m-1}, T_{m}\right)}-1\right)
$$

resulting in the scheme

$$
\begin{equation*}
P\left(t, T, T_{m}\right)=\left(\frac{T-T_{m-1}}{T_{m}-T_{m-1}}+\frac{T_{m}-T}{T_{m}-T_{m-1}} \frac{1}{P\left(t, T_{m-1}, T_{m}\right)}\right)^{-1} \tag{15.4}
\end{equation*}
$$
^eq15-4

While the interpolation schemes (15.1), (15.2), (15.3), and (15.4) are all simple to understand and to apply, they are ultimately flawed as they violate the basic no-arbitrage conditions

$$
\begin{align*}
P\left(0, T_{m-1}, T\right) & =\mathrm{E}^{T_{m-1}}\left(P\left(t, T_{m-1}, T\right)\right)  \tag{15.5}\\
P\left(0, T_{m-1}, T_{m}\right) & =\mathrm{E}^{T_{m-1}}\left(P\left(t, T_{m-1}, T_{m}\right)\right) \tag{15.6}
\end{align*}
$$
^eq15-5
^eq15-6

Any interpolation scheme, once we apply the expected value operator, imposes a certain relationship on discount bonds at time 0 , a relationship that, in general, will not be satisfied by actual market prices. Taking as an example (15.2), applying the expectation operator $\mathrm{E}^{T_{m-1}}$ and using (15.5), we obtain

$$
P\left(0, T_{m-1}, T\right)=\frac{T_{m}-T}{T_{m}-T_{m-1}}+\frac{T-T_{m-1}}{T_{m}-T_{m-1}} P\left(0, T_{m-1}, T_{m}\right)
$$

a relationship between time 0 discount bond prices that is unlikely to be satisfied a priori.

### 15.1.2 Back Stub, Arbitrage-Free Interpolation

To ensure that observable relationships between time 0 discount bond prices are respected, consider choosing an arbitrary constant $\alpha(T)$ and setting

$$
\begin{equation*}
P\left(t, T_{m-1}, T\right)=P\left(0, T_{m-1}, T\right)+\alpha(T)\left(P\left(t, T_{m-1}, T_{m}\right)-P\left(0, T_{m-1}, T_{m}\right)\right) . \tag{15.7}
\end{equation*}
$$
^eq15-7

Clearly, the additive scheme (15.7) will preserve (15.5) as long as (15.6) is satisfied - and (15.6) is essentially a no-arbitrage condition for discount bonds maturing on the tenor structure and is guaranteed by the LM model construction itself.

We can regard $\alpha(\cdot)$ as a function of maturity time ${ }^{1}$; for consistency $\alpha\left(T_{m}\right)$ must equal 1 and $\alpha\left(T_{m-1}\right)=0$, but beyond this there are few restrictions on $\alpha(T)$. Yet it is not advisable to specify $\alpha(T)$ arbitrarily, as this may affect model dynamics in unintended ways. To devise a reasonable approach to the definition of $\alpha(T)$, we note from (15.7) that

$$
d P\left(t, T_{m-1}, T\right)=O(d t)+\alpha(T) d P\left(t, T_{m-1}, T_{m}\right)
$$

[^69]On the other hand, in an HJM model we have

$$
\begin{aligned}
& d P\left(t, T_{m-1}, T\right) / P\left(t, T_{m-1}, T\right) \\
& \quad=O(d t)+\left(\sigma_{P}\left(t, T_{m-1}\right)-\sigma_{P}(t, T)\right)^{\top} d W(t)
\end{aligned}
$$

and

$$
\begin{aligned}
& d P\left(t, T_{m-1}, T_{m}\right) / P\left(t, T_{m-1}, T_{m}\right) \\
& \quad=O(d t)+\left(\sigma_{P}\left(t, T_{m-1}\right)-\sigma_{P}\left(t, T_{m}\right)\right)^{\top} d W(t)
\end{aligned}
$$

from which we conclude that $\alpha(T)$ is linked to the ratio of forward discount bond volatilities. Exploiting this link, we may define $\alpha(T)$ from, for instance. the equation

$$
\begin{aligned}
& P\left(t, T_{m-1}, T\right)\left\|\sigma_{P}\left(t, T_{m-1}\right)-\sigma_{P}(t, T)\right\| \\
& \quad=\alpha(T) P\left(t, T_{m-1}, T_{m}\right)\left\|\sigma_{P}\left(t, T_{m-1}\right)-\sigma_{P}\left(t, T_{m}\right)\right\|
\end{aligned}
$$

Then ${ }^{2}$

$$
\begin{align*}
\alpha(T) & =\frac{P\left(t, T_{m-1}, T\right)}{P\left(t, T_{m-1}, T_{m}\right)} \frac{\left\|\sigma_{P}\left(t, T_{m-1}\right)-\sigma_{P}(t, T)\right\|}{\left\|\sigma_{P}\left(t, T_{m-1}\right)-\sigma_{P}\left(t, T_{m}\right)\right\|} \\
& \approx \frac{P\left(0, T_{m-1}, T\right)}{P\left(0, T_{m-1}, T_{m}\right)} \frac{\left\|\sigma_{P}\left(t, T_{m-1}\right)-\sigma_{P}(t, T)\right\|}{\left\|\sigma_{P}\left(t, T_{m-1}\right)-\sigma_{P}\left(t, T_{m}\right)\right\|} \tag{15.8}
\end{align*}
$$
^eq15-8

As we have seem in Chapter 14, the LM model does not uniquely define all the bond volatilities in (15.8) and we would need to interpolate available volatilities to compute (15.8). Note that (15.8) turns the problem of interpolating bond prices into a problem of interpolating bond volatilities instead. This point of view is advantageous as there are few, if any, arbitrage restrictions on interpolating the volatilities of bonds, as opposed to the bonds themselves. One can, for example, choose a linear interpolation to obtain $\sigma_{P}(t, T)$ from $\sigma_{P}\left(t, T_{m-1}\right)$ and $\sigma_{P}\left(t, T_{m}\right)$ or, for a more sophisticated scheme, draw inspiration from the shape of forward volatilities in a meanreverting one-factor Gaussian model. To explore the latter idea, recall that in the one-dimensional Gaussian model with constant mean reversion (Section 10.1.2),

$$
\begin{equation*}
\sigma_{P}(t, T)=\sigma(t) \frac{1-e^{-\varkappa(T-t)}}{\varkappa} \tag{15.9}
\end{equation*}
$$
^eq15-9

where $\varkappa$ is the mean reversion and $\sigma(t)$ is the short rate volatility. Plugging into (15.8) and rearranging we get[^70]

$$
\alpha(T) \approx \frac{P\left(0, T_{m-1}, T\right)}{P\left(0, T_{m-1}, T_{m}\right)} \frac{1-e^{-\varkappa\left(T-T_{m-1}\right)}}{1-e^{-\varkappa\left(T_{m}-T_{m-1}\right)}}
$$

The mean reversion parameter $\varkappa$ could be either set as part of the user input, or obtained by best-fitting the Gaussian parametric form (15.9) to the volatility structure of the LM model.

### 15.1.3 Back Stub, Interpolation Inspired by the Gaussian Model

Above, we used a volatility parameterization inspired by the Gaussian model to construct an interpolation scheme for bond prices. An alternative, more direct, approach to extracting information from a Gaussian model in the interpolation exercise is to let the bond reconstitution formulas from the Gaussian model form the basis for interpolation. To demonstrate, recall that in the one-factor Gaussian model with constant mean reversion $\varkappa$ and short rate volatility $\sigma(t)$ (see Section 10.1.2.2),

$$
\begin{align*}
P\left(t, T_{m-1}, T\right)= & P\left(0, T_{m-1}, T\right) \exp \left(-G\left(T_{m-1}, T\right) e^{-\varkappa\left(T_{m-1}-t\right)} x(t)\right) \\
& \times \exp \left(-\frac{1}{2}\left(G(t, T)^{2}-G\left(t, T_{m-1}\right)^{2}\right) y(t)\right), \tag{15.10}
\end{align*}
$$
^eq15-10

and

$$
\begin{align*}
P\left(t, T_{m-1}, T_{m}\right)= & P\left(0, T_{m-1}, T_{m}\right) \exp \left(-G\left(T_{m-1}, T_{m}\right) e^{-\varkappa\left(T_{m-1}-t\right)} x(t)\right) \\
& \times \exp \left(-\frac{1}{2}\left(G\left(t, T_{m}\right)^{2}-G\left(t, T_{m-1}\right)^{2}\right) y(t)\right), \tag{15.11}
\end{align*}
$$
^eq15-11

where

$$
\begin{align*}
G(t, T) & =\frac{1-e^{-\varkappa(T-t)}}{\varkappa}  \tag{15.12}\\
y(t) & =e^{-2 \varkappa t} \int_{0}^{t} e^{2 \varkappa s} \sigma(s)^{2} d s \tag{15.13}
\end{align*}
$$
^eq15-12
^eq15-13

and $x(t)$ is the short rate state. Eliminating $x(t)$ in (15.10) and (15.11) defines a relationship between bond prices,

$$
\begin{aligned}
& \frac{1}{G\left(T_{m-1}, T\right)}\left(\ln \frac{P\left(t, T_{m-1}, T\right)}{P\left(0, T_{m-1}, T\right)}+\left(G(t, T)^{2}-G\left(t, T_{m-1}\right)^{2}\right) \frac{y(t)}{2}\right) \\
= & \frac{1}{G\left(T_{m-1}, T_{m}\right)}\left(\ln \frac{P\left(t, T_{m-1}, T_{m}\right)}{P\left(0, T_{m-1}, T_{m}\right)}+\left(G\left(t, T_{m}\right)^{2}-G\left(t, T_{m-1}\right)^{2}\right) \frac{y(t)}{2}\right)
\end{aligned}
$$

or, after a few additional manipulations,

$$
\begin{align*}
P\left(t, T_{m-1}, T\right)= & P\left(0, T_{m-1}, T\right)\left(\frac{P\left(t, T_{m-1}, T_{m}\right)}{P\left(0, T_{m-1}, T_{m}\right)}\right)^{\frac{G\left(T_{m-1} \cdot T\right)}{G\left(T_{m-1} \cdot T_{m}\right)}} \\
& \times \exp \left(\frac{1}{2} \frac{G\left(T_{m-1}, T\right)}{G\left(T_{m-1}, T_{m}\right)}\left(G\left(t, T_{m}\right)^{2}-G\left(t, T_{m-1}\right)^{2}\right) y(t)\right) \\
& \times \exp \left(-\frac{1}{2}\left(G(t, T)^{2}-G\left(t, T_{m-1}\right)^{2}\right) y(t)\right) \tag{15.14}
\end{align*}
$$
^eq15-14

In (15.14), the volatility $\sigma(t)$ and the mean reversion $\varkappa$ can be obtained, for example, by fitting the Gaussian volatility structure to the volatilities of Libor rates generated by the LM model itself. High level of precision is not required here; we can take, say, $\varkappa=0$ and $\sigma(t)=\left\|\sigma_{m}(t)\right\|$ where $\sigma_{m}(t)$ is the vector of volatilities for the $m$-th Libor rate $L_{m}(t)$ (the comment of footnote 2 applies here as well, to Libor volatilities).

The scheme (15.14) would be arbitrage-free in the context of a Gaussian model, i.e. if the expected value in (15.5) were computed in the Gaussian model. In the LM model the equality (15.5) would not hold exactly, but would be a good approximation as long as the choice of the volatility/mean reversion in the scheme is reasonably consistent with the actual LM model volatility structure. While we have no strong opinions on the matter, on the whole (15.14) tends to be our preferred choice for the back-stub interpolation as it is nearly arbitrage-free and is perhaps somewhat easier to implement and maintain than other schemes from Section 15.1.2.

### 15.1.4 Front Stub, Zero Volatility

Having considered various options for the back stub, let us now focus on the front stub. The simplest way to "complete" the LM model volatility specification is to specify that the front stub bond has no volatility, i.e.

$$
\begin{equation*}
\sigma_{P}\left(t, T_{q(t)}\right) \equiv 0 \tag{15.15}
\end{equation*}
$$
^eq15-15

in the notation of Section 14.2.3. This automatically specifies the front stub interpolation scheme, as

$$
P\left(t, T_{q(t)}\right)=P\left(T_{q(t)-1}, t, T_{q(t)}\right)
$$

where the right-hand side can be computed from the previous results on back stub interpolation.

The choice (15.15) was first proposed in Brace et al. [1997] and implies that the continuously rolling money market account $\beta(t)$ coincides with the discrete numeraire $B(t)$, whereby the risk-neutral measure is identical to the spot Libor measure. The bond volatilities for "core" bonds - that is, the bonds paying on the dates in the tenor structure - are explicitly given by

$$
\sigma_{P}\left(t, T_{n}\right)=\sum_{j=q(t)}^{n-1} \frac{\tau_{j} \sigma_{j}(t)}{1+\tau_{j} L_{j}(t)}
$$

While certainly technically convenient, (15.15) leads to unrealistic dynamics of the yield curve. To give an example, assume that the LM model is specified with a 6 month tenor; i.e. $\tau_{n}$ 's are all approximately 0.5 and $T_{1}=0.5$. For a 3 month option on a 3 month rate, i.e. a caplet with the payoff

$$
\begin{equation*}
(L(0.25,0.25,0.5)-K)^{+} \tag{15.16}
\end{equation*}
$$
^eq15-16

(15.15) implies that the rate $L(t, 0.25,0.5)$ has no volatility for $t \in[0,0.25]$, i.e. the option with the payoff (15.16) will be priced at its intrinsic value. Clearly, this is unrealistic and will not be consistent with an actual market price. This example is not as contrived as one might think, as many complex derivatives have at least some exposure to short-expiry, short-tenor volatility. The specification (15.15) forces such volatility to zero, and cannot be recommended.

The analysis of the shortcomings of (15.15) puts into focus the requirements that one would wish to impose on all "good" front-stub interpolation schemes. In particular, we will be looking for schemes that recover both market-implied forward rates, as well as volatilities for rates with nonstandard expiries and tenors. With "non-standard", we here refer to rates that are not aligned with the tenor structure. For example, for a 6 month tenor LM model, we can examine the dynamics and values of 3 month tenor forward rates maturing in 1 month, 2 months, ..., or, in 3 months time, we can look at rates with 3 month, 6 month, 9 month, . . . tenors. Note that such an analysis will also help uncover problems with the back-stub interpolation, should there be any.

### 15.1.5 Front Stub, Exogenous Volatility

At a fundamental level, the volatilities $\sigma_{P}\left(t, T_{q(t)}\right)$ for all $t$ constitute extra information that is required to specify the model behavior between tenor dates. This information would define the dynamics of the stub bonds $P\left(t, T_{q(t)}\right)$ (for all $t$ ) which, in principle, is all that is required to define the values of all bonds, as for any $0 \leq s<t$,

$$
\begin{equation*}
P(s, t)=P\left(s, T_{q(s)}\right) P\left(s, T_{q(s)}, T_{q(t)}\right) \mathrm{E}_{s}^{T_{q(t)}}\left(\frac{1}{P\left(t, T_{q(t)}\right)}\right) \tag{15.17}
\end{equation*}
$$
^eq15-17

As this formula specifies the values of all discount bonds, in principle it makes back-stub interpolation methods from earlier in the section redundant. In practice, however, calculation of the expected value in (15.17) is non-trivial, especially for $s \ll t$, so the application of (15.17) is best left for the case of $t-s$ being small, when accurate approximations are easier to obtain.

The information on the stub bond volatilities $\sigma_{P}\left(t, T_{q(t)}\right)$ should be supplied in addition to the basic Libor rate volatility structure of the LM model. To make such exogenous volatility specification easier on the model user, we recommend that some additional structure is added through
simplifying assumptions. For instance, an obvious simplification is to assume that the stub volatilities are both time-homogeneous and identical for all periods in the tenor structure; in other words, $\sigma_{P}(t, u)=\sigma_{\text {stub }}(u-t)$ for some chosen $\sigma_{\text {stub }}(\tau), t \leq u<T_{q(t)}$. This reduces the problem to specifying a stub volatility function of only one argument; this function could be calibrated to short-dated options on short-tenor rates.

In the scheme above, we notice that $\sigma_{\text {stub }}$ is vector-valued. To reduce the specification burden further, we can inherit the correlation structure from that of the core Libor forward rates, and set

$$
\sigma_{\text {stub }}(u-t)=\frac{\left\|\sigma_{\text {stub }}(u-t)\right\|}{\left\|\sigma_{1}(t)\right\|} \sigma_{1}(t),
$$

where $\sigma_{1}(t)$ is the vector of volatilities for the front Libor rate $L_{1}(t)$. With this scheme, we only need to specify a scalar function $\left\|\sigma_{\text {stub }}(\tau)\right\|$.

The dimensionality of the additional model inputs could be further reduced by using a particular parametric form. Drawing inspiration - yet again - from the one-factor Gaussian model, we could for instance specify-

$$
\begin{equation*}
\left\|\sigma_{\text {stub }}(\tau)\right\|=\bar{\sigma}_{\text {stub }} \frac{1-e^{-\varkappa_{\text {stub }}, \tau}}{\varkappa_{\text {stub }}} \tag{15.18}
\end{equation*}
$$
^eq15-18

for given constants $\bar{\sigma}_{\text {stub }}$ and $\varkappa_{\text {stub }}$. With this, we have reduced the specification of the front-stub volatility function to just two constants.

Once the stub volatility function is specified, it needs to be incorporated into a Monte Carlo simulation. Suppose we need to perform a time step from $t$ to $t^{\prime}$ where $t$ is one of the tenor dates, $t=T_{m-1}$ and $t^{\prime}<T_{m}$. Then, in addition to the standard MC step (see (14.64))

$$
\widehat{L}_{n}\left(t^{\prime}\right)=\widehat{L}_{n}(t)+\sigma_{n}(t)^{\top} \mu_{n}(t)\left(t^{\prime}-t\right)+\widehat{M}_{n}\left(t, t^{\prime}\right), n=m, \ldots, N-1
$$

we also need an update equation for $P\left(t^{\prime}, T_{m}\right)$. Over the period $s \in \left[T_{m-1}, T_{m}\right]$, the spot numeraire coincides, up to a constant, with the bond price $P\left(t, T_{m}\right)$ (see (14.8)). Hence, in spot measure, the forward bond price $P\left(s, t^{\prime}, T_{m}\right), t \leq s \leq t^{\prime}$, satisfies

$$
\begin{aligned}
& d P\left(s, t^{\prime}, T_{m}\right) / P\left(s, t^{\prime}, T_{m}\right)=\left\|\sigma_{\text {stub }}\left(s, t^{\prime}, T_{m}\right)\right\|^{2} d s \\
& \quad+\sigma_{\text {stub }}\left(s, t^{\prime}, T_{m}\right)^{\top} d W^{B}(s)
\end{aligned}
$$

where we have defined

$$
\sigma_{\text {stub }}\left(s, t^{\prime}, T_{m}\right) \triangleq \sigma_{\text {stub }}\left(T_{m}-s\right)-\sigma_{\text {stub }}\left(t^{\prime}-s\right)
$$

A simple log-Euler scheme for the bond is given by

$$
\begin{align*}
\widehat{P}\left(t^{\prime}, T_{m}\right)=\widehat{P}\left(t, t^{\prime}, T_{m}\right) & \exp \left(\sigma_{\text {stub }}\left(t, t^{\prime}, T_{m}\right)^{\top} \sqrt{t^{\prime}-t} Z\right) \\
& \times \exp \left(\frac{1}{2}\left\|\sigma_{\text {stub }}\left(t, t^{\prime}, T_{m}\right)\right\|^{2}\left(t^{\prime}-t\right)\right) \tag{15.19}
\end{align*}
$$
^eq15-19

where $Z$ is a draw from the standard one-dimensional Gaussian distribution. This scheme, together with the update equations for $\left\{\widehat{L}_{n}(t)\right\}_{n=m}^{N-1}$, defines the Monte Carlo step. Extensions to more sophisticated discretization schemes for $\widehat{P}\left(t^{\prime}, T_{m}\right)$ are straightforward, but rarely needed.

Since we assumed that $t=T_{m-1}$, the term $\widehat{P}\left(t, t^{\prime}, T_{m}\right)$ in (15.19) is available at time $t$ from Libor rates using back-stub interpolation methods discussed earlier in the section. If we, in addition, need to calculate $\widehat{P}\left(t^{\prime \prime}, T_{m}\right)$ for some $t^{\prime \prime}, t^{\prime}<t^{\prime \prime}<T_{m}$, we can either apply (15.19) to step from $t$ to $t^{\prime \prime}$ directly, or to suitably modify (15.19) to step from $t^{\prime}$ to $t^{\prime \prime}$. In the latter case we will, however, need to be able to simulate other short-dated discount bond values, i.e. $P\left(t^{\prime}, u\right)$ for $u<T_{m}$ (something that may be required for other purposes as well); this can be incorporated into the time-stepping algorithm in the same way as for $P\left(t^{\prime}, T_{m}\right)$.

More pragmatically, instead of simulating a whole collection of shorttenor bonds $P\left(\cdot, u_{1}\right), \ldots, P\left(\cdot, u_{k}\right), t^{\prime}<u_{1}<\ldots<u_{k} \leq T_{m}$, that may be required at time $t^{\prime}$, we may choose to propagate a single discount bond with the shortest time to maturity, i.e. $P\left(\cdot, u_{1}\right)$, with other bond prices subsequently obtained by one of the back-stub interpolation schemes discussed earlier. Alternatively we could propagate the stub bond $P\left(\cdot, T_{m}\right)$ and use (15.17) to approximate all $P\left(t^{\prime}, u_{1}\right), \ldots, P\left(t^{\prime}, u_{k-1}\right)$ (with, perhaps, some approximation to calculate the expected values):

$$
P\left(t^{\prime}, u_{l}\right)=P\left(t^{\prime}, T_{m}\right) \mathrm{E}_{t^{\prime}}^{T_{m}}\left(\frac{1}{P\left(u_{l}, T_{m}\right)}\right) .
$$

Proceeding by either of these methods, the yield curve at time $t^{\prime}$, as well as the spot numeraire, are then fully defined by the augmented state vector $\left\{P\left(t^{\prime}, T_{q\left(t^{\prime}\right)}\right), L_{q\left(t^{\prime}\right)}\left(t^{\prime}\right), \ldots, L_{N-1}\left(t^{\prime}\right)\right\}$. Moreover, by construction, the Monte Carlo scheme recovers the expected values of short-tenor rates and their volatilities, as specified by the stub volatility function (within Monte Carlo error and up to discretization bias, of course). The disadvantage of the scheme is the need to carry an extra state variable for each time step.

Remark 15.1.1. Parameterizations such as (15.18) introduce volatility functions (of stub bonds) that are not constant between tenor dates. A similar idea could be applied to core Libor rate volatilities as well, as an extension of the simple piecewise constant interpolation we discussed in Section 14.5.3. This could be required to be able to match shorter-dated options on Libor rates, e.g. a 3 month option on a 6 month Libor rate in an LM model with 6 months between dates in the tenor structure. For example, drawing inspiration from a Gaussian model we could amend (14.42) by a final-period Libor volatility specification of the form
$\left\|\lambda_{k}(t)\right\|=e^{\varkappa_{\text {shirrt }}\left(t-T_{k-1}\right)}\left\|\lambda_{k}\left(T_{k-1}\right)\right\|=e^{\varkappa_{\text {hihort }}\left(t-T_{k-1}\right)}\left\|\lambda_{k, k}\right\|, t \in\left[T_{k-1}, T_{k}\right)$, where $\varkappa_{\text {short }}$ is user-specified or calibrated to short-dated options.

### 15.1.6 Front Stub, Simple Interpolation

As an alternative to the scheme in the previous section, one could contemplate avoiding simulation of the stub bond altogether by interpolating it from the original state vector of Libor rates. To demonstrate, let us assume the same setup as in the previous section. Then, at time $t^{\prime}$, simulated forward bond prices $P\left(t^{\prime}, T_{n}, T_{n+1}\right), n \geq m=q\left(t^{\prime}\right)$, are available; the goal is to construct $P\left(t^{\prime}, T_{m}\right)$ from these. A simple interpolation scheme could, for example, specify that the simulated instantaneous forward rates $f\left(t^{\prime}, u\right)$ are constant over $u \in\left[t^{\prime}, T_{m+1}\right]$. This would link $P\left(t^{\prime}, T_{m}\right)$ to $P\left(t^{\prime}, T_{m}, T_{m+1}\right)$ in the following way:

$$
\begin{equation*}
P\left(t^{\prime}, T_{m}\right)=P\left(t^{\prime}, T_{m}, T_{m+1}\right)^{\frac{T_{m}-t^{\prime}}{T_{m+1}-T_{m}}} . \tag{15.20}
\end{equation*}
$$
^eq15-20

This scheme, while simple, does not satisfy our notion of a "good" one. In particular, it does not recover the time 0 forward discount bond value, as generally

$$
P\left(0, t^{\prime}, T_{m}\right) \neq \mathrm{E}^{t^{\prime}}\left(P\left(t^{\prime}, T_{m}, T_{m+1}\right)^{\frac{T_{m-} \prime^{\prime}}{T_{m+1}-T_{m}}}\right)
$$

Nor does the scheme recover the market volatility of $P\left(t^{\prime}, T_{m}\right)$, as now

$$
\ln P\left(t^{\prime}, T_{m}\right)=\frac{T_{m}-t^{\prime}}{T_{m+1}-T_{m}} \ln P\left(t^{\prime}, T_{m}, T_{m+1}\right)
$$

and the volatilities of $\ln P\left(t^{\prime}, T_{m}\right)$ and $\ln P\left(t^{\prime}, T_{m}, T_{m+1}\right)$ are then linked in a specific way that does not necessarily hold in the market.

Both of the shortcomings above can be addressed in a more sophisticated interpolation scheme, to be discussed in the next section. Before proceeding, however, we note that a similar, albeit somewhat more general, scheme has been proposed by Schlögl [2002]:

$$
\frac{1}{P\left(t^{\prime}, T_{m}\right)}=1+\left(T_{m}-t^{\prime}\right)\left(\xi\left(t^{\prime}\right) L_{m-1}\left(T_{m-1}\right)+\left(1-\xi\left(t^{\prime}\right)\right) L_{m}\left(t^{\prime}\right)\right)
$$

where an essentially arbitrary deterministic function $\xi(t)$ is chosen so that $1=\xi\left(T_{m-1}\right) \geq \xi\left(t^{\prime}\right) \geq \xi\left(T_{m}\right)=0$. While the volatility of the stub bond $P\left(t^{\prime}, T_{m}\right)$ could be manipulated by the choice of the function $\xi(t)$, the forward discount bond value is still not recovered by the model.

### 15.1.7 Front Stub, Interpolation Inspired by the Gaussian Model

The idea of employing a bond reconstitution formula from a Gaussian model, already used for the back stub, can be applied to the front stub as well. Here we assume that for short tenors, the LM model can be locally approximated by a one-factor Gaussian model. Recall that in the latter,

$$
P\left(t^{\prime}, T_{m}\right)=P\left(0, t^{\prime}, T_{m}\right) \exp \left(-G\left(t^{\prime}, T_{m}\right) x\left(t^{\prime}\right)-\frac{1}{2} G\left(t^{\prime}, T_{m}\right)^{2} y\left(t^{\prime}\right)\right)
$$

and

$$
\begin{aligned}
P\left(t^{\prime}, T_{m}, T_{m+1}\right)= & P\left(0, T_{m}, T_{m+1}\right) \\
& \times \exp \left(-\left(G\left(t^{\prime}, T_{m+1}\right)-G\left(t^{\prime}, T_{m}\right)\right) x\left(t^{\prime}\right)\right) \\
& \times \exp \left(-\frac{1}{2}\left(G\left(t^{\prime}, T_{m+1}\right)^{2}-G\left(t^{\prime}, T_{m}\right)^{2}\right) y\left(t^{\prime}\right)\right)
\end{aligned}
$$

Solving the second equation for $x\left(t^{\prime}\right)$ and substituting into the first one, we obtain

$$
\begin{aligned}
& \ln \frac{P\left(t^{\prime}, T_{m}\right)}{P\left(0, t^{\prime}, T_{m}\right)}+\frac{1}{2} G\left(t^{\prime}, T_{m}\right)^{2} y\left(t^{\prime}\right)=\frac{G\left(t^{\prime}, T_{m}\right)}{G\left(t^{\prime}, T_{m+1}\right)-G\left(t^{\prime}, T_{m}\right)} \\
& \quad \times\left(\ln \frac{P\left(t^{\prime}, T_{m}, T_{m+1}\right)}{P\left(0, T_{m}, T_{m+1}\right)}+\frac{1}{2}\left(G\left(t^{\prime}, T_{m+1}\right)^{2}-G\left(t^{\prime}, T_{m}\right)^{2}\right) y\left(t^{\prime}\right)\right)
\end{aligned}
$$

so that

$$
\begin{align*}
P\left(t^{\prime}, T_{m}\right)=P\left(0, t^{\prime}, T_{m}\right) & \left(\frac{P\left(t^{\prime}, T_{m}, T_{m+1}\right)}{P\left(0, T_{m}, T_{m+1}\right)}\right)^{\frac{G\left(t^{\prime} \cdot T_{m}\right)}{G\left(t^{\prime} \cdot T_{m+1}\right)-c\left(t^{\prime} \cdot T_{m}\right)}} \\
& \times \exp \left(\frac{1}{2} G\left(t^{\prime}, T_{m}\right) G\left(t^{\prime}, T_{m+1}\right) y\left(t^{\prime}\right)\right) \tag{15.21}
\end{align*}
$$
^eq15-21

The interpolation scheme depends on two parameters, the mean reversion $\varkappa$ in (15.12) and the short rate volatility $\sigma(t)$ in (15.13) of the Gaussian model. The latter could be approximated by the volatilities of the front Libor rate, ensuring that the forward price of the front-stub discount bond is approximately recovered by the model. The mean reversion plays an important role here as it defines the relative magnitude of the volatility of $(\log ) P\left(t^{\prime}, T_{m}\right)$ in relationship to $(\log ) P\left(t^{\prime}, T_{m}, T_{m+1}\right)$. As such, $\varkappa$ can be used to set the stub bond volatility to, or near, its market-implied value.

Empirical evidence shows a close fit between time 0 market-observed short-tenor rates and those computed from the model in the manner described above. Also, with the right choice of mean reversion, the market-implied front-stub volatility is recovered as well. We recommend this scheme for most applications.

### 15.2 Advanced Approximations for Swaption Pricing via Markovian Projection

Having wrapped up the topic of interpolation, let us now go back to the problem of approximate pricing of vanilla options - such as swaptions

- in LM models. While the pricing approximation for swaptions derived in Section 14.4.2 has proved to be remarkably successful for calibration purposes, it does have limitations, especially for longer-dated options on longer-dated swaps and for swaptions that are not at-the-money. There are several improvements that could be made to the basic approximation. For instance, in (14.33) the skew of a swap rate is taken to be the same as the (common) skew of all Libor rates, yet numerical simulation shows that this is not the case. For instance, in a log-normal LM model, a long-expiry, long-tenor swap rate would have a "super log-normal" skew, i.e. implied log-normal swaption volatilities would trend up with the strike. Hence, a more accurate estimate of the swap rate skew from Libor rate skews would be useful. Moreover, the accuracy of the swap rate volatility calculations can be improved by a more careful analysis than that of Section 14.4.2.

Before proceeding, let us first enlarge the model setup from Section 14.2.5 somewhat. Specifically, we wish to address the fact that the specification of the model in Section 14.2.4 (or Section 14.2.5) uses the same timehomogeneous local volatility function for all Libor rates. Such a setup implies that swaptions of all tenors and expiries have essentially the same volatility skew, a model feature that is inconsistent with current market reality. For the Libor market model to be able to match the swaption volatilities for all expiries, tenors and strikes, it is necessary to assume that different Libor rates have different local volatility functions, and that those functions explicitly depend on time. In the stochastic volatility model, it may also be necessary to assume that the volatility of variance is time-dependent, so that the curvatures of smiles of different swaptions are allowed to differ. More advanced methods are then required to derive approximations to swap rate volatilities in such a, more generic, specification. In total, we therefore consider the following generalization ${ }^{3}$ of the specification (14.15)-(14.16):

$$
\begin{align*}
d z(t) & =\theta\left(z_{0}-z(t)\right) d t+\eta(t) \psi(z(t)) d Z(t)  \tag{15.22}\\
d L_{n}(t) & =\sqrt{z(t)} \varphi_{n}\left(t, L_{n}(t)\right) \lambda_{n}(t)^{\top} d W^{T_{n+1}}(t), \quad n=1, \ldots, N-1 \tag{15.23}
\end{align*}
$$
^eq15-22
^eq15-23

In particular, the volatility of variance parameter $\eta(t)$ depends on $t$, and the DVFs $\varphi_{n}(t, x)$ are now specific to each Libor rate (i.e. depend on $n$ ) and may also depend on time $t$. Without loss of generality, we assume

$$
\varphi_{n}\left(t, L_{n}(0)\right)=1
$$

As in other applications of DVF modeling, we assume that the functions $\varphi_{n}(t, x)$ are well-approximated by their first-order expansions,[^71]

$$
\begin{align*}
\varphi_{n}(t, x) & \approx 1+b_{n}(t)\left(x-L_{n}(0)\right)  \tag{15.24}\\
b_{n}(t) & \triangleq \frac{\partial}{\partial x} \varphi_{n}\left(t, L_{n}(0)\right)
\end{align*}
$$
^eq15-24

In practical applications, this usually means that $\varphi_{n}(t, x)$ 's are either linear or power functions, see Table 14.1.

As we did in Section 14.6.2.1, let us denote by $\mathbf{L}(t)$ the vector of all Libor rates, i.e. $\mathbf{L}(t)=\left(L_{1}(t), \ldots, L_{N-1}(t)\right)^{\top}$, with the convention that $L_{i}(t) \equiv L_{i}\left(T_{i}\right)$ for $i<q(t)$. Continuing with the notations of Section 14.4.2, let $S(t)=S_{j, k-j}(t)$ be a particular swap rate. The dynamics of $S(t)$ in the model (15.22)-(15.23) are easy to write down,

$$
\begin{equation*}
d S(t)=\sqrt{z(t)} \lambda_{S}(t, \mathbf{L}(t))^{\top} d W^{A}(t) \tag{15.25}
\end{equation*}
$$
^eq15-25

where

$$
\begin{equation*}
\lambda_{S}(t, \mathbf{L}(t))=\sum_{n=j}^{k-1} \frac{\partial S(t)}{\partial L_{n}(t)} \varphi_{n}\left(t, L_{n}(t)\right) \lambda_{n}(t) \tag{15.26}
\end{equation*}
$$
^eq15-26

and $W^{A}(t)$ is a Brownian motion in the annuity measure $\mathrm{Q}^{A}$ for $S(t)$. Moreover, the dynamics can be written in a one-dimensional form,

$$
\begin{equation*}
d S(t)=\sqrt{z(t)}\left\|\lambda_{S}(t, \mathbf{L}(t))\right\| d Y^{A}(t) \tag{15.27}
\end{equation*}
$$
^eq15-27

where $Y^{A}(t)$ is a one-dimensional Brownian motion in $\mathrm{Q}^{A}$.
The following result serves as a starting point for various useful approximations.

Proposition 15.2.1. For the purposes of European swaption valuation, the dynamics of the swap rate $S(t)$ in $\mathrm{Q}^{A}$ are approximately given by the following displaced log-normal stochastic volatility SDE

$$
\begin{equation*}
d S(t) \approx \sqrt{z(t)} p_{S}(t)\left(1+b_{S}(t)(S(t)-S(0))\right) d Y^{A}(t) \tag{15.28}
\end{equation*}
$$
^eq15-28

with

$$
\begin{align*}
& p_{S}(t)=\left\|\lambda_{S}\left(t, \mathrm{E}^{A}(\mathbf{L}(t))\right)\right\|  \tag{15.29}\\
& b_{S}(t)=\frac{1}{p_{S}(t)} \sum_{n=j}^{k-1} \frac{\partial\left\|\lambda_{S}\left(t, \mathrm{E}^{A}(\mathbf{L}(t))\right)\right\|}{\partial L_{n}(t)} \frac{\int_{0}^{t}\left(\lambda_{n}(u)^{\top} \lambda_{S}(u, \mathbf{L}(0))\right) d u}{\int_{0}^{t}\left\|\lambda_{S}(u, \mathbf{L}(0))\right\|^{2} d u} \tag{15.30}
\end{align*}
$$
^eq15-29
^eq15-30

Proof. The proof relies on standard results on Markovian projection, see Appendix A. From (A.18), the European options on $S(t)$ in the model (15.27) have the same values as in the Markovian model

$$
d S(t)=\sqrt{z(t)}\left(\mathrm{E}^{A}\left(\left\|\lambda_{S}(t, \mathbf{L}(t))\right\|^{2} \mid S(t)\right)\right)^{1 / 2} d Y^{A}(t)
$$

First, we approximate

$$
\left(\mathrm{E}^{A}\left(\left\|\lambda_{S}(t, \mathbf{L}(t))\right\|^{2} \mid S(t)\right)\right)^{1 / 2} \approx \mathrm{E}^{A}\left(\left\|\lambda_{S}(t, \mathbf{L}(t))\right\| \mid S(t)\right)
$$

and then linearize $\left\|\lambda_{S}(t, \mathbf{L}(t))\right\|$ around $\mathrm{E}^{A}(\mathbf{L}(t))$,

$$
\begin{aligned}
\left\|\lambda_{S}(t, \mathbf{L}(t))\right\| \approx \| \lambda_{S}(t, & \left.\mathrm{E}^{A}(\mathbf{L}(t))\right) \| \\
& +\left(\nabla\left\|\lambda_{S}\left(t, \mathrm{E}^{A}(\mathbf{L}(t))\right)\right\|\right)\left(\mathbf{L}(t)-\mathrm{E}^{A}(\mathbf{L}(t))\right)
\end{aligned}
$$

where $\nabla=\left(\frac{\partial}{\partial L_{1}(t)}, \ldots, \frac{\partial}{\partial L_{N-1}(t)}\right)$ is the (row-vector) gradient.
Let us introduce the Gaussian approximation

$$
\begin{equation*}
d \widehat{L}_{n}(t)=\lambda_{n}(t)^{\top} d W^{A}(t), \quad d \widehat{S}(t)=\lambda_{S}(t, \mathbf{L}(0))^{\top} d W^{A}(t), \tag{15.31}
\end{equation*}
$$
^eq15-31

so that we can approximate

$$
\begin{aligned}
\mathrm{E}^{A}\left(\mathbf{L}(t)-\mathrm{E}^{A}(\mathbf{L}(t)) \mid S(t)=s\right) & \approx \mathrm{E}^{A}(\widehat{\mathbf{L}}(t)-\mathbf{L}(0) \mid \widehat{S}(t)=s) \\
& =\frac{\mathrm{E}^{A}((\widehat{\mathbf{L}}(t)-\mathbf{L}(0))(\widehat{S}(t)-S(0)))}{\mathrm{E}^{A}((\widehat{S}(t)-S(0))(\widehat{S}(t)-S(0)))}(s-S(0))
\end{aligned}
$$

where we have defined $\widehat{\mathbf{L}}(t)=\left(\widehat{L}_{1}, \ldots, \widehat{L}_{N-1}\right)^{\top}$. The result follows.
The price of a European swaption is given by the value of the option on $S\left(T_{j}\right)$ (times the annuity). The model (15.22) and (15.28) is a stochastic volatility model with time-dependent parameters. Using the methods of Chapter 9, effective, time-constant parameters can be derived, to facilitate fast pricing of European swaptions, as well as calibration of the model parameters to their market values. We do not repeat the relevant formulas here, but simply note that the total volatility, skew and volatility of variance of any swap rate are available as functions of the model parameters.

### 15.2.1 Advanced Formula for Swap Rate Volatility

In this section, we use the results of Proposition 15.2.1 to derive useful formulas for the swaption volatility. Recall (15.28), the approximate SDE for the swap rate used for European option pricing. The function $p_{S}(t)$ in (15.29) is the (time-dependent) swap rate volatility,

$$
\begin{aligned}
p_{S}(t) & =\left\|\lambda_{S}\left(t, \mathrm{E}^{A}(\mathbf{L}(t))\right)\right\| \\
& =\left\|\left.\sum_{n=j}^{k-1} \frac{\partial S(t)}{\partial L_{n}(t)}\right|_{\mathbf{L}(t)=\mathrm{E}^{A}(\mathbf{L}(t))} \varphi_{n}\left(t, \mathrm{E}^{A}\left(L_{n}(t)\right)\right) \lambda_{n}(t)\right\|,
\end{aligned}
$$

where we have used (15.26). The equivalent quantity in the standard approximation of Section 14.4.2 is given by (compare to (14.33))

$$
p_{S, \text { standard }}(t)=\left\|\left.\sum_{n=j}^{k-1} \frac{\partial S(t)}{\partial L_{n}(t)}\right|_{\mathbf{L}(t)=\mathbf{L}(0)} \varphi_{n}\left(t, L_{n}(0)\right) \lambda_{n}(t)\right\| .
$$

Hence, the improvements over the standard approximation come from evaluating the actual volatility function of the swap rate $\left(\lambda_{S}(t, \mathbf{L}(t))\right)$ at $\mathbf{L}(t)=\mathrm{E}^{A}(\mathbf{L}(t))$ rather than at $\mathbf{L}(t)=\mathbf{L}(0)$; this is similar to the improvements obtained for the quasi-Gaussian model in Section 13.1.4.2. Clearly $\mathrm{E}^{A}(\mathbf{L}(t)) \neq \mathbf{L}(0)$; the difference can be approximated with the help of the following proposition.

Proposition 15.2.2. For $j \leq n \leq k-1$, the expected value of the $n$-th Libor rate in the annuity measure is approximately given by

$$
\mathrm{E}^{A}\left(L_{n}(t)\right)=L_{n}(0)\left(1+c_{n}(t)\right)
$$

where

$$
c_{n}(t)=\frac{1}{L_{n}(0) Q_{n}(0)} \sum_{i=j}^{k-1} \frac{\partial Q_{n}(0)}{\partial L_{i}(0)} \int_{0}^{t}\left(\lambda_{i}^{\top}(s) \lambda_{n}(s)\right) d s, Q_{n}(t)=\frac{A(t)}{P\left(t, T_{n+1}\right)},
$$

with $A(t)$ defined in (14.29).
Proof. We have,

$$
\mathrm{E}^{A}\left(L_{n}(t)\right)=Q_{n}^{-1}(0) \mathrm{E}^{T_{n+1}}\left(Q_{n}(t) L_{n}(t)\right)
$$

Both $Q_{n}(t)$ and $L_{n}(t)$ are martingales in $\mathrm{Q}^{T_{n+1}}$. Applying Gaussian approximations,

$$
d \widehat{Q}_{n}(t) \approx \lambda_{Q_{n}}(t)^{\top} d W^{T_{n+1}}(t), \quad d \widehat{L}_{n}(t) \approx \lambda_{n}^{\top}(t) d W^{T_{n+1}}(t)
$$

where

$$
\lambda_{Q_{n}}(t)=\sum_{i=j}^{k-1} \frac{\partial Q_{n}(0)}{\partial L_{i}(0)} \lambda_{i}(t)
$$

we obtain

$$
\mathrm{E}^{T_{n+1}}\left(Q_{n}(t) L_{n}(t)\right) \approx Q_{n}(0) L_{n}(0)+\int_{0}^{t}\left(\lambda_{Q_{n}}(s)^{\top} \lambda_{n}(s)\right) d s
$$

and the statement follows.
The idea of employing $\mathrm{E}^{A}(\mathbf{L}(t))$ instead of $\mathbf{L}(0)$ in the standard swap rate volatility approximation can be used independently of any considerations involving time-dependent skews. In particular, the more accurate approximation can be applied directly in the statement of Proposition 14.4.3 when defining the value of $\lambda_{S}(t)$ in (14.33). The differences between the two formulas are small, but become noticeable for swaptions of longer-dated expiries and tenors.

### 15.2.2 Advanced Formula for Swap Rate Skew

In the approximate SDE (15.28) for the swap rate used for European option pricing, the parameter $b_{S}(t)$ controls the skew of the volatility smile. Define

$$
v_{n}(t)=\frac{\partial S(t)}{\partial L_{n}(t)}, \quad n=j, \ldots, k-1
$$

(compare to (14.31)), and

$$
v_{n, n^{\prime}}(t)=\frac{\partial^{2} S(t)}{\partial L_{n}(t) \partial L_{n^{\prime}}(t)}, \quad n, n^{\prime}=j, \ldots, k-1
$$

Proposition 15.2.3. The time-dependent swaption skew $b_{S}(t)$ is approximately given by

$$
b_{S}(t)=\sum_{i, n=j}^{k-1} \frac{r_{S, n}(t)}{r_{S, i}(t)} v_{\imath, n}(0) \xi_{i}(t)+\sum_{i=j}^{k-1} b_{i}(t) v_{i}(0) \xi_{\imath}(t)
$$

where

$$
\begin{gathered}
r_{S, i}(t)=\lambda_{i}(t)^{\top} \lambda_{S}(t, \mathbf{L}(0)), \quad r_{S}(t)=\left\|\lambda_{S}(t, \mathbf{L}(0))\right\|^{2} \\
\xi_{i}(t)=\frac{r_{S, i}(t) \int_{0}^{t} r_{S, 2}(u) d u}{r_{S}(t) \int_{0}^{t} r_{S}(u) d u}
\end{gathered}
$$

Proof. Recall from Proposition 15.2.1,

$$
\begin{aligned}
b_{S}(t) & =\left.\frac{1}{p_{S}(t)} \sum_{i=j}^{k-1} \frac{\partial\left\|\lambda_{S}(t, \mathbf{L}(t))\right\|}{\partial L_{\imath}(t)}\right|_{\mathbf{L}(t)=\mathrm{E}^{\wedge}(\mathbf{L}(t))} \frac{\int_{0}^{t} r_{S, i}(u) d u}{\int_{0}^{t} r_{S}(u) d u} \\
& =\left.\sum_{i=j}^{k-1} \frac{\partial \ln \left\|\lambda_{S}(t, \mathbf{L}(t))\right\|}{\partial L_{i}(t)}\right|_{\mathbf{L}(t)=\mathrm{E}^{\wedge}(\mathbf{L}(t))} \frac{\int_{0}^{t} r_{S, i}(u) d u}{\int_{0}^{t} r_{S}(u) d u}
\end{aligned}
$$

We have

$$
\begin{aligned}
& \frac{\partial \ln \left\|\lambda_{S}(t, \mathbf{L}(t))\right\|}{\partial L_{i}(t)}=\frac{1}{\left\|\lambda_{S}(t, \mathbf{L}(t))\right\|^{2}} \\
& \times \sum_{n, n^{\prime}=j}^{k-1}\left(\lambda_{n}(t)^{\top} \lambda_{n^{\prime}}(t)\right)\left(v_{n}(t) \varphi_{n}\left(t, L_{n}(t)\right)\right) \frac{\partial}{\partial L_{i}(t)}\left(v_{n^{\prime}}(t) \varphi_{n^{\prime}}\left(t, L_{n^{\prime}}(t)\right)\right)
\end{aligned}
$$

While using $\mathrm{E}^{A}(\mathbf{L}(t))$ instead of $\mathbf{L}(0)$ in the calculations of the swaption volatility (see the previous section) results in noticeable improvements in the approximation quality, this turns out to not be the case for approximations to the skew. Furthermore, as $\varphi_{n}\left(t, L_{n}(0)\right)=1$ for any $n$, the formulas resulting
from evaluating the expression for $b_{S}(t)$ above at $L_{n}(0)$ are compact and convenient, which in practice will justify any slight deterioration in precision. With this in mind, we note that

$$
\left.\frac{\partial}{\partial L_{i}(t)}\left(v_{n^{\prime}}(t) \varphi_{n^{\prime}}\left(t, L_{n^{\prime}}(t)\right)\right)\right|_{\mathbf{L}(t)=\mathbf{L}(0)}=v_{i, n^{\prime}}(0)+1_{\left\{i=n^{\prime}\right\}} v_{n^{\prime}}(0) b_{n^{\prime}}(t)
$$

(recall (15.24) for the definition of $b_{n}$ 's). Hence,

$$
\begin{aligned}
\left.\frac{\partial \ln \left\|\lambda_{S}(t, \mathbf{L}(t))\right\|}{\partial L_{i}(t)}\right|_{\mathbf{L}(t)=\mathbf{L}(0)} & =\frac{1}{r_{S}(t)} \sum_{n, n^{\prime}=j}^{k-1}\left(\lambda_{n}(t)^{\top} \lambda_{n^{\prime}}(t)\right) v_{n}(0) v_{i, n^{\prime}}(0) \\
& +\frac{1}{r_{S}(t)} b_{i}(t) v_{i}(0) \sum_{n=j}^{k-1}\left(\lambda_{n}(t)^{\top} \lambda_{i}(t)\right) v_{n}(0)
\end{aligned}
$$

We recognize

$$
\sum_{n=j}^{k-1}\left(\lambda_{n}(t)^{\top} \lambda_{n^{\prime}}(t)\right) v_{n}(0)=\lambda_{n^{\prime}}(t)^{\top} \lambda_{S}(t, \mathbf{L}(0))=r_{S, n^{\prime}}(t)
$$

so that

$$
\begin{aligned}
\left.\frac{\partial \ln \left\|\lambda_{S}(t, \mathbf{L}(t))\right\|}{\partial L_{i}(t)}\right|_{\mathbf{L}(t)=\mathbf{L}(0)} & \\
& =\left(\sum_{n^{\prime}=j}^{k-1} v_{i, n^{\prime}}(0) \frac{r_{S, n^{\prime}}(t)}{r_{S}(t)}\right)+b_{i}(t) v_{i}(0) \frac{r_{S, i}(t)}{r_{S}(t)} .
\end{aligned}
$$

The result follows.
Remark 15.2.4. The swaption skew consists of two parts, one that involves second-order derivatives of the swap rate with respect to Libor rates and captures overall convexity of a swap rate with respect to Libor rates, and the other being a weighted average of Libor skews.

Remark 15.2.5. Even with Libor rates sharing a common skew, i.e. in the "classic" LM model of Sections 14.2.4, 14.2.5, the swaption skew is not exactly equal to the (shared) Libor skews. If $b_{i}(t)=b$ for all $i$ and $t$, then

$$
b_{S}(t)=\sum_{i, n=j}^{k-1} \frac{r_{S, n}(t)}{r_{S, i}(t)} v_{i, n}(0) \xi_{i}(t)+b \sum_{i=j}^{k-1} v_{i}(0) \xi_{i}(t)
$$

Even if the convexity term is ignored, we have

$$
b_{S}(t)=b \sum_{i=j}^{k-1} v_{i}(0) \xi_{i}(t) \neq b
$$

### 15.2.3 Skew and Smile Calibration in LM Models

With pricing issues out of the way, let us now turn our attention to calibration, and see what modifications to the calibration algorithm of Section 14.5 are required by the more general model specification (15.22)-(15.23). We assume a typical set of swaption volatilities is given, and we have available a collection of market-implied volatility smiles across strikes for swaptions of different expiries and tenors. The model (15.22)-(15.23) has enough parameter flexibility to match

- At-the-money volatilities of all European swaptions, using the volatility structure $\left\|\lambda_{n}(t)\right\|$.
- Slopes of volatility smiles (skews) of all European swaptions, using the skew structure $b_{n}(t)$ (see (15.24)).
- Curvatures of smiles for swaptions of different expiries, for a given tenor, using term structure of volatilities of variance $\eta(t)$.

Of course, all these parameters are in addition to the correlation structure of Libor rates, as in the standard LM model case. Note that, for the last point, technically there is no flexibility in the model to change the volatility smile curvature for swaptions of the same expiry but different tenors. This is not really a serious limitation as the curvature of the smile tends to be constant across tenors for a given expiry.

Assume, as in Section 14.5, that we have chosen calibration targets that include $N_{S}$ swaptions, $V_{\text {swaption, } 1}, V_{\text {swaption, } 2}, \ldots, V_{\text {swaption, } N_{S}}$, and let us ignore caps; the considerations below extend trivially to cover them as well. Unlike previously, however, let us assume that each target includes not one swaption, but a collection of them of different strikes. Hence, we redefine the calibration as the goal to match volatility smiles of $N_{S}$ swaptions.

Having chosen a (vanilla) stochastic volatility model for European swaptions such as (8.3)-(8.4), the target volatility smiles can be summarized by a collection of market-implied SV parameters, namely volatilities $\widehat{\lambda}_{S_{i}}$, skews $\widehat{b}_{S_{i}}$ and volatilities of variance $\widehat{\eta}_{S_{i}}$ for $i=1, \ldots, N_{S}$ (a common mean reversion of volatility parameter is assumed). Recall that in Section 14.5.2 we denoted by $G$ a grid of the Libor volatilities $\left\|\lambda_{n}(t)\right\|$. To be able to generalize, redefine $G_{\lambda}=G$ and, in the same spirit, define $G_{b}$ to be the grid of Libor skews $b_{n}(t)$, and $G_{\eta}$ (a vector) to be the discretized term structure of volatilities of variance $\eta(t)$. The formulas from the previous section allow us to compute term volatilities, skews and volatilities of variance from the model. We denote them by

$$
\bar{\lambda}_{S_{i}}\left(G_{\lambda}, G_{b}, G_{\eta}\right), \quad \bar{b}_{S_{i}}\left(G_{\lambda}, G_{b}, G_{\eta}\right), \quad \bar{\eta}_{S_{i}}\left(G_{\lambda}, G_{b}, G_{\eta}\right) .
$$

One can incorporate the skew and smile calibration in the algorithm of Section 14.5 by adding extra terms to the calibration norm, replacing (14.54) with

$$
\begin{align*}
\mathcal{I}\left(G_{\lambda}, G_{b}, G_{\eta}\right)= & \frac{w_{S, \lambda}}{N_{S}} \sum_{i=1}^{N_{S}}\left(\bar{\lambda}_{S_{i}}\left(G_{\lambda}, G_{b}, G_{\eta}\right)-\widehat{\lambda}_{S_{i}}\right)^{2}  \tag{15.32}\\
& +\frac{w_{S, b}}{N_{S}} \sum_{i=1}^{N_{S}}\left(\bar{b}_{S_{i}}\left(G_{\lambda}, G_{b}, G_{\eta}\right)-\widehat{b}_{S_{i}}\right)^{2} \\
& +\frac{w_{S, \eta}}{N_{S}} \sum_{i=1}^{N_{S}}\left(\bar{\eta}_{S_{i}}\left(G_{\lambda}, G_{b}, G_{\eta}\right)-\widehat{\eta}_{S_{i}}\right)^{2} \\
& +\ldots
\end{align*}
$$
^eq15-32

with dots denoting various regularization terms for $G_{\lambda}, G_{b}$, and $G_{\eta}$. This, however, is not necessarily the best approach, as it increases the number of degrees of freedom in the non-linear optimization problem quite substantially, thus potentially significantly reducing the speed at which it could be solved. It is much better to take advantage of the structure of the problem and solve for volatilities, skews and volatilities of variance separately rather than all at the same time.

To motivate the method, we note that the impact of changes in the Libor skews $G_{b}$ on term swaption volatilities, or their volatilities of variance, is rather small. Likewise, changes in Libor volatilities $G_{\lambda}$ have only a small impact on term swaption skews, and so on. This near-orthogonality allows us to solve for various parameters sequentially. To facilitate this approach, we define three norms

$$
\begin{aligned}
& \mathcal{I}_{\lambda}\left(G_{\lambda}, G_{b}, G_{\eta}\right)=\frac{w_{S, \lambda}}{N_{S}} \sum_{i=1}^{N_{S}}\left(\bar{\lambda}_{S_{i}}\left(G_{\lambda}, G_{b}, G_{\eta}\right)-\widehat{\lambda}_{S_{i}}\right)^{2}+\ldots, \\
& \mathcal{I}_{b}\left(G_{\lambda}, G_{b}, G_{\eta}\right)=\frac{w_{S, b}}{N_{S}} \sum_{i=1}^{N_{S}}\left(\bar{b}_{S_{i}}\left(G_{\lambda}, G_{b}, G_{\eta}\right)-\widehat{b}_{S_{i}}\right)^{2}+\ldots, \\
& \mathcal{I}_{\eta}\left(G_{\lambda}, G_{b}, G_{\eta}\right)=\frac{w_{S, \eta}}{N_{S}} \sum_{i=1}^{N_{S}}\left(\bar{\eta}_{S_{i}}\left(G_{\lambda}, G_{b}, G_{\eta}\right)-\widehat{\eta}_{S_{i}}\right)^{2}+\ldots,
\end{aligned}
$$

and modify the algorithm from Section 14.5.7 as such.

1. First, make a guess for $G_{b}$ and $G_{\eta}$, denoted by $G_{b}^{0}$ and $G_{\eta}^{0}$. The guesses could be quite approximate. For example the grid $G_{b}^{0}$ could be set to the same value, an average of swaption term skews $N_{S}^{-1} \sum_{i=1}^{N_{S}} \hat{b}_{S_{i}}$, and the same for $G_{\eta}^{0}$.
2. Perform steps 1-5 from Section 14.5.7 with $\mathcal{I}=\mathcal{I}_{\lambda}\left(G_{\lambda}, G_{b}^{0}, G_{\eta}^{0}\right)$ until $G_{\lambda}^{1}$, the solution, is found. Note that we keep the skew and volatility of variance grids constant throughout
3. Minimize $\mathcal{I}_{b}\left(G_{\lambda}^{1}, G_{b}, G_{\eta}^{0}\right)$ by iterating over $G_{b}$; denote the solution by $G_{b}^{1}$.
4. Minimize $\mathcal{I}_{\eta}\left(G_{\lambda}^{1}, G_{b}^{1}, G_{\eta}\right)$ by iterating over $G_{\eta}$; denote the solution by $G_{\eta}^{1}$.

Typically, the triple of parameters ( $G_{\lambda}^{1}, G_{b}^{1}, G_{\eta}^{1}$ ) provides a good overall solution. If a better fit is desired, the steps could be iterated, starting with ( $G_{b}^{1}, G_{\eta}^{1}$ ) on Step 1. If the number of iterations is more than 1, it could be beneficial to stop after Step 2 (of the second or subsequent iterations) in order to have the best possible fit to the volatilities, usually the most important target.

### 15.3 Near-Markov LM Models

The LM model, even if driven by a single Brownian motion, is not Markovian in a low number of state variables. Rather, the state vector comprises all Libor rates and whatever state variables are needed for the stochastic volatility, if present. This is readily seen from the expression for the drift of each Libor rate (see e.g. (14.7)) which involves multiple other rates.

This state of affairs, however, has not stopped some researchers from attempting to approximate an LM model with a low-dimensional Markovian one. These attempts mostly involve i) restricting the volatility structure to a "separable" form, similar to that used to specify Markovian models in Chapters 12 and 13; and ii) removing path dependence in the Libor drifts by either "freezing" them at time 0 values or by employing various tricks not unlike those discussed previously in the context of drift estimation for large time steps in Monte Carlo.

The practical usefulness of such approximations for pricing and risk managing of derivatives is limited. The imposed restrictions on the volatility structure remove many of the main LM advantages in terms of calibration flexibility, and the necessity of approximating the Libor rate drifts makes the model arbitrageable and problematic to use for anything other that short-dated derivatives. Fundamentally, there is also something misplaced about trying to force LM models into a low-dimensional Markovian setting: if such a setting is desired, there are really no appealing reasons to use an LM model in the first place. Instead, one should from the outset pick one of the many perfectly good low-dimensional Markovian models we have covered earlier in the book. The models in Chapter 13 are particularly appealing, we think.

That said, a low-dimensional Markov approximation to an LM model can, however, still find use as a type of a model-based control variate. We discuss this application in more detail in Chapter 25; this chapter also covers the mechanics of the various steps involved in creating near-Markov LM models.

### 15.4 Swap Market Models

In a nutshell, the LM modeling principle revolves around using models for simple forward rates (Libor) that become tractable in properly selected
martingale measures. Longer-dated swap rates can be constructed iteratively by, in effect, adding up the individual forward rates that constitute the LM model primitives. As we have seen in Section 14.4.2, the swap rates constructed in this manner rarely have tractable dynamics, and swaption pricing formulas will nearly always involve approximations and/or numerical methods.

An alternative modeling approach due to Jamshidian [1997] turns the LM philosophy on its head by using forward swap rates as the fundamental model primitives, and constructing individual Libor rates as differences of such swap rates. This is similar to the idea behind swap Markov-functional model construction that we briefly discussed in Appendix 11.A. 3 of Chapter 11.

By specifying tractable dynamics for forward swap rates, swaption pricing formulas now often can be done exactly, whereas cap pricing must be done by approximation. While this so-called swap market model framework has some proponents (e.g. Galluccio and Hunter [2003], Galluccio et al. [2005]), it is fair to say that the approach has attracted much less interest from practitioners and academics than has the LM model. Our treatment of swap market models shall therefore be brief.

Given a tenor structure $0=T_{0}<T_{1}<\ldots<T_{N}$, our market primitives are now the swap rates

$$
S_{j}(t)=S_{j, N-j}(t), \quad 1 \leq j<N,
$$

where the notation is similar to that used in Section 14.4.2:

$$
S_{j}(t)=\frac{P\left(t, T_{j}\right)-P\left(t, T_{N}\right)}{A_{j}(t)}, \quad A_{j}(t)=A_{j, N-j}(t)=\sum_{n=j}^{N-1} P\left(t, T_{n+1}\right) \tau_{n}
$$

We emphasize that all forward swaps here have identical terminal maturity, namely $T_{N}$. We also emphasize that the information content of the primitives of a swap market model is identical to that of an LM model, as the forward Libor rates curve can be uniquely constructed from knowledge of $S_{i}(t)$, for all $i \in[q(t), N-1]$. To show this, write for any $j$

$$
\begin{aligned}
L_{j-1}(t) & =\frac{S_{j-1}(t) A_{j-1}(t)-S_{j}(t) A_{j}(t)}{P\left(t, T_{j}\right) \tau_{j-1}} \\
& =S_{j-1}(t) \frac{A_{j-1}(t)}{P\left(t, T_{j}\right) \tau_{j-1}}-S_{j}(t) \frac{A_{j-1}(t)-P\left(t, T_{j}\right) \tau_{j-1}}{P\left(t, T_{j}\right) \tau_{j-1}} \\
& =S_{j}(t)+\left(S_{j-1}(t)-S_{j}(t)\right) \frac{A_{j-1}(t)}{P\left(t, T_{j}\right) \tau_{j-1}}
\end{aligned}
$$

As

$$
\frac{P\left(t, T_{n+1}\right)}{P\left(t, T_{j}\right)}=\frac{P\left(t, T_{j+1}\right)}{P\left(t, T_{j}\right)} \frac{P\left(t, T_{j+2}\right)}{P\left(t, T_{j+1}\right)} \ldots \frac{P\left(t, T_{n+1}\right)}{P\left(t, T_{n}\right)}=\prod_{k=j}^{n}\left(1+L_{k}(t) \tau_{k}\right)^{-1},
$$

it follows from the definition of $A_{j-1}(t)$ that

$$
\begin{equation*}
L_{j-1}(t)=S_{j}(t)+\left(S_{j-1}(t)-S_{j}(t)\right)\left(1+\sum_{n=j}^{N-1} \prod_{k=j}^{n}\left(1+L_{k}(t) \tau_{k}\right)^{-1} \frac{\tau_{n}}{\tau_{j-1}}\right) . \tag{15.33}
\end{equation*}
$$
^eq15-33

Starting from $L_{N-1}(t)=S_{N-1}(t)$, equation (15.33) gives us an iterative formula to construct the Libor forward curve at time $t$, as claimed.

In the swap market model framework, the modeler specifies dynamics on the swap forward rates in their respective annuity measures. Let $W^{A_{j}}(t)$ be an $m$-dimensional Brownian motion in the annuity measure induced by $A_{j}(t)$. As shown earlier, $S_{j}(t)$ is a martingale in this measure, and we can write

$$
\begin{equation*}
d S_{j}(t)=\sigma_{S_{j}}(t)^{\top} d W^{A_{j}}(t), \quad j=q(t), \ldots, N-1, \tag{15.34}
\end{equation*}
$$
^eq15-34

for some adapted volatility function $\sigma_{S_{j}}(t)$ specific to the $j$-th swap rate. As for an LM model, we could now proceed to specialize the model by using DVF or SV type specifications for $\sigma_{S_{1}}(t)$, generally keeping an eye out for models that allow construction of easy-to-compute expressions for payer swaption prices

$$
V_{\text {swaption }, j}(t)=A_{j}(t) \mathrm{E}_{t}^{A_{j}}\left(\left(S_{j}\left(T_{j}\right)-c\right)^{+}\right) .
$$

As we mentioned earlier, no exact pricing formulas for caplets will normally exist (as should be obvious from the complicated form of (15.33)), a situation that also holds for any swaption that does not involve a swap maturing at time $T_{N}$. For these instruments, approximative formulas must be devised if a quick calibration of the model is desired. See Galluccio et al. [2005] for some details on this.

In the context of the LM models, Section 14.2.2 derived relations between the different forward martingale measures, allowing us to describe the dynamics of all forward rates in a single measure, as required in Monte Carlo simulations, say. For the swap market models, starting with the specification (15.34) we can derive similar relations between the different annuity measures, ultimately allowing for simulation of all swap rates $S_{1}, S_{2}, \ldots, S_{N-1}$ in a common measure. For details, the reader is referred to Jamshidian [1997] and Section 14.4 of Musiela and Rutkowski [1997].

### 15.5 Evolving Separate Discount and Forward Rate Curves

A single yield curve for discounting and calculating Libor rates is not always compatible with no-arbitrage constraints of cross-currency markets, nor is it
particularly realistic in stressed marked conditions, such as those experienced during the subprime crisis of 2007-2009. As explained in Section 6.5.2.2, separating the discounting curve from the forward, or index, curve will ensure that linear instruments (i.e. swaps and bonds) are correctly priced at time 0 . In this section we consider how to incorporate the idea of curve separation into a dynamic model of interest rates. While we use an LM setting for some of the material in this section, the basic ideas are generic and can be applied to virtually all models in this book.

### 15.5.1 Basic Ideas

Suppose two yield curves are given at time 0 , the discount curve $P(0, T)$, $T \geq 0$, and the index curve $\widetilde{P}(0, T), T \geq 0$ (note the change of notation from $P^{(L)}$ in Section 6.5.2.2 to $\widetilde{P}$ for convenience). The index curve corresponds to a particular Libor tenor $\tau$, and is defined through the requirement that forward Libor rates of tenor $\tau$ must equal the conditional expected values of the future spot Libor rates,

$$
\begin{equation*}
\widetilde{L}(t, T, T+\tau)=\mathrm{E}_{t}^{T+\tau}(L(T, T, T+\tau)), \tag{15.35}
\end{equation*}
$$
^eq15-35

where $\mathrm{E}_{t}^{T+\tau}$ denotes expectation in the $(T+\tau)$-forward measure and, by definition,

$$
\begin{equation*}
\widetilde{L}(t, T, T+\tau)=\frac{1}{\tau}\left(\frac{1}{\widetilde{P}(t, T, T+\tau)}-1\right) \tag{15.36}
\end{equation*}
$$
^eq15-36

For emphasis, we have added a tilde to the regular Libor rate notation to highlight the link to the index curve. Note that we take (15.35) to hold for arbitrary values of $T \geq t$.

It is important to point out that using different Libor tenors $\tau$ in the equation above would lead to different models of two-curve evolution, and our first choice shall focus on what tenor we actually want to use. In the next section we will be looking at extending the HJM instantaneous forward rate formalism, so our choice will be influenced by that fact. Later, in Section 15.5.3, we will consider Libor market models that, not surprisingly, lead to a somewhat different choice.

As we initially focus on adding a second curve evolution in the HJM setting, it is convenient to consider a specific choice of $\tau=0$ that corresponds to instantaneous forward rates. In particular, we note that for $\tau \rightarrow 0$, $\widetilde{L}(t, T, T+\tau)$ converges to $\widetilde{f}(t, T)$, so according to (15.35), for the next section we shall require that

$$
\begin{equation*}
\widetilde{f}(t, T)=\mathrm{E}_{t}^{T}(\widetilde{f}(T, T)) \tag{15.37}
\end{equation*}
$$
^eq15-37

As discussed in Chapter 6 (see (6.44)-(6.45)), it is convenient to represent the index curve through an additive spread $\varepsilon$ in continuously compounded forward rates. Specifically, at time 0 we write

$$
\begin{equation*}
\widetilde{P}(0, T)=P(0, T) e^{\int_{0}^{T} \varepsilon(0, u) d u}, \quad T \geq 0 \tag{15.38}
\end{equation*}
$$
^eq15-38

We already hinted in Section 6.5.2.4 that a standard way of including the spread in a dynamic model is to assume that it evolves deterministically. We will ultimately make the same recommendation here, but it is still instructive to develop a generic dynamic model of two-curve evolution first. In developing such a model, we have several possible choices for the model primitives. For example, we could impose dynamics on the index and discount curves in separation. Alternatively, we could impose dynamics on the spread and just one of the two curves, deriving the remaining curve from these two primitives. It is often desirable to have direct control of the spread process - e.g. to ensure that its dynamics and range are in line with historical observations - so we adopt the latter approach here. As for the choice of which yield curve to use as a direct model ingredient, this is largely a matter of taste and convenience. We shall demonstrate both approaches below, using the discount curve as the primary curve in the HJM model setting of Section 15.5.2; and the index curve as the primary curve in the LM model setting of Section 15.5.3.

### 15.5.2 HJM Extension

The framework for discount bonds is unchanged by the presence of an index curve, and we can therefore start with the standard HJM dynamics for the discount factors (see Section 4.4.1) in the risk-neutral measure Q,

$$
d P(t, T) / P(t, T)=r(t) d t-\sigma_{P}(t, T)^{\top} d W(t)
$$

where $W(t)$ is a $d$-dimensional Brownian motion and $\sigma_{P}(t, T)$ is a $d$ dimensional discount bond volatility function. Importantly, the $T$-forward measure $\mathrm{Q}^{T}$ is determined solely by the dynamics of the discount curve because the numeraire, the zero-coupon discount bond $P(t, T)$, is just a particular point on the time $t$ discount curve. In particular, the $T$-forward measure is defined by the familiar requirement that

$$
d W^{T}(t)=d W(t)+\sigma_{P}(t, T) d t
$$

be a driftless Brownian motion.
As for the instantaneous forward rates $f(t, T)$, the standard HJM result (see Lemma 4.4.1) shows that their dynamics are given by

$$
\begin{equation*}
d f(t, T)=\sigma_{f}(t, T)^{\top} \sigma_{P}(t, T) d t+\sigma_{f}(t, T)^{\top} d W(t) \tag{15.39}
\end{equation*}
$$
^eq15-39

Let $\varepsilon(t, u)$ be the forward rate spread defined by extending (15.38) to arbitrary $t$,

$$
\begin{equation*}
\widetilde{P}(t, T)=P(t, T) e^{\int_{t}^{T} \varepsilon(t, u) d u}, \quad T \geq t \tag{15.40}
\end{equation*}
$$
^eq15-40

Treating $\varepsilon$ as a forward rate spread is justified by considering $\widetilde{f}(t, T)$, the instantaneous forward rates calculated from the index curve $\tilde{f}(t, T)= -\partial \ln \widetilde{P}(t, T) / \partial T$, and observing that

$$
\begin{equation*}
\widetilde{f}(t, T)=f(t, T)-\varepsilon(t, T), \quad 0 \leq t \leq T . \tag{15.41}
\end{equation*}
$$
^eq15-41

Let us endow $\varepsilon(t, T)$ with the following dynamics ${ }^{4}$,

$$
\begin{equation*}
d \varepsilon(t, T)=\alpha_{\varepsilon}(t, T) d t+\sigma_{\varepsilon}(t, T)^{\top} d W(t) \tag{15.42}
\end{equation*}
$$
^eq15-42

for some, yet unspecified, adapted processes $\alpha_{\varepsilon}(t, T)$ and $\sigma_{\varepsilon}(t, T)$, the latter $d$-dimensional. Not surprisingly, the drift term in (15.42) cannot be set arbitrarily.

Proposition 15.5.1. To satisfy the martingale restrictions of (15.37), the drift of $\varepsilon(t, T)$ in (15.42) must obey

$$
\begin{equation*}
\alpha_{\varepsilon}(t, T)=\sigma_{\varepsilon}(t, T)^{\top} \sigma_{P}(t, T) \tag{15.43}
\end{equation*}
$$
^eq15-43

Proof. According to (15.37), $\widetilde{f}(t, T)$ must be a martingale in the $T$-forward measure. Its dynamics in that measure follow from (15.39)-(15.42),

$$
\begin{aligned}
& d \tilde{f}(t, T)=\left(\sigma_{f}(t, T)^{\top} \sigma_{P}(t, T)-\alpha_{\varepsilon}(t, T)\right) d t \\
& \quad+\left(\sigma_{f}(t, T)-\sigma_{\varepsilon}(t, T)\right)^{\top} d W^{T}(t)-\left(\sigma_{f}(t, T)-\sigma_{\varepsilon}(t, T)\right)^{\top} \sigma_{P}(t, T) d t
\end{aligned}
$$

and the result in the proposition follows by setting the $d t$ term to zero.
The SDEs for various quantities related to the index curve can now be derived from the discount and spread parameters.

Proposition 15.5.2 Define

$$
\tilde{\alpha}_{f}(t, T) \triangleq \tilde{\sigma}_{f}(t, T)^{\top} \sigma_{P}(t, T), \quad \tilde{\sigma}_{f}(t, T) \triangleq \sigma_{f}(t, T)-\sigma_{\varepsilon}(t, T)
$$

The dynamics of $\widetilde{f}(t, T)$ in the risk-neutral measure are then

$$
d \widetilde{f}(t, T)=\widetilde{\alpha}_{f}(t, T) d t+\widetilde{\sigma}_{f}(t, T)^{\top} d W(t)
$$

Further, set $\widetilde{r}(t)=\widetilde{f}(t, t)$ and[^72]

$$
\begin{aligned}
& \widetilde{\alpha}_{P}(t, T) \triangleq \int_{t}^{T} \widetilde{\sigma}_{f}(t, u)^{\top}\left(\widetilde{\sigma}_{P}(t, u)-\sigma_{P}(t, u)\right) d u \\
& \widetilde{\sigma}_{P}(t, T) \triangleq \int_{t}^{T} \widetilde{\sigma}_{f}(t, u) d u
\end{aligned}
$$

Then

$$
d \widetilde{P}(t, T) / \widetilde{P}(t, T)=\widetilde{r}(t) d t+\widetilde{\alpha}_{P}(t, T) d t-\widetilde{\sigma}_{P}(t, T)^{\top} d W(t)
$$

Proof. The result for $d \widetilde{f}(t, T)$ follows directly from (15.39) and (15.41)(15.43). From the equation

$$
\widetilde{P}(t, T)=e^{-\int_{t}^{T} \tilde{f}(t, u) d u}
$$

it follows that $Y(t, T)=\ln (\widetilde{P}(t, T))$ satisfies

$$
\begin{aligned}
d Y(t, T) & =\widetilde{f}(t, t) d t-\int_{t}^{T} \widetilde{\alpha}_{f}(t, u) d u d t-\widetilde{\sigma}_{P}(t, T)^{\top} d W(t) \\
& =\widetilde{r}(t) d t-\widetilde{\sigma}_{P}(t, T)^{\top} d W(t)-\int_{t}^{T} \widetilde{\sigma}_{f}(t, u)^{\top} \sigma_{P}(t, u) d u d t
\end{aligned}
$$

An application of Ito's lemma (to $e^{Y}$ ) then shows that

$$
d \widetilde{P}(t, T) / \widetilde{P}(t, T)=\widetilde{r}(t) d t+\widetilde{\alpha}_{P}(t, T) d t-\widetilde{\sigma}_{P}(t, T)^{\top} d W(t)
$$

where

$$
\widetilde{\alpha}_{P}(t, T)=-\int_{t}^{T} \widetilde{\sigma}_{f}(t, u)^{\top} \sigma_{P}(t, u) d u+\frac{1}{2} \widetilde{\sigma}_{P}(t, T)^{\top} \widetilde{\sigma}_{P}(t, T)
$$

Integration by parts shows that

$$
\tilde{\sigma}_{P}(t, T)^{\top} \tilde{\sigma}_{P}(t, T)=2 \int_{t}^{T} \tilde{\sigma}_{f}(t, u)^{\top} \tilde{\sigma}_{P}(t, u) d u
$$

and the result follows.
Remark 15.5.3. Note the presence of an "extra" drift term $\widetilde{\alpha}_{P}(t, T)$ in the dynamics for the pseudo-bonds $\widetilde{P}(t, T)$; sometimes, in an analogy with foreign exchange markets, this term is called a quanto correction.

For future use, let us define

$$
\begin{equation*}
Z(t, T) \triangleq e^{-\int_{1}^{T} \varepsilon(t, u) d u} \tag{15.44}
\end{equation*}
$$
^eq15-44

so that

$$
\begin{equation*}
\widetilde{P}(t, T)=P(t, T) / Z(t, T), \quad 0 \leq t \leq T \tag{15.45}
\end{equation*}
$$
^eq15-45

The drift and diffusion terms in the process for $Z(t)$,

$$
\begin{equation*}
d Z(t, T) / Z(t, T)=(r(t)-\widetilde{r}(t)) d t+\alpha_{Z}(t, T) d t-\sigma_{Z}(t, T)^{\top} d W(t) \tag{15.46}
\end{equation*}
$$
^eq15-46

are given by

$$
\begin{equation*}
\alpha_{Z}(t, T)=-\widetilde{\alpha}_{P}(t, T)-\sigma_{Z}(t, T)^{\top} \widetilde{\sigma}_{P}(t, T), \quad \sigma_{Z}(t, T)=\int_{t}^{T} \sigma_{\varepsilon}(t, u) d u \tag{15.47}
\end{equation*}
$$
^eq15-47

Options linked directly to the spread between the index and discount curves are rarely, if ever, traded, so a high level of sophistication in basis spread dynamics is seldom required. A simple one-factor Gaussian model for the spread, say, is more than sufficient for most applications. Nevertheless, richer dynamics are possible. Recall that if the forward rate volatility function is of separable form,

$$
\sigma_{f}(t, T)=g(t) h(T)
$$

where $g(\cdot)$ is a $d \times d$ matrix-valued process and $h(\cdot)$ is a $d$-dimensional deterministic vector-valued function (see e.g. (4.44), (12.2) or (13.70)), then the discount curve admits a $d$-dimensional Markovian representation, see, for example, Proposition 13.3.1. If, in addition, the forward rate spread volatility function $\sigma_{\varepsilon}(t, T)$ is also of separable form with the same $g(t)$ but a different $h_{\varepsilon}(T)$,

$$
\sigma_{\varepsilon}(t, T)=g(t) h_{\varepsilon}(T)
$$

then the forward rate spread curve $\varepsilon(t, T)$ and, by extension, the index curve also admit Markovian representations with the same Markovian state variables. This fact opens up a possibility of building efficient Markovian models of joint evolution of the discount and the index curve with nondeterministic spread. We leave this line of inquiry for the reader to explore, and instead move on to the changes required to extend the LM model framework to support two-curve dynamics.

### 15.5.3 Applications to LM Models

We now return to the LM framework, and continue with the notations of Section 14.2. In particular, we assume that a tenor structure (14.1) has been specified. Clearly, in the LM setup, it is natural to assume that (15.36) is satisfied for the Libor tenor used in the definition of the LM model. This would correspond to extending the LM model to drive both the discounting curve and the index curve that corresponds to the Libor tenor used in LM model definition. We remind the reader that if multiple index curve dynamics are required, each will have to be driven by its own LM model (of a given Libor tenor).

As discussed earlier, we here choose as our primary model variables (in addition to the spread) the set of "regular" Libor rates $\widetilde{L}_{n}(t) \triangleq \widetilde{L}\left(t, T_{n}, T_{n+1}\right)$,
$n=0, \ldots, N-1$, as defined off the index curve in (15.36). Apart from giving ourselves a chance to present a slightly different twist on the material of the previous section, the choice of the index Libors as building blocks allows us to use the volatilities of $\left\{\widetilde{L}_{n}(t)\right\}_{n=0}^{N-1}$ that are more directly related to the market-observable volatility information, i.e. the volatilities of swap rates of various maturities and tenors. Conveniently, for each $n, \widetilde{L}_{n}(t)$ is a martingale in $T_{n+1}$-forward measure by construction, see (15.35). Therefore, in close analogy to (14.5), we can define the dynamics by

$$
\begin{equation*}
d \widetilde{L}_{n}(t)=\widetilde{\sigma}_{n}(t)^{\top} d W^{n+1}(t) \tag{15.48}
\end{equation*}
$$
^eq15-48

where $W^{n+1}(t)$ is an $m$-dimensional Brownian motion in $\mathrm{Q}^{T_{n+1}}$ and $\widetilde{\sigma}_{n}(t)$ is an $m$-dimensional adapted process.

The discrete money market account $B(t)$, see (14.8), induces a useful common measure for all Libor rates in the single-curve case. The extension to the two-curve framework is straightforward: we apply (14.8) literally, i.e. use simply compounded forward rates computed off the discount curve in the definition of $B(t)$. This corresponds exactly to the actual trading strategy of re-investing into deposits of a given tenor, with the value of the strategy solely determined by the discount curve:

$$
\begin{equation*}
B(t)=P\left(t, T_{q(t)}\right) \prod_{n=0}^{q(t)-1} \frac{1}{P\left(T_{n}, T_{n+1}\right)} . \tag{15.49}
\end{equation*}
$$
^eq15-49

We define the spot measure $Q^{B}$ accordingly. Interestingly, a measure change from any $T$-forward measure to the spot measure is determined by the dynamics of discount factors $P(t, T)$ only, and we see that if $\widetilde{L}_{n}$ 's satisfy (15.48), then in the spot measure $\mathrm{Q}^{B}$, the process for $\widetilde{L}_{n}$ is given by

$$
\begin{equation*}
d \widetilde{L}_{n}(t)=\widetilde{\sigma}_{n}(t)^{\top}\left(\sum_{j=q(t)}^{n} \frac{\tau_{j} \sigma_{j}(t)}{1+\tau_{j} L_{j}(t)} d t+d W^{B}(t)\right) \tag{15.50}
\end{equation*}
$$
^eq15-50

where $W^{B}(t)$ is an $m$-dimensional Brownian motion in measure $\mathrm{Q}^{B}$. Here $L_{j}(t)$ are simply compounded forward rates calculated off the discount curve, and $\sigma_{j}(t)$ are their volatilities. It would, however, be more convenient to have the drift in terms of the Libor rate model primitives $\widetilde{L}_{j}(t)$ and the process for the spread.

Before deriving the result for the drifts, let us decide what variables to use to represent the spread evolution. Using the instantaneous forward spread curve $\varepsilon(t, T)$ is not convenient in the LM framework, and we settle on forward bond ratios $Z^{n}(t), n=0, \ldots, N-1$, as state variables, defined by

$$
\begin{equation*}
Z^{n}(t) \triangleq Z\left(t, T_{n+1}\right) / Z\left(t, T_{n}\right) \tag{15.51}
\end{equation*}
$$
^eq15-51

with the $Z(t, T)$ 's given in (15.44). This choice of state variables keeps the LM framework in line with (but not exactly equivalent to!) the HJM extension
from Section 15.5.2, although other choices are possible; we comment on that at the end of the section.

Let us assume that forward bond ratios follow the dynamics

$$
\begin{equation*}
d Z^{n}(t) / Z^{n}(t)=O(d t)-\sigma_{Z^{n}}(t)^{\top} d W^{B}(t) \tag{15.52}
\end{equation*}
$$
^eq15-52

in the spot measure, with $\sigma_{Z^{n}}(t)$ 's being their volatilities. Note that $\sigma_{Z^{\prime \prime}}(t)= \sigma_{Z}\left(t, T_{n+1}\right)-\sigma_{Z}\left(t, T_{n}\right)$, where $\sigma_{Z}(t, T)$ are defined in (15.46). The drifts of forward bond ratios, not surprisingly, are not free parameters. We have the following result on the dynamics of Libor rates and forward bond ratios.

Proposition 15.5.4. If the Libor rates $\widetilde{L}_{n}(t)$ satisfy (15.48), and forward bond ratios $Z^{n}(t)$ follow the dynamics of (15.52), then in the spot measure $\mathrm{Q}^{B}$, the process for the state variables of the model is given by

$$
\begin{align*}
d \widetilde{L}_{n}(t)= & \widetilde{\sigma}_{n}(t)^{\top}\left(\widetilde{\mu}_{n}(t) d t+d W^{B}(t)\right)  \tag{15.53}\\
d Z^{n}(t) / Z^{n}(t)= & \sigma_{Z^{n}}(t)^{\top} \\
& \times\left(\left(\frac{\tau_{n} \widetilde{\sigma}_{n}(t)}{1+\tau_{n} \widetilde{L}_{n}(t)}+\sigma_{Z^{n}}(t)-\widetilde{\mu}_{n}(t)\right) d t-d W^{B}(t)\right) \tag{15.54}
\end{align*}
$$
^eq15-53
^eq15-54

where

$$
\widetilde{\mu}_{n}(t)=\sum_{j=q(t)}^{n}\left(\frac{\tau_{j} \tilde{\sigma}_{j}(t)}{1+\tau_{j} \widetilde{L}_{j}(t)}+\sigma_{Z^{j}}(t)\right)
$$

and $W^{B}(t)$ is an $m$-dimensional Brownian motion in measure $\mathrm{Q}^{B}$.
Proof. Clearly

$$
\begin{equation*}
\frac{1}{1+\tau_{j} L_{j}(t)}=\frac{1}{1+\tau_{j} \widetilde{L}_{j}(t)} Z^{j}(t) \tag{15.55}
\end{equation*}
$$
^eq15-55

Differentiating (15.55) and matching the $d W$ terms, we obtain that

$$
\frac{\tau_{j} \sigma_{j}(t)}{1+\tau_{j} L_{j}(t)}=\frac{\tau_{j} \widetilde{\sigma}_{j}(t)}{1+\tau_{j} \widetilde{L}_{j}(t)}+\sigma_{Z^{\prime}}(t)
$$

and (15.53) follows.
To derive the drift in (15.54), we note that $\left(1+\tau_{n} L_{n}(t)\right) Z^{n}(t)=1+ \tau_{n} \widetilde{L}_{n}(t)$ is a $\mathrm{Q}^{n+1}$-martingale. Hence, in $T^{n+1}$-forward measure we have the following dynamics,

$$
d Z^{n}(t) / Z^{n}(t)=\sigma_{Z^{\prime \prime}}(t)^{\top}\left(\frac{\tau_{n} \sigma_{n}(t)}{1+\tau_{n} L_{n}(t)} d t-d W^{n+1}(t)\right)
$$

Switching to the spot measure, we obtain (15.54).

The simulation scheme for the two-curve LM model is similar to the one-curve case, with obvious modifications. The initial values of the model primitives, the Libor rates $\widetilde{L}_{n}(0)$ and the forward bond ratios $Z^{n}(t), n= 0, \ldots, N-1$, are derived from the initial values of the discount and index curves. Having specified the Libor rate volatilities $\widetilde{\sigma}_{n}(t)$ and forward bond spread volatilities $\sigma_{Z^{n}}(t), n=1, \ldots, N-1$, we simulate the Libor rates and the bond spreads using SDEs from Proposition 15.5.4. At each time $t$, we may construct an index curve $\widetilde{P}(t, T), T \geq t$, from the simulated Libor rates $\widetilde{L}_{n}(t), n=0, \ldots, N-1$, as in the one-curve case in Chapter 14. From the simulated index curve and the simulated bond spread curve $Z(t, T), T \geq t$, the discount curve $P(t, T), T \geq t$, is calculated via (15.40). The discount curve is used to update the numeraire by (15.49) and, together with the index curve, to calculate market rates needed for evaluating derivative payoffs. In particular, forward swap rates are calculated by projecting Libor rates off the index curve and discounting them off the discount curve (compare to (6.39)), so a rate fixing at $T_{j}$ covering $k-j$ periods is calculated by

$$
\begin{equation*}
S_{j, k-j}(t)=\frac{\sum_{i=j}^{k-1} \tau_{i} P\left(t, T_{i+1}\right) \widetilde{L}_{i}(t)}{\sum_{i=j}^{k-1} \tau_{i} P\left(t, T_{i+1}\right)} \tag{15.56}
\end{equation*}
$$
^eq15-56

In the two-curve setup, the derivation of swaption pricing approximations for model calibration can proceed as for the single-curve case. For this purpose, first rewrite (15.56) in terms of model primitives, to obtain

$$
S_{j, k-j}(t)=\frac{\sum_{i=j}^{k-1} \tau_{i} Z\left(t, T_{i+1}\right) \widetilde{P}\left(t, T_{i+1}\right) \widetilde{L}_{i}(t)}{\sum_{i=j}^{k-1} \tau_{i} Z\left(t, T_{i+1}\right) \widetilde{P}\left(t, T_{i+1}\right)}
$$

Proceeding as in the standard single-curve case and ignoring contributions from the $Z\left(t, T_{i+1}\right)$ terms in the swap rate dynamics - a respectable approximation even with stochastic spreads - we obtain the following dynamics in the appropriate annuity measure

$$
\begin{equation*}
d S_{j, k-j}(t)=\left.\sum_{n=1}^{N-1} \frac{\partial S_{j, k-j}(t)}{\partial \widetilde{L}_{n}(t)}\right|_{t=0} \tilde{\sigma}_{n}(t)^{\top} d W^{A_{j, k-\jmath}}(t) \tag{15.57}
\end{equation*}
$$
^eq15-57

This formula can be used as the starting point for the usual European swaption approximations. Compared to earlier results in Section 14.4.2, the presence of the basis spread leads to slightly altered expressions for the weights $\partial S(t) /\left.\partial \widetilde{L}_{n}(t)\right|_{t=0}$; the exact form of these weights is left for the reader to derive.

To conclude our discussion of dynamic two-curve modeling in the LM framework, we note that our decision to use forward bond spreads $Z^{j}(t)$ as variables driving the spread between the index and discount curves is not the only choice. For example, we could have used the spreads between the simply compounded forward rates computed off the index and discount curves,
$\widetilde{L_{j}}(t)-L_{j}(t), j=0, \ldots, N-1$, instead. This would lead to a somewhat different, but equally tractable, extension. As this area of interest rate modeling is still rapidly evolving, consensus on what are the right variables to use has not emerged yet.

### 15.5.4 Deterministic Spread

It is not hard to see that the values of derivative securities that do not directly reference the spread between the index and the discount curve which is the majority of them - respond approximately linearly to this spread. For one, the spread is usually confined to a rather tight range, ensuring that a linear approximation may be adequate. Moreover, changes in the spread will normally affect the discount curve substantially more than the index curve, as the latter is predominantly calibrated to linear market instruments such as FRAs and swaps. Values of most securities, even exotic ones, tend to be approximately linear to changes in discounting.

A fairly self-evident rule states that if the dependence of a value of a security on a given market factor is approximately linear, there is little reason to model that factor with a stochastic process for the purposes of valuation and hedging of that security. As a consequence, it is not uncommon to assume deterministic evolution for the spread between the discount and index curves. In the framework we developed, this may be achieved by setting the volatility of the spread $\sigma_{\varepsilon}(t, T)$ to zero. The two-curve LM model simplifies accordingly. Obviously, one does not need to simulate $\varepsilon(t, T)$ anymore, as

$$
\varepsilon(t, T)=\varepsilon(0, T), \quad Z(t, T)=Z(0, T) / Z(0, t)
$$

for all $t, T$ such that $0 \leq t \leq T$. Furthermore, the drift term in (15.53) reduces to the standard single-curve expression, and discount factors $P(t, T)$ are obtained from pseudo-discount ones by

$$
P(t, T)=\frac{P(0, t, T)}{\widetilde{P}(0, t, T)} \widetilde{P}(t, T), \quad T \geq t
$$

Derivation of swaption pricing expressions proceeds as in Section 15.5.3 where the randomness in the spreads was already ignored in (15.57).

### 15.6 SV Models with Non-Zero Correlation

Let us return to the SV-type LM model we considered in Section 14.2.5. In the spot measure, we recall that forward rate dynamics are postulated to be of the form

$$
\begin{aligned}
d L_{n}(t) & =\sqrt{z(t)} \varphi\left(L_{n}(t)\right) \lambda_{n}(t)^{\top}\left(\sqrt{z(t)} \mu_{n}(t) d t+d W^{B}(t)\right) \\
\mu_{n}(t) & =\sum_{j=q(t)}^{n} \frac{\tau_{j} \varphi\left(L_{j}(t)\right) \lambda_{j}(t)}{1+\tau_{j} L_{j}(t)}
\end{aligned}
$$

with

$$
d z(t)=\theta\left(z_{0}-z(t)\right) d t+\eta \psi(z(t)) d Z(t)
$$

See Section 14.2.5 for further details on the notation.
In previous work, we assumed that the scalar Brownian motion $Z(t)$ was independent of all components of the m-dimensional Brownian motion $W(t)$. This assumption allowed us to switch from the spot measure into more convenient forward measures, without altering the form of the process for $z(t)$. We shall now briefly consider relaxing this assumption. In particular, let us assume a non-zero deterministic correlation vector $\rho(t)$ between $Z(t)$ and $W^{B}(t)$. In this case, Lemma 14.2.6 tells us the dynamics of $z(t)$ in the forward measure $\mathrm{Q}^{T_{n+1}}$. As we defined $\left\langle d W^{B}(t), d Z(t)\right\rangle=\rho(t) d t$, then

$$
d z(t)=\theta\left(z_{0}-z(t)-\frac{\eta}{\theta} \psi(z(t)) \sqrt{z(t)} \rho(t)^{\top} \mu_{n}(t)\right) d t+\eta \psi(z(t)) d Z^{n+1}(t)
$$

where $Z^{n+1}(t)$ in a $Q^{T_{n+1}}$-Brownian motion. In certain important cases, this expression can be simplified and approximated somewhat, as demonstrated in the following corollary.
Corollary 15.6.1. When $\psi(z(t))=\sqrt{z(t)}$, we have

$$
\begin{equation*}
d z(t)=\theta\left(z_{0}-\alpha(t) z(t)\right) d t+\eta \sqrt{z(t)} d Z^{n+1}(t) \tag{15.58}
\end{equation*}
$$
^eq15-58

where

$$
\begin{equation*}
\alpha(t)=1+\frac{\eta}{\theta} \rho(t)^{\top} \mu_{n}(t) \approx 1+\frac{\eta}{\theta} \rho(t)^{\top} \sum_{j=q(t)}^{n} \frac{\tau_{j} \varphi\left(L_{j}(0)\right) \lambda_{j}(t)}{1+\tau_{j} L_{j}(0)} . \tag{15.59}
\end{equation*}
$$
^eq15-59

Notice that in Corollary 15.6.1 we have used an approximation under which the multiplier $\alpha(t)$ in the drift term of (15.58) is approximately deterministic. The approximation is based on the same "freeze along forwards" idea as was used in Section 14.4.2 and elsewhere, and makes the SDE for $z(t)$ affine. As such, from the results in Chapter 8 we would expect that a Fourier-based approach would allow for (approximate) caplet pricing if the dynamics for $L_{n}(t)$ are themselves affine, e.g. if $\varphi(x)=a x+b$. We verify this below for the log-normal case ( $\varphi(x)=x$ ); for displaced log-normal dynamics the result follows along similar lines.

Proposition 15.6.2. Define $X_{n}(t)=\ln L_{n}(t)$ where

$$
\begin{gathered}
d L_{n}(t) / L_{n}(t)=\sqrt{z(t)} \lambda_{n}(t)^{\top} d W^{n+1}(t) \\
d z(t)=\theta\left(z_{0}-\alpha(t) z(t)\right) d t+\eta \sqrt{z(t)} d Z^{n+1}(t)
\end{gathered}
$$

for a deterministic function $\alpha(t)$ and $\left\langle d W^{n+1}(t), d Z^{n+1}(t)\right\rangle=\rho(t) d t$. The moment-generating function

$$
\Psi_{X_{n}}(u)=\mathrm{E}^{T_{n+1}}\left(e^{u X_{n}\left(T_{n}\right)}\right)
$$

can be written as

$$
\Psi_{X_{n}}(u)=e^{A(0, T)+X_{n}(0) B(0, T)+z(0) C(0, T)},
$$

where $A, B$, and $C$ solve a system of Riccati ODEs

$$
\begin{aligned}
0= & A^{\prime}(t, T)+\theta z_{0} C(t, T) \\
0= & B^{\prime}(t, T) \\
0= & C^{\prime}(t, T)-\frac{1}{2}\left\|\lambda_{n}(t)\right\|^{2} B(t, T)+\frac{1}{2}\left\|\lambda_{n}(t)\right\|^{2} B(t, T)^{2} \\
& -\theta \alpha(t) C(t, T)+\frac{1}{2} \eta^{2} C(t, T)^{2}+\eta \lambda_{n}(t)^{\top} \rho(t) B(t, T) C(t, T),
\end{aligned}
$$

with the terminal conditions

$$
A(T, T)=0, \quad B(T, T)=u, \quad C(T, T)=0 .
$$

Proof. Positing the form

$$
\mathrm{E}^{T_{n+1}}\left(e^{u X_{n}\left(T_{n}\right)} \mid X_{n}(t)=x, z(t)=v\right)=e^{A(t, T)+x B(t, T)+v C(t, T)}
$$

and substituting it into the Feynman-Kac PDE that corresponds to the dynamics of $\left(X_{n}(t), z_{n}(t)\right)$, we obtain the following equation

$$
\begin{aligned}
0=A^{\prime} & (t, T)+x B^{\prime}(t, T)+v C^{\prime}(t, T) \\
& -\frac{1}{2} v\left\|\lambda_{n}(t)\right\|^{2} B(t, T)+\frac{1}{2} v\left\|\lambda_{n}(t)\right\|^{2} B(t, T)^{2} \\
& +\theta\left(z_{0}-\alpha(t) v\right) C(t, T)+\frac{1}{2} \eta^{2} v C(t, T)^{2} \\
& +v \eta \lambda_{n}(t)^{\top} \rho(t) B(t, T) C(t, T)
\end{aligned}
$$

Combining the terms in $x$ and $v$ together, the result follows.
The result of the proposition above can be combined with Theorem 8.4.1 to allow for analytical pricing of caplets by Fourier methods. Pricing of swaptions follows a similar line of attack. To start, we write down the drift of $z(t)$ under the corresponding annuity measure; while somewhat more complicated in appearance, it will still be linear in $z$ after application of the "freezing" technique. Then, the moment-generating function for the logarithm of the swap rate is available, allowing for analytic pricing of swaptions via Fourier methods. We trust the reader to fill in missing details.

### 15.7 Multi-Stochastic Volatility Extensions

### 15.7.1 Introduction

In the specification (15.22)-(15.23) of the LM model, a single stochastic volatility process $\sqrt{z(t)}$ is used to scale the diffusion coefficients of all forward rates. As such, the volatility structure of the model is only allowed (nearly) parallel moves up and down. While sufficiently rich to introduce the volatility smile for all European swaptions, one has to wonder about the limitations of this one-factor specification.

The value of many exotic interest rate derivatives, sometimes called "first generation" exotics, is primarily linked to the overall level of the interest rate curve. This class includes, for example, Bermudan swaptions, callable inverse floaters, or callable range accruals on a Libor or CMS rate. For such instruments, having a single stochastic volatility factor applied to all rates is typically adequate.

On the other hand, interest rate exotics linked to the spread of two CMS rates, such as CMS spread callable swaps (see Section 5.13.3) or CMS spread TARNs (see Section 5.15.2) derive their value from the distribution of the slope of the interest rate curve. Just like a single-factor model is unsuitable for pricing such exotics - being unable to represent the changes in the slope of the curve - a common stochastic volatility factor applied to all rates does not always allow for sufficient control over the distribution of the slope of the interest rate curve. In particular, such a specification does not typically allow for much control over the finer features of the volatility smile of the spread, e.g. its slope or curvature.

We will have quite a bit more to say about modeling the smile of a CMS spread later on in Chapter 17, but let us take these observations as a rationale (or excuse) for sketching an extension of the LM model that allows for some de-correlation in stochastic volatility factors applied to different rates. While many roads could be taken, the route we shall suggest here

- Incorporates the features necessary for realistic CMS spread modeling.
- Remains relatively parsimonious.
- Can still be calibrated using analytical approximations to prices of caps and swaptions.

We must point out that the whole area of multi-dimensional stochastic volatility interest rate modeling is quite new, so we keep the discussion suitably brief, with details to be filled by future research.

### 15.7.2 Setup

One can take a view, admittedly not inconsistent with the philosophy of LM modeling, that each Libor rate should be driven by its own stochastic variance process. However, it should be obvious that such specification
would be quite unwieldy, and would likely not lend itself easily to closedform approximations for swaption prices. With the point of our proposed extension focused primarily on controlling the smiles of CMS spreads, we instead consider a parsimonious extension that involves only two stochastic variance processes.

We define

$$
\begin{equation*}
d z^{i}(t)=\theta^{i}\left(z_{0}^{i}-z^{i}(t)\right) d t+\eta^{i} \sqrt{z^{2}(t)} d Z^{i}(t), \quad i=1,2 . \tag{15.60}
\end{equation*}
$$
^eq15-60

Note the time-independence of parameters; while a more general specification is certainly possible, we keep our focus on more important details. Moreover, we require

$$
\left\langle d Z^{1}(t), d Z^{2}(t)\right\rangle=0
$$

For the applications we have in mind, we need to allow for non-zero correlation between Brownian motions driving the Libor rates and those driving the stochastic variances. Hence, (15.60) is understood to hold under $Q^{B}$, the spot Libor measure only. Under the same measure, we assume that the Libor rates follow

$$
\begin{align*}
& d L_{n}(t) / \varphi\left(L_{n}(t)\right)=\sqrt{z^{1}(t)} \lambda_{n}^{1}(t)^{\top}\left(d W^{1}(t)+\mu_{n}^{1}(t) d t\right) \\
& \quad+\sqrt{z^{2}(t)} \lambda_{n}^{2}(t)^{\top}\left(d W^{2}(t)+\mu_{n}^{2}(t) d t\right), \quad n=1, \ldots, N-1 \tag{15.61}
\end{align*}
$$
^eq15-61

Here, $W^{1}(t)$ and $W^{2}(t)$ are two independent copies of a $d$-factor Brownian motion. Moreover, we assume that

$$
\begin{equation*}
\left\langle d W^{i}(t), d Z^{i}(t)\right\rangle=\chi^{i}, \quad i=1,2 \tag{15.62}
\end{equation*}
$$
^eq15-62

but also

$$
\begin{equation*}
\left\langle d W^{1}(t), d Z^{2}(t)\right\rangle=\left\langle d W^{2}(t), d Z^{1}(t)\right\rangle=0 \tag{15.63}
\end{equation*}
$$
^eq15-63

Under the $T_{n+1}$-forward measure, the Libor rate $L_{n}(t)$ is a martingale,

$$
\begin{align*}
d L_{n}(t) / \varphi\left(L_{n}(t)\right) & =\sqrt{z^{1}(t)} \lambda_{n}^{1}(t)^{\top} d W^{1, T_{n+1}}(t) \\
+ & \sqrt{z^{2}(t)} \lambda_{n}^{2}(t) d W^{2, T_{n+1}}(t), \quad n=1, \ldots, N-1 \tag{15.64}
\end{align*}
$$
^eq15-64

### 15.7.3 Pricing Caplets and Swaptions

Following the same techniques used in Section 15.6 above, it is straightforward to show that in the $T_{n+1}$-forward measure, the dynamics of the stochastic variance processes are
$d z^{i}(t)=\theta^{2}\left(z_{0}^{i}-z^{i}(t)\right) d t-\eta^{i} \nu^{i, n+1}(t, \mathbf{L}(t)) z^{i}(t) d t+\eta^{2} \sqrt{z^{i}(t)} d Z^{i, T_{n}+1}(t)$,
$i=1,2$, where

$$
\nu^{i, n+1}(t, \mathbf{L}(t))=\left(\chi^{i}\right)^{\top} \sum_{j=q(t)}^{n} \frac{\tau_{j} \varphi\left(L_{j}(t)\right) \lambda_{j}^{i}(t)}{1+\tau_{j} L_{j}(t)}
$$

Freezing the drifts, in the same manner as in Corollary 15.6.1, we get
$d z^{i}(t)=\theta^{i}\left(z_{0}^{i}-z^{i}(t)\right) d t-\eta^{i} \nu^{i, n+1}(t, \mathbf{L}(0)) z^{i}(t) d t+\eta^{i} \sqrt{z^{i}(t)} d Z^{i, T_{n+1}}(t)$,
and we obtain that (15.64), (15.65) constitute an affine specification (thanks to (15.63)). Hence, the moment-generating function could be represented in a quasi-closed form of an exponential of coefficients computable from Riccati equations, and prices of caplets could be obtained by Fourier inversion of the moment-generating function. We omit straightforward details.

As far as swaptions are concerned, as is the case in many other LM model specifications, we could derive the dynamics of a swap rate of essentially the same form as that of the Libor rates. Then, the same arguments as above can be applied to compute European swaption prices.

### 15.7.4 Spread Options

The intention of introducing two stochastic variance processes above is to be able to control the volatility smile ${ }^{5}$ of the spread option, while keeping the smiles of individual swap rates fixed and in calibration with observed swaption values. As shown later in Chapter 17, such control can be achieved if we have a mechanism to affect i) the correlation of the stochastic variance processes affecting the two rates in the spread; and ii) the "cross" correlations between a forward rate and the variance process of the other forward rate in the spread. The specification outlined above allows for this, as we shall now demonstrate.

Consider two forward Libor rates, $L_{n}(t)$ and $L_{m}(t), n \neq m$. For simplicity we assume that in $(15.61) \varphi(x) \equiv 1$ and $\lambda_{k}^{i}(t) \equiv \lambda_{k}^{i}, i=1,2, k=n, m$. Focusing on the diffusion terms only, we have (drifts and probability measure are irrelevant),
$d L_{k}(t)=O(d t)+\sqrt{z^{1}(t)}\left(\lambda_{k}^{1}\right)^{\top} d W^{1}(t)+\sqrt{z^{2}(t)}\left(\lambda_{k}^{2}\right)^{\top} d W^{2}(t), k=n, m$,
$d z^{i}(t)=O(d t)+\eta^{i} \sqrt{z^{i}(t)} d Z^{i}(t), \quad i=1,2$.

In the following easily proven result, we rewrite the dynamics in the form convenient for correlation analysis.

Proposition 15.7.1. Assume that $L_{n}$ and $L_{m}$ satisfy (15.66)-(15.67) above, with the correlation structure as in (15.62) and (15.63). The joint dynamics

[^73]of the two Libor rates and their stochastic variance processes can then be written as

$$
d L_{k}(t)=O(d t)+\sqrt{u^{k}(t)} d U^{k}(t), d u^{k}(t)=O(d t)+\sqrt{\eta^{k}(t)} d X^{k}(t), k=n, m
$$
where

$$
\begin{aligned}
u^{k}(t) & =z^{1}(t)\left\|\lambda_{k}^{1}\right\|^{2}+z^{2}(t)\left\|\lambda_{k}^{2}\right\|^{2} \\
\eta^{k}(t) & =z^{1}(t)\left\|\lambda_{k}^{1}\right\|^{4}\left(\eta^{1}\right)^{2}+z^{2}(t)\left\|\lambda_{k}^{2}\right\|^{4}\left(\eta^{2}\right)^{2}
\end{aligned}
$$
and

$$
\begin{aligned}
& d U^{k}(t)=\frac{1}{\sqrt{u^{k}(t)}}\left(\sqrt{z^{1}(t)}\left(\lambda_{k}^{1}\right)^{\top} d W^{1}(t)+\sqrt{z^{2}(t)}\left(\lambda_{k}^{2}\right)^{\top} d W^{2}(t)\right) \\
& d X^{k}(t)=\frac{1}{\sqrt{\eta^{k}(t)}}\left(\sqrt{z^{1}(t)}\left\|\lambda_{k}^{1}\right\|^{2} \eta^{1} d Z^{1}(t)+\sqrt{z^{2}(t)}\left\|\lambda_{k}^{2}\right\|^{2} \eta^{2} d Z^{2}(t)\right)
\end{aligned}
$$

It follows from Proposition 15.7.1 that our model setup allows for essentially independent control of the correlations between the variance processes $u^{n}$ and $u^{m}$, as well as the correlations between $L_{m}$ and $u^{n}$ and $L_{n}$ and $u^{m}$ as should be clear from the expression for one of the correlations (others are similar):

$$
\begin{aligned}
\operatorname{Corr}\left(d L^{n}(t), d u^{m}(t)\right) & =\frac{1}{\sqrt{u^{n}(t) \eta^{m}(t)}} \\
\times & \left(z^{1}(t)\left\|\lambda_{m}^{1}\right\|^{2} \eta^{1}\left(\lambda_{n}^{1}\right)^{\top} \chi^{1}+z^{2}(t)\left\|\lambda_{m}^{2}\right\|^{2} \eta^{2}\left(\lambda_{n}^{2}\right)^{\top} \chi^{2}\right)
\end{aligned}
$$

The results extend to swap rate spreads in a predictable, and largely mechanical, fashion, allowing us to set up an LM calibration that targets quantities linked to the shape of smile of various CMS spread options of interest. As the reader might expect, the formulas become rather cumbersome and we do not list them here.

At this point we have gathered enough results for our later discussion on spread option pricing, so we conclude the analysis here. We return to spread options in stochastic volatility models in Chapter 17.

### 15.7.5 Another Use of Multi-Dimensional Stochastic Volatility

In models for foreign exchange or for equity prices, multi-dimensional stochastic volatility is typically used not as a way to refine spread option pricing, but rather as i) a mechanism to induce multiple ${ }^{6}$ time-scales in the meanreverting behavior of volatility; or ii) to control the evolution of the volatility smile through time. A discussion of multiple time-scales in empirical data

[^74]can be found in Perello et al. [2004], and Kainth and Saravanamuttu [2007] (among others) discuss applications specific to foreign exchange, with an emphasis on smile dynamics. We also note that Andersen and BrothertonRatcliffe [2005] introduce a tractable alternative to (15.61) that has perfect correlation between the variances of all Libor forwards, but still introduces multiple time-scales in the variance dynamics; this setup (the details of which we omit) is mainly useful for applications where the correlation between variances of different forwards is thought to be high.


[^0]:    ${ }^{1}$ The original paper by Ho and Lee was set exclusively in discrete time. The continuous-time version of the model developed here is, we feel, significantly more transparent.

[^1]:    ${ }^{2}$ Due to the presence of the $T$-dependent "convexity adjustment" term $\sigma_{r}^{2}(T-t)$ in the drift of the forward rate process, net forward curve movements are not perfectly parallel. Were this the case, it is well-known that the model would be arbitrageable.

[^2]:    ${ }^{3}$ An exception to this observation is the humped volatility term structure that can often be observed in caplet markets. We return to this issue in Section 10.1.2.3.

[^3]:    ${ }^{4}$ On occasion there is significant uncertainty in the market about the intentions of monetary authorities, in which case the caplet hump may disappear temporarily.

[^4]:    ${ }^{5}$ As argued above, normally we would pick $\varkappa(t)$ to be a constant. A more detailed examination of the estimation of mean reversions - and the role it plays in Bermudan swaption pricing - can be found in Chapter 13.

[^5]:    ${ }^{6}$ Note that it is common to set $T_{0}=0$ when defining swaption strips, a convention that slightly clashes with the notation used above when deriving swaption formulas (where the swaption maturity $T_{0}>0$ ). We shall later, in Chapter 14, develop more formal notation for indexation, but for now trust the reader's ability to adapt generic swaption formulas to the swaption strip convention.

[^6]:    ${ }^{7}$ An alternative approach that is advocated by some is to set the mean reversion $\%$ to zero when determining $x_{\max }$ and $x_{\min }$, in which case $\alpha$ can be reduced.

[^7]:    ${ }^{8}$ The results of the next sections do, however, often generalize to timedependence in $\alpha$ and $\beta$. See Remark 10.2.3, for example.

[^8]:    ${ }^{9}$ As pointed out in Section 9.1, depending on the level of accuracy required, the Runge-Kutta numerical solution of the ODEs can sometimes have higher computational efficiency.

[^9]:    ${ }^{10}$ For time-homogeneous models, closed-form pricing formulas for options on discount bonds exist for some models, including the CIR model (see Cox et al. (385).

[^10]:    ${ }^{11}$ Recall that we solve $b(t, T+\Delta)$ backward in time from the known boundary condition at $t=T+\Delta$.
    ${ }^{12}$ For instance, if $\sigma_{g}(t)$ is the volatility function in the Gaussian model, then we can extract an estimate for $\sigma(t)$ from the relation $\sigma(t) \sqrt{\alpha+f(0, t) \beta} \approx \sigma_{g}(t)$.

[^11]:    ${ }^{13}$ Alternatively, we can also set $\varkappa_{i}=\left(\varkappa\left(t_{i}\right)+\varkappa\left(t_{i+1}\right)\right) / 2$, and so forth.

[^12]:    ${ }^{1}$ A computationally efficient recursive procedure for bond pricing can be found in Hansen and Jørgensen [1998].

[^13]:    ${ }^{2}$ Even though the log-normal distribution has finite moments of all orders, the moment-generating function is infinite at any positive argument.
    ${ }^{3}$ For details on Jensen's incquality, see the proof of Proposition 11.1.3.
    ${ }^{4}$ Recall from Section 4.1.2 that the risk-neutral expectation of a random variable must, in the absence of arbitrage, equal its traded futures price.

[^14]:    ${ }^{5}$ The proof of this statement is straightforward, and follows from a standard comparison theorem for SDEs (see p. 293 of Karatzas and Shreve [1997]).

[^15]:    ${ }^{6}$ The authors also develop eigenfunction expansions for cases where the shadow rate $r^{*}(t)$ follow time-homogeneous diffusions more complicated than the Vasicek model.

[^16]:    ${ }^{7}$ The ODEs are easily derived by substituting the right-hand side of (11.21) into the no-arbitrage PDE for $P\left(t, T^{\prime}\right)$ in the FV model.

[^17]:    ${ }^{8}$ The same holds for exercise values of callable securities.
    ${ }^{9}$ An initial guess for $\vartheta_{i-1}$ could be $\vartheta_{\imath-1}=\vartheta_{\imath-2}$. Subsequent guesses would be performed by a root-search algorithm.

[^18]:    ${ }^{10}$ For some models, the density can grow to infinity at the boundary, notably in the CIR model for $x \rightarrow 0+$ when the Feller condition is violated. Should that be the case, a more careful analysis of boundary conditions is required, see e.g. Lucic [2008].

[^19]:    ${ }^{11}$ This model has been advocated by Brigo and Mercurio [2001] as an easy-toimplement alternative to a true time-dependent CIR process. For reasons explained above, this model has certain drawbacks that require careful evaluation.

[^20]:    ${ }^{12}$ Multi-dimensional extensions are possible.

[^21]:    ${ }^{13}$ Instead of discretizing $Y(t)$, we could also discretize $I(t)=\ln Y(t)$, as in Section 10.1.6.1. For variation we use $Y(t)$ in this section.

[^22]:    ${ }^{14}$ With this approach, the numeraire-mapping function is undefined for times that are not in the tenor structure. The "discrete" nature of the resulting model is one of the common criticisms of the MF approach. Pragmatically, it means that all dates of interest for a particular derivative security should be added to the tenor structure, or interpolation schemes not unlike those considered in Section 15.1 need to be designed.

[^23]:    ${ }^{15}$ It is straightforward to extend the arguments to the case of displaced $\log$ normal Libor rates, say. We leave this to the reader.

[^24]:    ${ }^{16}$ We could use market prices for this, or we could use option prices computed from a vanilla model that we wish for the MF model to emulate.

[^25]:    ${ }^{17}$ Indeed, we are unaware of the existence of any published MF models with stochastic volatility.
    ${ }^{18}$ As we shall see in Chapter 13, one component in the state vector is locally deterministic (i.e. it involves no Brownian motion term), so in a sense the quasiGaussian model has a state process dimension of "one-and-a-half", allowing for significant speed-ups in the numerical implementation.

[^26]:    ${ }^{1}$ The judgment of whether a security is convex in forward rate twists and tilts can often be quite difficult. Some securities that one might guess should be sensitive to forward rate correlation in fact only display material sensitivity to forward rate auto-correlation. Bermudan swaptions are a good example; see Chapter 19 for more details.

[^27]:    ${ }^{2}$ Some readers may recognize $J_{Q}(T)$ as the product integral of $-Q(t)$ on $[0, T]$; see Dollard and Friedman [1979]. In the probability literature, the product integral is often referred to as the fundamental matrix, see Arnold [1974] or Karatzas and Shreve [1997].
    ${ }^{3}$ Recall that the exponential of a square matrix $A$ is defined as $e^{A}= \sum_{k=0}^{\infty} A^{k} / k!$.

[^28]:    ${ }^{4}$ If $g(t)$ is not lower diagonal, we can always change variables (via a Cholesky decomposition, say) to make it so.

[^29]:    ${ }^{5}$ The material in this section is largely inspired by Balasanov [1996].

[^30]:    ${ }^{6}$ PC analysis was introduced in Section 3.1.3; its application to statistical analysis of interest rate curve movements will be described in details later, in Chapter 14.

[^31]:    ${ }^{7}$ This choice is made largely for reasons of familiarity. We indicate later (at the very end of this section) how the choice of different state variables may streamline the method.
    ${ }^{8}$ We here use $\sigma(t)=\sigma_{x}(t)^{\top}$.

[^32]:    ${ }^{9}$ Reflecting a rather intuitive fact that in a Gaussian model a swap rate is approximately Gaussian.
    ${ }^{10}$ Or other derivatives with first-order correlation dependence. See Chapter 17 for details on yield curve spread options.

[^33]:    ${ }^{11}$ As is standard practice in much of the literature on multi-factor affine models, for notational simplicity we here assume that $\mu$ and $\sigma$ do not depend on time $t$. See the comments at the end of Section 12.2.4 for extensions to time-dependent parameters.

[^34]:    ${ }^{12}$ Throughout this section, E denotes expectation in the risk-neutral measure Q.

[^35]:    ${ }^{13}$ If $a+b^{\top} x=c+d^{\top} x$ holds for a non-empty open set, then $a=c$ and $b=d$.

[^36]:    ${ }^{14}$ Longstaff and Schwartz [1993] suggest making $\varkappa_{2}$ a function of time.

[^37]:    ${ }^{15}$ In the language of Section 11.2.3, the stochastic volatility driver is here evidently of the spanned type. See also Section 12.3.6.

[^38]:    ${ }^{16}$ Contrary to some claims in the literature, see e.g. Boyarchenko and Levendorski [2007].

[^39]:    ${ }^{1}$ Also known as pseudo-Gaussian or Cheyette models.

[^40]:    ${ }^{2}$ We use a small $q$ in the abbreviation $q G$ to avoid notational conflict with the quadratic Gaussian (QG) models of Chapter 12.

[^41]:    ${ }^{3}$ Constant Elasticity of Variance, see Section 7.2.

[^42]:    ${ }^{4}$ We normally use $y(t, T, M)$ for continuously compounded forward yield, see Section 4.1.1, and $F(t, T, M)$ for a futures rate, see Section 4.1.2, but for the remainder of this chapter only we allow ourselves a slight notational inconsistency to avoid the possibility of confusion between $y$, the state variable, and $y$, the forward yield.

[^43]:    ${ }^{5}$ Note that it is also possible to calibrate mean reversion to a whole collection of European swaptions (i.e. more than two) sharing the same expiry; the motivation for such calibration and an outline of the algorithm are presented in Section 19.4.4.

[^44]:    ${ }^{6}$ Removing deterministic drift components from variables before discretizing a PDE is a useful trick for any model.

[^45]:    ${ }^{7}$ For readers not fully familiar with Libor market models, we recommend reading Chapters 14 and 15 before proceeding with this section.

[^46]:    ${ }^{1}$ As we have seen in earlier chapters, for special choices of the forward rate volatility we can sometimes identify a finite-dimensional Markovian representation of the forward curve that eliminates the need to store the entire curve. This is not possible in general, however.

[^47]:    ${ }^{2}$ In fact, as discussed in Jamshidian [1997], one does not need to assume that a short rate process exists when constructing an LM model.

[^48]:    ${ }^{3}$ We briefly return to the general case in Section 15.6.

[^49]:    ${ }^{4 r}$ The use of sliding forward rates, i.e. forward rates with a fixed time to maturity rather than a fixed time of maturity, is often known as the Musiela parameterization.

[^50]:    ${ }^{5}$ For each date in the time grid $t_{j}$ we construct the forward curve from market observable swaps, futures, and deposits, using the techniques from Chapter 6.

[^51]:    ${ }^{6}$ The introduction of user-specified weights into this norm is a straightforward extension.

[^52]:    ${ }^{7}$ In our experience, the majorization method in Pietersz and Groenen [2004] is faster than the method in Proposition 14.3.2 but less robust, at least for large and irregular correlation matrices.

[^53]:    ${ }^{8}$ An alternative would be to rely on extrapolation.

[^54]:    ${ }^{9} \mathrm{~A}$ more refined approach, especially for low values of time-to-maturity, is advisable for some applications where the fine structure of short-term volatilities is important. See the discussion in Remark 15.1.1.

[^55]:    ${ }^{10}$ Earlier results show that the true instantaneous covariance matrix for forward rates may involve DVF- or SV-type scales on the elements of $c$. For the purposes of calibration of $\lambda_{k}$, we omit these scales.

[^56]:    ${ }^{11}$ Since the matrix $D$ shrinks in $T_{n}$, we need to repeat the PCA analysis at each tenor date. Alternatively, but suboptimally, we can do the PCA analysis only once at time 0 , pruning the results as needed for other values of $T_{n}$.

[^57]:    ${ }^{12}$ We also revisit the subject in the context of callable Libor exotics in Section 18.1.
    ${ }^{13}$ One would here limit $T_{s}$ to be no larger than $T_{N}$, so the total number of swaptions would be less that $N_{t} \cdot N_{x}$. See our discussion of redundant grid entries in Section 14.5.2.

[^58]:    ${ }^{14}$ Depending on how table boundary elements are treated, notice that the range for $i$ and $j$ may not always be as stated in (14.51).

[^59]:    ${ }^{15}$ Another approach to producing more equitable scaling involves using relative (=percentage) price errors, rather than absolute price errors.

[^60]:    ${ }^{16}$ donlp2 SQP/ECQP algorithm, available on www.mathematik. tu-darmstadt.de:8080/ags/ag8/Mitglieder/spellucci_de.html.
    ${ }^{17}$ www.imsl.com.
    ${ }^{18}$ www.nag.com.

[^61]:    ${ }^{19}$ As many of the algorithms we have in mind compute, at the very least, the sensitivity of each calibration target to each input variable, the computational complexity is at least $O\left(m^{2}\right)$; if the order of complexity is higher, the case for problem splitting is even more compelling.

[^62]:    ${ }^{20}$ Note that our original norm $\mathcal{I}$ now also is a function of $\xi$, since cap and swaption prices depend on the correlation matrix $R$.

[^63]:    ${ }^{21}$ By representing spread options through implied term correlations, the information extracted from spread options is more "orthogonal" to that extracted from caps and swaptions, something that can help improve the numerical properties of the calibration algorithm, particularly if split calibration approach is used.
    ${ }^{22}$ Assuming that there are only a few parameters that define the shape of the function. We generally recommend using simple skew functions that can be described by a single-parameter family, such as linear or power functions.

[^64]:    ${ }^{23}$ In addition to these time-stepping schemes for the forward rates, it may be necessary to simultaneously evolve stochastic volatility variables if one works with models such as those in Section 14.2.5.

[^65]:    ${ }^{24}$ This precludes the existence of other large gaps in the simulation time line prior to time $T_{j}$. When using a hybrid measure such as $\bar{P}_{j}$, we would need to

[^66]:    ${ }^{25}$ See Section 3.2.9 for an introduction to the Brownian bridge, albeit for a somewhat different application.

[^67]:    ${ }^{26}$ But not error-free, of course - there will still be a statistical zero-mean error on the simulation results. See Section 14.6.4 below.

[^68]:    ${ }^{27}$ Suppose that we estimate $\mathrm{E}(X) \approx \mathrm{E}\left(X^{\prime}+Y^{\prime}-\mu_{Y}\right)$, where $\mu_{Y}=\mathrm{E}\left(\mathrm{Y}^{\prime}\right)+O\left(\Delta^{p}\right)$ and $\mathrm{E}\left(X^{\prime}\right)=\mathrm{E}(X)+O\left(\Delta^{p}\right)$. Then clearly also $\mathrm{E}\left(X^{\prime}+Y^{\prime}-\mu_{Y}\right)=\mathrm{E}(X)+O\left(\Delta^{p}\right)$.

[^69]:    ${ }^{1}$ And, implicitly, calendar time, which we ignore here as we work with a fixed $t$.

[^70]:    ${ }^{2}$ As discount bond volatility vectors $\sigma_{P}\left(t, T^{\prime}\right)$ are generally non-deterministic, a suitable approximation is required to make sense of this formula. This is most easily done by freezing any state variables appearing in bond volatilities, such as Libor rates, at their time 0 values.

[^71]:    ${ }^{3}$ As before, we here assume zero correlation between $Z(t)$ and $W(t)$, but relax this assumption in Section 15.6.

[^72]:    ${ }^{4}$ Arguably, diffusive dynamics do not reflect the observed movements of the spread particularly well, as the spread between discount and index curves tends to remain stable for extended period of times, followed by sometimes violent dislocations. On the other hand, our needs in capturing the distribution of the spread within a model are usually quite modest, rarely exceeding the requirement to capture an overall level of its dispersion. As long as the volatility of the spread is reasonable, this requirement will be satisfied by a diffusion model.

[^73]:    ${ }^{5}$ One way to define such a smile is in terms of the implied Bachelier volatility of the spread itself; see Sections 14.4.3 and 17.4.1 for more details.

[^74]:    ${ }^{6}$ Basically this means that $\kappa_{1}$ is either much larger or much smaller than $\kappa_{2}$.

