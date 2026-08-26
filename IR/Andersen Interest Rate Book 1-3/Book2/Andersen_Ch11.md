---
tags: [quant/knowledge, quant/books, interest-rates]
---

## One-Factor Short Rate Models II

While the affine specification (including the Gaussian case) of Chapter 10 is, without doubt, the most popular one-factor short rate model in practice, quite a few other models have been proposed in the literature. In this chapter we cover the most important of these models, paying special attention to the case where the short rate is log-normal. We also briefly discuss some issues in the econometric estimation of short rate models, and introduce the important concept of unspanned stochastic volatility.

As most of the models introduced in this chapter have no analytical bond reconstitution formulas, their calibration to the initial term structure requires numerical work. Accordingly, the second half of this chapter is dedicated to numerical methods for pricing and, especially, calibration of models based on generic short rate SDE. Particularly useful in this regard is the discussion in Section 11.3.2 on efficient finite difference schemes based on the important concept of forward induction.

### 11.1 Log-Normal Short Rate Models

Given the pervasiveness of the log-normal Black-Scholes model in derivatives pricing theory, it should come as no surprise that many authors have attempted to specify one-factor short rate models where the dynamics of $r(t)$ are of the form $d r(t)=O(d t)+\sigma_{r}(t) r(t) d W(t)$ with deterministic $\sigma_{r}(t)$. This section reviews this class of models which, somewhat surprisingly, turns out to have a number of rather severe drawbacks.

### 11.1.1 The Black-Derman-Toy Model

The Black-Derman-Toy (BDT) model was originally specified in a discretetime binomial setting in Black et al. [1990], but subsequent research has shown the continuous-time limit of the model to be

$$
\begin{equation*}
r(t)=U(t) e^{\sigma_{r}(t) W(t)} \tag{11.1}
\end{equation*}
$$
^eq11-1

where $U(t)$ and $\sigma_{r}(t)$ are deterministic functions, and $W(t)$ is a scalar Brownian motion in the risk-neutral measure. Notice that $r(t)$ is here an outright function of $W(t)$, a property sometimes known as path independence in the short rate dynamics. The following lemma rewrites the BDT model in more familiar terms.

Lemma 11.1.1. Let $r(t)$ be given by (11.1), and let the prime denote a time derivative. Then

$$
\begin{align*}
d \ln r(t) & =\left(\vartheta(t)+\frac{\sigma_{r}^{\prime}(t)}{\sigma_{r}(t)} \ln r(t)\right) d t+\sigma_{r}(t) d W(t)  \tag{11.2}\\
d r(t) / r(t) & =\left(\vartheta(t)+\frac{1}{2} \sigma_{r}(t)^{2}+\frac{\sigma_{r}^{\prime}(t)}{\sigma_{r}(t)} \ln r(t)\right) d t+\sigma_{r}(t) d W(t)
\end{align*}
$$
^eq11-2

where

$$
\vartheta(t)=U^{\prime}(t) / U(t)-\ln (U(t)) \sigma_{r}^{\prime}(t) / \sigma_{r}(t)
$$

Proof. Set $y(t)=\ln r(t)=\ln U(t)+\sigma(t) W(t)$ and apply Ito's lemma to get

$$
\begin{aligned}
d y(t) & =\left(\frac{d \ln U(t)}{d t}+\sigma_{r}^{\prime}(t) W(t)\right) d t+\sigma_{r}(t) d W(t) \\
& =\left(U^{\prime}(t) / U(t)+(y(t)-\ln U(t)) \sigma_{r}^{\prime}(t) / \sigma_{r}(t)\right) d t+\sigma_{r}(t) d W(t)
\end{aligned}
$$

so that (11.2) follows. A second application of Ito's lemma to $r(t)=e^{y(t)}$ then gives the result for $d r(t)$.

Of the various specifications of the BDT model, the most convenient is probably (11.2) which describes the logarithm of $r(t)$ as a mean-reverting Gaussian process with a mean reversion speed

$$
\varkappa(t)=-\frac{\sigma_{r}^{\prime}(t)}{\sigma_{r}(t)}
$$

As $\ln r(t)$ is Gaussian, $r(t)$ is log-normal. In the formulation (11.2), $\vartheta(t)$ can be considered a free parameter, the value of which can - for a fixed volatility function $\sigma_{r}(t)$ - be determined by calibrating the model to the initial yield curve. As the BDT model has no known bond reconstitution formula, this fit must be done numerically. The original presentation in Black et al. [1990] outlines such a routine, based on brute-force backward induction in a binomial tree; this algorithm is, however, computationally inefficient and should not be used. For a much faster approach, see Sections 11.3.2.1 and 11.3.2.2.

Besides the lack of the bond reconstitution formula, the BDT model is plagued by a number of issues that makes it unsuitable for practical applications. One of the issues is explained in Section 11.1.3. In addition, it
is problematic that the mean reversion speed of the model is beyond user control and is fully determined from the short rate volatility and its time derivative. In particular, for those values of $t$ where $\sigma_{r}(t)$ grows in $t$, the mean reversion will be negative, i.e. the model will imply "mean-fleeing" behavior. It should be obvious that this feature of the model is undesirable.

### 11.1.2 Black-Karasinski Model

In order to rectify the problems surrounding the mean reversion in the BDT model, Black and Karasinski [1991] (BK) took the obvious step of introducing a model where mean reversion for $\ln r(t)$ is exogenously specified. In other words, we write

$$
\begin{equation*}
d \ln r(t)=\varkappa(t)(\vartheta(t)-\ln r(t)) d t+\sigma_{r}(t) d W(t) \tag{11.3}
\end{equation*}
$$
^eq11-3

Equivalently, we may write

$$
r(t)=e^{x(t)}
$$

where $x(t)$ is a standard mean-reverting Gaussian process; accordingly, the BK short rate dynamics are straightforward to simulate in the Monte Carlo method. The BK dynamics generalize and improve those of the BDT model but still do not allow for a closed-form discount bond reconstitution formula.

### 11.1.3 Issues in Log-Normal Models

The BK model - and its special case, the BDT model - have short rates that are log-normally distributed. A similar distribution of rates would arise for risk-neutral dynamics of geometric Brownian motion type

$$
\begin{equation*}
d r(t)=\mu_{r}(t) r(t) d t+\sigma_{r}(t) r(t) d W(t) \tag{11.4}
\end{equation*}
$$
^eq11-4

A time-homogeneous version of this model was considered in Rendleman and Bartter [1980]; for the time-homogeneous case (very complicated) formulas ${ }^{1}$ for discount bond prices exist, see Dothan [1978] and Hogan and Weintraub [1993]. The model (11.4) has no mean reversion, and cannot be recommended for practical applications, however.

A common problem shared by all the log-normal short rate models (11.2), (11.3), (11.4) is the fact that the expected value of the inverse of future discount bond price is infinite, i.e.

$$
\begin{equation*}
\mathrm{E}_{t}\left(\frac{1}{P\left(t^{\prime}, T\right)}\right)=\infty, \quad t<t^{\prime}<T \tag{11.5}
\end{equation*}
$$
^eq11-5

[^12]This result is formally shown in Hogan and Weintraub [1993], but is hardly surprising since the expectation of $e^{c X}, c>0$, is well-known to be infinite when $X$ is a log-normal random variable ${ }^{2}$. An important corollary of this result is listed below, originally due to Sandmann and Sondermann [1997].

Corollary 11.1.2. Define a forward Libor rate

$$
L(t, T)=\frac{1}{\tau}\left(\frac{P(t, T)}{P(t, T+\tau)}-1\right),
$$

where $\tau>0$ is some accrual factor. Assuming that (11.5) holds, then also

$$
\begin{equation*}
\mathrm{E}_{t}(L(T, T))=\infty, \quad T>t \tag{11.6}
\end{equation*}
$$
^eq11-6

and

$$
\begin{equation*}
\mathrm{E}_{t}\left(e^{\int_{t^{\prime}}^{\prime T} r(u) d u}\right)=\infty \tag{11.7}
\end{equation*}
$$
^eq11-7

Proof. As $1+\tau L(T, T)=1 / P(T, T+\tau)$, equation (11.6) follows directly from (11.5). To show (11.7), we use Jensen's inequality ${ }^{3}$ to show

$$
\begin{equation*}
\frac{1}{P\left(t^{\prime}, T\right)}=\frac{1}{\mathrm{E}_{l^{\prime}}\left(e^{-\int_{t^{\prime}}^{\prime T^{\prime}} r(u) d u}\right)} \leq \mathrm{E}_{t^{\prime}}\left(e^{\int_{t^{\prime}}^{\prime \prime} r(u) d u}\right) \tag{11.8}
\end{equation*}
$$
^eq11-8

Taking expectations conditional on the time $t$ filtration yields (11.7).
Both (11.6) and (11.7) have unfortunate economic consequences. Formula (11.7) predicts that the expected return of investing in the continuously compounded money market account for a finite period of time is infinite; and (11.6) predicts that all Libor futures rates should be infinite ${ }^{4}$. A related problem was discussed in the context of log-normal HJM models in Section 4.5.3.

### 11.1.4 Sandmann-Sondermann Transformation

The problems outlined in Corollary 11.1.2 are a significant drawback of log-normal short rate models, and one that should disqualify their use in many applications. On the other hand, market data may dictate that interest rates should, in fact, be log-normal. This is not as big a dilemma as it may appear, as there are ways to build models with a strong log-normal flair, yet avoiding (11.5). One way is to use HJM models of the type

$$
d f(t, T)=O(d t)+\sigma(t) r(t) d W(t)
$$

[^13]where $f$ is, as always, the instantaneous forward rate. We return to this type of models in Chapter 13. An interesting alternative was proposed by Sandmann and Sondermann [1997] and in essence involves writing the lognormal dynamics for a discretely compounded rate $r_{d}(t)$, rather than the infinitesimal (continuously compounded) rate $r(t)$. Specifically, we relate $r_{d}$ and $r$ through the expression

$$
\begin{equation*}
e^{r(t) \delta}=1+r_{d}(t) \delta \quad \Longrightarrow \quad r(t)=\delta^{-1} \ln \left(1+r_{d}(t) \delta\right), \tag{11.9}
\end{equation*}
$$
^eq11-9
where $\delta>0$ is some finite compounding interval. Sandmann and Sondermann [1997] set $\delta=1$, but any finite positive value will, in fact, do.

The effect of working with discretely compounded rates is summarized in the following result.

Proposition 11.1.3. Let $r_{d}(t)$ be log-normal for all $t$. Then, for $t<t^{\prime}<T$,

$$
\begin{align*}
\mathrm{E}_{t}(L(T, T)) & <\infty,  \tag{11.10}\\
\mathrm{E}_{t}\left(e^{\int_{t^{\prime}}^{T} r(u) d u}\right) & <\infty \tag{11.11}
\end{align*}
$$
^eq11-10
^eq11-11

Proof. As the proof of Proposition 11.1.3 is quite instructive, we give full details. From (11.8), to prove both (11.10) and (11.11) it suffices to show that the expectation of $e^{\int^{\prime} \prime^{\prime} r(u) d u}$ is finite. For this, let us recall that Jensen's inequality for integrals states that for a real-valued function $g$ and a concave function $\varphi$

$$
\begin{equation*}
\varphi\left(\int_{t^{\prime}}^{T} g(u) f(u) d u\right) \geq \int_{t^{\prime}}^{T} \varphi(g(u)) f(u) d u \tag{11.12}
\end{equation*}
$$
^eq11-12

provided that $f(u)$ is non-negative and $\int_{t^{\prime}}^{T} f(u) d u=1$. If $\varphi$ is convex, the inequality is reversed. Now write

$$
\delta^{-1} \int_{t^{\prime}}^{T} \ln \left(1+r_{d}(u) \delta\right) d u=\delta^{-1} \int_{t^{\prime}}^{T} \frac{1}{T-t^{\prime}} \ln \left(\left(1+r_{d}(u) \delta\right)^{T-t^{\prime}}\right) d u
$$

and apply (11.12) to the right-hand side with $f(u)=1 /\left(T-t^{\prime}\right), \varphi(u)=\ln (u)$, and $g(u)=\left(1+r_{d}(u) \delta\right)^{T-t^{\prime}}$ to show that

$$
\begin{equation*}
\delta^{-1} \ln \left(\int_{t}^{T} \frac{1}{T-t^{\prime}}\left(1+r_{d}(u) \delta\right)^{T-t^{\prime}} d u\right) \geq \delta^{-1} \int_{t^{\prime}}^{T} \ln \left(1+r_{d}(u) \delta\right) d u \tag{11.13}
\end{equation*}
$$
^eq11-13

From (11.13) it then follows immediately that

$$
\begin{aligned}
& \mathrm{E}_{t}\left(\exp \left(\int_{t^{\prime}}^{T} r(u) d u\right)\right) \\
& \quad=\mathrm{E}_{t}\left(\exp \left(\int_{t^{\prime}}^{T} \delta^{-1} \ln \left(1+r_{d}(u) \delta\right) d u\right)\right) \\
& \quad \leq \mathrm{E}_{t}\left(\exp \left(\delta^{-1} \ln \left(\int_{t^{\prime}}^{T} \frac{1}{T-t^{\prime}}\left(1+r_{d}(u) \delta\right)^{T-t^{\prime}} d u\right)\right)\right) \\
& \quad=\mathrm{E}_{t}\left(\left(\frac{1}{T-t^{\prime}} \int_{t^{\prime}}^{T}\left(1+r_{d}(u) \delta\right)^{T-t^{\prime}} d u\right)^{1 / \delta}\right)
\end{aligned}
$$

Assume that $0<\delta<1$, and set $f(u)=1 /\left(T-t^{\prime}\right), \varphi(u)=u^{1 / \delta}$, and $g(u)=\left(1+r_{d}(u) \delta\right)^{T-t^{\prime}}$. By Jensen's inequality (here (11.12) is reversed, since $\varphi$ is now convex) we must have

$$
\begin{equation*}
\left(\frac{1}{T-t^{\prime}} \int_{t^{\prime}}^{T}\left(1+r_{d}(u) \delta\right)^{T-t^{\prime}} d u\right)^{1 / \delta} \leq \frac{1}{T-t^{\prime}} \int_{t^{\prime}}^{T}\left(1+r_{d}(u) \delta\right)^{\left(T-t^{\prime}\right) / \delta} d u \tag{11.14}
\end{equation*}
$$
^eq11-14

Therefore

$$
\begin{equation*}
\mathrm{E}_{t}\left(e^{\int_{t^{\prime}}^{T} r(u) d u}\right) \leq \frac{1}{T-t^{\prime}} \mathrm{E}_{t}\left(\int_{t^{\prime}}^{T}\left(1+r_{d}(u) \delta\right)^{\left(T-t^{\prime}\right) / \delta} d u\right), \quad 0<\delta<1 \tag{11.15}
\end{equation*}
$$
^eq11-15

Since finite powers of log-normal random variables have finite expected values, (11.11) has been shown for $0<\delta<1$. For $\delta \geq 1$ we have $1 / \delta \leq 1$, and

$$
\begin{aligned}
& \mathrm{E}_{t}\left(\left(\frac{1}{T-t^{\prime}} \int_{t^{\prime}}^{T}\left(1+r_{d}(u) \delta\right)^{T-t^{\prime}} d u\right)^{1 / \delta}\right) \\
& \quad \leq\left(\mathrm{E}_{t}\left(\frac{1}{T-t^{\prime}} \int_{t^{\prime}}^{T}\left(1+r_{d}(u) \delta\right)^{T-t^{\prime}} d u\right)\right)^{1 / \delta}
\end{aligned}
$$

and (11.11) follows from the same arguments.
Comparison of Corollary 11.1.2 and Proposition 11.1.3 shows that models of the BK type

$$
\begin{equation*}
d \ln r_{d}(t)=\varkappa(t)\left(\vartheta(t)-\ln r_{d}(t)\right) d t+\sigma(t) d W(t) \tag{11.16}
\end{equation*}
$$
^eq11-16

and geometric Brownian motion models

$$
\begin{equation*}
d r_{d}(t) / r_{d}(t)=\mu(t) d t+\sigma(t) d W(t) \tag{11.17}
\end{equation*}
$$
^eq11-17

become significantly more reasonable when the dynamics are written in $r_{d}(t)$, rather than $r(t)$. We invite the reader to apply Ito's lemma to the
transformation (11.9) to uncover the $r$-dynamics for the models (11.16)(11.17).

Remark 11.1.4. When applied to the HJM model class, the "trick" of shifting from continuously compounded to discretely compounded rates lays the foundation for the class of so-called Libor market models. We return to these models in Chapters 14 and 15.

### 11.2 Other Short Rate Models

### 11.2.1 Power-Type Models and Empirical Model Estimation

A natural extension of the Gaussian and affine short rate SDEs involves retaining the linear mean-reverting drift term of these models, but using a general power function in the diffusion term. That is, we write

$$
\begin{equation*}
d r(t)=\varkappa(t)(\vartheta(t)-r(t)) d t+\sigma(t) r(t)^{p} d W(t), \quad p>0 \tag{11.18}
\end{equation*}
$$
^eq11-18

The time-inhomogeneous Gaussian and CIR models correspond to the choices $p=0$ and $p=1 / 2$, respectively. The special case of (11.18) where $p=1$ was suggested in Brennan and Schwartz [1980] and Courtadon [1982], and is quite similar to the BK model - to the extent that the case $p=1$ shares ${ }^{5}$ the unfortunate properties of the BK model listed in Corollary 11.1.2.

The general model (11.18) is similar to the CEV model described in Chapter 7, and manipulation of the parameter $p$ may allow for a better fit of the model to observed volatility smiles in interest rate options. Due to its intractability - no bond reconstitution formula exists for $p \notin\{0,1 / 2\}$ the model is, however, rarely used in derivatives pricing applications. (For related, but significantly more tractable, models with power-type diffusion terms, see Chapter 13.) Starting with the influential article by Chan et al. [1992] (CKLS), the specification (11.18) has, however, been quite popular in econometric work. As the CKLS paper is one of the most frequently cited references on estimation of one-factor short rate models, let us make a (very) brief foray into econometrics, to review the CKLS conclusions and some of the criticism their work has subsequently drawn.

The CKLS estimation procedure is based on US Treasury bond data from the period 1964-1989. Assuming that $r(t)$ can be approximated by the one-month yield on US Treasury bonds, they estimate eight models with parameter restrictions, and one model with no restrictions. Generally speaking, the empirical results indicate that the value of the parameter $p$ is the most important in determining whether a model is accepted or rejected.

[^14]The unrestricted estimate of $p$ is close to 1.5 , and values less than around $p=1$ are rejected in their tests. The Vasicek and CIR models are, for instance, rejected, whereas the Brennan-Schwarz/Courtadon model (with $p=1$ ) is accepted.

The fact that the CKLS estimates suggest that $p \geq 1$ is surprising in light of the generally downward-sloping volatility skews in fixed income derivatives, and also raises considerable questions about model regularity, as one would expect from Corollary 11.1.2. Indeed, as shown in Honore [1998b], $r(t)$ will a.s. explode (i.e. reach $\infty$ ) in finite time if $p>1$. Beyond this, we notice that the CKLS analysis has received criticism on a number of procedural and data-related fronts. For instance, Honore [1998a] (and quite a few others) point out that the 1 month Treasury yield may be an unreliable proxy for the short rate. Repeating the analysis with a carefully computed value of $r(t)$ (obtained by exploiting the fact that a one-factor model implies a one-to-one correspondence between any discount yield and the short rate), Honore revises the CKLS estimate for $p$ significantly downward, to around $p=0.8$. Bliss and Smith [1997] also point out that the data used by CKLS covers the period October 1979 - September 1982, when the US Federal Reserve followed unusual monetary policies ("The Fed Experiment"). Properly accounting for this, Bliss and Smith revise the CKLS estimate for $p$ down to around 1.0. Applying different estimation methods and different data sets, Andersen and Lund [1997] and Christensen et al. [2001] estimate $p$ to around 0.7 and 0.8 , respectively.

Moving away from observations of only a short rate proxy, Gibbons and Ramaswamy [1993] test the ability of the CIR model to simultaneously describe the evolution in four zero-coupon rates; with data covering the same period as the CKLS study, they accept the hypothesis $p=1 / 2$. Finally, to muddy waters even further, Ait-Shalia [1996] points out (in an analysis that has subsequently been criticized as lacking robustness in Chapman and Pearson [2000]) that the simple linear drift term in (11.18) is fundamentally misspecified and should be adjusted to include non-linear terms such as $1 / r(t)$ and $r(t)^{2}$.

By now, it should be clear to the reader that the problem of estimating short rate models is not close to being conclusively solved, despite an impressively long list of papers associated to it. In much contemporary empirical research, the importance of the choice of $p$ is generally downplayed, with the affine class ( $p=1 / 2$ ) enjoying considerable current popularity due to its analytical tractability.

### 11.2.2 The Black Shadow Rate Model

As described in Chapter 10, one drawback of the Gaussian short rate model class is the implication that interest rates can become negative with positive probability. As long as investors can stuff their mattresses with currency, (nominal) interest rates must, however, always remain non-negative
to preclude arbitrage. In practice, the probability of negative rates may be small enough to ignore, but as argued in Rogers [1996] prices of certain contingent claims may be highly sensitive to even a remote probability of negative rates, in which case the Gaussian model should obviously be avoided. Possible model alternatives with non-negative rates include the log-normal models in Section 11.1, as well as the affine models in Section 10.2.

Rather than altering the underlying model, an alternative "fix" of the Gaussian model involves simply taking the positive part of the Gaussian short rate process, i.e. writing

$$
\begin{equation*}
r(t)=\left(r^{*}(t)\right)^{+} \tag{11.19}
\end{equation*}
$$
^eq11-19

where $r^{*}(t)$ is a Gaussian process

$$
\begin{equation*}
d r^{*}(t)=\varkappa(t)\left(\vartheta(t)-r^{*}(t)\right) d t+\sigma_{r}(t) d W(t) \tag{11.20}
\end{equation*}
$$
^eq11-20

This approach was first proposed in Black [1995] and Rogers [1995]. The form of (11.19) suggests an analogy where the interest rate $r(t)$ is an option, granting a choice between an underlying shadow short rate $r^{*}(t)$ and zero. In other words, whenever an interest rate product has a negative rate, we invest our money in currency instead.

The truncation in (11.19) may at first glance appear rather crude. For instance, the process for $r(t)$ retains full volatility $\sigma_{r}(t)$ as long as $r(t)>0$, and can then suddenly get extinguished completely for potentially long stretches of time. In contrast, alternative models with non-negative rates such as $\mathrm{BDT}, \mathrm{BK}$, and CIR generally all imply that the interest rate volatility will gradually vanish (linearly for BDT/BK, as a square-root for CIR) as the short rate tends to zero. Interestingly, there is some empirical evidence (from the US in the 1930's and Japan in the 1990's) that suggests that very low interest rates are often accompanied by higher absolute rate volatility than standard models would predict, see Goldstein and Keirstead [1997]. Such evidence may lend some credibility to models such as (11.19).

The process (11.19) is not analytically tractable, and numerical methods (such as those of Section 11.3) must be applied to price discount bonds and other fixed income securities. For the case where all parameters in (11.20) are constants - i.e. $r^{*}(t)$ follows a Vasicek model - Gorovoi and Linetsky [2004] list ${ }^{6}$ a complicated eigenfunction expansion series for discount bond prices. Of course, as the constant-parameter model will not be able to match the current yield curve, the result in Gorovoi and Linetsky [2004] has limited uses in practical applications.

Finally, the reader may very well ask whether perhaps (11.19) could be replaced by a reflecting or absorbing barrier at zero, or by the application of

[^15]a suitable transformation such as $r(t)=r^{*}(t)^{2}$. The latter idea was discussed in Section 10.2.6 and the former is investigated in Goldstein and Keirstead [1997] where eigenfunction expansions are derived for the time-homogeneous case. In Black [1995], the author objects to reflecting barriers on economic grounds.

### 11.2.3 Spanned and Unspanned Stochastic Volatility: the Fong and Vasicek Model

In previous chapters, we demonstrated how to incorporate stochastic volatility as a mechanism to model the volatility smile in vanilla models. One wonders how to proceed with such a construction for term structure models. We postpone much of this discussion to later chapters, but shall here take the opportunity to discuss what constitutes a true stochastic volatility model for interest rate evolution, and why the short rate framework is not particularly amendable to stochastic volatility extensions. More specifically, we shall introduce the notion of spanned and unspanned stochastic volatility. For concreteness, our discussion focuses on a model proposed by Fong and Vasicek [1991] which, despite initial appearances, is in fact not a true stochastic volatility model.

The Fong-Vasicek (FV) model is characterized by risk-neutral SDEs

$$
\begin{aligned}
& d r(t)=\varkappa_{r}\left(\vartheta_{r}-r(t)\right) d t+\sqrt{z(t)} d W_{1}(t) \\
& d z(t)=\varkappa_{z}\left(\vartheta_{z}-z(t)\right) d t+\eta \sqrt{z(t)} d W_{2}(t)
\end{aligned}
$$

where $W_{1}$ and $W_{2}$ are correlated Brownian motions, $\left\langle d W_{1}(t), d W_{2}(t)\right\rangle=\rho d t$, and $\varkappa_{r}, \vartheta_{r}, \varkappa_{z}, \vartheta_{z}, \eta$ are positive constants. We recognize the FV model as being essentially a time-homogeneous GSR model augmented by a stochastic variance process $z(t)$; the process for $z(t)$ is identical to that of the Heston model (see Chapter 8). Bond prices in the FV model can be shown to satisfy

$$
\begin{equation*}
P(t, t+\delta)=e^{A(\delta)+r(t) B(\delta)+z(t) C(\delta)} \tag{11.21}
\end{equation*}
$$
^eq11-21

for deterministic functions $A, B, C$ satisfying a coupled system of ODEs. Details about these ODEs ${ }^{7}$ and their rather complicated analytical solution can be found in Selby and Strickland [1995]. For our purposes here, what matters is not the precise form of $A, B$, and $C$, but rather the fact that $P(t, T)$ in the FV model is a deterministic function of the two state variables $r(t)$ and $z(t)$. As a consequence, one can theoretically hedge out exposure to both $r(t)$ and $z(t)$ by simply taking positions in two discount bonds with different maturities. Equivalently, given observations at time $t$ of the prices of two discount bonds with different maturities, we can invert (11.21) to uncover the current values of the two variables $r(t)$ and $z(t)$.

[^16]When a "stochastic volatility" variable $z(t)$ can be hedged by positions in discount bonds, we say that $z(t)$ is spanned by the discount curve. If all stochastic volatility variables are spanned by the discount curve, moves in, say, implied volatilities of caps and swaptions would always be accompanied by moves in the yield curve, making vega hedging theoretically unnecessary. In reality, however, there is much evidence that interest rate option volatilities cannot be perfectly hedged by trading only discount bonds; see Casassus et al. [2005] for a review of the literature. Formally, this implies that the volatilities of discount bond prices depend on a vector of random state variables $\left(z_{1}(t), \ldots, z_{n}(t)\right)$ that are not included in the state variables used in reconstitution formulas for the discount curve. That is,

$$
\begin{equation*}
\frac{\partial P(t, T)}{\partial z_{i}(t)}=0, \quad i=1, \ldots, n, \tag{11.22}
\end{equation*}
$$
^eq11-22

vet, for $Q(t, T) \equiv \mathrm{E}_{t}\left(d(P(t, T))^{2}\right) / d t$,

$$
\begin{equation*}
\frac{\partial Q(t, T)}{\partial z_{i}(t)} \neq 0, \quad i=1, \ldots, n \tag{11.23}
\end{equation*}
$$
^eq11-23

Random variables satisfying (11.22)-(11.23) are said to represent unspanned stochastic volatility (USV). Whenever we talk about true stochastic volatility in this book, we always refer to models with USV, i.e. to models that prescribe moves in rate volatility that cannot be inferred from moves in the level and the shape of the discount curve.

A detailed, and highly recommended, account of USV can be found in Collin-Dufresne and Goldstein [2002b]. Among many results, the authors prove that in a time-homogeneous setting, bivariate affine models for the term structure of interest rates cannot exhibit USV. This explains our results for the FV model, and also demonstrates that several other classical stochastic-volatility models (e.g. Longstaff and Schwartz [1992]) are not in the USV class. Further analysis of a number of stochastic volatility models proposed in the literature (a surprising number of which do not, in fact, allow for USV) can be found in Collin-Dufresne and Goldstein [2002b] and Casassus et al. [2005].

### 11.3 Numerical Methods for General One-Factor Short Rate Models

The models described in Sections 11.1 and 11.2 all (with the exception of the Fong-Vasicek model) involve risk-neutral SDEs for the short rate of the type

$$
\begin{equation*}
d r(t)=\mu_{r}(t, r(t)) d t+\sigma_{r}(t, r(t)) d W(t) \tag{11.24}
\end{equation*}
$$
^eq11-24

for certain user-specified deterministic functions $\mu_{r}(t, r)$ and $\sigma_{r}(t, r)$. As discussed, many of these models do not allow for closed-form expressions
relating discount bonds to the state of $r(t)$, necessitating the application of numerical methods for even simple tasks such as calibrating the model to the initial yield curve. While this issue ultimately should make one pause when it comes to deciding whether a model is suited for practical applications, we still want to cover some techniques that are useful in handling generic models such as (11.24). Some of the techniques we shall discuss, e.g. calibration through forward induction, have broad applicability.

### 11.3.1 Finite Difference Methods

As always, we set $x(t)=r(t)-f(0, t)$, such that the generic pricing PDE for a derivative $V=V(t, x)$ becomes

$$
\begin{equation*}
\frac{\partial V}{\partial t}+\mu_{x}(t, x) \frac{\partial V}{\partial x}+\frac{1}{2} \sigma_{x}(t, x)^{2} \frac{\partial^{2} V}{\partial x^{2}}=(x+f(0, t)) V, \tag{11.25}
\end{equation*}
$$
^eq11-25

where we have defined

$$
\begin{equation*}
\mu_{x}(t, x)=\mu_{r}(t, x+f(0, t))-\frac{\partial f(0, t)}{\partial t}, \quad \sigma_{x}(t, x)=\sigma_{r}(t, x+f(0, t)) \tag{11.26}
\end{equation*}
$$
^eq11-26

Given terminal conditions for $V(T, x)$, as well as suitable boundary conditions in the $x$-domain, we can solve this equation by the generic finite difference methods of Chapter 2. For later use, let us quickly recall that a standard $\theta$-method finite difference scheme on an equidistant $x$-grid $\left\{x_{j}\right\}_{j=0}^{m+1}$ would result in a matrix scheme of the type

$$
\begin{align*}
& {\left[\mathbf{I}-\theta \Delta_{t} \mathbf{A}\left((1-\theta) t_{i+1}+\theta t_{i}\right)\right] \widehat{\mathbf{V}}\left(t_{i}\right)=} \\
& \quad\left[\mathbf{I}+(1-\theta) \Delta_{t} \mathbf{A}\left((1-\theta) t_{i+1}+\theta t_{i}\right)\right] \widehat{\mathbf{V}}\left(t_{i+1}\right)+\mathbf{B}\left(t_{i}, t_{i+1}\right) \tag{11.27}
\end{align*}
$$
^eq11-27

where $\widehat{\mathrm{V}}(t)=\left(\widehat{V}\left(t, x_{1}\right), \ldots, \widehat{V}\left(t, x_{m}\right)\right)^{\top}$ with $\widehat{V}(t, x)$ denoting the approximation to the true solution $V(t, x), \mathbf{A}(t)$ is a tri-diagonal matrix, and $\mathrm{B}\left(t_{i}, t_{i+1}\right)$ is a vector representing any boundary conditions that cannot be folded into the matrix $\mathbf{A}$. We solve the matrix system (11.27) backward in time, starting from a given value of $\hat{\mathrm{V}}(T)$. Determination of the spatial boundaries in the $x$-domain - that is, the values of $x_{0}$ and $x_{m+1}$ - can, as always, be set by probabilistic arguments, through estimation of the first and second moment of $x(T)$. While these estimates should be targeted to the specifics of the model at hand, if all else fails we can always rely on Gaussian estimate, e.g. something like

$$
x_{0}=\mu_{x}(0,0) T-\alpha \sigma_{x}(0,0) \sqrt{T}, \quad x_{m+1}=\mu_{x}(0,0) T+\alpha \sigma_{x}(0,0) \sqrt{T},
$$

where $\alpha$ is some confidence interval multiplier (e.g. 4 or 5 ).
The discussion in Section 10.1.5.2 about $x$-domain boundary conditions apply to (11.27) as well, but determining the terminal condition $\widehat{\mathrm{V}}(T)$ can
be problematic ${ }^{8}$, as option payouts will often involve several discount factors at time $T$ (e.g. to price a swap or to compute a Libor rate). As there is generally no way of computing such discount bonds analytically for the model (11.24), we are forced to compute the discount bond prices themselves by finite difference methods. To compute the value of the discount bond $P\left(T, T^{*}, x\right), T^{*}>T$, at time $T$ we:

1. Extend the finite difference grid to time $T^{*}$.
2. Set the boundary condition $P\left(T^{*}, T^{*}, x_{j}\right)=1, j=0, \ldots, m+1$.
3. On some suitable time grid, use (11.27) on $P$ to step backward to time $T$.
4. Use the finite difference estimates $\widehat{P}\left(T, T^{*}, x_{j}\right), j=0, \ldots, m+1$, to fill in $\widehat{\mathrm{V}}(T)$.

To the extent that $\mathrm{V}(T)$ involves multiple discount bonds, we perform the algorithm above on all of the required discount bonds; the grid must then be extended to the maturity of the longest-dated discount bond needed in the payout computation. In some cases this can dramatically increase computation times relative to models where closed-form discount bond reconstitution formulas exist. For instance, for a 3 month option on a 30 year swap, a model with an analytical formula for discount bonds (e.g. the affine short rate model) would require us only to build the finite difference grid out to 3 months; when such a formula does not exist, we are forced to use a 30 year finite difference grid.

### 11.3.2 Calibration to Initial Yield Curve

Assume that the volatility function $\sigma_{r}(t, r)$ has been fixed, but $\mu_{r}(t, r)= \mu_{r}(t, r ; \vartheta(t))$ has a free time-dependent parameter $\vartheta(t)$, the value of which we wish to set in such a way that the initial discount curve $P(0, T), T>0$, is correctly recovered by the model. Suppose that we assume, as is common, that $\vartheta(t)$ is piecewise constant on some time grid $0=t_{0}<t_{1}<\ldots<t_{N}$, with $\vartheta_{i-1}$ denoting the value of $\vartheta$ that applies on the interval $\left[t_{i-1}, t_{i}\right)$. A brute force approach to the calibration of $\vartheta(t)$ could proceed as follows.

1. Assuming that $\vartheta_{0}, \ldots, \vartheta_{i-2}$ are known, make a guess ${ }^{9}$ for $\vartheta_{i-1}$.
2. Setting the terminal boundary value to $V\left(t_{i}, x\right)=1$ for all $x$, use the backward finite difference grid algorithm (11.27) to compute the value of $P\left(0, t_{i}\right)$.
3. If the computed value of $P\left(0, t_{i}\right)$ equals that quoted in the market stop; otherwise return to Step 1.
[^17]This approach is similar to an algorithm proposed in Black et al. [1990] (albeit the authors worked only with binomial trees) and involves very high computational costs as the numerical search for each of the parameters $\vartheta_{0}, \vartheta_{1}, \ldots$ involves solving a full finite difference grid in each loop. Specifically, assume that on average each search for $\vartheta_{i}$ involves $M$ iterations over steps $1-3$ above. With $N \vartheta$-values to find, the computational effort of a finite difference grid with $m$ spatial steps (say) is $O(N m)$; it follows that the total computational cost for the calibration of $\vartheta(t)$ is

$$
O\left(M N^{2} m\right),
$$

which is often prohibitively expensive in practice.

### 11.3.2.1 Forward Induction

In the setting of binomial and trinomial trees, the brute-force BDT algorithm was markedly improved upon in Jamshidian [1991a], using a technique known as forward induction. The basic idea is to work with a forward equation, rather than with the backward equation (11.27). Two varieties of this approach are feasible in a finite difference setting, depending on whether the forward equation is developed by direct discretization of the continuous-time Fokker-Planck equation, or by rearrangement of a discretized backward equation. We cover the former approach here, and the latter in Section 11.3.2.2 below.

Let $G(t, x ; s, y), s \geq t$, denote the time $t$ value of a security that pays out a Dirac delta amount iff $x(s)=y$, given that $x(t)=x$. Clearly then

$$
\begin{equation*}
P(0, T)=\int_{-\infty}^{\infty} G(0,0 ; T, y) d y \tag{11.28}
\end{equation*}
$$
^eq11-28

We already saw this financial contract, the so-called Arrow-Debreu security, at the end of Section 1.8. Its price $G$ - being the value of a perfectly valid derivative contract - certainly satisfies a backward Kolmogorov equation

$$
\begin{equation*}
\frac{\partial G}{\partial t}+\mu_{x}(t, x ; \vartheta(t)) \frac{\partial G}{\partial x}+\frac{1}{2} \sigma_{x}(t, x)^{2} \frac{\partial^{2} G}{\partial x^{2}}=(x+f(0, t)) G, \quad(s, y) \text { fixed } \tag{11.29}
\end{equation*}
$$
^eq11-29

with the terminal value condition $G(s, y ; s, y)=\delta(y)$. Clearly $G(t, z ; s, y)$ is closely related to the transition density $p(t, z ; s, y)$ defined in Section 1.8, and can be expected to satisfy a forward Kolmogorov equation, too. The correct equation is

$$
\begin{align*}
-\frac{\partial G}{\partial s}-\frac{\partial}{\partial y}\left(\mu_{x}(s, y ; \vartheta(s)) G\right)+ & \frac{1}{2} \frac{\partial^{2}}{\partial y^{2}}\left(\sigma_{x}(t, y)^{2} G\right) \\
& =(y+f(0, t)) G, \quad(t, x) \text { fixed } \tag{11.30}
\end{align*}
$$
^eq11-30

subject to the initial condition $G(t, x ; t, x)=\delta(x)$. This PDE is identical to the Fokker-Planck equation listed in Section 1.8, except for the fact that the right-hand side is not zero, but contains a discounting term.

Fixing $(t, x)=(0,0)$, the PDE (11.30) can be discretized by finite difference methods in standard fashion, although we keep in mind that the PDE is to be solved forward in time from its initial (Dirac) condition, rather than backward. In most cases ${ }^{10}$, the appropriate spatial boundary conditions for the finite difference solution $\widehat{G}$ are

$$
\widehat{G}\left(0,0 ; s, y_{0}\right)=\widehat{G}\left(0,0 ; s, y_{m+1}\right)=0,
$$

which assumes that $y_{0}$ and $y_{m+1}$ have been set such that the probability density at these levels is negligible. In a discrete setting, the initial condition $G(t, x ; t, x)=\delta(x)$ is translated to $\widehat{G}(0,0 ; 0, y)=(\Delta y)^{-1} 1_{\{y=0\}}$ for $\Delta y$ suitably defined, in agreement with the averaging principles of Section 2.5.2. Due to the strongly discontinuous initial condition, Rannacher stepping (see Section 2.5) should always be used.

We are now ready to state our revised algorithm for calibration of $\vartheta(t)$, working again with the assumption that $\vartheta(t)$ is piecewise constant on a time grid $0=t_{0}<t_{1}<\ldots<t_{N}$. We assume that $\vartheta_{0}, \ldots, \vartheta_{i-2}$ have been found, as has been $\widehat{G}\left(0,0 ; t_{i-1}, y\right)$ for all $y \in\left\{y_{j}\right\}_{j=1}^{m}$ in the finite difference grid.

1. Make a guess for $\vartheta_{i-1}$.
2. Solve (11.30) one time step forward, to time $t_{\imath}$, and save $\widehat{G}\left(0,0 ; t_{i}, y_{j}\right)$, $j=1, \ldots, m$.
3. Compute the discount bond price $P\left(0, t_{i}\right)=\sum_{j} \widehat{G}\left(0,0 ; t_{i}, y_{j}\right)\left(y_{j}-y_{j-1}\right)$.
4. If the computed value of $P\left(0, t_{i}\right)$ equals that quoted in the market stop; otherwise return to Step 1.

The cost of Steps 2 and 3 in this algorithm is $O(m)$, where $m$ is the number of points in the $y$-direction of the finite difference grid. The total computational effort of this algorithm is therefore

$$
O(M m N),
$$

where $M$ is the average number of root search iterations over Steps 1-4 above. We recall that the effort of the brute-force backward equation approach was $O\left(M m N^{2}\right)$, so the use of forward induction saves us a factor of $N$. As $N$ is often of the order of $N=100$, these savings can be considerable. In typical applications, $M$ is often in the order $M=2$ to 4 , so calibration of the model (11.24) to the initial discount curve should only be a few times slower than pricing by finite difference methods a single option maturing at time $t_{N}$ (the cost of which we recall to be $O(m N))$.

[^18]
### 11.3.2.2 Forward-from-Backward Induction

The backward and forward Kolmogorov equations (11.29)-(11.30) are consistent in the continuous-time limit, but not necessarily so when discretized, finite difference style. As a result, the function $\vartheta(t)$ uncovered from the algorithm in Section 11.3.2.1 will generally not allow a finite difference grid based on the backward equation (11.25) to recover the initial term structure of discount bonds without errors, even if discretized on $t$ - and $x$-grids that are identical to those used for the calibration of $\vartheta(t)$. As long as the time line is sufficiently finely spaced, the errors tend to be very small, however, and rarely a cause for concern. Nevertheless, it should be noted that it is, in fact, possible to restate the forward induction algorithm in such a way that the algorithm becomes precisely compatible with the brute-force backward equation approach we discussed earlier.

To develop this approach, we start out with the discretized backward equation (11.27) and rearrange it to yield

$$
\begin{equation*}
\widehat{\mathrm{V}}\left(t_{i}\right)=\mathbf{T}_{i}^{2+1} \widehat{\mathrm{~V}}\left(t_{i+1}\right)+\mathbf{G}_{i}^{2+1}, \tag{11.31}
\end{equation*}
$$
^eq11-31

where, with $\mathbf{A}_{i}^{i+1} \triangleq \mathbf{A}\left((1-\theta) t_{i+1}+\theta t_{i}\right)$,

$$
\begin{aligned}
\mathbf{T}_{i}^{i+1} & \triangleq\left[\mathbf{I}-\theta \Delta_{t} \mathbf{A}_{i}^{i+1}\right]^{-1}\left[\mathbf{I}+(1-\theta) \Delta_{t} \mathbf{A}_{i}^{i+1}\right] \\
& =\mathbf{I}+\left[\mathbf{I}-\theta \Delta_{t} \mathbf{A}_{i}^{i+1}\right]^{-1} \mathbf{A}_{i}^{i+1} \Delta_{t}
\end{aligned}
$$

and

$$
\mathbf{G}_{i}^{i+1}=\left[\mathbf{I}-\theta \Delta_{t} \mathbf{A}_{i}^{i+1}\right]^{-1} \mathbf{B}\left(t_{i}, t_{i+1}\right) .
$$

Repeated application of (11.31) yields, for some $l$,

$$
\begin{equation*}
\widehat{\mathrm{V}}(0)=\mathbf{T}_{0}^{l} \widehat{\mathrm{~V}}\left(t_{l}\right)+\mathrm{G}_{0}^{l} \tag{11.32}
\end{equation*}
$$
^eq11-32

where $\mathrm{T}_{0}^{l}$ and $\mathrm{G}_{0}^{l}$ can be found iteratively from the equations

$$
\begin{align*}
& \mathbf{T}_{0}^{i+1}=\mathbf{T}_{0}^{i} \mathbf{T}_{\imath}^{i+1}, \quad \mathbf{T}_{0}^{0}=\mathbf{I}  \tag{11.33}\\
& \mathbf{G}_{0}^{i+1}=\mathbf{T}_{0}^{i} \mathbf{G}_{i}^{i+1}+\mathbf{G}_{0}^{i}, \quad \mathbf{G}_{0}^{0}=0 . \tag{11.34}
\end{align*}
$$
^eq11-33
^eq11-34

Before we proceed, let us introduce some notation. First, we let $1_{j}$ denote an $m$-dimensional column vector with $j$-th element equal to 1 and all other elements equal to zero. Also we set 1 to mean a column vector with all elements equal to 1 . Consider now the zero-coupon bond maturing at time $t_{2}$, and assume that for the grid $\left\{x_{j}\right\}_{j=0}^{m+1}$ the initial value of $x=0$ sits in position $\beta$, i.e. $x_{\beta}=0$. By the definition of a discount bond, (11.32) allows us to write

$$
\begin{equation*}
P\left(0, t_{i}\right)=1_{\beta}^{\top}\left(\mathbf{T}_{0}^{i} 1+\mathbf{G}_{0}^{i}\right)=\mathbf{1}^{\top} \mathbf{D}_{0}^{i}+g_{0}^{i}, \tag{11.35}
\end{equation*}
$$
^eq11-35

where $\mathrm{D}_{0}^{i}$ is an $m$-dimensional vector and $g_{0}^{i}$ is a scalar:

$$
\mathbf{D}_{0}^{i}=\left(\mathbf{T}_{0}^{2}\right)^{\top} \mathbf{1}_{\beta}, \quad g_{0}^{2}=\left(\mathbf{G}_{0}^{i}\right)^{\top} \mathbf{1}_{\beta} .
$$

Assuming that the influence from the finite difference grid boundary is small, we evidently have

$$
P\left(0, t_{i}\right) \approx \sum_{j=1}^{m}\left(\mathbf{D}_{0}^{i}\right)_{j}
$$

where $\left(\mathbf{D}_{0}^{i}\right)_{j}$ denotes the $j$-th element of $\mathbf{D}_{0}^{i}$. Comparison with (11.28) shows that $\mathrm{D}_{0}^{i}$ can be interpreted as the discrete-time Arrow-Debreu security vector for maturity $t_{i}$ (up to the scaling $\Delta x$ ). From (11.33)-(11.34) we have

$$
\begin{align*}
\mathbf{D}_{0}^{2+1} & =\left(\mathbf{T}_{0}^{i} \mathbf{T}_{i}^{i+1}\right)^{\top} \mathbf{1}_{\beta}=\left(\mathbf{T}_{\imath}^{i+1}\right)^{\top} \mathbf{D}_{0}^{i}, \quad \mathbf{D}_{0}^{0}=\mathbf{1}_{\beta}  \tag{11.36}\\
g_{0}^{i+1} & =\left(\mathbf{T}_{0}^{i} \mathbf{G}_{i}^{i+1}+\mathbf{G}_{0}^{i}\right)^{\top} \mathbf{1}_{\beta}=\left(\mathbf{G}_{i}^{i+1}\right)^{\top} \mathbf{D}_{0}^{i}+g_{0}^{i}, \quad g_{0}^{0}=0 \tag{11.37}
\end{align*}
$$
^eq11-36
^eq11-37

Recalling the definition of $\mathbf{T}_{2}^{i+1}$, it follows that

$$
\begin{align*}
\mathrm{D}_{0}^{2+1} & =\left(\mathbf{T}_{i}^{i+1}\right)^{\top} \mathrm{D}_{0}^{i} \\
& =\left[\mathbf{I}+(1-\theta) \Delta_{t} \mathbf{A}_{i}^{i+1}\right]^{\top}\left[\mathbf{I}-\theta \Delta_{t}\left(\mathbf{A}_{i}^{2+1}\right)^{\top}\right]^{-1} \mathbf{D}_{0}^{i} \\
& \triangleq\left[\mathbf{I}+(1-\theta) \Delta_{t}\left(\mathbf{A}_{i}^{i+1}\right)^{\top}\right] \mathbf{Y}_{0}^{i} \tag{11.38}
\end{align*}
$$
^eq11-38

where $\mathbf{Y}_{0}^{i}$ is a vector satisfying a tri-diagonal matrix equation

$$
\begin{equation*}
\left[\mathbf{I}-\theta \Delta_{t}\left(\mathbf{A}_{i}^{i+1}\right)^{\top}\right] \mathbf{Y}_{0}^{i}=\mathbf{D}_{0}^{i} \tag{11.39}
\end{equation*}
$$
^eq11-39

With this, we are ready to state our revised fitting algorithm for $\vartheta(t)$. We assume that $\vartheta_{0}, \ldots, \vartheta_{i-2}$ have been found, as has $\mathbf{D}_{0}^{i-1}$ and $g_{0}^{i-1}$.

1. Make a guess for $\vartheta_{i-1}$.
2. Solve (11.39) for $\mathbf{Y}_{0}^{i-1}$.
3. Compute $\mathrm{D}_{0}^{i}$ from (11.38); and $g_{0}^{i}$ from (11.37).
4. Compute the discount bond price $P\left(0, t_{\imath}\right)$ from (11.35).
5. If the computed value of $P\left(0, t_{i}\right)$ equals that quoted in the market, then stop; otherwise return to Step 1.

The computational efforts for Steps 2, 3, and 4 are all $O(m)$, so the algorithm above is of the same computational complexity as the algorithm in Section 11.3.2.1.

### 11.3.2.3 Yield Curve and Volatility Calibration

Volatility calibration of a general short rate model of the type (11.24) is a rather involved affair. The typical scheme - moving model volatilities around until the prices of calibration targets match the market - is beset with complications such as

- Bond reconstitution formulas are unavailable so the model needs to be numerically recalibrated to the initial yield curve after each update of the model volatilities.
- All bonds used in the payoffs of calibration targets need to be computed numerically for each volatility update.
- Bond and calibration target values depend on the entire volatility function, making decoupling of individual target calibrations difficult.

These difficulties, however, have not prevented some major investment banks from risk-managing large derivatives portfolios with such a setup. As we cannot possibly do justice to all the tricks that would be required to make this operational, we content ourselves with presenting a mere outline of a possible algorithm.

We consider the model of the type (11.24) but, for notational convenience, write the volatility term in a separable form

$$
\begin{equation*}
d r(t)=\mu_{r}(t, r(t) ; \vartheta(t)) d t+\sigma_{r}(t) \psi(r(t)) d W(t) \tag{11.40}
\end{equation*}
$$
^eq11-40

Here, the purely-time-dependent function $\sigma_{r}(t)$ is used to calibrate to swaptions, and $\vartheta(t)$ is used to match the initial yield curve. Implicitly, $\vartheta(t)$ depends on $\sigma_{r}(t)$.

As in Section 10.1.4, we assume that we are given a collection of $N-1$ swaptions defined on a maturity grid $0=T_{0}<T_{1}<\ldots<T_{N}$ such that the $i$-th swaption expires at times $T_{i}, i=1, \ldots, N-1$. For concreteness assume that all underlying swaps mature on $T_{N}$. We further assume that $\sigma_{r}(t)$ is discretized in a piecewise constant manner on the maturity grid, with $\sigma_{i}$ denoting the flat value on $\left[T_{i}, T_{i+1}\right)$.

Before discussing calibration, let us outline an efficient algorithm for pricing all swaptions in the calibration set given a collection of volatilities $\sigma_{0}, \ldots, \sigma_{N-1}$. We implicitly assume that the model is rewritten using a state variable $x$ as in (11.25), and the $x$-domain is discretized with $\left\{x_{j}\right\}_{j=0}^{m+1}$. As the algorithm involves both forward and backward induction, it is important to follow the approach of Section 11.3.2.2 and use a forward algorithm that is fully compatible with the backward one.

1. Update the volatility function of the model with the new values $\sigma_{0}, \ldots, \sigma_{N-1}$.
2. Using the forward induction algorithm from Section 11.3.2.2, calibrate $\vartheta(t)$ for all $t \in\left[0, T_{N}\right]$.
3. On Step 2, Arrow-Debreu prices $\widehat{G}\left(0,0 ; T_{i}, \cdot\right)$ are calculated; save them for all $i=1, \ldots, N-1$.
4. For each $n=N, \ldots, 2$ :
a) Create a new copy of the finite difference grid.
b) Populate a payoff $P\left(T_{n}, T_{n}\right)=1$ at time $T_{n}$.
c) Calculate $P\left(T_{i}, T_{n}\right)$ from $P\left(T_{i+1}, T_{n}\right)$ by backward induction for $i=n-1, \ldots, 1$.
5. For each $i=1, \ldots, N-1$ :
a) Create the $T_{i}$-expiry swaption payoff from $P\left(T_{i}, T_{n}\right), n=i+1, \ldots, N$, calculated on Step 4.
b) Integrate the payoff against $\widehat{G}\left(0,0 ; T_{i}, \cdot\right)$ stored on Step 3.
c) This gives the value of the $T_{i}$-expiry swaption.

This algorithm describes how to map a set of model volatilities $\sigma_{0}, \ldots, \sigma_{N-1}$ into model prices of calibration targets. In principle, one can now perform a multi-dimensional optimization to match swaption prices from the model to the market. Given that the number of swaptions could be large - 30 or 40 would not be uncommon - and each valuation is rather expensive, the resulting algorithm, while not necessarily completely impractical, would require significant computational resources.

A further improvement entails adopting the iterative bootstrap algorithm outlined in Section 10.2.5. We recall that the main premise of this algorithm approach was that the value of the $T_{i}$-expiry swaption depended on $\sigma_{r}(s)$ for $s \in\left[0, T_{i}\right]$ in a much stronger way than on $\sigma_{r}(s)$ for $s \in\left[T_{i}, T_{N}\right]$; this tends to also be true for the model of the type (11.40) for a wide range of volatility specifications. To incorporate this observation into the algorithm above, we would work our way forward from $i=0$ and determine $\sigma_{i}$ by one-dimensional root-search to match the market value of the $T_{i+1}$-expiry swaption; all values $\sigma_{j}, j \neq i$ would be kept constant in the root search. The bootstrap loop would work its way from $i=0$ to $i=N-1$ and would be repeated a few times until convergence, in the manner described in Section 10.2.5. We trust that the reader gets the idea and can fill in remaining details.

### 11.3.2.4 The Dybvig Parameterization

In some models, it may be the case that the SDE (11.24) can be reformulated as

$$
\begin{equation*}
r(t)=s(t)+\vartheta(t) \tag{11.41}
\end{equation*}
$$
^eq11-41

where $\vartheta(t)$ is a free time-dependent function to be fitted against discount bond prices, and $s(t)$ satisfies an SDE

$$
\begin{equation*}
d s(t)=\mu_{s}(t, s(t)) d t+\sigma_{s}(t, s(t)) d W(t) \tag{11.42}
\end{equation*}
$$
^eq11-42

We notice that

$$
P(0, T)=\mathrm{E}\left(e^{-\int_{0}^{T} r(u) d u}\right)=e^{-\int_{0}^{T} \vartheta(u) d u} \mathrm{E}\left(e^{-\int_{0}^{T} s(u) d u}\right)
$$

such that

$$
\begin{equation*}
\int_{0}^{T} \vartheta(u) d u=\ln \left(\frac{\mathrm{E}\left(e^{-\int_{0}^{T} s(u) d u}\right)}{P(0, T)}\right) \tag{11.43}
\end{equation*}
$$
^eq11-43

To the extent that the numerator in the right-hand side of (11.43) is easy to compute - e.g. if the SDE for $s(t)$ permits a closed-form solution calibration of $\vartheta(t)$ can conveniently be found by direct differentiation, see (11.44).

The specification in (11.41) was proposed in Dybvig [1997] and is, as we have already seen in Section 10.1.1.2, quite natural in the context of Gaussian models. For non-Gaussian models, a "fudge" approach in (11.41) may be less desirable, as the domain of $r(t)$ is hard to control. For instance, suppose that $s(t)$ is a time-homogeneous CIR process

$$
d s(t)=\varkappa\left(s_{0}-s(t)\right) d t+\sigma \sqrt{s(t)} d W(t)
$$

which is guaranteed to produce only non-negative values of $s(t)$. The combined process ${ }^{11} r(t)=s(t)+\vartheta(t)$, however, will have domain $r(t) \in[\vartheta(t), \infty)$ which is rather awkward as $\vartheta(t)$ is largely out of the user's control. This is reflected in the SDE for $r(t)$, where now $\vartheta(t)$ enters into the volatility term:

$$
d r(t)=\varkappa\left(s_{0}+\vartheta(t)+\vartheta^{\prime}(t) / \varkappa-r(t)\right) d t+\sigma \sqrt{r(t)-\vartheta(t)} d W(t)
$$

where $\vartheta^{\prime}(t)=d \vartheta(t) / d t$. By affecting the short rate volatility, the interpretation of $\vartheta(t)$ as only serving to fit the yield curve can no longer be maintained. This conclusion holds not only for CIR models, but for all models where $\sigma_{s}$ depends on $s(t)$, since

$$
d r(t)=\mu_{s}(t, r(t)-\vartheta(t)) d t+\vartheta^{\prime}(t) d t+\sigma_{s}(t, r(t)-\vartheta(t)) d W(t)
$$

While sometimes very convenient, the Dybvig parameterization should consequently be approached with considerable care.

### 11.3.2.5 Link to HJM Models

By construction the Dybvig parameterization in Section 11.3.2.4 ensures that the model is calibrated to the initial forward curve $f(0, t), t>0$. As the resulting model is driven by Brownian motions, we know that it must be in the HJM class. An interesting question arises: what is the type of an HJM model that is defined through the Dybvig procedure. To answer this, let $s(t)$ satisfy (11.42) and define

$$
Q(t, T, s)=\mathrm{E}\left(e^{-\int_{t}^{T} s(u) d u} \mid s(t)=s\right)
$$

From (11.43), we then get[^19]

$$
\begin{equation*}
\vartheta(t)=\frac{\partial}{\partial t} \ln \left(\frac{Q(0, t, s(0)}{P(0, t)}\right) \tag{11.44}
\end{equation*}
$$
^eq11-44
and

$$
P(t, T)=e^{-\int_{t}^{T} \vartheta(u) d u} Q(t, T, s(t))=\frac{Q(0, t, s(0))}{Q(0, T, s(0))} \frac{P(0, T)}{P(0, t)} Q(t, T, s(t)) .
$$

As

$$
f(t, T)=-\frac{\partial}{\partial T} \ln P(t, T)=-\frac{\partial}{\partial T} Q(t, T, s(t))+\vartheta(T)
$$

it follows that

$$
d f(t, T)=O(d t)+\sigma_{f}(t, T) d W(t)
$$

where

$$
\begin{aligned}
\sigma_{f}(t, T) & =-\frac{\partial^{2}}{\partial T \partial s} Q(t, T, s(t)) \sigma_{s}(t, s(t)) \\
& =-\frac{\partial^{2}}{\partial T \partial s} Q(t, T, f(t, t)-\vartheta(t)) \sigma_{s}(t, f(t, t)-\vartheta(t))
\end{aligned}
$$

At time $t$, the forward rate volatility structure generated by the short rate model evidently depends on the forward curve at time $t$ (through $f(t, t)$ ) as well as the function $\vartheta(t)$. As $\vartheta(t)$ depends (through (11.44)) on the forward curve at time 0 , it follows that the HJM dynamics here generally have "memory" of the initial condition at time 0 . If one were to alter these initial conditions, the form of HJM dynamics would fundamentally change.

Looking back to Section 10.1.2.2 where the Gaussian short rate model was developed from a separable HJM model, no dependence in the HJM dynamics on initial conditions arose. Hence, it is clear that not all short rate models generate "memory" in the HJM dynamics. This raises the obvious question: under which circumstances will a finite-dimensional Markov HJM model have no dependence in its dynamics of the initial forward curve at time 0? The answer to this is listed in Filipovic and Teichmann [2004] which shows that "essentially" all such models must be time-inhomogeneous affine models. There are considerable technical details involved in the exact statement of the result, all of which can be found in Filipovic and Teichmann [2004].

### 11.3.2.6 The Hagan and Woodward Parameterization

Many of the short rate models considered so far have a free time-dependent parameter in the drift of the state variable dynamics; the forward induction algorithm can then be used to set this parameter to match the initial yield curve. Hagan and Woodward [1999a] propose an interesting twist on the idea, where the free parameter is introduced into the numeraire.

Hagan and Woodward [1999a] start with an observation that only two ingredients are required to create an interest rate model:

- A set of stochastic processes that drive the evolution of interest rates.
- A positive-valued process that is used as a numeraire.

Once the numeraire is specified, the values of all instruments are recovered by the standard pricing formula, see Chapter 1. Critically, the numeraire does not need to be the money market account, or any other "identifiable" security such as a discount bond or an annuity - a positive process is all that is required (so in fact we need not a numeraire but a deflator as defined in Section 1.3, but in this section we use the two terms interchangeably).

We define the stochastic process that drives interest rates by a general one-dimensional ${ }^{12}$ process

$$
\begin{equation*}
d x(t)=\mu_{x}(t, x(t)) d t+\sigma_{x}(t, x(t)) d W(t), \quad x(0)=0 \tag{11.45}
\end{equation*}
$$
^eq11-45

Furthermore, we choose the deflator to be a function of the state variable $x(t)$. Without loss of generality, we specify

$$
\begin{equation*}
N(t)=\frac{1}{P(0, t)} e^{h(t, x(t))+a(t)} \tag{11.46}
\end{equation*}
$$
^eq11-46

Here $h(t, x)$ is user-specified, and $a(t)$ is used to fit the model to the initial yield curve. It is often natural to normalize the parameters such that

$$
h(t, 0)=0, \quad a(0)=0
$$

Once the deflator is specified, we assume that (11.45) is, in fact, given under the measure $\mathrm{Q}^{N}$ associated with this deflator.

Let $\mathrm{E}^{N}$ denote expectations in measure $\mathrm{Q}^{N}$. The time $t$ price of a $T$ maturity discount bond in this model, as a function of the state variable $x$, is given by

$$
\begin{align*}
P(t, T, x) & =N(t) \mathrm{E}^{N}\left(N(T)^{-1} \mid x(t)=x\right) \\
& =\frac{P(0, T)}{P(0, t)} e^{h(t, x)+a(t)-a(T)} \mathrm{E}_{t}^{N}\left(e^{-h(T, x(T))}\right) . \tag{11.47}
\end{align*}
$$
^eq11-47

Consistency with the initial yield curve requires

$$
P(0, T, x(0))=P(0, T), \quad T \geq 0
$$

and we obtain the following condition on $a(T)$,

$$
\begin{equation*}
a(T)=\ln \mathrm{E}^{N}\left(e^{-h(T, x(T))}\right), \quad T \geq 0 . \tag{11.48}
\end{equation*}
$$
^eq11-48

Hagan and Woodward [1999a] show that if the model (11.45)-(11.46) is consistent with the initial yield curve, i.e. the condition (11.48) is satisfied, the model is in fact arbitrage free.

[^20]To obtain $a(t)$ from (11.48), one should not solve for $E^{N}\left(e^{-h(T, x(T))}\right)$ with a backward PDE. Instead, similarly to Section 11.3.2, one should use forward PDE to obtain $p(t, x)$, the density $\mathrm{Q}^{N}(x(t) \in d x) / d x$, for $(t, x) \in \mathbb{R}^{+} \times \mathbb{R}$. The forward Kolmogorov equation states that

$$
-\frac{\partial p}{\partial t}(t, x)-\frac{\partial}{\partial x}\left(\mu_{x}(t, x) p(t, x)\right)+\frac{1}{2} \frac{\partial^{2}}{\partial x^{2}}\left(\sigma_{x}(t, x)^{2} p(t, x)\right)=0, p(0, x)=\delta(x) .
$$

Once $p(t, x)$ is determined, we obtain

$$
a(T)=\ln \int e^{-h(T, x)} p(T, x) d x, \quad T \geq 0
$$

where the integral is taken over the range of the random variable $x(T)$. Interestingly, the calibration of $a(t)$ is independent of the initial yield curve $P(0, T), T \geq 0$. The calibration is somewhat faster than the forward induction algorithm of Section 11.3.2.1 for general short rate models as it requires only a single forward pass of the finite difference scheme.

Zero-coupon discount bonds are obtained via (11.47), i.e., generally, numerically. For special choices of $\mu_{x}(t, x), \sigma_{x}(t, x)$ and $h(t, x)$, closed-form formulas could be available.

Let us define $\gamma(t, x ; T)$ by

$$
\gamma(t, x(t) ; T)=\mathrm{E}_{t}^{N}\left(e^{-h(T, x(T))}\right)
$$

so that

$$
P(t, T, x)=\frac{P(0, T)}{P(0, t)} e^{h(t, x)+a(t)-a(T)} \gamma(t, x ; T)
$$

Then instantaneous forward rates are given by

$$
f(t, T, x)=f(0, T)+a^{\prime}(T)-\frac{\partial \ln \gamma(t, x ; T)}{\partial T}
$$

Applying Ito's lemma to $e^{-h(T, x(T))}$ and setting $T=t$, we observe that the short rate $r(t)=f(t, t, x(t))=r(t, x(t))$ is given by

$$
\begin{align*}
r(t, x)= & f(0, t)+a^{\prime}(t)-\left.\frac{\partial \ln \gamma(t, x ; T)}{\partial T}\right|_{T=t}  \tag{11.49}\\
= & f(0, t)+a^{\prime}(t)+\frac{\partial h(t, x)}{\partial t}+\mu_{x}(t, x) \frac{\partial h(t, x)}{\partial x} \\
& +\frac{1}{2} \sigma_{x}(t, x)^{2}\left(\frac{\partial^{2} h(t, x)}{\partial x^{2}}-\left(\frac{\partial h(t, x)}{\partial x}\right)^{2}\right) \tag{11.50}
\end{align*}
$$
^eq11-49
^eq11-50

To make matters more concrete, let us now specialize to the case where $h(t, x)=h(t) x$ and $\mu_{x}(t, x)=0$ in the general framework (11.45)-(11.46). One can show that this class of models includes the one-factor Gaussian
short rate model and the affine models. If one relaxes the requirement that $\mu_{x}(t, x) \equiv 0$, then the BK model is also in the class. With this restricted parameterization (11.50) yields, after ignoring small convexity terms,

$$
\begin{equation*}
r(t) \approx f(0, t)+a^{\prime}(t)+h^{\prime}(t) x(t) \tag{11.51}
\end{equation*}
$$
^eq11-51

and, approximately,

$$
d r(t) \approx\left(\vartheta_{r}(t)-\varkappa_{r}(t) r(t)\right) d t+\sigma_{r}(t, r(t)) d W(t)
$$

where

$$
\begin{equation*}
\varkappa_{r}(t)=-\frac{h^{\prime \prime}(t)}{h^{\prime}(t)}, \tag{11.52}
\end{equation*}
$$
^eq11-52

with $\vartheta_{r}(t), \sigma_{r}(t, r)$ appropriately defined. Hence, the numeraire scaling $h(t)$ can be conceptually linked to the mean reversion parameter for the short rate.

As a practical example of the general approach, Hagan and Woodward [1999a] propose the following class of " $\beta-\eta$ " models:

$$
\begin{aligned}
d x(t) & =\lambda(t)(1+\beta x(t))^{\eta} d W(t) \\
N(t) & =\frac{1}{P(0, t)} e^{h(t) x(t)+a(t)}
\end{aligned}
$$

For this specification, the transition density of $x(t)$ is known in closed form, allowing for (more or less) analytical calibration to the initial yield curve. The parameters $\beta, \eta$ are used to match the skew of the volatility smile. Note the resemblance between the volatility term in this model and that of a vanilla displaced-CEV model in Section 7.2.4. As pointed out in that section, adding a displacement to the CEV function does not significantly alter the range of available volatility smiles. Hence, should this approach be pursued, we recommend the specialization of the $\beta-\eta$ model with $\eta=1$, and perhaps an extension of the skew parameter $\beta$ to be time dependent. Specifying a constant mean reversion $\varkappa_{r}(t) \equiv \varkappa_{r}$ and solving (11.52) for $h(t)$ (which we normalize, conveniently and without loss of generality, by $h(0)=0$ and $h^{\prime}(0)=1$ ), we obtain

$$
\begin{aligned}
d x(t) & =\lambda(t)(1+\beta(t) x(t)) d W(t) \\
N(t) & =\frac{1}{P(0, t)} \exp \left(\frac{1-e^{-\varkappa_{r} t}}{\varkappa_{r}} x(t)+a(t)\right)
\end{aligned}
$$

In conclusion, we note a fairly strong resemblance between the Dybvig parameterization and the approach of Hagan and Woodward. Indeed, comparing (11.51) and (11.41) we see that $a^{\prime}(t)$ in the former plays pretty much the same role as $\vartheta(t)$ in the latter. The initial yield curve fit conditions, (11.48) and (11.43), are also rather similar. Hence, our words of caution with regards to the Dybvig parameterization apply here as well.

### 11.3.3 Monte Carlo Simulation

### 11.3.3.1 SDE Discretization

For the purposes of securities pricing by Monte Carlo methods, we are generally interested in advancing not only $r(t)$ through time, but also the inverse of the money market numeraire

$$
\begin{equation*}
\exp \left(-\int_{0}^{t} r(u) d u\right)=P(0, t) \exp \left(-\int_{0}^{t} x(u) d u\right) \triangleq P(0, t) Y(t) \tag{11.53}
\end{equation*}
$$
^eq11-53

where we recall that $x(t) \equiv r(t)-f(0, t)$. Our starting point can be ${ }^{13}$ the vector SDE

$$
\begin{equation*}
d\binom{x(t)}{Y(t)}=\binom{\mu_{x}(t, x(t))}{-x(t) Y(t)} d t+\binom{\sigma_{x}(t, x(t))}{0} d W(t) \tag{11.54}
\end{equation*}
$$
^eq11-54

where the functions $\mu_{x}$ and $\sigma_{x}$ were defined in (11.26) above.
In general, simulation of (11.54) requires usage of discretization methods, several of which were introduced in Chapter 3. For instance, the Euler scheme for (11.54) would advance the SDE for time $t_{i}$ to $t_{i+1}$ according to

$$
\binom{\widehat{x}_{i+1}}{\widehat{Y}_{i+1}}=\binom{\widehat{x}_{i}}{\widehat{Y}_{i}}+\binom{\mu_{x}\left(t_{i}, \widehat{x}_{i}\right)}{-\widehat{Y}_{i} \widehat{x}_{i}} \Delta_{i}+\binom{\sigma_{x}\left(t_{i}, \widehat{x}_{i}\right)}{0} Z_{i} \sqrt{\Delta_{i}}, \Delta_{i}=t_{i+1}-t_{i},
$$

where $\widehat{x}_{i}=\widehat{x}\left(t_{i}\right), \widehat{Y}_{i}=\widehat{Y}\left(t_{i}\right)$, and $Z_{i} \sim \mathcal{N}(0,1)$ is a sample from a standard Gaussian distribution. The Euler scheme is of (weak) convergence order one in the time step. To improve this, a second-order Milstein scheme can be constructed by Ito-Taylor expanding (11.54) to second order, using the technique in Section 3.2.6.3. The construction is tedious but straightforward (see Chapter IV in Andersen [1996] for the details), so we skip it and only show the final result

$$
\begin{align*}
\widehat{x}_{i+1} & =\widehat{x}_{i}+\left(\mu_{x}\left(t_{i}, \widehat{x}_{i}\right)-\frac{1}{2} \mathcal{L}_{1} \sigma_{x}\left(t_{i}, \widehat{x}_{i}\right)\right) \Delta_{i}+\sigma_{x}\left(t_{i}, \widehat{x}_{i}\right) Z_{i} \sqrt{\Delta_{i}} \\
& +\frac{1}{2}\left(\mathcal{L}_{1} \mu_{x}\left(t_{i}, \widehat{x}_{i}\right)+\mathcal{L}_{0} \sigma_{x}\left(t_{i}, \widehat{x}_{i}\right)\right) Z_{\imath} \Delta_{i} \sqrt{\Delta_{i}} \\
& +\frac{1}{2} \mathcal{L}_{0} \mu_{x}\left(t_{i}, \widehat{x}_{i}\right) \Delta_{i}^{2}+\frac{1}{2} \mathcal{L}_{1} \sigma_{x}\left(t_{i}, \widehat{x}_{i}\right) Z_{i}^{2} \Delta_{i}  \tag{11.55}\\
\widehat{Y}_{i+1} & =\widehat{Y}_{i}\left(1-\widehat{x}_{i} \Delta_{i}+\frac{1}{2}\left(\widehat{x}_{i}^{2}-\mu_{x}\left(t_{i}, \widehat{x}_{i}\right)\right) \Delta_{i}^{2}-\frac{1}{2} \sigma_{x}\left(t_{i}, \widehat{x}_{i}\right) Z_{i} \Delta_{i} \sqrt{\Delta_{i}}\right) \tag{11.56}
\end{align*}
$$
^eq11-55
^eq11-56

where we have introduced differential operators[^21]

$$
\mathcal{L}_{0}=\frac{\partial}{\partial t}+\mu_{x}(t, x) \frac{\partial}{\partial x}+\frac{1}{2} \sigma_{x}(t, x)^{2} \frac{\partial^{2}}{\partial x^{2}}, \quad \mathcal{L}_{1}=\sigma_{x}(t, x) \frac{\partial}{\partial x}
$$

The Milstein scheme (11.55)-(11.56) is rather formidable-looking, and its practical efficiency tends to be quite model-dependent. Still, using an affine model as a test case Andersen [1996] shows that the Milstein scheme outperforms the Euler scheme handily, even after taking the additional computational burden of (11.55)-(11.56) into consideration. Schemes with order higher than two can be constructed along similar principles, but will, in our experience, rarely be worth the hassle. We also remind the reader that higher-order schemes can be constructed by Richardson extrapolation, as discussed in Section 3.2.7. Andersen [1996] reports modest gains for a third-order scheme constructed by Richardson extrapolation of the Milstein scheme above.

At this point, let us note that for European-style securities paying a single cash flow at time $T$, the ideas of Section 10.1.6.3 can be applied here, and the burden of simulating $Y(t)$ could be avoided by a change to the $T$-forward measure $\mathrm{Q}^{T}$. For securities that pay intermediate cash flows, however, matters are more complicated as these flows must effectively be future-valued to time $T$. For instance, a random coupon $c$ paid at time $T^{\prime}<T$ will require us to compute the numeraire-deflated value $c / P\left(T^{\prime}, T\right)$. But here, unfortunately, the quantity $P\left(T^{\prime}, T\right)$ is generally not known analytically at time $T^{\prime}$ as a function of the model state variables. Of course, without affecting the economics of the trade one could invest the proceeds $c$ into a money market account $\beta$ at time $T^{\prime}$, yielding the payout $c \beta(T) / \beta\left(T^{\prime}\right)$ at time $T$. Evaluating this payout, however, would again require us to keep track of $Y(t)$, at least on the interval [ $T^{\prime}, T$ ]. This problem, however, can be avoided by using the spot measure instead of the forward measure, as outlined in Section 10.1.6.3. Much more material about numeraire simulation strategies can be found in Chapter 14.

Finally, a note on variance reduction for short rate model simulation. A systematic discussion of variance reduction techniques for short rate models can be found in Chapter IV of Andersen [1996] and Andersen and Boyle [2000]. Most of the methods discussed in these sources can be found in the survey of Section 3.4 and shall not be repeated. We do, however, highlight here the particularly useful idea of applying importance sampling based on information extracted from a tractable (e.g. Gaussian or affine) approximation to the short rate SDE. We postpone the discussion of this technique, which relies on the material in Section 3.4.4.3, to Chapter 25.

### 11.3.3.2 Practical Issues with Monte Carlo Methods

As was the case for finite difference methods (see discussion in Section 11.3.1), whenever an explicit bond reconstitution formula is lacking, the effort and complexity required to price derivatives by Monte Carlo methods
increase significantly. For instance, consider applying Monte Carlo methods to the pricing of short-dated expiry $T$ call option on a long-dated maturity $T^{*}$ discount bond. This price ( $V$ ) is computed as a risk-neutral expectation

$$
\begin{aligned}
V(0) & =P(0, T) \mathrm{E}\left(Y(T)^{-1}\left(P\left(T, T^{*}\right)-K\right)^{+}\right) \\
& =P(0, T) \mathrm{E}\left(Y(T)^{-1}\left(\mathrm{E}_{T}\left(e^{-\int_{T}^{T^{*}} x(u) d u-\int_{T}^{T^{*}} f(0, u) d u}\right)-K\right)^{+}\right) \\
& =P\left(0, T^{*}\right) \mathrm{E}\left(Y(T)^{-1}\left(\mathrm{E}_{T}\left(e^{-\int_{T}^{T^{*}} x(u) d u}\right)-K \frac{P(0, T)}{P\left(0, T^{*}\right)}\right)^{+}\right)
\end{aligned}
$$

An immediate problem is here the fact that the inner time $T$ expectation

$$
\mathrm{E}_{T}\left(e^{-\int_{T}^{T^{*}} x(u) d u}\right)
$$

is not explicitly known as a function of $x(T)$, but must itself be computed by numerical methods. A brute-force approach involves estimating the expectation by Monte Carlo methods, launching a "simulation-within-asimulation" at time $T$. The computational expense involved in such a scheme would most likely be prohibitive. Alternatives involve using a regression on a space of basis functions to estimate the function

$$
Q\left(T, T^{*}, x\right)=\mathrm{E}\left(e^{-\int_{T}^{T^{*}} x(u) d u} \mid x(T)=x\right)
$$

we discuss this approach in some detail in Chapter 18.
Alternatively, we can always estimate $Q\left(T, T^{*}, x\right)$ by finite difference methods, as in Section 11.3.1. Combining finite difference methods and Monte Carlo methods for the purposes of pricing a European option on a discount bond makes little practical sense, of course, as we would always prefer finite difference methods for this payout. For path-dependent options, however, this idea may in fact be the best way of computing option prices. Loosely, such a scheme would use a finite difference grid to pre-compute zero-coupon bond prices $P\left(t_{i}, \cdot, x\right)$ on a grid $x \in\left\{x_{j}\right\}_{j=0}^{m+1}$, at all dates $t_{i}$, $i=1,2, \ldots, N$, required by the path-dependent payout function considered. When paths for $x(t)$ in a subsequent step are generated by Monte Carlo simulation, interpolation of the $N$ discount bond price vectors of dimension $m+2$ available in the finite difference grid would allow us to compute rapidly discount bond prices $Q\left(t_{i}, \cdot, x\left(t_{i}\right)\right), i=1,2, \ldots$, at all relevant dates. We can use the schemes in Section 11.3.3.1 for the purpose of drawing paths of $x(t)$. If, however, we wish to make the dynamics for $x(t)$ perfectly consistent with the finite difference grid, we can use the forward induction techniques of Sections 11.3.2.1 and 11.3.2.2 to work out the (discrete) transition probabilities for $x(t)$ implied by the finite difference grid. Paths for $x(t)$ can then be generated directly from these probabilities. With this approach, we only draw values of $x(t)$ on the spatial grid $\left\{x_{j}\right\}_{j=0}^{m+1}$ of the finite difference grid, and therefore never have to apply interpolation methods when looking up discount bond prices.

## 11.A Appendix: Markov-Functional Models

The purpose of this appendix is to give a brief account of the class of Markovfunctional (MF) models. We only consider the one-factor case. Extensions to higher dimensions are possible, but practical implementation challenges tend to increase substantially for dimensions higher than one. MF models were introduced in Kennedy et al. [2000] and while their popularity is generally waning, they are still used in some banks.

## 11.A.1 State Process and Numeraire Mapping

We have already observed in Section 11.3.2.6 that to define an arbitrage-free interest rate model we really only need two ingredients: a stochastic process that drives the evolution of interest rates, and a functional form that maps that process into a numeraire. The development of Markov-functional models normally starts with specializing this setup to a numeraire taken to be the discount bond $P\left(\cdot, T^{*}\right)$ to the final maturity of interest $T^{*}$, and a Markov stochastic process that is Gaussian in the corresponding terminal measure $\mathrm{Q}^{*}$ (see Section 4.2.4). Assuming arbitrarily that $x(t)$ is a $\mathrm{Q}^{*}$-martingale, we write

$$
\begin{equation*}
d x(t)=\sigma(t) e^{\varkappa t} d W^{*}(t), \quad x(0)=0 \tag{11.57}
\end{equation*}
$$
^eq11-57

where $W^{*}$ is a one-dimensional Brownian motion in the terminal measure, and where we for simplicity have assumed that the mean reversion $\varkappa$ is constant. The role of $\varkappa$ is to control inter-temporal correlations in the model; see Section 13.1.8.1 for the importance of this. The transition density of $x(t)$ in $\mathrm{Q}^{*}$ is, trivially,

$$
\begin{aligned}
p(y, t ; z, s) & \triangleq Q^{*}(x(s) \in[z, z+d z] \mid x(t)=y) / d z \\
& =\frac{1}{\sqrt{2 \pi} v(t, s)} \exp \left(\frac{-(z-y)^{2}}{2 v(t, s)^{2}}\right), \quad s>t
\end{aligned}
$$

where

$$
v(t, s)^{2}=\int_{t}^{s} \sigma(u)^{2} e^{2 \approx u} d u
$$

In the spirit of Section 11.3.2.6, we define $P\left(t, T^{*}\right)$ to be a deterministic function of the state variable process $x(t)$,

$$
\begin{equation*}
P\left(t, T^{*}\right)=P\left(t, T^{*}, x(t)\right), \quad P\left(t, T^{*}, x\right)=H(t, x), \tag{11.58}
\end{equation*}
$$
^eq11-58

for some exogenously given function $H: \mathbb{R}^{2} \rightarrow[0,1]$. As we recall, this is sufficient to define all discount bonds in the model since, for any $0 \leq t< T \leq T^{*}$,

$$
\begin{equation*}
P(t, T)=P\left(t, T^{*}\right) \mathrm{E}_{t}^{*}\left(\frac{1}{P\left(T, T^{*}\right)}\right) \tag{11.59}
\end{equation*}
$$
^eq11-59

where $\mathrm{E}^{*}$ denotes expectation in measure $\mathrm{Q}^{*}$. This allows us to express all discount bonds as functions of $x(t)$ :

$$
\begin{align*}
P(t, T) & =P(t, T, x(t)) \\
P(t, T, x) & =P\left(t, T^{*}, x\right) \mathrm{E}^{*}\left(\left.\frac{1}{P\left(T, T^{*}, x(T)\right)} \right\rvert\, x(t)=x\right) \\
& =H(t, x) \int_{-\infty}^{\infty} \frac{p(x, t ; z, T)}{H(T, z)} d z \tag{11.60}
\end{align*}
$$
^eq11-60

The formula (11.59) can be specialized to $t=0$, yielding

$$
\begin{equation*}
P(0, T)=P\left(0, T^{*}\right) \mathrm{E}^{*}\left(\frac{1}{P\left(T, T^{*}\right)}\right)=P\left(0, T^{*}\right) \mathrm{E}^{*}\left(\frac{1}{H(T, x(T))}\right) \tag{11.61}
\end{equation*}
$$
^eq11-61

which constitutes a no-arbitrage condition on the mapping function $H(\cdot, \cdot)$. This condition is often used to choose a particular function $H(\cdot, \cdot)$ from a given parametric family; compare this to condition (11.48) in Section 11.3.2.6.

In practice, the numeraire mapping function $H(t, x)$ in (11.58) is often specified only indirectly, through definition of functional forms for market rates, such as Libor or swap rates. The following two sections explore variations on this idea.

## 11.A. 2 Libor MF Parameterization

Let us assume that a tenor structure

$$
0=T_{0}<T_{1}<\ldots<T_{N}=T^{*}, \quad \tau_{n}=T_{n+1}-T_{n}
$$

is given, and define spanning forward Libor rates by

$$
\begin{equation*}
L_{n}(t) \triangleq L\left(t, T_{n}, T_{n+1}\right)=\left(\frac{P\left(t, T_{n}\right)}{P\left(t, T_{n+1}\right)}-1\right) \tau_{n}^{-1}, \quad n=0, \ldots, N-1 \tag{11.62}
\end{equation*}
$$
^eq11-62

(see (4.2)). It turns out that, if we can specify the mapping of the state process $x(\cdot)$ into Libor rates on their fixing dates, $L_{n}\left(T_{n}\right)$, for all $n=1, \ldots, N-1$ :

$$
L_{n}\left(T_{n}\right)=l_{n}\left(x\left(T_{n}\right)\right), \quad n=1, \ldots, N-1
$$

then this is sufficient to recover the numeraire-mapping function $H\left(T_{n}, \cdot\right)$, $n=1, \ldots, N$, on tenor dates ${ }^{14}$ and consequently define the MF model by

[^22](11.59). We show this by induction on the fixing time $T_{n}$, for $n=N-1, \ldots, 1$. The starting point of the induction follows directly from (11.62) as we have

$$
\begin{equation*}
H\left(T_{N-1}, x\right)=P\left(T_{N-1}, T^{*}, x\right)=P\left(T_{N-1}, T_{N}, x\right)=\left(1+\tau_{N-1} l_{N-1}(x)\right)^{-1} \tag{11.63}
\end{equation*}
$$
^eq11-63

For the induction step, let us assume that $H\left(T_{i}, x\right)$ are known for $i= n+1, \ldots, N-1$. By (11.59) we have

$$
\frac{P\left(T_{n}, T_{n+1}\right)}{P\left(T_{n}, T^{*}\right)}=\mathrm{E}_{T_{n}}^{*}\left(\frac{1}{P\left(T_{n+1}, T^{*}\right)}\right)
$$

which implies that

$$
\begin{equation*}
\frac{1}{H\left(T_{n}, x\right)}=\left(1+\tau_{n} l_{n}(x)\right) \mathrm{E}_{T_{n}}^{*}\left(\left.\frac{1}{H\left(T_{n+1}, x\left(T_{n+1}\right)\right)} \right\rvert\, x\left(T_{n}\right)=x\right), \tag{11.64}
\end{equation*}
$$
^eq11-64

and the statement follows.
The consistency condition (11.61) is often used to select a particular function $l_{n}(\cdot)$ from a parametric family, for each $n$. To explain this, let us first consider what functional forms for $l_{n}(\cdot)$ are typically used. Suppose we desire to build a model where Libor rates on the tenor structure are close to log-normal ${ }^{15}$. Then, with $v_{n}=v\left(0, T_{n}\right)$, we fundamentally would want something like

$$
\begin{equation*}
L_{n}\left(T_{n}\right) \approx L_{n}(0) \exp \left(k_{n} x\left(T_{n}\right)-\frac{1}{2} k_{n}^{2} v_{n}^{2}\right) \tag{11.65}
\end{equation*}
$$
^eq11-65

where

$$
k_{n}=\frac{e^{-\varkappa T_{n}}-e^{-\varkappa T_{n+1}}}{\varkappa \tau_{n}},
$$

to hold for all $n$. The particular form of $k_{n}$ 's, as well as the volatility parameterization (11.57), are strongly inspired by the Gaussian short rate model, see Proposition 10.1.7. (To see the connection more clearly the reader should note that the state variable $x(t)$ here is related to $x(t)$ in Proposition 10.1.7 by a multiplicative scaling of $e^{\varkappa t}$; compare (11.57) to (10.16) and disregard $y(t)$ in the latter.) Note that the quantity

$$
T_{n}^{-1 / 2} k_{n} v_{n}
$$

after calibration should be close to the implied Black volatility of a caplet maturing at time $T_{n}$. To preclude arbitrage, (11.65) cannot be used as is for all $n$ (since only $L_{N-1}(t)$ is a martingale in $Q^{*}$ ), so we could, for instance, add a "convexity multiplier" $c_{n}$ and write

$$
l_{n}\left(x ; c_{n}\right)=c_{n} L_{n}(0) \exp \left(k_{n} x-\frac{1}{2} k_{n}^{2} v_{n}^{2}\right) .
$$

[^23]While $v_{n}$ 's (or, equivalently, the mean reversion $\varkappa$ and the model volatility $\sigma(t))$ can be treated as free constants to be calibrated to option prices, we would use (11.61) to set $c_{n}$ 's such that the initial yield curve is replicated by the model. It is trivial to see that $c_{N-1}=1$ and

$$
H_{N-1}\left(T_{N-1}, x\right)=\left(1+\tau_{N-1} l_{N-1}(x ; 1)\right)^{-1}
$$

wherefore other $c_{n}$ 's may be obtained as solutions to

$$
\mathrm{E}_{T_{n}}^{*}\left(\left(1+\tau_{n} l_{n}\left(x\left(T_{n}\right) ; c_{n}\right)\right) \mathrm{E}_{T_{n}}^{*}\left(\left.\frac{1}{H\left(T_{n+1}, x\left(T_{n+1}\right)\right)} \right\rvert\, x\left(T_{n}\right)\right)\right)=P\left(0, T_{n}\right)
$$

for $n=N-2, \ldots, 1$.

## 11.A.3 Swap MF Parameterization

Defining an MF model in terms of Libor rates is especially convenient if the model is meant to price a security that depends primarily on Libor rates (on their fixing dates), e.g. a TARN (see Section 5.15 and Chapter 20). In particular, the Libor MF parameterization allows one direct control over volatilities and other distributional characteristics of Libor rates, which makes it fairly straightforward to set up a calibration scheme that is suitable for the security (see e.g. Section 20.1.3). Libor rates, however, are not always the primary driving factors; for example, prices of Bermudan swaptions are arguably more directly linked to distributions of swap rates (see Section 19.2). Fortunately, MF models can be formulated in terms of swap rates as well, as we shall now demonstrate. It is worth noting that the relationship between Libor and swap MF models is similar to that between Libor and swap market models which we explore in Section 15.4.

For concreteness let us consider a set of so-called "core" swap rates,

$$
\begin{align*}
& S_{n}(t) \triangleq S_{n, N-n}(t)=\frac{P\left(t, T_{n}\right)-P\left(t, T^{*}\right)}{A_{n}(t)}  \tag{11.66}\\
& A_{n}(t) \triangleq A_{n, N-n}(t)=\sum_{i=n}^{N-1} \tau_{i} P\left(t, T_{i+1}\right)
\end{align*}
$$
^eq11-66

$n=1, \ldots, N-1$, where we used the notations (4.8), (4.10) for $A_{k, m}$ and $S_{k, m}$. We assume that the core swap rates on their fixing dates are specified as deterministic functions $s_{n}(x)$ of the state process,

$$
S_{n}\left(T_{n}\right)=s_{n}\left(x\left(T_{n}\right)\right), \quad n=1, \ldots, N-1 .
$$

As for the Libor-based specification above, we claim that the knowledge of the functions $\left\{s_{n}(\cdot)\right\}_{n=1}^{N-1}$ (together with the dynamics of the state process $x(t)$ ) is sufficient to define the numeraire mapping $H(\cdot, \cdot)$ and, therefore, an arbitrage-free model of interest rates.

The proof also proceeds by induction. As we have that $S_{N-1}\left(T_{N-1}\right)= L_{N-1}\left(T_{N-1}\right)$, the starting point of induction is given by (11.63), i.e.

$$
H\left(T_{N-1}, x\right)=\left(1+\tau_{N-1} s_{N-1}(x)\right)^{-1}
$$

For the induction step $n+1 \rightarrow n$, we note from (11.66) that

$$
\begin{aligned}
\frac{1}{P\left(T_{n}, T^{*}\right)} & =1+S_{n}\left(T_{n}\right) \frac{A_{n}\left(T_{n}\right)}{P\left(T_{n}, T^{*}\right)} \\
& =1+S_{n}\left(T_{n}\right) \sum_{i=n}^{N-1} \tau_{i} \frac{P\left(T_{n}, T_{i+1}\right)}{P\left(T_{n}, T^{*}\right)}
\end{aligned}
$$

and so we have (compare to (11.64))

$$
\frac{1}{H\left(T_{n}, x\right)}=1+s_{n}(x) \sum_{i=n}^{N-1} \tau_{i} \mathrm{E}_{T_{n}}^{*}\left(\left.\frac{1}{H\left(T_{i+1}, x\left(T_{i+1}\right)\right)} \right\rvert\, x\left(T_{n}\right)=x\right) .
$$

As in the Libor specification, we can choose functions $\left\{s_{n}(\cdot)\right\}$ to approximate log-normal (or displaced log-normal) distribution of swap rates. Also in direct analogy to the Libor case, we typically have some no-arbitrage conditions to satisfy, usually by means of setting some parameters in the specific functional form of $\left\{s_{n}(\cdot)\right\}$. We leave these details for the reader to explore.

## 11.A.4 Non-Parametric Calibration

So far we defined Markov-functional models by specific parametric mappings of the state process into Libor and swap rates. Originally, however, the class of models was introduced in a non-parametric way (see Hunt and Kennedy [2000] for a typical treatment), where mapping functions are deduced from market prices of caplets or swaptions across all strikes. While we typically prefer the parametric approach (for reasons we touch upon below), let us nevertheless quickly review the non-parametric method for completeness.

Through equation (11.60), we can turn the payout of any $T$-maturity security that depends on the state of the yield curve into a function of $x(T)$, $g(x(T) ; K, T)$ say, where $K$ is some payout parameter (virtually always a strike). The time 0 price of this security is

$$
\begin{align*}
V(0 ; K, T) & =P\left(0, T^{*}\right) \mathrm{E}^{*}\left(\frac{g(x(T))}{P\left(T, T^{*}\right)}\right) \\
& =P\left(0, T^{*}\right) \int_{-\infty}^{\infty} p(0,0 ; z, T) \frac{g(z ; K)}{H(T, z)} d z \tag{11.67}
\end{align*}
$$
^eq11-67

Assuming that $H(T, x)$ is invertible in $x$, we may write (11.67) as

$$
\begin{equation*}
V(0 ; K, T)=P\left(0, T^{*}\right) \int_{-\infty}^{\infty} p(0,0 ; z, T) q(H(T, z) ; K, T) d z \tag{11.68}
\end{equation*}
$$
^eq11-68

where

$$
q(H(T, z) ; K, T) \triangleq \frac{g(z ; K, T)}{H(T, z)}
$$

If $V\left(0 ; K, T_{n}\right)$ is known ${ }^{16}$ for a continuum of parameters (strikes) $K$, (11.68) defines an integral equation that may allow one to uncover the function $H\left(T_{n}, \cdot\right)$ (or, often more conveniently, $l_{n}(\cdot)$ or $s_{n}(\cdot)$ ).

Solution of (11.68) is typically done for a fixed number of $M$ strikes, with $H(T, x)$ solved for on a grid $\left\{x_{j}\right\}_{j=1}^{M}$. In practice, this procedure is difficult to make fully robust, and the numerical solution is often prone to instabilities at long maturities, even if sophisticated special-purpose numerical techniques are employed (see Hunt and Kennedy [2000] for such techniques, many of which rely on the fact that polynomials can be integrated exactly against the Gaussian density). Even when numerically stable, a non-parametric solution for $H(\cdot, \cdot)$ may imply unrealistic evolution of the volatility smile through time, a general feature of local volatility models as explained in Section 7.1.3. To avoid these issues we may either pre-smooth the option prices used for calibration purposes (e.g. by best-fitting a CEV or a displaced log-normal model to the market smile), or, preferably in our opinion, we may use a low-dimensional parametric form for $H(t, x)$ as in Sections 11.A. 2 and 11.A.3.

## 11.A. 5 Numerical Implementation

Numerical securities valuation in an MF model is typically quite simple as the state process is Gaussian. Let us assume that the function $H(t, x)$ has been established, and consider, say, implementation of the model in a finite difference grid. Let the derivative value function be $V(t, x(t))$, and set $V^{*}(t, x)=V(t, x) / H(t, x)$ (such that $V(0,0)=P\left(0, T^{*}\right) V^{*}(0,0)$ ). As $V^{*}(t, x)$ must be a $\mathrm{Q}^{*}$-martingale, we can write

$$
\begin{equation*}
\frac{\partial V^{*}}{\partial t}+\frac{1}{2} \sigma(t)^{2} e^{2 \varkappa t} \frac{\partial^{2} V^{*}}{\partial x^{2}}=0 \tag{11.69}
\end{equation*}
$$
^eq11-69

subject to appropriate terminal and intermediate jump payout conditions. In evaluating terminal and intermediate payout conditions, we would typically need to apply the numerical expression (11.60) to establish the state of the yield curve. We should note that the MF literature generally prefers to use Gaussian integration methods (rather than standard PDE solvers) to evaluate the PDE (11.69), see Hunt and Kennedy [2000] for details.

[^24]
## 11.A. 6 Comments and Comparisons

The one-factor MF model competes with a number of models in this book, especially the quasi-Gaussian class in Chapter 13. The quasi-Gaussian model allows for arbitrary local volatility (as does the MF model), but has closed-form formulas that allow for reconstituting the term structure of discount bonds analytically from the underlying state variable, rather than through numerical integration (see (11.60)). In addition, the quasiGaussian model is substantially more "direct" in its modeling of the forward curve and has an easy-to-state term structure of instantaneous forward rate volatilities. This, in turn, makes the model more transparent in its causality structure - especially when it comes to the evolution of the volatility term structure and smile - and often makes it easy to devise good closed-form approximations for swaption and cap prices. As a consequence, calibration of quasi-Gaussian models to option and bond prices is virtually always much faster than for MF models. In addition, quasi-Gaussian models are quite straightforward to extend to high dimensions and to stochastic volatility dynamics; these extensions are far more difficult ${ }^{17}$ for MF models. On the flip side, a quasi-Gaussian model involving a single Brownian motion involves a two-dimensional state vector process, which makes derivatives pricing by finite difference methods slower than for MF models ${ }^{18}$. For most applications, total computation time of calibration and valuation is, however, less for the quasi-Gaussian model.

Due to its flexibility, extensibility, transparency, and ease of numerical implementation (no integration tricks are required), we generally prefer the quasi-Gaussian model over the MF model, and consequently dedicate an entire chapter to the former - and only this appendix to the latter. For those interested in learning more about MF models than what we offered here, Kennedy et al. [2000] and Hunt and Kennedy [2000] are good starting points.

[^25]
