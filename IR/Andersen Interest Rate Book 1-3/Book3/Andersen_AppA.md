---
tags: [quant/knowledge, quant/books, interest-rates]
---

## Markovian Projection

Markovian projection is a powerful method for simplifying complex process dynamics to a form that enables rapid calibration of model parameters to quoted option prices. We use the method several times in this book, see for instance Chapters 13 and 15. The usefulness of Markovian projection, however, extends beyond interest rate modeling applications. In this appendix, we develop the relevant theory behind the method and present additional examples and applications.

## A. 1 Marginal Distributions of Ito Processes

Models used in quantitative finance generally serve to define dynamics of market observables. Some models impose such dynamics directly on the observables; this is, for instance, the case for vanilla models (see Chapters 7, 8, 16 and 17) where the evolution of swap rates is modeled explicitly. Outside of interest rate modeling, equity and FX models typically fall in this category as well. In other cases, the dynamics of market observables are specified indirectly, through modeling of abstract Markovian state variables that drive the market observables through functional relationships (i.e., reconstitution formulas). This style of modeling is common in term structure models for commodities and interest rates (see Chapter 13 for a typical model).

Regardless of type, all models ultimately need to be calibrated to liquidly traded options, most often European call/put options on market observables. To facilitate efficient model calibration, it is generally helpful if exact or approximate analytical expressions exist for European options. Often, the derivation of such results is significantly aided by an initial simplification of the underlying dynamic processes, either because these processes are outright too complex to handle analytically, or because non-linear reconstitution formulas translate simple state variable dynamics into intractable dynamics for market observables. The viability of such simplifications stems from
the simple structure of European options, which only depend on the onedimensional marginal distributions of the market observable process. As it turns out, irrespective of how complicated a process for a particular marker observable is, it is often possible to find much simpler process that preserves the marginal distributions.

A systematic way of finding process simplifications is based on the following fundamental result, see Theorem 4.6 in Gyöngy [1986].

Theorem A.1.1 (Gyöngy). Let $X(t)$ be given by an $S D E$

$$
\begin{equation*}
d X(t)=\lambda(t) d W(t) \tag{(A. 1 !}
\end{equation*}
$$
^eqa-1

where $W(t)$ is a one-dimensional Brownian motion under some probability measure P . Assume that the process $\lambda(t)$ is adapted, bounded, and uniformly bounded away from 0 , such that (A.1) admits a unique strong solution_ Define $b(t, x)$ by

$$
\begin{equation*}
b(t, x)^{2}=\mathrm{E}\left(\lambda(t)^{2} \mid X(t)=x\right) \tag{A.2}
\end{equation*}
$$
^eqa-2

where E is the expected value operator for the pricing measure P . Then the SDE

$$
\begin{equation*}
d Y(t)=b(t, Y(t)) d W(t), \quad Y(0)=X(0) \tag{A.3}
\end{equation*}
$$
^eqa-3

admits a weak solution $Y$ that has the same one-dimensional distributions as $X$.

Remark A.1.2. The original result by Gyöngy also includes a drift in the dynamics of $X$, considering

$$
d X(t)=\mu(t) d t+\lambda(t) d W(t)
$$

instead of (A.1). The theorem then still holds with (A.3) replaced by

$$
d Y(t)=a(t, Y(t)) d t+b(t, Y(t)) d W(t), \quad Y(0)=X(0)
$$

where

$$
a(t, x)=\mathrm{E}(\mu(t) \mid X(t)=x) .
$$

In financial applications we nearly always have $\mu(t)=0$ as we tend to consider the dynamics of $X$ in its own martingale measure. For this reason we do not consider drifts in what follows.

Proof. The original proof in Gyöngy [1986] is fairly involved. A rigorous proof under much weaker assumptions than we stated (see Proposition A.1.4 below) is given in Brunick [2008] and is also highly technical. We do not reproduce either of these proofs, but instead we present a somewhat informal argument ${ }^{1}$ originally due to Dupire, see Dupire [1994], Dupire [1997], who independently discovered essentially the same results as in Gyöngy [1986].

[^116]The function $b(t, x)$ is often called the Dupire local volatility function for the process $X$.

Let us denote

$$
c(t, K) \triangleq c(0, S(0) ; t, K)=\mathrm{E}\left((X(t)-K)^{+}\right)
$$

to be the values of European call options on $X$ for expiries $t$ and strikes $K$. It follows from Proposition 7.4.2 that, if we define $Y$ by

$$
\begin{equation*}
d Y(t)=b(t, Y(t)) d W(t) \tag{A.4}
\end{equation*}
$$
^eqa-4

with

$$
\begin{equation*}
b(t, K)^{2}=\frac{2 \frac{\partial}{\partial t} c(t, K)}{\frac{\partial^{2}}{\partial K^{2}} c(t, K)}, \tag{A.5}
\end{equation*}
$$
^eqa-5

then the values of European call options in the model (A.4) will be equal to $c(t, K)$ for all expiries $t$ and strikes $K$, i.e. will be the same as in the model (A.1). To compute the right-hand side, we first write (the use of Dirac delta functions in the integrands can be justified by Tanaka's formula, see Section 1.9.2 and Karatzas and Shreve [1997])

$$
d(X(t)-K)^{+}=1_{\{X(t)>K\}} d X(t)+\frac{1}{2} \delta(X(t)-K) \lambda(t)^{2} d t
$$

and, since $X(t)$ is a P -martingale,

$$
\mathrm{E}(X(t)-K)^{+}-(X(0)-K)^{+}=\frac{1}{2} \int_{0}^{t} \mathrm{E}\left(\delta(X(t)-K) \lambda(t)^{2}\right) d t
$$

Clearly

$$
\mathrm{E}\left(\delta(X(t)-K) \lambda(t)^{2}\right)=\mathrm{E}(\delta(X(t)-K)) \times \mathrm{E}\left(\lambda(t)^{2} \mid X(t)=K\right)
$$

and

$$
\mathrm{E}(\delta(X(t)-K))=\frac{\partial^{2}}{\partial K^{2}} \mathrm{E}(X(t)-K)^{+}=\frac{\partial^{2}}{\partial K^{2}} c(t, K) .
$$

In particular,

$$
\begin{aligned}
\frac{\partial}{\partial t} c(t, K) & =\frac{\partial}{\partial t}\left(\mathrm{E}(X(t)-K)^{+}-(X(0)-K)^{+}\right) \\
& =\frac{1}{2} \frac{\partial^{2}}{\partial K^{2}} c(t, K) \times \mathrm{E}\left(\lambda(t)^{2} \mid X(t)=K\right)
\end{aligned}
$$

Substituting this equality into (A.5) we obtain

$$
b(t, K)^{2}=\mathrm{E}\left(\lambda(t)^{2} \mid X(t)=K\right),
$$

consistent with (A.2). $\square$

Since $X$ and $Y$ have the same one-dimensional marginal distributions, the prices of European options on $X$ and $Y$ for all strikes $K$ and expiries $T$ will be identical (a result that is also implicit in our proof of the theorem). Thus. for the purposes of European option valuation, a potentially complicated process $X$ can be replaced with a simpler Markov process $Y$; we call $Y$ the Markovian projection of $X$. Notice that the process $Y$ conveniently is of the local volatility type considered in Chapter 7, for which we have developed many exact or approximate methods for valuation of European options.

Theorem A.1.1 can be extended in a number of ways. Possibly the simplest extension involves relaxing the assumption that the Brownian motion $W(t)$ in (A.1) is one-dimensional, in which case the following trivial corollary holds.

Corollary A.1.3. Suppose that $X$ follows multi-factor dynamics

$$
\begin{equation*}
d X(t)=\lambda(t)^{\top} d W(t) \tag{A.6}
\end{equation*}
$$
^eqa-6

with $W(t)$ a $d$-dimensional Brownian motion and $\lambda(t)$ a $d$-dimensional adapted process whose norm is bounded and uniformly bounded away from 0. Define the SDE

$$
\begin{equation*}
d Y(t)=b(t, Y(t)) d \widetilde{W}(t), \quad Y(0)=X(0) \tag{A.7}
\end{equation*}
$$
^eqa-7

where $\widetilde{W}(t)$ a one-dimensional Brownian motion, and

$$
b(t, x)^{2}=\mathrm{E}\left(\lambda(t)^{\top} \lambda(t) \mid X(t)=x\right)
$$

Then (A.7) admits a weak solution $Y(t)$ that has the same one-dimensional distributions as $X(t)$.

Proof. Clearly $X$ can be written in one-dimensional form

$$
d X(t)=\left(\lambda(t)^{\top} \lambda(t)\right)^{1 / 2} d \widetilde{W}(t)
$$

where

$$
d \widetilde{W}(t)=\frac{1}{\left(\lambda(t)^{\top} \lambda(t)\right)^{1 / 2}} \lambda(t)^{\top} d W(t)
$$

Simple quadratic variance calculations show that $\widetilde{W}(t)$ is a one-dimensional Brownian motion, and the corollary then is a direct consequence of Theorem A.1.1.

The original result by Gyöngy required that the variance process $\lambda(t)^{\top} \lambda(t)$ in (A.6) be both bounded and uniformly bounded away from 0 . These are rather severe limitations that are violated in some standard models of mathematical finance, including the Heston model of Chapter 8. Brunick [2008] has proved the same result under much milder regulatory conditions and also extended it to the case where the asset process $X$ itself is multi-dimensional.

Proposition A.1.4. Let $X(t)$ be a $p$-dimensional stochastic process given by the strong solution of the SDE

$$
d X(t)=\lambda(t)^{\top} d W(t)
$$

where now $\lambda(t)$ is a $(d \times p)$-matrix-valued adapted process and $W(t)$ is a $d$-dimensional Brownian motion. We assume that

$$
\mathrm{E}\left(\int_{0}^{t}\left\|\lambda(s)^{\top} \lambda(s)\right\| d s\right)<\infty
$$

for all $t \geq 0$. Then there exists a $(d \times p)$-matrix-valued function $b$ such that

$$
b(t, x)^{\top} b(t, x)=\mathrm{E}\left(\lambda(t)^{\top} \lambda(t) \mid X(t)=x\right)
$$

for any $t \geq 0$, the SDE

$$
d Y(t)=b(t, Y(t))^{\top} d W(t), \quad Y(0)=X(0)
$$

admits a weak solution $Y$, and the random vector $Y(t)$ has the same distribution as the random vector $X(t)$ for any $t \geq 0$.

More generally, Brunick [2008] also proves that we can construct a "mimicking" process $Y$ such that marginal distributions of some functional of $X$, rather than of $X$ itself, are matched by a suitable functional of $Y$. The definition of functionals for which this result works is fairly technical, but the allowed set includes such financially relevant cases as functions of the running average of $X$, functions of the running maximum (or minimum) of $X$, and many others. To avoid technicalities, let us consider only the case of a running maximum of a one-dimensional asset process. Specifically, let us define $M(t, \xi)$ to be the running maximum for a given process $\xi$, i.e.

$$
M(t, \xi)=\sup _{0 \leq s \leq t} \xi(s), \quad t \geq 0
$$

Proposition A.1.5. Let $X$ follow a one-dimensional diffusion

$$
d X(t)=\lambda(t) d W(t)
$$

where $W(t)$ is a one-dimensional Brownian motion, and $\lambda(t)$ is a scalar adapted process that satisfies

$$
\mathrm{E}\left(\int_{0}^{t} \lambda(s)^{2} d s\right)<\infty, \quad t \geq 0
$$

Then i) there exists a deterministic function $b(t, x, m)$ such that

$$
b(t, x, m)^{2}=\mathrm{E}\left(\lambda(t)^{2} \mid X(t)=x, M(t, X)=m\right)
$$

holds for all $t \geq 0$; ii) the SDE

$$
d Y(t)=b(t, Y(t), M(t, Y)) d W(t), \quad Y(0)=X(0)
$$

admits a weak solution $Y$; and iii) the pair $(Y(t), M(t, Y))$ has the same distribution as the pair $(X(t), M(t, X))$ for any $t \geq 0$.

Example A.1.6. A European up-and-out barrier call option (see Section 2.1) with expiry $T$, strike $K$ and barrier level $B$ is an option with the payoff

$$
1_{\{M(T, X)<B\}}(X(T)-K)^{+} .
$$

The values of such options on $X$ for all expiries, strikes and barrier levels match those on the mimicking process $Y$ defined by Proposition A.1.5.

While the projection defined by Proposition A.1.5 matches all (up-andout) barrier option prices, the standard projection result in Theorem A.1.1 does not, as it does not preserve joint distributions of the process observed at multiple times. This is an important point that should be kept in mind in a calibration setting: the Markovian projection of Theorem A.1.1 should be used solely as the means to calibrate to European option prices and not to more complicated derivatives ${ }^{2}$.

## A. 2 Approximations for Conditional Expected Values

According to Theorem A.1.1, the coefficients for the SDE of the Markovian projection are obtained by calculating conditional expected values as in (A.2). This, in the majority of interesting cases, is a non-trivial task. Below, we consider several possible approximations.

## A.2.1 Gaussian Approximation

Of the few probability distributions that allow us to calculate conditional expected values in closed form, the most important is, of course, the Gaussian distribution (a fact we use extensively in many places of the book, see for instance Chapter 17). Not surprisingly, we can get good mileage out of the idea of approximating the original distributions of $X(t)$ and $\lambda(t)^{2}$ with Gaussian distributions, in order to calculate the conditional expected value in (A.6). Many variations are possible here; we present one approach in Proposition A.2.1 below. To fix our setup, we assume that $X$ follows the SDE (A.1) with a process $\lambda(t)$ given by the SDE

$$
d \lambda(t)^{2}=\nu(t) d t+\varepsilon(t) d Z(t)
$$

with two adapted stochastic processes $\nu(t)$ and $\varepsilon(t)$, and a Brownian motion $Z(t)$.

Proposition A.2.1. The conditional expected value in (A.2) can be approximated by[^117]

$$
\begin{align*}
b(t, x)^{2} & \approx \overline{\lambda(t)^{2}}+s(t)(x-X(0))  \tag{A.8}\\
s(t) & =\frac{\int_{0}^{t} \bar{\varepsilon}(s) \bar{\lambda}(s) \bar{\rho}(s) d s}{\int_{0}^{t} \overline{\lambda(s)^{2}} d s}
\end{align*}
$$
^eqa-8
where $\bar{\varepsilon}(t), \overline{\lambda(t)^{2}}, \bar{\rho}(t), \bar{\lambda}(t)$ are deterministic approximations to $\varepsilon(t), \lambda(t)^{2}$, $\rho(t) \triangleq\langle d W(t), d Z(t)\rangle / d t$, and $\lambda(t)$, respectively. In particular, we can take

$$
\begin{aligned}
& \bar{\varepsilon}(t)=\mathrm{E}(\varepsilon(t)), \quad \overline{\lambda(t)^{2}}=\mathrm{E}\left(\lambda(t)^{2}\right)=\mathrm{E}\left(\int_{0}^{t} \nu(s) d s\right) \\
& \bar{\rho}(t)=\mathrm{E}(\langle d W(t), d B(t)\rangle / d t), \quad \bar{\lambda}(t)=\sqrt{\overline{\lambda(t)^{2}}}
\end{aligned}
$$

Proof. First, we approximate the dynamics of $\left(X(t), \lambda(t)^{2}\right)$ with the Gaussian processes

$$
\begin{aligned}
d X(t) & \approx \bar{\lambda}(t) d W(t) \\
d \lambda(t)^{2} & \approx \bar{\nu}(t) d t+\bar{\varepsilon}(t) d \bar{Z}(t)
\end{aligned}
$$

where $\bar{\nu}(t)=\mathrm{E}(\nu(t))$ and $\bar{Z}(t)$ is a Brownian motion such that $\langle d W(t), d \bar{Z}(t)\rangle=\bar{\rho}(t) d t$. The result then follows from the standard conditioning formula for Gaussian variables $U, V$ :

$$
\begin{equation*}
\mathrm{E}(U \mid V)=\mathrm{E}(U)+\frac{\operatorname{Cov}(U, V)}{\operatorname{Var}(V)}(V-\mathrm{E}(V)) \tag{A.9}
\end{equation*}
$$
^eqa-9

Rather than approximating $X(t)$ and $\lambda(t)^{2}$ directly by Gaussian processes, we can instead use deterministic functions of Gaussian processes. For instance, if we use exponential "mapping functions", we would then arrive at a lognormal (rather than Gaussian) approximation. Furthermore, instead of approximating the drift of $\lambda(t)^{2}$ as a deterministic function, we could instead approximate it by a linear function of $\lambda(t)^{2}$ itself, which would retain a Gaussian distribution for $\lambda(t)^{2}$. Ultimately, the original form of the SDEs for the asset and variance in a given model would typically suggest the most proper usage of the Gaussian approximation principle.

As evident from (A.8), the local variance function that emerges from the Gaussian approximation method is linear, whereby the approximating model (A.3) will always generate a monotonic implied volatility smile. On the other hand, we may know $a$-priori that the true volatility smile of the original model (A.1) is not close to linear - it could be U-shaped, say. In such cases, a wholesale replacement of the original dynamics by the approximating linear local variance model is unlikely to be satisfactory. As it turns out, it is possible to apply the Gaussian approximation in a more sophisticated manner, leading to a better approximating model. We discuss this in Section A.3.1 below, but first we outline an alternative approach to estimating conditional expected values.

## A.2.2 Least-Squares Projection

Section 16.6.2 develops the least-squares projection method for conditional expected values, using the insight that a conditional expected value can be defined as a projection onto a suitable functional space. If we project onto some subspace of the full functional space, we obtain an approximation to the conditional expected value, as stated formally in Proposition 16.6.2. For our purposes here, we focus only on the particularly tractable case of linear subspaces, utilized in Section 16.6.4 to produce a linear least-squares projection. Restating the linear projection in the notations of this appendix. we obtain the following result.

Proposition A.2.2. The linear least-squares approximation to the conditional expected value in ( $A . \mathscr{Q}$ ) is given by

$$
\begin{equation*}
b(t, x)^{2} \approx \overline{\lambda(t)^{2}}+s(t)(x-X(0)) \tag{A.10}
\end{equation*}
$$
^eqa-10

where

$$
\begin{equation*}
s(t)=\frac{\mathrm{E}\left(\lambda(t)^{2}(X(t)-X(0))\right)}{\mathrm{E}\left((X(t)-X(0))^{2}\right)}, \quad \overline{\lambda(t)^{2}}=\mathrm{E}\left(\lambda(t)^{2}\right) \tag{A.11}
\end{equation*}
$$
^eqa-11

or, more compactly,

$$
s(t)=\frac{\operatorname{Cov}\left(\lambda(t)^{2}, X(t)\right)}{\operatorname{Var}(X(t))}
$$

We notice a strong similarity between the expressions for the Dupire local volatility approximations in Propositions A.2.1 and A.2.2. By design (as we projected the variance only on linear functions of $X(t)$ ), the local variance function approximation in (A.10) is still linear in $x$, and the expression (A.11) reduces to the formulas of Proposition A.2.1 if we apply a Gaussian approximation. This is not surprising, as the linear least-squares projection is known to produce an exact result for conditional expectations of Gaussian variables.

As Propositions A.2.1 and A.2.2 approximate local variance with a linear function of spot, both suggest that the SDE for the approximating process is of the displaced square-root type:

$$
d Y(t)=\left(\overline{\lambda(t)^{2}}+s(t)(Y(t)-Y(0))\right)^{1 / 2} d W(t)
$$

While such processes are analytically tractable (see Sections 7.2 .4 and 10.2), it is often more convenient to work with processes of the displaced log-normal type (see Proposition 7.2.12) where the local volatility function is linear:

$$
\begin{equation*}
d Y(t)=\sigma(t)(1+b(t)(Y(t)-Y(0))) d W(t) \tag{A.12}
\end{equation*}
$$
^eqa-12

To obtain an approximating process of this type, we can expand the square root function to the first order around $x=Y(0)$, yielding the following result.

Proposition A.2.3. The displaced log-normal approximation to the process $X(t)$ in (A.1) is given by (A.12), where

$$
\begin{aligned}
\sigma(t) & =\sqrt{\overline{\lambda(t)^{2}}}=\sqrt{\mathrm{E}\left(\lambda(t)^{2}\right)} \\
b(t) & =\frac{s(t)}{2 \sigma(t)^{2}}=\frac{\operatorname{Cov}\left(\lambda(t)^{2}, X(t)\right)}{2 \mathrm{E}\left(\lambda(t)^{2}\right) \operatorname{Var}(X(t))}
\end{aligned}
$$

The same result was obtained in Antonov and Misirpashaev [2009a] by a direct application of the least-squares method, i.e. by solving the minimization problem

$$
\mathrm{E}\left(\left(\lambda(t)^{2}-\sigma(t)^{2}(1+b(t)(X(t)-X(0)))^{2}\right)^{2}\right) \rightarrow \min
$$

in $\sigma(t)$ and $b(t)$ and keeping only the leading $O\left(\lambda(t)^{2}\right)$ terms.

## A. 3 Applications to Local Stochastic Volatility Models

## A.3.1 Markovian Projection onto a Stochastic Volatility Model

In applying the Markovian projection method, we are limited by the accuracy of approximations to conditional expected values. As reviewed in Section A.2, the methods that are generally available approximate local volatility or variance functions with linear ${ }^{3}$ functions which, as discussed in Section A.2.1, is insufficient for the case where $X(t)$ has complex dynamics. Fortunately, Theorem A.1.1 provides us with means to approximate a given model by a model of essentially any type, not just local volatility. The following, borderline trivial, corollary to Theorem A.1.1 is the key ${ }^{4}$.

Corollary A.3.1. If two processes $X_{1}$ and $X_{2}$ have the same Markovian projections (i.e., they imply identical Dupire local volatility functions), then European put/call option prices on $X_{1}$ and $X_{2}$ both are identical for all strikes and expiries.

Let us demonstrate the usefulness of Corollary A.3.1 by applying it to a stochastic volatility model. Let $X_{1}(t)$ follow a stochastic volatility SDE (recall from Section 8.1 that for non-linear functions $b(t, x)$ such models are sometimes called local stochastic volatility, or LSV, models)

$$
d X_{1}(t)=b_{1}\left(t, X_{1}(t)\right) \sqrt{z_{1}(t)} d W(t)
$$

[^118]where $z_{1}(t)$ is some variance process. Suppose we would like to derive approximations for European options on $X_{1}$. One possibility is to approximate $X_{1}$ with a local volatility model, using Theorem A.1.1 directly. As we discussed, this is unlikely to work well, at least if we compute conditional expected values with the approximations developed in Section A.2. Instead, we can use Corollary A.3.1 and approximate $X_{1}$ with a stochastic volatility process that employs a more tractable process for stochastic variance. Let us call this variance process $z_{2}$, and consider a model of the form

$$
\begin{equation*}
d X_{2}(t)=b_{2}\left(t, X_{2}(t)\right) \sqrt{z_{2}(t)} d W(t) . \tag{A.13}
\end{equation*}
$$
^eqa-13

Then Corollary A.3.1 and Theorem A.1.1 imply that to match European option prices in the two models for all strikes and expiries, we need to set $b_{2}(t, x)$ such that

$$
\begin{equation*}
b_{2}(t, x)^{2}=b_{1}(t, x)^{2} \frac{\mathrm{E}\left(z_{1}(t) \mid X_{1}(t)=x\right)}{\mathrm{E}\left(z_{2}(t) \mid X_{2}(t)=x\right)} . \tag{A.14}
\end{equation*}
$$
^eqa-14

While we still need to apply formulas from Section A. 2 to approximate conditional expected values in (A.14), the fact that we calculate the ratio of two expected values gives us some hope for error cancellation - i.e. even if each individual approximation is not particularly accurate, they are inaccurate "in the same way" and the overall error diminishes when the ratio is formed. To maximize the error cancellation effect, it is obviously beneficial to choose $z_{2}$ as similar to $z_{1}$ as possible, while still retaining analytical tractability.

Using the SV model (A.13) as the target model for Markovian projection exercise will benefit from the fact that the model (A.13) is quite rich, even for linear local volatility functions $b_{2}(t, x)$. If $z_{2}$ is a square-root process and $b_{2}(t, x)$ is linear (in $x$ ), then the resulting model reduces to the displaced Heston model (8.3)-(8.4) which, as we saw in Chapter 8, is both tractable and capable of generating a wide variety of implied volatility smiles.

In general, limitations of available conditional expected value approximations impose certain restrictions on designing approximating models. In particular, as we want the "output" local volatility $b_{2}(t, x)$ to be as close to linear as possible - so that the inevitable linear approximation is not far off - we should choose the stochastic variance process $z_{2}$ in such a way that the characteristics of this process, and not the shape of the local volatility, explain as much of the curvature of the implied volatility smile of the model for $X_{1}$ as possible (note that this also holds true for $X_{1}$ processes of the non-SV type).

Let us turn our attention to another common application of Corollary A.3.1. Suppose that $X_{1}$ follows the SDE

$$
\begin{equation*}
d X_{1}(t)=\lambda(t) \sqrt{z(t)} d W(t), \tag{A.15}
\end{equation*}
$$
^eqa-15

where $z(t)$ is a stochastic variance process and $\lambda(t)$ is now a stochastic process in its own right. For instance, $\lambda(t)$ could be a complicated function
of state variables in a term structure model of interest rates, which is a relevant example when $S(t)$ represents a swap rate. We would like to replace the SDE (A.15) with a local stochastic volatility model,

$$
d X_{2}(t)=b\left(t, X_{2}(t)\right) \sqrt{z(t)} d W(t)
$$

where we use the same stochastic variance process $z(t)$ as in (A.15). Then, according to Corollary A.3.1 we need to set

$$
\begin{equation*}
b(t, x)^{2}=\frac{\mathrm{E}\left(\lambda(t)^{2} z(t) \mid X_{1}(t)=x\right)}{\mathrm{E}\left(z(t) \mid X_{2}(t)=x\right)} \tag{A.16}
\end{equation*}
$$
^eqa-16

This formula can be simplified when $\lambda(t)$ and $\Sigma(t)$ are (approximately) conditionally independent given $X_{1}(t)$, in which case we get

$$
\begin{equation*}
b(t, x)^{2} \approx \mathrm{E}\left(\lambda(t)^{2} \mid X_{1}(t)=x\right) \frac{\mathrm{E}\left(z(t) \mid X_{1}(t)=x\right)}{\mathrm{E}\left(z(t) \mid X_{2}(t)=x\right)} \tag{A.17}
\end{equation*}
$$
^eqa-17

In many situations, it can be safely assumed that

$$
\mathrm{E}\left(z(t) \mid X_{1}(t)=x\right) \approx \mathrm{E}\left(z(t) \mid X_{2}(t)=x\right)
$$

in which case the formula simplifies further,

$$
\begin{equation*}
b(t, x)^{2} \approx \mathrm{E}\left(\lambda(t)^{2} \mid X_{1}(t)=x\right) \tag{A.18}
\end{equation*}
$$
^eqa-18

This formula forms the basis for European swaption approximations in term structure models with stochastic volatility; we use it for both quasi-Gaussian models (see Section 13.3.3) and Libor market models (see Proposition 15.2.1).

## A.3.2 Fitting the Market with a Local Stochastic Volatility Model

While the focus of this appendix is on approximating more complicated models with simpler ones, direct calibration of local stochastic volatility (LSV) models to the market is another possible application of the techniques we consider.

Let $S(t)$ be the value of a given market variable (for example, a swap rate or an equity price), with initial value $S(0)=S_{0}$. Suppose market prices of European call or put options are known for all expiries $T$ and strikes $K$. These can be easily converted (see (A.5)) into a market-implied Dupire local volatility $b_{\mathrm{mkt}}(t, x)$, such that the market European option prices are reproduced by the model

$$
d Y(t)=b_{\mathrm{mkt}}(t, Y(t)) d W(t), \quad Y(0)=S_{0}
$$

Suppose we postulate a stochastic variance process $z(t)$ (such as the squareroot process of the Heston model, or whatever multi-factor variance process
is in favor in equity modeling circles at the moment), and aim to construct a stochastic volatility model

$$
\begin{equation*}
d S(t)=b(t, S(t)) \sqrt{z(t)} d W(t), \quad S(0)=S_{0} \tag{A.19}
\end{equation*}
$$
^eqa-19

consistent with market European option prices. As follows from Theorem A.1.1, the local volatility function $b(t, x)$ is then given by

$$
\begin{equation*}
b(t, x)^{2}=\frac{b_{\mathrm{mkt}}(t, x)^{2}}{\mathrm{E}(z(t) \mid S(t)=x)} \tag{A.20}
\end{equation*}
$$
^eqa-20

As mentioned before, the challenge of computing $\mathrm{E}(z(t) \mid S(t)=x)$ makes the method difficult to implement in practice. One choice is to apply finite difference methods to compute the conditional expected value numerically in a forward Kolmogorov PDE for ( $S(t), z(t)$ ), see Ren et al. [2007]. Here, we instead look for analytic approximations.

Linear projections of the type explored in Section A. 2 are possible, but are likely to be inaccurate in this case if $b_{\mathrm{mkt}}(t, x)$ has a high degree of convexity. We instead wish to explore methods based on comparing (A.19) to another stochastic volatility model in which European option prices can be cheaply computed. Suppose we have identified such a "proxy" model, defined through a known local volatility function $\widetilde{b}(t, x)$,

$$
\begin{equation*}
d X(t)=\widetilde{b}(t, X(t)) \sqrt{z(t)} d W(t), \quad X(0)=S_{0} \tag{A.21}
\end{equation*}
$$
^eqa-21

where $z(t)$ is the same process as in (A.19). For tractability, we often choose $\widetilde{b}(t, x)$ to be a linear function of $x$. We will have more to say about the parameterization of the proxy model in Section A.3.3, but for now it suffices to assume that this model allows us to quickly and efficiently calculate European call (or put) option prices. These prices can be turned into a "proxy" Dupire local volatility function $b_{\text {proxy }}(t, x)$ by means of (A.5). Rewriting (A.20), we then have

$$
\begin{equation*}
\mathrm{E}(z(t) \mid X(t)=x)=\frac{b_{\mathrm{proxy}}(t, x)^{2}}{\widetilde{b}(t, x)^{2}} \tag{A.22}
\end{equation*}
$$
^eqa-22

In other words, from a proxy stochastic volatility model which easily computed European option prices, we can efficiently compute the conditional expected values $\mathrm{E}(z(t) \mid X(t)=x)$. One way to take advantage of this observation is to combine (A.20) and (A.22) as follows.

Proposition A.3.2. The local volatility function $b(t, x)$ that makes the model (A.19) consistent with the market is given by

$$
b(t, x)=\widetilde{b}(t, x) \frac{b_{\mathrm{mkt}}(t, x)}{b_{\mathrm{proxy}}(t, x)} \frac{\mathrm{E}(z(t) \mid X(t)=x)}{\mathrm{E}(z(t) \mid S(t)=x)}
$$

where $X(t)$ follows the "proxy" model (A.21) with a known local volatility function $\widetilde{b}(t, x)$.

The ratio $b_{\text {mkt }}(t, x) / b_{\text {proxy }}(t, x)$ can, as discussed, usually be computed efficiently. Approximating

$$
\begin{equation*}
\frac{\mathrm{E}(z(t) \mid X(t)=x)}{\mathrm{E}(z(t) \mid S(t)=x)} \approx 1 \tag{A.23}
\end{equation*}
$$
^eqa-23

we obtain the following useful corollary.
Corollary A.3.3. Under the approximation (A.23) we have that

$$
b(t, x) \approx \tilde{b}(t, x) \frac{b_{\mathrm{mkt}}(t, x)}{b_{\mathrm{proxy}}(t, x)}
$$

To obtain a more sophisticated approximation than that of Corollary A.3.3, we can attempt to improve on (A.23) by looking for an (approximate) functional relationship between $X(t)$ and $S(t)$. Denote

$$
\begin{gather*}
h(t, x)=\int_{x_{0}}^{x} \frac{d y}{b(t, y)}, \quad \widetilde{h}(t, x)=\int_{x_{0}}^{x} \frac{d y}{\widetilde{b}(t, y)}  \tag{A.24}\\
H(t, x)=\widetilde{h}^{-1}(t, h(t, x))
\end{gather*}
$$
^eqa-24

where $\widetilde{h}^{-1}(t, x)$ is the inverse of $\widetilde{h}(t, x)$ in the second (i.e., $\left.x\right)$ argument. Furthermore, denote

$$
\tilde{X}(t)=H(t, S(t)) .
$$

Then

$$
\begin{aligned}
d \widetilde{X}(t) & =\left.\frac{\partial}{\partial x} H(t, x)\right|_{x=S(t)} d S(t)+O(d t) \\
& =\widetilde{b}(t, \widetilde{X}(t)) \sqrt{z(t)} d W(t)+O(d t)
\end{aligned}
$$

We see that $\widetilde{X}(t)$ and $X(t)$ have identical diffusion coefficients, which suggests the approximation

$$
\begin{equation*}
X(t) \approx H(t, S(t)) . \tag{A.25}
\end{equation*}
$$
^eqa-25

This leads to the following result.
Proposition A.3.4 The local volatility function $b(t, x)$ that makes the model (A.19) consistent with the market is approximately given by

$$
\begin{equation*}
b(t, x) \approx \tilde{b}(t, H(t, x)) \frac{b_{\mathrm{mkt}}(t, x)}{b_{\mathrm{proxy}}(t, H(t, x))} \tag{A.26}
\end{equation*}
$$
^eqa-26

with $H(t, x)$ given by (A.24).

Proof. By (A.25),

$$
\begin{aligned}
\mathrm{E}(z(t) \mid S(t)=x) & =\mathrm{E}(z(t) \mid H(t, S(t))=H(t, x)) \\
& \approx \mathrm{E}(z(t) \mid X(t)=H(t, x))
\end{aligned}
$$

From (A.22),

$$
\mathrm{E}(z(t) \mid X(t)=H(t, x))=\frac{b_{\mathrm{proxy}}(t, H(t, x))^{2}}{\tilde{b}(t, H(t, x))^{2}}
$$

and the result follows from (A.20).
We emphasize that $H(t, x)$ depends on the (unknown) function $b(t, x)$, hence (A.26) is, in fact, an equation for $b(t, x)$. This equation can be solved. Let us first denote

$$
\begin{equation*}
h_{\mathrm{proxy}}(t, x)=\int_{x_{1}}^{x} \frac{d y}{b_{\mathrm{proxy}}(t, y)}, \quad h_{\mathrm{mkt}}(t, x)=\int_{x_{1}}^{x} \frac{d y}{h_{\mathrm{mkt}}(t, y)} \tag{A.27}
\end{equation*}
$$
^eqa-27

The following then holds (see Henry-Labordére [2009]).
Proposition A.3.5. The mapping function $H(t, x)$ is given by

$$
\begin{equation*}
H(t, x)=h_{\text {proxy }}^{-1}\left(t, h_{\mathrm{mkt}}(t, x)\right) \tag{A.28}
\end{equation*}
$$
^eqa-28

where $h_{\text {proxy }}^{-1}(t, x)$ is the inverse of $h_{\text {proxy }}(t, x)$ defined by (A.27) in the second (x) argument. Furthermore,

$$
\begin{equation*}
b(t, x) \approx \tilde{b}\left(t, h_{\mathrm{proxy}}^{-1}\left(t, h_{\mathrm{mkt}}(t, x)\right)\right) \frac{b_{\mathrm{mkt}}(t, x)}{b_{\mathrm{proxy}}\left(t, h_{\mathrm{proxy}}^{-1}\left(t, h_{\mathrm{mkt}}(t, x)\right)\right)} \tag{A.29}
\end{equation*}
$$
^eqa-29

Proof. Differentiating $H(t, x)$ in (A.24) with respect to $x$ we obtain

$$
\frac{\partial H(t, x)}{\partial x}=\frac{\partial h(t, x)}{\partial x} /\left.\frac{\partial \widetilde{h}(t, f)}{\partial f}\right|_{f=H(t, x)}=\frac{\tilde{b}(t, H(t, x))}{b(t, x)} .
$$

Therefore, we can rewrite (A.26) as an (approximate) equation

$$
\begin{equation*}
\frac{\partial H(t, x)}{\partial x}=\frac{b_{\mathrm{proxy}}(t, H(t, x))}{b_{\mathrm{mkt}}(t, x)} \tag{A.30}
\end{equation*}
$$
^eqa-30

Treating this as an ODE in $x$ for fixed $t$, we solve it to find

$$
\begin{equation*}
\int_{x_{0}}^{H(t, x)} \frac{d y}{b_{\mathrm{proxy}}(t, y)}=\int_{x_{1}}^{x} \frac{d y}{b_{\mathrm{mkt}}(t, y)} \tag{A.31}
\end{equation*}
$$
^eqa-31

resulting in (A.28). Then (A.29) follows from (A.26) and (A.28).

Remark A.3.6. Henry-Labordére [2009] notices that in fact (A.30) implies a condition more general than (A.31), namely that

$$
\int_{H_{0}}^{H(t, x)} \frac{d y}{b_{\mathrm{proxy}}(t, y)}=\int_{x_{11}}^{x} \frac{d y}{b_{\mathrm{mkt}}(t, y)}
$$

for any $H_{0}$. He proposes to choose $H_{0}$ so that the difference in drifts of $S(t)$ and $X(t)$ is minimized (our approximation above matches diffusion terms only).

## A.3.3 On Calculating Proxy Local Volatility

The techniques in Section A.3.2 above all hinge on the critical assumption that we can pick a proxy stochastic volatility model that allows for efficient computation of the Dupire local volatility function $b_{\text {proxy }}(t, x)$. Provided that $z(t)$ follows the standard mean-reverting square-root process (as in the Heston model), then an obvious choice for the proxy model is a displaced Heston model, with

$$
\begin{equation*}
\tilde{b}(t, x)=\tilde{b}_{1} x+\tilde{b}_{2} \tag{A.32}
\end{equation*}
$$
^eqa-32

Special cases include the $\widetilde{b}(t, x)=\widetilde{b}_{1} x$ (the original Heston model) and $\widetilde{b}(t, x)=\widetilde{b}_{2}$ (the "Gaussian" Heston model), but any choice of constants $\widetilde{b}_{1}$ and $\tilde{b}_{2}$ will allow for quick pricing of European put and call options; see Chapter 8 for details. Using the averaging techniques of Chapter 9, we can, in fact, extend (A.32) to linear local volatility function with time-dependent coefficients,

$$
\begin{equation*}
\widetilde{b}(t, x)=\widetilde{b}_{1}(t) x+\widetilde{b}_{2}(t) \tag{A.33}
\end{equation*}
$$
^eqa-33

The time-dependent coefficients in (A.33) can be chosen to make the proxy model resemble as much as possible the true model for $S$, thereby improving the quality of various approximations made. Ideally we should use

$$
\tilde{b}_{2}(t)=b\left(t, S_{0}\right), \quad \tilde{b}_{1}(t)=\frac{\partial}{\partial x} b\left(t, S_{0}\right)
$$

which is the first-order approximation to the Dupire local volatility $b(t, x)$ along the forward value of $S(t)$. Of course, the value and derivative of $b(t, x)$ are unknown a-priori, but one can easily envision various approximations or, perhaps, an iterative procedure where an approximation for the Dupire local volatility in step $n$ is used to define the proxy local volatility $\widetilde{b}$ for step $n+1$.

The choice of a mean-reverting square-root process for $z(t)$ for variance leads to a specification that is quite amenable to the methods of Section A.3.2 and is often sufficient. However, more advanced applications such as those considered in Section 15.7 (and some popular models in equity modeling, see e.g. Bergomi [2009]) involve multi-factor stochastic variance
dynamics and require extra effort as European options are then not always easy to calculate or approximate. As should be clear from Section A.3.2, we actually do not need to be able to calculate European option prices in the proxy model; all we really need is the Dupire local volatility for the proxymodel, $b_{\text {proxy }}(t, x)$. This function is defined by (see (A.22))

$$
\begin{equation*}
b_{\text {proxy }}(t, x)^{2}=\widetilde{b}(t, x)^{2} \mathrm{E}(z(t) \mid X(t)=x) \tag{A.34}
\end{equation*}
$$
^eqa-34

It turns out that there exists a reasonably efficient algorithm for calculating the right-hand side of this question for a large selection of proxy models. To make this statement precise, let us define the proxy stochastic volatility. model by

$$
\begin{align*}
d X(t) / X(t) & =\lambda \sqrt{z(t)}\left(\rho d W(t)+\left(1-\rho^{2}\right)^{1 / 2} d W_{X}(t)\right)  \tag{A.35}\\
d z(t) & =\nu(t) d t+\varepsilon_{1}(t) d W(t)+\varepsilon_{2}(t) d W_{z}(t) \tag{A.36}
\end{align*}
$$
^eqa-35
^eqa-36

for a deterministic $\lambda>0$ (in the notation of (A.21) therefore $\widetilde{b}(t, x)=\lambda x$ ). Here $W(t), W_{X}(t)$ and $W_{z}(t)$ are independent Brownian motions and $\nu(t)$, $\varepsilon_{1}(t)$ and $\varepsilon_{2}(t)$ are sufficiently regular adapted processes. The algorithm is based on the following result (see Romano and Touzi [1997] and Lee [2001]).

Proposition A.3.7. For the model (A.35)-(A.36), we have

$$
\begin{equation*}
\mathrm{E}(z(t) \mid X(t)=x)=\frac{\mathrm{E}(z(t) \xi(t, x))}{\mathrm{E}(\xi(t, x))}, \quad \xi(t, x)=\frac{1}{D(t)^{1 / 2}} \phi\left(\frac{\ln (x)-m(t)}{D(t)^{1 / 2}}\right), \tag{A.37}
\end{equation*}
$$
^eqa-37

where

$$
\begin{aligned}
& m(t)=\ln (X(0))+\lambda \rho \int_{0}^{t} \sqrt{z(s)} d W(s)-\frac{\lambda^{2}}{2} \int_{0}^{t} z(s) d s \\
& D(t)=\lambda^{2}\left(1-\rho^{2}\right) \int_{0}^{t} z(s) d s
\end{aligned}
$$

and $\phi(z)=(2 \pi)^{-1 / 2} e^{-z^{2} / 2}$ is the standard Gaussian PDF.
Proof. Proceeding informally, we observe that

$$
\begin{equation*}
\mathrm{E}(z(t) \mid X(t)=x)=\frac{\mathrm{E}(z(t) \delta(X(t)-x))}{\mathrm{E}(\delta(X(t)-x))}, \tag{A.38}
\end{equation*}
$$
^eqa-38

where $\delta(x)$ is the Dirac delta function. Let $\mathcal{F}_{t}^{z}$ be the filtration generated by $W(s)$ and $W_{z}(s), 0 \leq s \leq t$. Clearly $z(t)$ is adapted to $\mathcal{F}_{t}^{z}$. We can write

$$
\begin{aligned}
X(t)=X(0) \exp \left(\lambda \rho \int_{0}^{t} \sqrt{z(s)} d W(s)-\frac{\lambda^{2}}{2} \int_{0}^{t} z(s) d s\right) & \\
& \quad \times \exp \left(\lambda\left(1-\rho^{2}\right)^{1 / 2} \int_{0}^{t} \sqrt{z(s)} d W_{X}(s)\right)
\end{aligned}
$$

where the first exponential is adapted to $\mathcal{F}_{t}^{z}$ and the second is driven by a Brownian motion $W_{\mathrm{X}}$ that is independent of this filtration. Conditioned on the filtration $\mathcal{F}_{t}^{z}, \ln (X(t))$ is Gaussian with known moments,

$$
\begin{equation*}
\ln (X(t)) \mid \mathcal{F}_{t}^{z} \sim \mathcal{N}(m(t), D(t)) \tag{A.39}
\end{equation*}
$$
^eqa-39

Notice that

$$
\mathrm{E}(\delta(X(t)-x))=\mathrm{E}\left(\mathrm{E}\left(\delta(X(t)-x) \mid \mathcal{F}_{t}^{z}\right)\right)
$$

where $\mathrm{E}\left(\delta(X(t)-x) \mid \mathcal{F}_{t}^{2}\right)$ is the $\log$-normal density for the conditional distribution of $X(t)$ defined by (A.39). We therefore have

$$
\mathrm{E}\left(\mathrm{E}\left(\delta(X(t)-x) \mid \mathcal{F}_{t}^{z}\right)\right)=\mathrm{E}\left(\frac{1}{x D(t)^{1 / 2}} \phi\left(\frac{\ln (x)-m(t)}{D(t)^{1 / 2}}\right)\right)
$$

and

$$
\begin{aligned}
\mathrm{E}(z(t) \delta(X(t)-x)) & =\mathrm{E}\left(z(t) \mathrm{E}\left(\delta(X(t)-x) \mid \mathcal{F}_{t}^{z}\right)\right) \\
& =\mathrm{E}\left(\frac{z(t)}{x D(t)^{1 / 2}} \phi\left(\frac{\ln (x)-m(t)}{D(t)^{1 / 2}}\right)\right)
\end{aligned}
$$

The result of the proposition then follows from (A.38).
To compute $b_{\text {proxy }}(t, x)$ in a proxy model of the type (A.35)-(A.36), Henry-Labordére [2009] suggests simulating paths of the variance process $z(t)$ and then tabulating the values of $\mathrm{E}(z(t) \mid X(t)=x)$ by calculating the formula (A.37) for a selection of values of $x$ and $t$. To obtain the values of $b_{\text {proxy }}(t, x)$ for all $t$ and $x$, we would then use (A.34) with some sort of interpolation to fill in values of $\mathrm{E}(z(t) \mid X(t)=x)$ between the tabulated ones.

## A. 4 Basket Options in Local Volatility Models

So far, our primary application for Markovian projection has been the pricing of European options on scalar processes. However, Markovian projection can easily be extended to options on multiple underlyings, e.g. basket options. Such options appear naturally in interest rate modeling - for example, we can often think of swap rates as baskets of Libor rates - a representation that we use in Section 15.2 together with a Markovian projection to approximate the swap rate process in a Libor market model. CMS spread options could also be thought of as options on baskets of two assets with weights $\pm 1$. Furthermore, baskets serve as a good example of the practical usage of the Markovian projection method and the various "tricks" that go along with it.

Let us consider a collection of $N$ assets $\mathbf{S}(t)=\left(S_{1}(t), \ldots, S_{N}(t)\right)^{\top}$, each driven by its own local volatility model

$$
\begin{equation*}
d S_{n}(t)=\sigma_{n}(t) \varphi_{n}\left(t, S_{n}(t)\right) d W_{n}(t), \quad n=1, \ldots, N \tag{A.40}
\end{equation*}
$$
^eqa-40

$$
\begin{equation*}
\bar{\lambda}(t)^{2} \triangleq \sum_{n, m=1}^{N} w_{n} w_{m} \sigma_{n}(t) \sigma_{m}(t) \rho_{n, m} . \tag{A.48}
\end{equation*}
$$
^eqa-48

In the same spirit, we approximate $W(t)$ in (A.45) by $\bar{W}(t)$, where

$$
d \bar{W}(t) \triangleq \frac{1}{\bar{\lambda}(t)} \sum_{n=1}^{N} w_{n} \sigma_{n}(t) d W_{n}(t)
$$

is also a Brownian motion. Finally, let us denote

$$
\begin{equation*}
\rho_{n}(t) \triangleq\left\langle d W_{n}(t), d \bar{W}(t)\right\rangle / d t=\frac{1}{\bar{\lambda}(t)} \sum_{m=1}^{N} w_{m} \sigma_{m}(t) \rho_{n, m} . \tag{A.49}
\end{equation*}
$$
^eqa-49

With all this notation in place, we then can state the following proposition.
Proposition A.4.1. Let the skew functions $\varphi_{n}$ be as stated in (A.41). A displaced log-normal approximation to the dynamics of the basket $S(t)$ in (A.43) is then given by (A.46) where

$$
\begin{align*}
\sigma(t)^{2} & \approx \sum_{n=1}^{N} \sum_{m=1}^{N} w_{n} w_{m} \sigma_{n}(t) \sigma_{m}(t) \rho_{n, m}  \tag{A.50}\\
b(t) & \approx \frac{\sum_{n=1}^{N} \bar{\lambda}(t) \sigma_{n}(t) \rho_{n}(t)\left(\int_{0}^{t} \bar{\lambda}(s) \sigma_{n}(s) \rho_{n}(s) d s\right) w_{n} b_{n}(t)}{\bar{\lambda}(t)^{2} \int_{0}^{t} \bar{\lambda}(s)^{2} d s}, \tag{A.51}
\end{align*}
$$
^eqa-50
^eqa-51

with the $\rho_{n}$ 's defined in (A.49) and $\bar{\lambda}(t)$ defined in (A.48).
Proof. The approximation (A.50) follows from

$$
\sigma(t)^{2}=\mathrm{E}\left(\lambda(t)^{2}\right) \approx \mathrm{E}\left(\bar{\lambda}(t)^{2}\right)=\bar{\lambda}(t)^{2}
$$

and (A.48). To find the covariance term in (A.47), we recall (A.41) and write $\operatorname{Cov}\left(\lambda(t)^{2}, S(t)\right)$

$$
\begin{aligned}
= & \sum_{n, m=1}^{N} w_{n} w_{m} \sigma_{n}(t) \sigma_{m}(t) \rho_{n, m} \\
& \times \operatorname{Cov}\left(\left(1+b_{n}(t)\left(S_{n}(t)-S_{n}(0)\right)\right)\left(1+b_{m}(t)\left(S_{m}(t)-S_{m}(0)\right)\right), S(t)\right) \\
\approx & 2 \bar{\lambda}(t) \sum_{n=1}^{N} w_{n} \sigma_{n}(t) \rho_{n}(t) b_{n}(t) \operatorname{Cov}\left(S_{n}(t), S(t)\right)
\end{aligned}
$$

where we disregarded the terms $\operatorname{Cov}\left(\left(S_{n}(t)-S_{n}(0)\right)\left(S_{m}(t)-S_{m}(0)\right), S(t)\right)$ as being of higher order in volatility. Furthermore,
$\operatorname{Cov}\left(S_{n}(t), S(t)\right)$

$$
\begin{align*}
& =\sum_{m=1}^{N} w_{m} \operatorname{Cov}\left(S_{n}(t), S_{m}(t)\right)  \tag{A.52}\\
& =\sum_{m=1}^{N} \int_{0}^{t} w_{m} \rho_{n . n} \sigma_{n}(s) \sigma_{m}(s) \mathrm{E}\left(\varphi_{n}\left(s, S_{n}(s)\right) \varphi_{m}\left(s, S_{m}(s)\right)\right) d s \\
& \approx \int_{0}^{t} \sigma_{n}(s) \bar{\lambda}(s) \rho_{n}(s) d s
\end{align*}
$$
^eqa-52

To summarize,

$$
\operatorname{Cov}\left(\lambda(t)^{2}, S(t)\right) \approx 2 \bar{\lambda}(t) \sum_{n=1}^{N} w_{n} \sigma_{n}(t) \rho_{n}(t) b_{n}(t) \int_{0}^{t} \sigma_{n}(s) \bar{\lambda}(s) \rho_{n}(s) d s
$$

In the same spirit, from (A.52),

$$
\begin{aligned}
\operatorname{Var}(S(t)) & =\sum_{n=1}^{N} w_{n} \operatorname{Var}\left(S_{n}(t), S(t)\right) \\
& \approx \sum_{n=1}^{N} w_{n} \int_{0}^{t} \sigma_{n}(s) \bar{\lambda}(s) \rho_{n}(s) d s=\int_{0}^{t} \bar{\lambda}(s)^{2} d s
\end{aligned}
$$

and the proposition follows from (A.47). $\square$

Remark A.4.2. Proposition A.4.1 depends on numerous ad-hoc approximations, the nature of which we did not characterize rigorously. A more detailed analysis can be found in Antonov and Misirpashaev [2009a] where it is shown that our expressions for $\sigma(t)$ and $b(t)$ are leading order terms in the smallvolatility limit.

The parameter $b(t)$ in (A.46) represents the slope of the local volatility function for the basket $S$, and the expression (A.51) relates this slope to a weighted average of the slopes $b_{n}(t)$ of the individual volatility functions for the basket components. This approximation works best when the skews $b_{n}(t)$ are of the same sign and the weights $w_{n}$ are all positive. For skews $b_{n}$ (or weights $w_{n}$ ) of mixed signs, the approximation is not very accurate, however. This should be intuitively clear, since the difference of two processes with positive skews (say) can easily have a $U$-shaped smile, which is obviously not well-approximated by a projection of the difference onto a displaced log-normal process. To handle this case, one possible solution is to use a projection on a stochastic volatility process (even though the components of the basket are local volatility processes), as described in Section A. 5 below. Alternatively, we can use a projection on a local volatility process with a quadratic local volatility as mentioned in footnote 3. Finally, we can always fall back on the copula-based methods of Chapter 17.

## A. 5 Basket Options in Stochastic Volatility Models

We continue investigating basket options, but now augment the model (A.40), (A.41) with stochastic volatility. Specifically, we replace (A.40) with

$$
d S_{n}(t)=\sigma_{n}(t) \varphi_{n}\left(t, S_{n}(t)\right) \sqrt{z_{n}(t)} d W_{n}(t), \quad n=1, \ldots, N
$$

where (A.41) still holds, and where individual stochastic variance processes $z_{n}(t)$ are defined by

$$
\begin{equation*}
d \tilde{z}_{n}(t)=\theta_{n}(t)\left(1-z_{n}(t)\right) d t+\eta_{n}(t) \sqrt{z_{n}(t)} d W_{N+n}(t), \quad z_{n}(0)=1 \tag{A.53}
\end{equation*}
$$
^eqa-53

$n=1, \ldots, N$. Here, $\left(W_{1}(t), \ldots, W_{2 N}\right)$ is a $2 N$-dimensional Brownian motion with correlations

$$
\left\langle d W_{i}, d W_{\jmath}\right\rangle=\rho_{\imath, j} d t, \quad i, j=1, \ldots, 2 N
$$

To simplify the already cumbersome notation, let us absorb the basket weights $w_{n}$ into a redefinition of the asset processes: $S_{n} \leftarrow w_{n} S_{n}$. The basket value $S(t)$ is now

$$
S(t)=\sum_{n=1}^{N} S_{n}(t)
$$

where we now have

$$
d S(t)=\sum_{n=1}^{N} \sigma_{n}(t) \varphi_{n}\left(t, S_{n}(t)\right) \sqrt{z_{n}(t)} d W_{n}(t)=\lambda(t) d W(t)
$$

with

$$
\begin{align*}
\lambda(t)^{2} & =\sum_{n, m=1}^{N} \sigma_{n}(t) \varphi_{n}\left(t, S_{n}(t)\right) \sigma_{m}(t) \varphi_{m}\left(t, S_{m}(t)\right) \rho_{n, m} \sqrt{z_{n}(t) z_{m}(t)} \\
d W(t) & =\frac{1}{\lambda(t)} \sum_{n=1}^{N} \sigma_{n}(t) \varphi_{n}\left(t, S_{n}(t)\right) \sqrt{z_{n}(t)} d W_{n}(t) \tag{A.54}
\end{align*}
$$
^eqa-54

$S$ above is driven by a multi-dimensional stochastic volatility process, and, as we discussed previously in Section A.3.1, projecting $S$ on a displaced log-normal local volatility process is unlikely to lead to accurate option approximations. Following Antonov et al. [2009], we instead investigate projections on a displaced Heston process.

Let us first assume that the skew parameter $b(t)$ of the target approximation is given exogenously (we will discuss its computation later in this section). With $b(t)$ given, we rewrite the dynamics of the process $S(t)$ in a way more suitable for approximations:

$$
\begin{equation*}
d S(t)=\sigma(t)(1+b(t)(S(t)-S(0))) \sqrt{z(t)} d W(t) \tag{A.55}
\end{equation*}
$$
^eqa-55

where

$$
\begin{equation*}
z(t) \triangleq \frac{\Lambda(t)^{2}}{\sigma(t)^{2}}, \quad \Lambda(t)^{2} \triangleq \frac{\lambda(t)^{2}}{(1+b(t)(S(t)-S(0)))^{2}} \tag{A.56}
\end{equation*}
$$
^eqa-56

and

$$
\begin{equation*}
\sigma(t)^{2} \triangleq \mathrm{E}\left(\Lambda(t)^{2}\right) . \tag{A.57}
\end{equation*}
$$
^eqa-57

For future reference, we apply Ito's lemma to $z(t)$ and write

$$
\begin{equation*}
d z(t)=\nu(t) d t+\varepsilon(t) \sqrt{z(t)} d Z(t) \tag{A.58}
\end{equation*}
$$
^eqa-58

where $Z(t)$ is a Brownian motion such that $\langle d W(t), d Z(t)\rangle=\chi(t) d t$, where the exact form of stochastic processes $\nu(t), \varepsilon(t), \chi(t)$ is not important for the moment.

By the multi-dimensional extension of Gyöngy's theorem in Proposition A.1.4, we replicate the exact distribution of the pair $(S(t), z(t))$ for each $t \geq 0$ with $(\widetilde{S}(t), \widetilde{\tilde{z}}(t))$, where

$$
\begin{align*}
d \widetilde{S}(t)= & \sigma(t)(1+b(t)(\widetilde{S}(t)-\widetilde{S}(0))) \sqrt{\widetilde{z}(t)} d W(t) \\
d \widetilde{z}(t)= & \mathrm{E}(\nu(t) \mid S(t)=\widetilde{S}(t), z(t)=\widetilde{z}(t)) d t \\
& +\left(\mathrm{E}\left(\varepsilon(t)^{2} \mid S(t)=\widetilde{S}(t), z(t)=\widetilde{z}(t)\right)\right)^{1 / 2} d Z(t) \tag{A.59}
\end{align*}
$$
^eqa-59

We cannot calculate the conditional expectations in (A.59) exactly, so we proceed to assume a particular parametric form for the process $(\widetilde{S}(t), \widetilde{z}(t))$,

$$
\begin{align*}
d \widetilde{S}(t) & \approx \sigma(t)(1+b(t)(\widetilde{S}(t)-\widetilde{S}(0))) \sqrt{\widetilde{z}(t)} d W(t)  \tag{A.60}\\
d \widetilde{z}(t) & \approx \theta(t)(1-\widetilde{z}(t)) d t+\eta(t) \sqrt{\widetilde{z}(t)} d Z(t) \tag{A.61}
\end{align*}
$$
^eqa-60
^eqa-61

with $\langle d W(t), d Z(t)\rangle=\rho(t) d t$. The unknown parameters $\theta(t), \eta(t)$ and $\rho(t)$ emerge as solutions to the following optimization problems,

$$
\begin{align*}
\mathrm{E}\left(\left(\varepsilon(t)^{2}-\eta(t)^{2} z(t)\right)^{2}\right) & \rightarrow \min \\
\mathrm{E}\left((\varepsilon(t) \chi(t) z(t)-\eta(t) \rho(t) z(t))^{2}\right) & \rightarrow \min \\
\mathrm{E}\left((\nu(t)-\theta(t)(1-z(t)))^{2}\right) & \rightarrow \min \tag{A.62}
\end{align*}
$$
^eqa-62

From the usual first-order optimality conditions, we then find

$$
\begin{align*}
\eta(t)^{2} & =\frac{\mathrm{E}\left(\varepsilon(t)^{2} z(t)\right)}{\mathrm{E}\left(z(t)^{2}\right)} \\
\rho(t) & =\frac{\mathrm{E}\left(\varepsilon(t) \chi(t) z(t)^{2}\right)}{\eta(t) \mathrm{E}\left(z(t)^{2}\right)} \\
\theta(t) & =\frac{\mathrm{E}(\nu(t)(1-z(t)))}{\mathrm{E}\left((1-z(t))^{2}\right)} \tag{A.63}
\end{align*}
$$
^eqa-63

The expectations required in these equations can be approximated from the definition of $z(t)$ in (A.56) and the coefficients $\nu(t), \varepsilon(t), \chi(t)$ in (A.58). The (laborious) calculations are performed in Antonov et al. [2009], and we do not reproduce them here.

We have not yet specified how to set the skew function $b(t)$ in (A.60), originally appearing in (A.55) and (A.56). One idea here is to use the results from the local volatility approximation in Section A. 4 and simply use (A.51). Another alternative suggested in Antonov et al. [2009] finds $b(t)$ by minimizing the defect $D_{\theta}(t)$ for the solution (A.63) of the problem (A.62),

$$
\begin{equation*}
D_{\theta}(t)=\mathrm{E}\left(\nu(t)^{2}\right)-\frac{(\mathrm{E}(\nu(t)(1-z(t))))^{2}}{\mathrm{E}\left((1-z(t))^{2}\right)} \tag{A.64}
\end{equation*}
$$
^eqa-64

The defect $D_{\theta}(t)$ measures the error in the objective function in (A.62) for the solution (A.63), and clearly is a function of $b(t)$. By minimizing the defect, an ODE for $b(t)$ arises that can be solved in closed form. Again, we refer the interested reader to Antonov et al. [2009] for details.

Our final topic in this section is the calculation of the effective volatility $\sigma(t)$ in (A.57). Among all parameters of the model (A.60)-(A.61), $\sigma(t)$ typically is the one that has the biggest impact on the quality of approximations. From (A.54) and (A.56), the standard "freezing" (where $S_{n}(t)=S_{n}(0)$ for all $n=1, \ldots, N$ ) approximation gives us

$$
\begin{equation*}
\sigma(t)^{2}=\sum_{n, m=1}^{N} \sigma_{n}(t) \sigma_{m}(t) \rho_{n, m} \mathrm{E}\left(\sqrt{z_{n}(t) z_{m}(t)}\right) \tag{A.65}
\end{equation*}
$$
^eqa-65

In principle, the same freezing idea could be applied to $z_{n}$ 's, leading us to an approximation that is identical to the one we already derived for the local volatility case, see (A.50). This is, indeed, the leading term of the expansion of $\mathrm{E}\left(\Lambda(t)^{2}\right)$ in small volatilities (see Antonov et al. [2009]), but for typical parameter settings encountered in interest rate modeling the quality of approximations based on this choice for $\sigma(t)$ is rather poor - volatility of variance parameters $\eta_{n}$ are often simply too far from being "small" (larger than 1 is typical).

To improve the accuracy of the overall Markovian projection onto a displaced Heston model, we need to find a way to calculate $\sigma(t)$ in (A.65) without the assumption of small variances of $z_{n}$ 's. Clearly, our ability to do so hinges on accurate approximations to $\mathrm{E}\left(\sqrt{z_{n}(t) \tilde{z}_{m}(t)}\right)$ in (A.65), which could be interesting for other purposes as well. We discuss such a "non-perturbative" approximation in Appendix A.A (see Proposition A.A.1) where we also consider a related problem of approximating $\mathrm{E}\left(\sqrt{z_{n}(t)}\right)$ (Lemma A.A.2).

With these approximations in place, Antonov et al. [2009] demonstrate excellent performance of Markovian projection onto a displaced Heston for basket and spread options.

To conclude, let us quickly summarize the entire algorithm. First, we calculate a non-perturbative approximation to $\sigma(t)$ given by (A.65), with
the square roots obtained by the methods of Appendix A.A. Second, we calculate the optimal skew function $b(t)$ by minimizing the function $D_{\theta}(t)$ in (A.64). As it turns out, this minimization problem (or, rather, a collection of minimization problems indexed by $t$ ) leads to a first-order ODE on $b(t)$ with coefficients that only depend on $\sigma(t)$ (in the small volatility limit), as derived in Antonov et al. [2009]. Finally, having established $\sigma(t)$ and $b(t)$ : we can now solve for optimal coefficients for the stochastic variance process (and its correlation with the asset process) by solving (A.62) for each $t>0$.

## A.A Appendix: Approximations for $\mathrm{E}\left(\sqrt{z_{n}(t) z_{m}(t)}\right)$ and $\mathrm{E}\left(\sqrt{z_{n}(t)}\right)$

As discussed in Section A.5, European option approximations in multidimensional SV models require calculations of certain expected values of square-root processes, a subject we consider in this section. Let us simplify the notations of (A.53) a little, and consider a two-dimensional square-root process

$$
\begin{equation*}
d z_{i}(t)=\theta_{i}\left(1-z_{i}(t)\right) d t+\eta_{i} \sqrt{z_{i}(t)} d W_{i}(t), \quad z_{i}(0)=1 \tag{A.66}
\end{equation*}
$$
^eqa-66

where $i=1,2$, and $\left\langle d W_{1}(t), d W_{2}(t)\right\rangle=\rho d t$. We first consider the problem of approximating $\mathrm{E}\left(\sqrt{z_{1}(t) z_{2}(t)}\right)$.

Proposition A.A.1. For a two-dimensional square-root process (A.66), let us define

$$
\begin{align*}
E\left(\theta_{1}, \theta_{2}, \eta_{1}, \eta_{2}, \rho\right)=\left(1+\rho \eta_{1} \eta_{2} q_{1} q_{2} \frac{1-e^{-\left(\theta_{1}+\theta_{2}-\rho \eta_{1} \eta_{2} q_{1} q_{2}\right) t}}{\theta_{1}+\theta_{2}-\rho \eta_{1} \eta_{2} q_{1} q_{2}}\right) & \\
& \times \mathrm{E}\left(\sqrt{z_{1}(t)}\right) \mathrm{E}\left(\sqrt{z_{2}(t)}\right) \tag{A.67}
\end{align*}
$$
^eqa-67

where $q_{i}, i=1,2$, are obtained by solving

$$
1=\frac{\left(\mathrm{E}\left(\sqrt{z_{i}(t)}\right)\right)^{2}}{2 \theta_{i}-\eta_{i}^{2} q_{i}^{2}}\left(2 \theta_{i}-\eta_{i}^{2} q_{i}^{2} e^{-\left(2 \theta_{i}-\eta_{i}^{2} q_{2}^{2}\right) t}\right) .
$$

This function gives an approximation to $\mathrm{E}\left(\sqrt{z_{1}(t) z_{2}(t)}\right)$, which is

1. Exact in the limit $\rho=0$.
2. Exact in the limit $\rho=1, \theta_{1}=\theta_{2}, \eta_{1}=\eta_{2}$.
3. Has correct leading behavior in the expansion in powers of $\eta_{1}$ and $\eta_{2}$.

The proposition is proved in Section A.A.1. We note that enforcement of the first two non-perturbative conditions substantially improves the accuracy of the approximation to $\sigma(t)$ in (A.65) when $\eta_{1}, \eta_{2} \geq 1$.

The approximation (A.67) relies on our ability to calculate $\mathrm{E}\left(\sqrt{z_{i}(t)}\right)$, the calculation that is also of independent interest sometimes. To state the result, we assume that $z(t)$ follows the square-root process (8.4):

$$
\begin{equation*}
d z(t)=\theta\left(z_{0}-z(t)\right) d t+\eta \sqrt{z(t)} d Z(t), \quad z(0)=z_{0} \tag{A.68}
\end{equation*}
$$
^eqa-68

and recall the definitions (8.6):

$$
d=4 \theta z_{0} / \eta^{2}, \quad n(t, T)=\frac{4 \theta e^{-\theta(T-t)}}{\eta^{2}\left(1-e^{-\theta(T-t)}\right)}
$$

The following result, proven in Section A.A.2, derives the required representation.

Lemma A.A.2. In the model (A.68) we have

$$
\begin{equation*}
\mathrm{E}(\sqrt{z(t)})=\left(\frac{2 e^{-\theta t}}{n(0, t)}\right)^{1 / 2} e^{-n(0, t) / 2} \sum_{j=0}^{\infty} \frac{(n(0, t) / 2)^{j}}{j!} \frac{\Gamma(d / 2+j+1 / 2)}{\Gamma(d / 2+j)} \tag{A.69}
\end{equation*}
$$
^eqa-69

where $\Gamma(z)$ is the Gamma function.
Remark A.A.3. The series in (A.69) converges rapidly, so only the first few terms need to be computed. In each term, the ratio of Gamma functions can be evaluated by standard algorithms present in most numerical software packages. We find that the following approximation is sufficient for our purposes:

$$
\Gamma(x+1 / 2) / \Gamma(x)= \begin{cases}x\left(0.5619 x^{2}-1.3353 x+1.6651\right), & x \in[0,0.9] \\ \sqrt{x}\left(1-\frac{1}{8 x}+\frac{1}{128 x^{2}}\right), & x>0.9\end{cases}
$$

where the first line is obtained by fitting $\Gamma(x+1 / 2) / \Gamma(x)$ with a third-degree polynomial over the interval $[0,1]$, and the second line is the truncation of the series expansion of $\Gamma(x+1 / 2) / \Gamma(x)$ valid for $x>1$, see Cevher et al. [2007]. The cut-off point of 0.9 is chosen to make the function continuous and (nearly) $C^{1}$.

## A.A. 1 Proof of Proposition A.A. 1

## A.A.1.1 Step 1. Reduction to Covariance

Taking expected values of both sides of (A.66) we see that $\mathrm{E}\left(z_{i}(t)\right)=1, \quad i=$ 1,2 . Furthermore, we note that

$$
\begin{aligned}
d\left(z_{1}(t) z_{2}(t)\right)= & z_{1}(t)\left(\theta_{2}\left(1-z_{2}(t)\right) d t+\eta_{2} \sqrt{z_{2}(t)} d W_{2}(t)\right) \\
& +z_{2}(t)\left(\theta_{1}\left(1-z_{1}(t)\right) d t+\eta_{1} \sqrt{z_{1}(t)} d W_{1}(t)\right) \\
& +\eta_{1} \eta_{2} \rho \sqrt{z_{1}(t) z_{2}(t)} d t
\end{aligned}
$$

so that

$$
\begin{equation*}
\frac{d}{d t} \mathrm{E}\left(z_{1}(t) z_{2}(t)\right)=\left(\theta_{1}+\theta_{2}\right)\left(1-\mathrm{E}\left(z_{1}(t) z_{2}(t)\right)\right)+\eta_{1} \eta_{2} \rho \mathrm{E}\left(\sqrt{z_{1}(t) z_{2}(t)}\right) \tag{A.70}
\end{equation*}
$$
^eqa-70

Hence, calculating $\mathrm{E}\left(\sqrt{z_{1}(t) z_{2}(t)}\right)$ is equivalent to calculating $\mathrm{E}\left(z_{1}(t) z_{2}(t)\right)$.

## A.A.1.2 Step 2. Linear Approximation to Volatility

To proceed, we consider a linear approximation to (A.66) obtained by substituting

$$
\begin{equation*}
\sqrt{z_{i}} \rightarrow p_{i}+q_{i} z_{i} \tag{A.71}
\end{equation*}
$$
^eqa-71

with the coefficients to be found. Define by $c(t)=c\left(t ; p_{1}, q_{1}, p_{2}, q_{2}\right)$ the value of $\mathrm{E}\left(\widetilde{\tilde{z}}_{1}(t) \widetilde{z}_{2}(t)\right)$ in the approximate model,

$$
\begin{equation*}
d \widetilde{z}_{i}(t)=\theta_{i}\left(1-\widetilde{z}_{i}(t)\right) d t+\eta_{i}\left(p_{i}+q_{i} \widetilde{z}_{i}(t)\right) d W_{i}(t), \quad \widetilde{z}_{i}(0)=1 \tag{A.72}
\end{equation*}
$$
^eqa-72

Simple calculations yield an ODE on $c(t)$,

$$
\begin{aligned}
& \frac{d}{d t} c(t)=-\left(\theta_{1}+\theta_{2}\right)(c(t)-1) \\
& \quad+\rho \eta_{1} \eta_{2}\left(p_{1}+q_{1}\right)\left(p_{2}+q_{2}\right)+\rho \eta_{1} \eta_{2} q_{1} q_{2}(c(t)-1)
\end{aligned}
$$

which can be solved,

$$
\begin{equation*}
c(t)=1+\rho \eta_{1} \eta_{2}\left(p_{1}+q_{1}\right)\left(p_{2}+q_{2}\right) \frac{1-e^{-\left(\theta_{1}+\theta_{2}-\rho \eta_{1} \eta_{2} q_{1} q_{2}\right) t}}{\theta_{1}+\theta_{2}-\rho \eta_{1} \eta_{2} q_{1} q_{2}} . \tag{A.73}
\end{equation*}
$$
^eqa-73

Then, from (A.70) and the equality (A.73) that we use to approximate $\mathrm{E}\left(\left(z_{1}(t) z_{2}(t)\right)\right.$ in the original model (A.66), we get

$$
\begin{align*}
\mathrm{E} & \left(\sqrt{z_{1}(t) z_{2}(t)}\right)  \tag{A.74}\\
\approx & \frac{\left(p_{1}+q_{1}\right)\left(p_{2}+q_{2}\right)}{\theta_{1}+\theta_{2}-\rho \eta_{1} \eta_{2} q_{1} q_{2}} \\
& \times\left(\theta_{1}+\theta_{2}-\rho \eta_{1} \eta_{2} q_{1} q_{2} e^{-\left(\theta_{1}+\theta_{2}-\rho \eta_{1} \eta_{2} q_{1} q_{2}\right) t}\right) \\
= & \left(p_{1}+q_{1}\right)\left(p_{2}+q_{2}\right)\left(1+\rho \eta_{1} \eta_{2} q_{1} q_{2} \frac{1-e^{-\left(\theta_{1}+\theta_{2}-\rho \eta_{1} \eta_{2} q_{1} q_{2}\right) t}}{\theta_{1}+\theta_{2}-\rho \eta_{1} \eta_{2} q_{1} q_{2}}\right)
\end{align*}
$$
^eqa-74

## A.A.1.3 Step 3. Coefficients of the Linear Approximation

Applying (A.74) to $z_{2}(t)=z_{1}(t)$ we obtain

$$
\begin{equation*}
1=\mathrm{E}\left(\left(z_{1}(t)\right)=\frac{\left(p_{1}+q_{1}\right)^{2}}{2 \theta_{1}-\eta_{1}^{2} q_{1}^{2}}\left(2 \theta_{1}-\eta_{1}^{2} q_{1}^{2} e^{-\left(20_{1}-\eta_{1}^{2} q_{1}^{2}\right) t}\right)\right. \tag{A.75}
\end{equation*}
$$
^eqa-75

which gives one equation on $p_{1}, q_{1}$. The other one is obtained by using as $z_{2}(t)$ an independent copy of $z_{1}(t)$,

$$
\begin{equation*}
\left(\mathrm{E}\left(\sqrt{z_{1}(t)}\right)\right)^{2}=\left(p_{1}+q_{1}\right)^{2} \tag{A.76}
\end{equation*}
$$
^eqa-76

The coefficients $p_{1}$ and $q_{1}$ are determined as a solution to the system (A.75), (A.76), provided $\mathrm{E}\left(\sqrt{z_{1}(t)}\right)$ is known. Similar system holds for $p_{2}, q_{2}$. This completes the derivation of the approximating function (A.67).

## A.A.1.4 Step 4. Order of Approximation

The first two features of the approximation listed in the Proposition A.A. 1 are valid by construction. The validity of the last feature is obvious because the usage of a linear approximation (A.71) is consistent with the first order of perturbative expansion in volatilities. This is also easy to verify directly. Indeed, a straightforward perturbative calculation gives ${ }^{5}$

$$
\begin{gathered}
\mathrm{E}\left(\sqrt{z_{i}(t)}\right)=1-\frac{\eta_{i}^{2}\left(1-e^{-2 \theta_{i} t}\right)}{16 \theta_{i}}+O\left(\eta_{i}^{3}\right), \\
\frac{\mathrm{E}\left(\sqrt{z_{1}(t) z_{2}(t)}\right)}{\mathrm{E}\left(\sqrt{z_{1}(t)}\right) \mathrm{E}\left(\sqrt{z_{2}(t)}\right)}=1+\frac{\rho \eta_{1} \eta_{2}\left(1-e^{-\left(\theta_{1}+\theta_{2}\right) t}\right)}{4\left(\theta_{1}+\theta_{2}\right)}+O\left(\left(\max \left(\eta_{1}, \eta_{2}\right)\right)^{3}\right),
\end{gathered}
$$

which is in agreement with the leading order of expansion of (A.67) in $\eta_{1}$ and $\eta_{2}$. Note that we used the zero order expansion for coefficients $q_{\imath}=1 / 2+O\left(\eta_{i}\right)$ to prove the agreement.

## A.A. 2 Proof of Lemma A.A. 2

By Proposition 8.3.2, $z(t)$ is distributed as $e^{-\theta t} / n(0, t)$ times a non-central chi-square distributed random variable $\xi$ with $\nu=d$ degrees of freedom and non-centrality parameter $\gamma=n(0, t)$. Thus we obtain that

$$
\mathrm{E}(\sqrt{z(t)})=\left(\frac{e^{-\theta t}}{n(0, t)}\right)^{1 / 2} \mathrm{E}(\sqrt{\xi})
$$

where $\xi$ has the density (see (8.5))

$$
\frac{\partial}{\partial z} \Upsilon\left(z ; \nu, \gamma^{\prime}\right)=e^{-\gamma / 2} \sum_{j=0}^{\infty} \frac{(\gamma / 2)^{j}}{j!2^{\nu / 2+j} \Gamma(\nu / 2+j)} z^{\nu / 2+j-1} e^{-z / 2}
$$

[^119]Then

$$
\begin{align*}
\mathrm{E}(\sqrt{\xi}) & =e^{-\gamma / 2} \sum_{j=0}^{\infty} \frac{(\gamma / 2)^{j}}{j!2^{\nu / 2+j} \Gamma(\nu / 2+j)} \int_{0}^{\infty} z^{1 / 2+\nu / 2+j-1} e^{-z / 2} d z \\
& =e^{-\gamma / 2} \sum_{j=0}^{\infty} \frac{(\gamma / 2)^{j}}{j!2^{\nu / 2+j} \Gamma(\nu / 2+j)} 2^{1 / 2+\nu / 2+j} \Gamma(\nu / 2+j+1 / 2)  \tag{A.77}\\
& =\sqrt{2} e^{-\gamma / 2} \sum_{j=0}^{\infty} \frac{(\gamma / 2)^{j}}{j!} \frac{\Gamma(\nu / 2+j+1 / 2)}{\Gamma(\nu / 2+j)}
\end{align*}
$$
^eqa-77


[^0]:    ${ }^{3}$ Vega is the volatility sensitivity, see Remark 8.9.3.

[^1]:    ${ }^{4}$ It can be shown that the sticky strike rule is arbitrageable.

[^2]:    ${ }^{5}$ These can be computed by differentiating swaption prices twice with respect to the strike, see Section 7.1.2.

[^3]:    ${ }^{6}$ In some applications of TSR models, it may suffice that this expression holds only for a single value of $M$ (namely the payment date of the security in question). In such cases, certain simplifications may be possible, as demonstrated in Section 16.6.4.

[^4]:    ${ }^{7}$ Of course establishing caplet prices may itself require some work, as only prices of full caps are quoted. See Section 16.2 for more on this.

[^5]:    ${ }^{8}$ Recall that $T_{1}$ is the first payment date of the swap underlying the swap rate $S$.
    ${ }^{9}$ It is nevertheless often useful to be able to calculate CMS convexity adjustments in a given term structure model in closed form, e.g. for assessing the loss of precision of the model when pricing exotics linked to CMS rates. We return to this task later in the chapter.

[^6]:    ${ }^{10}$ In theory; numerical implementation could lead to slight differences.

[^7]:    ${ }^{11}$ For $T^{\prime} \neq T_{p}^{\prime} \leq T_{1}$ the function $\rho$ to be defined momentarily would be multiplied by a term that expresses $P\left(T^{\prime}, T_{p}^{\prime}\right)$ as an (approximate) function of $x_{0}=L_{0}(T)$ as in, e.g., Section 16.5.1.

[^8]:    ${ }^{12}$ Here we essentially assume that the slope of the yield curve is independent of the "pure" FX component of the forward FX rate.

[^9]:    ${ }^{13}$ For very long-dated quanto contracts, the FX volatility smile may start to matter. The method we develop here can be extended to incorporate the FX smile, by techniques discussed in Chapter 17.
    ${ }^{14}$ This is not exact as the $T$-expiry FX option is written on $X_{T}(T)$, not $X_{T_{p}}(T)$. The difference is rarely material as $T_{p}-T$ is often small, but more elaborate schemes are not difficult to design, if desired.

[^10]:    ${ }^{15}$ This is similar to the replication method for computing CMS and Libor-inarrears convexity adjustments.

[^11]:    ${ }^{16}$ We use $\sigma$ instead of our customary $\lambda$ to avoid notational conflict with LM model volatilities used later on.

[^12]:    ${ }^{17}$ We note in passing that this technique can also be used to control errors in volatility expansion formulas, such as the one used in SABR, which may yield negative densities in the wings unless some kind of regularization is performed.

[^13]:    ${ }^{1}$ I.e. $\mathrm{E}\left(Z_{2}\right)=0, \operatorname{Var}\left(Z_{2}\right)=1, i=1, \ldots, d$.
    ${ }^{2}$ In this chapter we will often consider undiscounted values of derivatives, as we mostly work with cash flows that pay at a single payment time $T_{p}$. As it should always be clear whether discounting is applied or not, we do not introduce new notation.

[^14]:    ${ }^{3}$ If the marginal distribution functions $\Psi_{1}, \Psi_{2}, \ldots, \Psi_{d}$ are all continuous, then the copula function is unique.

[^15]:    ${ }^{4}$ To demonstrate the dependence of implied correlation on the marginal distribution of swap rates, Appendix 17.A uses a setting with displaced diffusion processes to quantify the effect on ATM spread option prices due to changes in the swap rate volatility skew.

[^16]:    ${ }^{5}$ This is most easily seen by assuming that the two rates are in fact the same. In this case the correlation between $S_{1}(T)$ and $S_{2}(T)$ is obviously 1 , whereas the correlation between $S_{1}(T)$ and $S_{2}(t)$ would be less than 1 as the increment $S_{2}(t)-S_{2}(T)$ would typically be only weakly dependent on $S_{1}(T)$.

[^17]:    ${ }^{6}$ Note we drop the superscript $T_{p}$ for notational convenience for the duration of this section.

[^18]:    ${ }^{7}$ A trapezoidal scheme is a possible candidate. Most Gaussian quadrature rules, however, are not ideal, since their integration grids are not directly user-specifiable, but emerge as roots of a particular function.
    ${ }^{8}$ Sometimes also called Gauss-Legendre quadrature after the name of the family of polynomials whose roots define the integration grid.

[^19]:    ${ }^{9}$ Note that in Section 17.6.3 we already derived a one-dimensional integral representation for this payoff, but we use it anyway to demonstrate the main idea of the method we develop here.
    ${ }^{10}$ Not annualized, i.e. not divided by $\sqrt{T}$.

[^20]:    ${ }^{11}$ Unless the marginal distributions of the rates are Gaussian, of course.
    ${ }^{12}$ At the time of writing, the Gaussian copula correlations implied from at-themoney spread options were largely independent of the expiry, an observation that is inconsistent with predictions of almost all multi-factor term structure models. Some observers attribute this market feature to the proliferation of copulas with time-independent parameters.

[^21]:    ${ }^{13}$ If the QE scheme is used, whenever a uniform random variable $U$ is required, we would write $U=\Phi(Z)$ with $Z$ being a standard Gaussian random variable. This way, the QE scheme can be driven solely by Gaussian random variables.

[^22]:    ${ }^{14}$ This is a consequence of using just a single stochastic variance scaling applied to all Libor rates.

[^23]:    ${ }^{15}$ A copula is easily shown to be invariant with respect to monotonic transformations of the underlying variables.

[^24]:    ${ }^{16}$ We can extend the definition to handle non-ATM strikes, but this shall not be needed for the purposes of this appendix.

[^25]:    ${ }^{1}$ A practice soinetimes known as mark-to-model.

[^26]:    ${ }^{2}$ As expressed h..say, Black volatilities.

[^27]:    ${ }^{3}$ We discuss the details of the implementation of the regression algorithm later, in Section 18.3.10.

[^28]:    ${ }^{4}$ In this basic scheme we implicitly assume that underlying value $U_{n}$ at time $T_{n}^{\prime}$ can be calculated in closed form from the simulated yield curve at time $T_{n}^{\prime}$ (or, more generally, from the model state variables observed on and before time $T_{n}^{\prime}$ ). This is a strong restriction which we relax later.

[^29]:    ${ }^{5} \mathrm{~A}$ naive numerical implementation of regression will not have this property, see Section 18.3.10.

[^30]:    ${ }^{6}$ The value of a coupon must obviously be known at the time it is paid.

[^31]:    ${ }^{7}$ A similar issue arises with the regressed hold values, as discussed in Section 18.3.4.

[^32]:    ${ }^{8}$ Jumping a bit ahead, we note that it is essential to use this optimization in conjunction with an independent post-simulation, see Section 18.3.6.

[^33]:    ${ }^{9}$ For more detailed discussion of nested Monte Carlo simulation, see Section 18.3.8.
    ${ }^{10}$ Which may include refinements such as that in Section 18.3.5, of course.

[^34]:    ${ }^{11}$ In addition to random Monte Carlo error, simulation of some models may also involve a systematic error stemming from the time-discretization of the model dynamics. Such discretization errors are not accounted for in (18.41), nor in any previous argument. In most cases, however, the discretization bias will be negligible relative to the random Monte Carlo error.

[^35]:    ${ }^{12}$ Note that in testing the viability of a class of exercise rules through an upper bound simulation, it is often acceptable to work with a reduced set of exercise opportunities - e.g. change a quarterly exercise schedule to an annual one, say in order to save computation time (see (18.40)).

[^36]:    ${ }^{13}$ This limit is essentially equivalent to imposing a so-called "carry" restriction on exercise, a notion that we explore in some detail in Section 18.3.10.2 and Proposition 19.7.1.

[^37]:    ${ }^{14}$ While for convenience we write $\zeta$ as a function of time $t$, the regression variables only need to be defined at times when we perform regression, i.e. at exercise dates $T_{1}, \ldots, T_{N-1}$.
    ${ }^{15}$ Beyond representing information available at time $T, \zeta(T)$ also serves to define the function space that is obtainable by least-squares regression.

[^38]:    ${ }^{16}$ Glasserman and Yu [2004] determine that the number of paths should grow exponentially in the number of regression variables. Under different assumptions, Moni [2005] shows that the number of paths should grow polynomially. Polynomial rate of growth is also derived in Egloff et al. [2007].
    ${ }^{17}$ Section 18.3.10 lists a number of regularization approaches that may help guard against some of the issues associated with spurious regression variables.

[^39]:    ${ }^{18}$ As discussed in Section 18.3.10, a simple pre-normalization of variables before the regression may still be useful in preventing numerical problems.

[^40]:    ${ }^{19}$ The interested reader could start with Campos et al. [2005].
    ${ }^{20} \mathrm{~A}$ potentially more accurate measure is "modified R -squared", see Campos et al. [2005]. For our applications where the dimension of the dataset ( $K$ ) is very large, modified R -squared and ordinary R -squared are almost identical, however.
    ${ }^{21}$ See e.g. www.pcgive.com/pcgets.

[^41]:    ${ }^{22} \mathrm{~A}$ more general result is derived in Section 19.7.2.

[^42]:    ${ }^{23}$ On a related note, Rasmussen [2005] suggests starting simulation in the past to ensure sufficient variability in the state space for small times, and shows that this step improves estimation of the exercise boundary.

[^43]:    ${ }^{24}$ Although we can always try to increase the range of applicability of the qG model with some of the techniques of Chapter 21.

[^44]:    ${ }^{25}$ For more detail on this topic, see Chapter 26.

[^45]:    ${ }^{1}$ But see Sections 19.4.7 and 19.4.8 below for exceptions.

[^46]:    ${ }^{2}$ At the time of writing this is done by Markit, see www.markit.com.

[^47]:    ${ }^{3}$ For more information on which volatility to calibrate to, see Section 19.3 below.

[^48]:    ${ }^{44}$ The optimal exercise strategy for a Bermudan swaption must be as least as good as simply picking at time 0 one of the exercise dates and never changing one's mind.

[^49]:    ${ }^{5}$ Although this is not true for the family of Markov-functional models, see Appendix 11.A in Chapter 11; in fact the difficulty of handling amortizing/accreting Bermudan swaptions is often cited as one of the problems with such models.

[^50]:    ${ }^{6}$ It is probably advisable to impose smoothness constraints on the timedependent mean reversion function while performing such calibration.

[^51]:    ${ }^{7}$ Term correlations between swap rates in one-factor models are not exactly $100 \%$ due to time-dependence in parameters and presence of mean reversion.
    ${ }^{8}$ Appendix 19.B gives a quick tour of the classical moment matching ideas.

[^52]:    ${ }^{9} \mathrm{Or}$, equivalently, pays 1 at time $T_{N}$ in addition to running Libor coupons on a unit notional.

[^53]:    ${ }^{10}$ We note that exercise in the middle of a fixed-rate period is most often accompanied by an exercise fee, a deterministic amount of money payable upon exercise that is agreed upon in advance. The fee is typically calculated to reflect the value of the part of the fixed-rate coupon accrued from the beginning of the period to the exercise date.

[^54]:    ${ }^{11}$ Recall from Section 18.3.8.2 that $\iota\left(T_{n}\right)=1_{\left\{U_{n}\left(T_{n}\right)>H_{n}\left(T_{n}\right)\right\}}$ defines the rule for exercising at $T_{n}$, assuming we have not exercised previously.

[^55]:    ${ }^{12}$ An alternative, and even less expensive, technique to apply policy improvement would rely on the carry argument developed in Section 19.7.2. The resulting bound requires no option price computations.
    ${ }^{13}$ For a discussion of exercise strategies for best-of options (also known as MAX-options), see Broadie and Detemple [1997].
    ${ }^{14}$ Indeed, we notice that Strategy I is, in fact, optimal for a 1 -factor Markov short rate model. For a 1 -factor LM model, however, Strategy I is not optimal (although, as we shall see later, it appears to perform very well).

[^56]:    ${ }^{15}$ Of course, when applying the portfolio replication method, we should explicitly link mean reversion to market inputs through the local projection method of Sections 19.2 and 18.4.

[^57]:    ${ }^{16}$ See Figures 1-3 in Andersen and Andreasen [2001] for visual confirmation of this as well as of (19.37).

[^58]:    ${ }^{1}$ There is some evidence that many participants in the volatility swap market tend to use fairly naive, low-dimensional models for valuation. As a result, if correlations for the LM model are extracted from spread options, say, the LM model may produce forward volatilities that are lower than the market consensus (see Appendix 19.A for the rationale). Any arbitrages induced by such "segmentation" between the markets for spread options and volatility swaps are hard to exploit in practice, so the market differences can be quite persistent.

[^59]:    ${ }^{2}$ Sometimes the coupon is linked to a swap rate that starts at shout time $\eta_{n}$ rather than at the end of the period $T_{n+1}$. This modification to our discussion is easy to incorporate, and we do not consider this case separately.
    ${ }^{3}$ In the $T_{n+1}$-forward measure, $S_{n+1}(t)$ will typically have a small convexityinduced drift, as $S_{n+1}(t)$ here represents some CMS rate. However, the period [ $T_{n}, T_{n+1}$ ] rarely exceeds one year, and over one year the drift of a CMS rate will normally be quite close to zero.

[^60]:    ${ }^{4}$ Here $s$ also stands for "strike setting" and $e$ for "expiry".

[^61]:    ${ }^{5}$ In the sense that there is no sensitivity to the yield curve, provided that all basis-point (Gaussian) volatilities are kept fixed. If we assume that the basis-point

[^62]:    ${ }^{1}$ For exotic swaps that require a model for their valuation, the "market" is often understood in a broad sense to represent the values of underlying exotic swaps in a vanilla model of choice.

[^63]:    ${ }^{2}$ Needless to say, this is not something that we generally recommend.

[^64]:    ${ }^{3}$ We briefly consider least-squares norms later in the section.

[^65]:    ${ }^{4}$ While we only consider one coupon and one swap rate, it is trivial to extend our discussion to the standard case of multiple coupons depending on different rates.

[^66]:    ${ }^{1}$ In addition, derivatives portfolios are exposed to credit risk, i.e. the risk that the counterparty to a derivatives transaction will default on its obligations. Management of credit risk is outside the scope of this book.

[^67]:    ${ }^{2}$ We use the term vector loosely, since some elements of $\Theta_{\text {prm }}$ (e.g. the discount curve) may be continuous functions.

[^68]:    ${ }^{3}$ Recall from Chapter 6 that some yield curve construction algorithms require control parameters (such as tension parameters and precision tolerances), so $\Theta_{\text {prm }}(t)$ may be required for the construction of the discount curve from market inputs.

[^69]:    ${ }^{4}$ For simple securities in simple models, it is possible to derive certain analytical results concerning the costs and risks of actual (discrete-time, costly) hedging strategies. A classical paper in this area is Leland [1985], although the approach has recently received some criticism (see, e.g., Kabanov and Safarian [1997]). Other relevant papers include, among many others, Soner et al. [1995], Barles and Soner [1998], and Derman and Kamal [1999]. As all derivatives traders manage risk at the portfolio (or book) level, the "net" security owned by traders is far more complicated than those covered by most of the literature. In addition, traders tend to rebalance their trading books according to rules much more complex than those assumed in academic papers. As a result, proper charging for transactions normally requires a heavy element of human judgment, and will depend strongly on the portfolio context.

[^70]:    ${ }^{5}$ If, say, a single vega is nevertheless required, it is often most reasonable to report the sensitivity to a parallel shift of the function $\sigma(t)$ at all values of $t$. We can think of this as roughly representing the sensitivity with respect to the first principal component of volatility curve moves. A similar principle can be applied to compute a single rho.

[^71]:    ${ }^{6}$ Sometimes the risk system is asked to perform the Taylor expansion around a series of different (perturbed) market data scenarios, not just the current market data. The resulting collections of reports are known as ladders.
    ${ }^{7}$ Most derivatives traders have views on the future market evolution and are allowed to express their views in proprietary ("prop") positions, meant to make money if the trader's views turn out to be correct. In this situation, the trader will purposely leave the portfolio open to certain market risk exposures. Banks typically enforce strict limits on the size of these exposures.

[^72]:    ${ }^{8}$ In a pinch, it is often reasonable to simply assume that $x^{*}(t)=T^{-1}(X(0)(T- t)+K t$ ).

[^73]:    ${ }^{9}$ Also known by the more grammatically sensible names $P \& L$ explanation or $P \& L$ attribution.

[^74]:    ${ }^{10}$ Notice that $\Lambda_{\alpha}$ virtually always is a negative number. Sometimes it is the absolute value of this number that is reported as the VaR.
    ${ }^{11}$ Specifically, VaR is not a coherent risk measure, in the sense defined in Artzner et al. [1999].

[^75]:    ${ }^{1}$ But not always - see Chapter 24.

[^76]:    ${ }^{2}$ The observant reader may have noticed a strong connection to material in Chapter 3 on pathwise and likelihood ratio methods for Monte Carlo applications.

[^77]:    ${ }^{3}$ Sometimes one gets better numerical properties if a spline is fit through all grid values at $t=0$; deltas and gammas can then be computed by differentiating this spline.

[^78]:    ${ }^{4}$ In fact, several commonly used risk measures are explicitly meant to be computed by finite-sized shifts. For instance, many swaption traders' definition of gamma is the change in delta, for a 10 basis point move in the yield curve.
    ${ }^{5}$ Here and elsewhere in this chapter, the term "singularity" refers to a point in the state space at which the payoff function, or one of its derivatives, is discontinuous. Examples include the barrier for a digital option and the strike of a European put or call.

[^79]:    ${ }^{6}$ If there is more than one singularity introduced in the event, the decomposition above holds locally around each singularity. We consider a single singularity case only, with trivial extension to multiple ones.

[^80]:    ${ }^{7}$ In other words, we assume that the discretization is fine enough to guarantee that there is only one singularity per interval.

[^81]:    ${ }^{8}$ In line with the footnote 7 we assume that the rectangles are small enough so that we can assume there is only two crossing points.

[^82]:    ${ }^{9} X_{\imath}(\omega)$ can, of course, be discontinuous, but this is not our focus at the moment.

[^83]:    ${ }^{10}$ The concept of "partial" membership in a set should be familiar to those schooled in fuzzy logic (see Zadeh [1965]), and tube Monte Carlo can, in fact, be considered a probabilistically motivated fuzzy logic algorithm. For more discussion of fuzzy logic applications to Monte Carlo sensitivity computations in finance, see Withington and Lucic [2009].

[^84]:    ${ }^{1}$ Of course, heeding advice from Chapter 23, we should use the same seed and the same number of paths in the base and shocked scenarios.

[^85]:    ${ }^{2}$ Naturally, the sensitivities to these rates can be projected into any other "basis", i.e. a set of rates used to define and aggregate curve sensitivities. For the relevant techniques, see Section 6.4.3.

[^86]:    ${ }^{3}$ Note that this measure density is not a probability density, as it does not integrate to 1 .

[^87]:    ${ }^{4}$ As we mentioned before in Section 6.4, once the deltas in a particular "basis" are computed, it is a matter of simple linear algebra to re-express them in any other basis, e.g. the one used by a risk management system.

[^88]:    ${ }^{5}$ We could also use the faster approximation (24.17) here; we leave relevant modifications for the reader to explore.

[^89]:    ${ }^{6}$ Also, as calibration and time-discretization is normally set up in a productspecific manner, it can be awkward (and even sub-optimal) to attempt to price many securities simultaneously in a single Monte Carlo loop.

[^90]:    ${ }^{7}$ Note that we are here comparing against a brute-force perturbation method where the exercise boundary - rather than the exercise time - is kept fixed under perturbations. Had we instead kept the exercise time fixed in perturbed scenarios, we would likely have obtained greeks of quality comparable to those produced by the pathwise method. Recall our comments at the end of Section 24.3.1.
    ${ }^{8}$ Another potential drawback is the need to know the transition density of the underlying process, although one always has the option of using a Gaussian approximation based on an Euler discretization of the true process. As discussed in Chen and Glasserman [2007b], the limit of this procedure for small time steps is deeply connected to the Malliavin calculus.

[^91]:    ${ }^{1}$ Sometimes a TARN is structured so that the first digital is virtually worthless, but the second one is important. The discussion that follows should then be modified accordingly.

[^92]:    ${ }^{2}$ We comment on the multi-factor case later.

[^93]:    ${ }^{3}$ See Section 3.4.1 for a discussion of efficiency measures for variance reduction schemes.
    ${ }^{4}$ Such as the quasi-Gaussian (qG) model, which normally involves simulation of many fewer state variables than the LM model, see Section 13.1.9.3.

[^94]:    ${ }^{5}$ Approximations to the drift could be improved by using the Brownian bridge techniques described in Section 14.6.2.5, but the impact of these improvements on the intended control variate applications is negligible.

[^95]:    ${ }^{6}$ We generally do not recommend this interpolation scheme, but it makes for a good example.

[^96]:    ${ }^{7}$ We can use the triangular form here because the square root of a variancecovariance matrix can always be written this way by application of the Cholesky decomposition, see Section 3.1.2.1. A more general form, however, could be beneficial for fitting as discussed later, see footnote 9 .

[^97]:    ${ }^{8}$ Note that, in line with (14.42), we assume that $\lambda$ 's and $\alpha$ 's are piecewise constant over $\left[T_{\imath}, T_{\imath+1}\right), i=0, \ldots, N-2$.

[^98]:    ${ }^{9}$ Denson and Joshi [2009] also extend the specification (25.47) by allowing an extra term for $\widehat{\lambda}_{n}^{1}(t)$ with $\widehat{\lambda}_{n}^{1}(t)=\sigma_{n}^{1} \alpha^{11}(t)+\sigma_{n}^{2} \alpha^{12}(t)$, which may improve the fit somewhat.

[^99]:    ${ }^{10}$ This is a good idea for any Markovian model, see Section 18.3.9.1.

[^100]:    ${ }^{11}$ The swap market models of Section 15.4 are, however, not affected by this issue.

[^101]:    ${ }^{12}$ While we use the fixed rate of the swap as a strike for a cap, a potentially better control variate could be constructed by using a strike at or near the exercise boundary of the Bermudan swaption.

[^102]:    ${ }^{13}$ Multiplicative duality is developed in Jamshidian [1995] and is closely related to the (additive) duality of Section 1.10.2. A comparison of multiplicative and additive duality for upper bound simulations can be found in Chen and Glasserman [2007a].

[^103]:    ${ }^{1}$ In multi-factor models, for each $t$ and $T^{\prime}, \sigma_{f}(t, T)$ is obviously a multidimensional vector, but it is not this dimensionality that we are interested in here.

[^104]:    ${ }^{2}$ We are ignoring the strike dimension for now.

[^105]:    ${ }^{3}$ See footnote 16 in Chapter 6 for a similar definition of sensitivity to a shock to a yield curve.

[^106]:    ${ }^{4}$ Between the two extremes of the flat and bucketed shocks lie row shocks $\delta_{n}$ of the form $\left(\delta_{n}\right)_{i, j}=1_{\{i=n\}}, n=1, \ldots, N_{t}$.

[^107]:    ${ }^{5}$ These are typically in the order of $0.1 \%$ in implied volatility terms (with typical market swaption volatilities being in the $10-50 \%$ range).

[^108]:    ${ }^{6}$ If the benchmark set does not include all coterminal swaptions for a given Bermudan swaption, non-zero vegas are still possible immediately below the coterminal diagonal due to interpolation effects.

[^109]:    ${ }^{7}$ See Section 14.4.2 for examples of such formulas. Notice that the Jacobian $\left\{\partial \Lambda_{i, j} / \partial \delta_{n, m}\right\}$ is often available for free as part of the initial calibration of the LM model, especially if calibration relies on a gradient-based optimization method.

[^110]:    ${ }^{8}$ This could be done arbitrarily, but for concreteness we do it in row-major order, i.e. rows of the matrix are stacked end-to-end to come up with a vector.
    ${ }^{9}$ Of course, the pathwise or likelihood differentiation methods of Chapter 24 could have been used here as well.

[^111]:    ${ }^{10} \mathrm{An}$ observation also made by Pietersz and Pelsser [2004], although in a somewhat different context.

[^112]:    ${ }^{11}$ We comment below on the choice of $\alpha_{n, m}$ as well as other required constants.

[^113]:    ${ }^{12}$ But not sensitivities to term correlations of swap rates, which would lead to another inconsistency, with double-counting of risk.

[^114]:    ${ }^{13}$ Note that it is impractical to include all swap rate correlations in the calibration set. Instead. one would typically choose a set (or perhaps a few sets) of correlations of two specific swap rate tenors, such as 10 year and 2 year, over a collection of time periods.

[^115]:    ${ }^{14}$ This situation would often also be reflected in the usage of different swaption weights in the calibration norm (14.51), with precision weights on illiquid swaptions set lower than on liquid ones.

[^116]:    ${ }^{1}$ A version of which we have already seen in Proposition 7.4.2.

[^117]:    ${ }^{2}$ Although some creative approximations for such securities can occasionally be derived from Markovian projection, see Section 13.1.9.4.

[^118]:    ${ }^{3}$ Markovian projection on processes with quadratic local volatility (see Section 7.3) was developed in Antonov and Misirpashaev [2009b], using Wiener chaos expansion techniques. Predictably it outperforms projections on linear local volatility processes but is significantly more complicated.
    ${ }^{4}$ Dupire [1997] dubs this corollary the universal law of volatility.

[^119]:    ${ }^{5}$ Due to the reflecting symmetry of the Brownian motion, an arbitrary average of the variables $z_{1}(t)$ and $z_{2}(t)$ is an even function of $\eta_{1}$ and $\eta_{2}$. Thus, the order of the approximations below is effectively higher.

