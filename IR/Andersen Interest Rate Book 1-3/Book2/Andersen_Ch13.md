---
tags: [quant/knowledge, quant/books, interest-rates]
---

## The Quasi-Gaussian Model with Local and Stochastic Volatility

In this chapter we consider extensions to one- and multi-factor Gaussian short rate models (Chapters 10, 11 and 12) with local and stochastic volatility. The extensions come at additional computational cost, as extra state variables are required to preserve the Markovian structure of the model. Following the pioneering work of Jamshidian [1991b], we use the term quasiGaussian ${ }^{1}$ for the models in this chapter; their development for practical applications was undertaken in Andreasen [2001], Andersen and Andreasen [2002] and Andreasen [2005], building on early work by Jamshidian [1991b], Babbs [1990], Cheyette [1991] and Ritchken and Sankarasubramanian [1995]. Low-dimensional versions of quasi-Gaussian models are, in our opinion, among the best - if not the best - low-factor short rate models, as they combine flexibility of volatility smile specification, relative ease of calibration, and efficient numerical implementation. Higher-dimensional quasi-Gaussian models, while not yet mainstream, provide an alluring alternative to the better-established Libor market models (see Chapter 14).

We start this chapter by developing a one-factor quasi-Gaussian model with a local volatility function. The problems of volatility and mean reversion calibration are given considerable attention, and are followed by a discussion of various numerical methods used for model implementation. A straightforward extension to stochastic volatility is presented next, followed by development of multi-factor quasi-Gaussian models.

### 13.1 One-Factor Quasi-Gaussian Model

### 13.1.1 Definition

Recall that any HJM model is defined by a volatility structure of instantaneous forward rates. In particular, for any "reasonable" random function $\sigma_{f}(t, T)=\sigma_{f}(t, T, \omega)$, the following SDE defines a valid HJM model,[^39]

$$
\begin{equation*}
d f(t, T)=\sigma_{f}(t, T)\left(\left(\int_{t}^{T} \sigma_{f}(t, u) d u\right) d t+d W(t)\right), \quad 0 \leq t \leq T<\infty \tag{13.1}
\end{equation*}
$$
^eq13-1

Here $W(t)$ is a one-dimensional Brownian motion in the risk-neutral measure, and $\{f(t, T)\}_{T \geq t}$ is a collection of instantaneous forward rates.

It is shown in Section 4.5.2 that a one-factor Markovian Gaussian model is obtained by imposing a separability condition on the deterministic volatility structure of instantaneous forward rates, see (4.44). A general class of onefactor quasi-Gaussian ( qG$)^{2}$ models is obtained by retaining the separability condition, but relaxing the deterministic requirement in a specific way. In particular, the component of the volatility structure that is a function of calendar time (the function $g$ ), is now allowed to be stochastic:

$$
\begin{equation*}
\sigma_{f}(t, T, \omega)=g(t, \omega) h(T) \tag{13.2}
\end{equation*}
$$
^eq13-2

In line with the notations of Section 10.1.2.2, we define

$$
\begin{align*}
\varkappa(t) & =-\frac{h^{\prime}(t)}{h(t)}  \tag{13.3}\\
G(t, T) & =\frac{\int_{t}^{T} h(s) d s}{h(t)} \\
\sigma_{r}(t, \omega) & =\sigma_{f}(t, t, \omega)=g(t, \omega) h(t)
\end{align*}
$$
^eq13-3

The proof of Proposition 10.1.7 carries through unchanged even for stochastic $g(t, \omega)$, and we obtain the following result.
Proposition 13.1.1. Consider a general one-factor $q G$ model, i.e. the HJM model (13.1) with the separable volatility condition (13.2). Define stochastic processes $x(t)$ and $y(t)$ by

$$
\begin{align*}
d x(t) & =(y(t)-\varkappa(t) x(t)) d t+\sigma_{r}(t, \omega) d W(t)  \tag{13.4}\\
d y(t) & =\left(\sigma_{r}(t, \omega)^{2}-2 \varkappa(t) y(t)\right) d t \\
x(0) & =y(0)=0
\end{align*}
$$
^eq13-4

In the general $q G$ model all zero-coupon discount bonds are deterministic functions of the processes $x(t)$ and $y(t)$,

$$
P(t, T)=P(t, T, x(t), y(t))
$$

where

$$
\begin{equation*}
P(t, T, x, y)=\frac{P(0, T)}{P(0, t)} \exp \left(-G(t, T) x-\frac{1}{2} G(t, T)^{2} y\right) \tag{13.5}
\end{equation*}
$$
^eq13-5

[^40]the instantaneous forward rates are given by

$$
\begin{equation*}
f(t, T)=f(0, T)+\frac{h(T)}{h(t)}(x(t)+G(t, T) y(t)) \tag{13.6}
\end{equation*}
$$
^eq13-6
and the short rate is

$$
\begin{equation*}
r(t)=f(t, t)=f(0, t)+x(t) \tag{13.7}
\end{equation*}
$$
^eq13-7

The proposition demonstrates that the evolution of the whole interest rate curve, as parameterized by either forward rates or discount bonds, in the model can be reduced to the evolution of just two state variables $x(t)$ and $y(t)$, with dynamics given by (13.4). Unlike many of the models in Chapter 11, the qG model has a closed-form bond reconstitution formula for arbitrary choices of $g(t, \omega)$; this tractability comes at the cost at requiring two state variables ( $x$ and $y$ ), even though the Brownian motion $W(t)$ is only one-dimensional. Observe that in general, the function $y(t)$ is not deterministic, except in the case of pure Gaussian dynamics, i.e. when $\sigma_{r}(t, \omega)$ is a deterministic function of $t$. However, even when it is not deterministic, $y(t)$ does not have the diffusion term; we call such processes locally deterministic.

The roles of the two state variables $x(t)$ and $y(t)$ in the qG model are rather different. The variable $x(t)$ constitutes the main yield curve driver, as evidenced in (13.7), whereas $y(t)$ is an auxiliary "convexity" variable required to uphold the no-arbitrage condition; in general, it is convenient to think of the model as having "one and a half" factors.

### 13.1.2 Local Volatility

A one-factor qG model with local volatility is obtained by requiring $g(\cdot)$ to be a deterministic, time-dependent function of the state variables,

$$
\begin{equation*}
g(t)=g(t, x(t), y(t)) \tag{13.8}
\end{equation*}
$$
^eq13-8

Then, the short rate volatility $\sigma_{r}(\cdot)$ is also a function of the state variables,

$$
\begin{equation*}
\sigma_{r}(t)=\sigma_{r}(t, x(t), y(t))=g(t, x(t), y(t)) h(t), \tag{13.9}
\end{equation*}
$$
^eq13-9

and the dynamics of the state variables in the local volatility qG model are given by

$$
\begin{align*}
& d x(t)=(y(t)-\varkappa(t) x(t)) d t+\sigma_{r}(t, x(t), y(t)) d W(t)  \tag{13.10}\\
& d y(t)=\left(\sigma_{r}(t, x(t), y(t))^{2}-2 \varkappa(t) y(t)\right) d t \tag{13.11}
\end{align*}
$$
^eq13-10
^eq13-11

Clearly, (13.10)-(13.11) define a two-dimensional Markovian process. As all zero-coupon discount bonds are functions of these two state variables by Proposition 13.1.1, the local volatility qG model is Markovian in two state variables.

For future use, we denote

$$
\begin{equation*}
\sigma_{r}^{0}(t) \triangleq \sigma_{r}(t, 0,0) \tag{13.12}
\end{equation*}
$$
^eq13-12

If $\sigma_{r}(t, x, y)$ is independent of $x, y$, the model reduces to a purely Gaussian model with the deterministic short rate volatility $\sigma_{r}^{0}(t)$.

### 13.1.3 Swap Rate Dynamics

For the purpose of European swaption pricing in the qG model, we shall need to establish swap rate dynamics in an annuity measure. For this purpose, let us fix a tenor structure

$$
0<T_{0}<T_{1}<T_{2}<\ldots<T_{N}
$$

with

$$
\tau_{n}=T_{n+1}-T_{n} .
$$

Consider a forward swap rate $S(t)$ with the first fixing $T_{0}$ and the last payment $T_{N}$ (see Section 4.1.3), i.e.

$$
\begin{align*}
& S(t) \triangleq S_{0, N}(t)=\frac{P\left(t, T_{0}\right)-P\left(t, T_{N}\right)}{A(t)}  \tag{13.13}\\
& A(t) \triangleq A_{0, N}(t)=\sum_{n=0}^{N-1} \tau_{n} P\left(t, T_{n+1}\right) \tag{13.14}
\end{align*}
$$
^eq13-13
^eq13-14

It follows from Proposition 13.1.1 that all zero-coupon bonds $P(t, \cdot)$ are deterministic functions of $x(t)$ and $y(t)$, and hence so is the forward swap rate; accordingly we define

$$
\begin{equation*}
S(t, x, y) \triangleq \frac{P\left(t, T_{0}, x, y\right)-P\left(t, T_{N}, x, y\right)}{\sum_{n=0}^{N-1} \tau_{n} P\left(t, T_{n+1}, x, y\right)} \tag{13.15}
\end{equation*}
$$
^eq13-15

The following proposition, a simple extension of the results from Section 10.1.3.2, determines the dynamics of the swap rate under its corresponding annuity measure, i.e. the measure $\mathrm{Q}^{A}$ for which $A(t)$ is a numeraire.

Proposition 13.1.2. We have

$$
\begin{equation*}
d S(t)=\left(\frac{\partial S}{\partial x}(t, x(t), y(t))\right) \sigma_{r}(t, x(t), y(t)) d W^{A}(t) \tag{13.16}
\end{equation*}
$$
^eq13-16

where $W^{A}(t)$ is a Brownian motion in measure $\mathrm{Q}^{A}$. Here

$$
\begin{align*}
\frac{\partial S}{\partial x}(t, x, y)= & -\frac{1}{A(t, x, y)}\left(P\left(t, T_{0}, x, y\right) G\left(t, T_{0}\right)-P\left(t, T_{N}, x, y\right) G\left(t, T_{N}\right)\right) \\
& +\frac{S(t, x, y)}{A(t, x, y)} \sum_{n=0}^{N-1} \tau_{n} P\left(t, T_{n+1}, x, y\right) G\left(t, T_{n+1}\right) \tag{13.17}
\end{align*}
$$
^eq13-17

Proof. By definition

$$
S(t)=S(t, x(t), y(t))
$$

The statement of the proposition follows by applying Ito's lemma and dropping $d t$ terms from the expression, as justified by the fact that the swap rate is a martingale in the annuity measure, per Lemma 4.2.4.

### 13.1.4 Approximate Local Volatility Dynamics for Swap Rate

The SDE (13.16) shows that a swap rate follows a local volatility process where the volatility is a function of the short rate state $x$ and the auxiliary variable $y$. Since the model is essentially one-factor, it is reasonable to assume that there is a strong linkage between a swap rate and the state of the short rate. Hence, it seems plausible that the dynamics for a swap rate could be written with a diffusion term that is just a function of the swap rate itself. Such a simplification would be convenient in many applications, as methods from Chapter 7 could be called upon to solve the resulting SDE or to price options. The following proposition proven by the methods of Markovian projection (see Appendix A) makes matters precise.

Lemma 13.1.3. The values of all European options on the swap rate $S(t)$ in the model (13.16) are identical to values computed in a vanilla model with time-dependent local volatility function:

$$
\begin{equation*}
d S(t)=\varphi(t, S(t)) d W^{A}(t) \tag{13.18}
\end{equation*}
$$
^eq13-18

where

$$
\begin{equation*}
\varphi(t, s)^{2}=\mathrm{E}^{A}\left(\left.\left(\frac{\partial S}{\partial x}(t, x(t), y(t)) \sigma_{r}(t, x(t), y(t))\right)^{2} \right\rvert\, S(t)=s\right) \tag{13.19}
\end{equation*}
$$
^eq13-19

While evaluating conditional expectations such as (13.19) is often rather difficult, here we are aided considerably by the essentially one-dimensional structure of the problem. If we assume that $y(t)$ is well approximated by a deterministic function $\bar{y}(t)$ - an approximation we shall use repeatedly in this chapter - then $S(t)$ would just be a deterministic function of $x(t)$ and time $t$. The opposite would also be true, i.e. $x(t)$ would be a deterministic function of $S(t), x(t)=X(t, S(t))$. If this function were available, the evaluation of the conditional expected value in (13.19) would boil down to evaluating $\frac{\partial S}{\partial x}(t, x(t), y(t)) \sigma_{r}(t, x(t), y(t))$ at $x(t)=X(t, s)$ and $y(t)=\bar{y}(t)$,

$$
\varphi(t, s) \approx \frac{\partial S}{\partial x}(t, X(t, s), \bar{y}(t)) \sigma_{r}(t, X(t, s), \bar{y}(t))
$$

where $\partial S / \partial x$ is given by (13.17).
Before presenting various methods for approximating $\bar{y}$ and the function $X$, let us emphasize again that once the volatility function $\varphi(t, s)$ is determined, swaption values can be computed from (13.18) by methods developed for local volatility vanilla models in Chapter 7.

### 13.1.4.1 Simple Approximation

A simple approximation for the function $\varphi$ is obtained by slightly extending the idea from Section 10.1.3.2. Setting

$$
\bar{y}(t)=0,
$$

and applying a linear approximation

$$
\begin{align*}
S(t, x, 0) & \approx S(t, 0,0)+\frac{\partial S}{\partial x}(t, 0,0) x,  \tag{13.20}\\
\frac{\partial S}{\partial x}(t, x, 0) & \approx \frac{\partial S}{\partial x}(t, 0,0),
\end{align*}
$$
^eq13-20

we obtain

$$
\begin{equation*}
x \approx \frac{S(t, x, 0)-S(t, 0,0)}{\partial S(t, 0,0) / \partial x} \tag{13.21}
\end{equation*}
$$
^eq13-21

Hence we arrive at the approximation

$$
\begin{equation*}
\varphi(t, s) \approx \frac{\partial S}{\partial x}(t, 0,0) \sigma_{r}\left(t, \frac{s-S(t, 0,0)}{\partial S(t, 0,0) / \partial x}, 0\right) \tag{13.22}
\end{equation*}
$$
^eq13-22

### 13.1.4.2 Advanced Approximation

The approximation (13.22) is quite accurate provided that volatility is low or moderate. For an approximation with a greater range of applicability, we can consider improving (13.20) by using a higher-order expansion and by taking greater care in selecting the expansion point. Starting with the latter, we note that the conditional expectation in (13.19) is taken in the annuity measure, suggesting that the expected values of $x(t)$ and $y(t)$ under $\mathrm{Q}^{A}$ provide a good expansion point.

Proposition 13.1.4. Let

$$
\bar{y}(t)=\mathrm{E}^{A}(y(t)) .
$$

Then, approximately,

$$
\begin{equation*}
\bar{y}(t) \approx h(t)^{2} \int_{0}^{t} \sigma_{r}^{0}(s)^{2} h(s)^{-2} d s, \quad t \in\left[0, T_{0}\right] \tag{13.23}
\end{equation*}
$$
^eq13-23

where, per (13.3),

$$
h(t)=\exp \left(-\int_{0}^{t} \varkappa(u) d u\right)
$$

Proof. Recall the dynamics of $y(t)$ in the quasi-Gaussian model (Proposition 13.1.1),

$$
d y(t)=\left(\sigma_{r}(t, x(t), y(t))^{2}-2 \varkappa(t) y(t)\right) d t, \quad y(0)=0
$$

Taking expected values, we obtain

$$
\begin{equation*}
d \mathrm{E}^{A}(y(t))=\left(\mathrm{E}^{A}\left(\sigma_{r}(t, x(t), y(t))^{2}\right)-2 \varkappa(t) \mathrm{E}^{A}(y(t))\right) d t \tag{13.24}
\end{equation*}
$$
^eq13-24

subject to $\mathrm{E}^{A}(y(0))=0$. Approximating, for the purposes of this calculation,

$$
\mathrm{E}^{A}\left(\sigma_{r}(t, x(t), y(t))^{2}\right) \approx \sigma_{r}^{0}(t)^{2}
$$

the equation (13.24) yields

$$
d \mathrm{E}^{A}(y(t)) \approx\left(\sigma_{r}^{0}(t)^{2}-2 \varkappa(t) \mathrm{E}^{A}(y(t))\right) d t, \quad \mathrm{E}^{A}(y(0))=0 .
$$

Solving this ODE leads to (13.23).
As above, let $X(t, s)$ be the function inverse, in $x$, to $S(t, x, \bar{y}(t))$, i.e.

$$
\begin{equation*}
S(t, X(t, s), \bar{y}(t)) \equiv s \tag{13.25}
\end{equation*}
$$
^eq13-25

and let $x_{0}(t)$ be given as the solution of

$$
\begin{equation*}
S\left(t, x_{0}(t), \bar{y}(t)\right)=S(0) \tag{13.26}
\end{equation*}
$$
^eq13-26

where $S(0)$ is the forward swap rate at time 0 .
Remark 13.1.5. The function $S(t, x, \bar{y}(t))$ is known in closed form from (13.15) and is smooth and monotonic in $x$. As such, (13.26) can be solved for $x_{0}(t)$ in just a few iterations of the Newton algorithm. A good starting point for the search is $x=0$.

It turns out that $x_{0}(t)$ in (13.26) is a good expansion point itself, and can also be used to calculate an even better one:

Lemma 13.1.6. The function $x_{0}(t)$ is a first-order approximation to $\mathrm{E}^{A}(x(t))$. An approximation to second order is given by

$$
\bar{x}(t)=x_{0}(t)+\frac{\partial^{2} X}{\partial s^{2}}(t, S(0)) \operatorname{Var}^{A}(S(t)), \quad t \in\left[0, T_{0}\right] .
$$

Proof. Expanding $X(t, s)$ around $s=S(0)$ to first order, we obtain

$$
\begin{equation*}
x(t)-x_{0}(t)=\left.\frac{\partial X}{\partial s}(t, s)\right|_{s=S(0)}(S(t)-S(0))+O\left((S(t)-S(0))^{2}\right) . \tag{13.27}
\end{equation*}
$$
^eq13-27

Taking expected values and using the fact that $S(t)$ is a martingale in measure $\mathrm{Q}^{A}$, we get

$$
\mathrm{E}^{A}(x(t))-x_{0}(t)=O\left(\mathrm{E}^{A}\left((S(t)-S(0))^{2}\right)\right)
$$

i.e. $x_{0}(t)$ is an approximation to $\mathrm{E}^{A}(x(t))$ to first order. The approximation $\bar{x}(t)$ is obtained by expanding (13.27) to second order,

$$
\begin{aligned}
x(t)-x_{0}(t)= & \frac{\partial X}{\partial s}(t, S(0))(S(t)-S(0)) \\
& +\frac{1}{2} \frac{\partial^{2} X}{\partial s^{2}}(t, S(0))(S(t)-S(0))^{2}+O\left((S(t)-S(0))^{3}\right)
\end{aligned}
$$

Then

$$
\mathrm{E}^{A}(x(t))=x_{0}(t)+\frac{\partial^{2} X}{\partial s^{2}}(t, S(0)) \operatorname{Var}^{A}(S(t))+O\left(\mathrm{E}^{A}\left((S(t)-S(0))^{3}\right)\right)
$$

Remark 13.1.7. As high precision is not required when calculating the variance $\operatorname{Var}^{A}(S(t))$ in Lemma 13.1.6, it can be evaluated by considering a simple Gaussian approximation to the dynamics of $S(t)$, i.e. using $x(t)=0$, $y(t)=0$ in (13.16),

$$
d S(t) \approx \frac{\partial S}{\partial x}(t, 0,0) \sigma_{r}^{0}(t) d W^{A}(t)
$$

which would yield

$$
\operatorname{Var}^{A}(S(t)) \approx \int_{0}^{t}\left(\frac{\partial S}{\partial x}(s, 0,0) \sigma_{r}^{0}(s)\right)^{2} d s
$$

The second derivative $\partial^{2} X(t, s) / \partial s^{2}$ can be computed by differentiating the implicit definition (13.25) twice.

Having now established an expansion point, let us proceed to determine an approximation to $X(t, s)$ with higher accuracy than the linear one in (13.21). Empirically, it can be observed that $S(t, x)$ is closely approximated as a quadratic function of $x$ across a wide range of the argument $x$. This suggests a second-order expansion of (13.25) around $\bar{x}(t)$, and approximating $X(t, s)$ with $\xi=\xi(t, s)$, the solution of the following quadratic equation in $\xi$ :

$$
\begin{align*}
S(t, \bar{x}(t), \bar{y}(t))+\frac{\partial S}{\partial x}(t, \bar{x}(t), \bar{y}(t))(\xi-\bar{x}(t)) & \\
& +\frac{1}{2} \frac{\partial^{2} S}{\partial x^{2}}(t, \bar{x}(t), \bar{y}(t))(\xi-\bar{x}(t))^{2}=s \tag{13.28}
\end{align*}
$$
^eq13-28

With $S(t, \bar{x}(t), \bar{y}(t)), \quad \partial S(t, \bar{x}(t), \bar{y}(t)) / \partial x$ and $\partial^{2} S(t, \bar{x}(t), \bar{y}(t)) / \partial x^{2}$ precomputed, the evaluation of $\xi(t, s)$ for any $s$ is essentially instantaneous, and we obtain the following efficient approximation for $\varphi(t, s)$.

Proposition 13.1.8. An approximation to $\varphi(t, s)$ in (13.19) is given by

$$
\varphi(t, s) \approx \frac{\partial S}{\partial x}(t, \xi(t, s), \bar{y}(t)) \sigma_{r}(t, \xi(t, s), \bar{y}(t))
$$

where $\xi(t, s)$ is the solution to the quadratic equation (13.28), with $\bar{x}(t)$ given by Lemma 13.1.6 and $\bar{y}(t)$ given by Proposition 13.1.4.

### 13.1.5 Linear Local Volatility

As demonstrated above, the local volatility qG model (13.9) has the flexibility to generate essentially arbitrary local volatility dynamics for swap rates. Using the results of Lemma 13.1.3 and Proposition 13.1.8, the function $\sigma_{r}(t, x, y)$ could therefore, in principle, be calibrated non-parametrically (see Dupire [1994] and the discussion in Section 7.1.3) to the implied volatilities of a collection of swaptions across all strikes. However, as explained in Section 7.1.3, we recommend the volatility function $\sigma_{r}(t, x, y)$ to be chosen from a parametric family of monotone, downward sloping functions of state variable(s). While power functions that give rise to models with CEV-type ${ }^{3}$ dynamics could be used, as explained in Remark 7.2.14 linear functions provide a less-involved alternative capable of producing essentially the same range of volatility smiles as CEV models.

With the above in mind, let us consider the following short rate local volatility function

$$
\begin{equation*}
\sigma_{r}(t, x, y)=\lambda_{r}(t)\left(\alpha_{r}(t)+b_{r}(t) x\right) \tag{13.29}
\end{equation*}
$$
^eq13-29

The scale function $\alpha_{r}(t)$ is redundant (as it can be absorbed in $\lambda_{r}(t)$ ) and may be set exogenously; Section 13.1.6 discusses a convenient choice. The functions $\lambda_{r}(t)$ (volatility) and $b_{r}(t)$ (skew) are calibrated to the market. Under (13.29), the local volatility of the swap rate $S$ is given, approximately, by

$$
\begin{equation*}
\varphi(t, s) \approx \lambda_{r}(t) \frac{\partial S}{\partial x}(t, \xi(t, s), \bar{y}(t))\left(\alpha_{r}(t)+b_{r}(t) \xi(t, s)\right) \tag{13.30}
\end{equation*}
$$
^eq13-30

with $\xi(t, s)$ as in Proposition 13.1.8. As the local volatility of the short rate is linear in $x$, it seems reasonable that the local volatility of the swap rate would be well approximated by a linear function as well. To exploit this, let $\varphi$ be given as in (13.30) and notice that

$$
\begin{aligned}
\varphi(t, S(0)) \approx & \lambda_{r}(t) \frac{\partial S}{\partial x}(t, \xi(t, S(0)), \bar{y}(t))\left(\alpha_{r}(t)+b_{r}(t) \xi(t, S(0))\right) \\
\frac{\partial \varphi}{\partial s}(t, S(0)) \approx & \lambda_{r}(t) \frac{\partial S}{\partial x}(t, \xi(t, S(0)), \bar{y}(t)) \frac{\partial \xi}{\partial s}(t, S(0)) \\
& \times\left[\frac{\frac{\partial^{2} S}{\partial x^{2}}(t, \xi(t, S(0)), \bar{y}(t))}{\frac{\partial S}{\partial x}(t, \xi(t, S(0)), \bar{y}(t))}\left(\alpha_{r}(t)+b_{r}(t) \xi(t, S(0))\right)+b_{r}(t)\right]
\end{aligned}
$$

[^41]Clearly

$$
\xi(t, S(0)) \approx \bar{x}(t)
$$

and, with $\xi(t, \cdot)$ being an approximate inverse to $S(t, \cdot, \bar{y}(t))$,

$$
\frac{\partial \xi}{\partial s}(t, S(0)) \approx \frac{1}{\frac{\partial S}{\partial x}(t, \bar{x}(t), \bar{y}(t))}
$$

It follows that

$$
\begin{align*}
\varphi(t, S(0)) & \approx \lambda_{r}(t) \frac{\partial S}{\partial x}(t, \bar{x}(t), \bar{y}(t))\left(\alpha_{r}(t)+b_{r}(t) \bar{x}(t)\right)  \tag{13.31}\\
\frac{\partial \varphi}{\partial s}(t, S(0)) & \approx \lambda_{r}(t)\left[\frac{\frac{\partial^{2} S}{\partial x^{2}}(t, \bar{x}(t), \bar{y}(t))}{\frac{\partial S}{\partial x}(t, \bar{x}(t), \bar{y}(t))}\left(\alpha_{r}(t)+b_{r}(t) \bar{x}(t)\right)+b_{r}(t)\right] \tag{13.32}
\end{align*}
$$
^eq13-31
^eq13-32

The following corollary to Proposition 13.1.8 holds.
Corollary 13.1.9. Under the assumption of linear local short rate volatility (13.29) for the quasi-Gaussian model (13.10), the dynamics of the swap rate $S(t)$ are approximated by

$$
\begin{equation*}
d S(t) \approx \lambda_{S}(t)\left(b_{S}(t) S(t)+\left(1-b_{S}(t)\right) S(0)\right) d W^{A}(t) \tag{13.33}
\end{equation*}
$$
^eq13-33

where

$$
\begin{align*}
\lambda_{S}(t) & =\lambda_{r}(t) \frac{1}{S(0)} \frac{\partial S}{\partial x}(t, \bar{x}(t), \bar{y}(t))\left(\alpha_{r}(t)+b_{r}(t) \bar{x}(t)\right)  \tag{13.34}\\
b_{S}(t) & =\frac{S(0)}{\left(\alpha_{r}(t)+b_{r}(t) \bar{x}(t)\right)} \frac{b_{r}(t)}{\frac{\partial S}{\partial x}(t, \bar{x}(t), \bar{y}(t))}+\frac{S(0) \frac{\partial^{2} S}{\partial x^{2}}(t, \bar{x}(t), \bar{y}(t))}{\left(\frac{\partial S}{\partial x}(t, \bar{x}(t), \bar{y}(t))\right)^{2}} \tag{13.35}
\end{align*}
$$
^eq13-34
^eq13-35

Proof. Under a linear approximation to the local volatility function of the swap rate we have, with $\varphi$ defined in (13.30),

$$
d S(t) \approx\left(\varphi(t, S(0))+\frac{\partial \varphi}{\partial s}(t, S(0))(S(t)-S(0))\right) d W^{A}(t)
$$

which, after rearranging the terms, yields

$$
d S(t) \approx \frac{\varphi(t, S(0))}{S(0)}\left(S(0)+S(0) \frac{\partial \varphi(t, S(0)) / \partial s}{\varphi(t, S(0))}(S(t)-S(0))\right) d W^{A}(t)
$$

Defining

$$
\lambda_{S}(t) \triangleq \frac{\varphi(t, S(0))}{S(0)}, \quad b_{S}(t) \triangleq S(0) \frac{\partial \varphi(t, S(0)) / \partial s}{\varphi(t, S(0))}
$$

the result follows from (13.31)-(13.32).
We recognize (13.33) as a displaced log-normal SDE with time-dependent volatility $\lambda_{S}(t)$ and skew $b_{S}(t)$. Using averaging techniques from Section 7.6.2, we can convert it into a displaced log-normal SDE with time-constant parameters $\bar{\lambda}_{S}$ and $\bar{b}_{S}$, see Proposition 7.2.12. For convenience, we list the resulting swaption pricing formula below.

Proposition 13.1.10. Consider a payer swaption with strike (i.e., coupon) $c$ and expiry $T_{0}$ on the swap rate $S(t)$ defined in (13.13). In the quasiGaussian model (13.10) with linear short rate volatility (13.29), the swaption price can be approximated by the displaced log-normal option formula

$$
\begin{aligned}
& V_{\text {swaption }}(0) \approx A(0)\left[\left(S(0)+S(0)\left(1-\bar{b}_{S}\right) / \bar{b}_{S}\right) \Phi\left(d_{+}\right)\right. \\
&\left.-\left(c+S(0)\left(1-\bar{b}_{S}\right) / \bar{b}_{S}\right) \Phi\left(d_{-}\right)\right]
\end{aligned},
$$

where

$$
\begin{align*}
\bar{\lambda}_{S} & =\left(\frac{1}{T_{0}} \int_{0}^{T_{0}} \lambda_{S}(t)^{2} d t\right)^{1 / 2}  \tag{13.36}\\
\bar{b}_{S} & =\int_{0}^{T_{0}} b_{S}(t) w_{S}(t) d t  \tag{13.37}\\
w_{S}(t) & =\frac{\lambda_{S}(t)^{2} \int_{0}^{t} \lambda_{S}(s)^{2} d s}{\int_{0}^{T_{0}}\left(\lambda_{S}(u)^{2} \int_{0}^{u} \lambda_{S}(s)^{2} d s\right) d u}
\end{align*}
$$
^eq13-36
^eq13-37

with $\lambda_{S}(t), b_{S}(t)$ given by Corollary 13.1.9.

### 13.1.6 Linear Local Volatility for a Swaption Strip

The quasi-Gaussian model (13.10) is typically calibrated to a swaption strip (recall the definition in Section 10.1.4) on a maturity grid $0=T_{0}<\ldots<T_{N}$, i.e. a collection of $N-1$ swaptions with the $n$-th swaption expiring on $T_{n}$, $n=1, \ldots, N-1$. Let us suppose a swaption strip is specified, and that the $n$-th swaption has an underlying swap with $\mu(n)$ periods. For each $n$, we denote the corresponding swap rate and the annuity by

$$
S_{n}(t) \triangleq S_{n, \mu(n)}(t), \quad A_{n}(t) \triangleq A_{n, \mu(n)}(t), \quad n=1, \ldots, N-1
$$

For the model with the linear volatility specification (13.29), Proposition 13.1.10 will, after proper adjustment of maturities (see footnote 6 in Chapter 10 ), allow us to value all swaptions of all strikes in the strip by using the
displaced log-normal model with the effective parameters computed from the local volatility function of the model.

As mentioned earlier, in the specification (13.29) the function $\alpha_{r}(t)$ is only included for convenient scaling. To find a good value for it, let us consider the relationship between the local swap rate skews $b_{S_{\prime \prime}}(t)$ and the local short rate skew $b_{r}(t)$ in (13.35). Ignoring small terms,

$$
\begin{equation*}
b_{S_{n}}(t) \approx \frac{S_{n}(0)}{\alpha_{r}(t)} \frac{b_{r}(t)}{\partial S_{n}(t, \bar{x}(t), \bar{y}(t)) / \partial x} \tag{13.38}
\end{equation*}
$$
^eq13-38

It is often convenient to parametrize the model in such a way that the values of model parameters (here $b_{r}(t)$ ) are roughly of the same order of magnitude as the output parameters (here $b_{S_{n}}(t)$ ). This allows one to quickly check whether model parameters are sensible, and may well lead to better numerical properties of the calibration algorithm. Based on (13.38), we elect to set $\alpha_{r}(t)$ equal to $S_{n}(0)$ (of course we need to account for different values of $n$ ), and rescale $b_{r}(t)$ to incorporate the term $\partial S_{n} / \partial x$ (again, for different $n)$. In summary, we specialize the definition (13.29) to be

$$
\begin{equation*}
\sigma_{r}(t, x, y)=\sum_{n=1}^{N-1} \lambda_{n}\left(S_{n}(0)+b_{n} D_{n} x\right) 1_{\left\{t \in\left(T_{n-1}, T_{n}\right]\right\}}, \tag{13.39}
\end{equation*}
$$
^eq13-39

where the skew scalings $D_{n}$ are given by

$$
D_{n}=\frac{\partial S_{n}}{\partial x}(t, 0,0)
$$

This definition recognizes the fact that the behavior of $\lambda_{r}(t)$ and $b_{r}(t)$ between the knot dates $\left\{T_{n}\right\}$ is of no consequence, wherefore these functions can be taken to be piecewise constant,

$$
\begin{aligned}
& \lambda_{r}(t)=\sum_{n=1}^{N-1} \lambda_{n} 1_{\left\{t \in\left(T_{n-1}, T_{n}\right)\right\}} \\
& \alpha_{r}(t)=\sum_{n=1}^{N-1} S_{n}(0) 1_{\left\{t \in\left(T_{n-1}, T_{n}\right]\right\}} \\
& b_{r}(t)=\sum_{n=1}^{N-1} b_{n} D_{n} 1_{\left\{t \in\left(T_{n-1}, T_{n}\right]\right\}}
\end{aligned}
$$

### 13.1.7 Volatility Calibration

Let us assume that a swaption strip is given, and the model is parameterized with the local volatility of the form (13.39). The model parameters $\lambda_{n}$ and $b_{n}, n=1, \ldots, N-1$, need to be determined by calibrating the model to
market prices of swaptions in the swaption strip. For now, let us suppose that the mean reversion function $\varkappa(t)$ in (13.10) is specified externally we will return to its calibration later in the chapter.

As each swap rate has an approximately displaced log-normal distribution in the model, the calibration objective could be expressed as the problem of matching displaced log-normal parameters, as given by the model for each swap rate, to a similar set of market-implied parameters. We already saw a similar approach in Section 9.3.4 and recall that performing the calibration in model parameter space, rather than in the space of calibration instrument values, avoids the expense of invoking option pricing formulas within the calibration loop.

Accordingly, assume that a collection of market parameters, i.e. displaced log-normal volatilities and skews $\left(\hat{\lambda}_{S_{n}}, \hat{b}_{S_{n}}\right), n=1, \ldots, N-1$, is given. In practice, these parameters are obtained by fitting a series of constantparameter displaced log-normal vanilla models to the observed swaption volatility smiles at all expiries $T_{1}, \ldots, T_{N-1}$. A best-fit model calibration across swaption strikes is possible (at the expense of using a numerical optimizer), or we may simply set the volatility $\widehat{\lambda}_{S}$, to match at-the-money swaption volatilities, while the skew $\widehat{b}_{S_{n}}$ is fit to the slope of the volatility smile at-the-money, or to the volatility at some relevant non-ATM strike.

It is clear from the swaption pricing formula in Proposition 13.1.10 that the value of a swaption with expiry $T_{n}$ depends on model parameters ( $\lambda_{i}, b_{i}$ ) for $i=1, \ldots, n$ only. Hence, the qG model can be calibrated by a bootstrap method, similarly to the pure Gaussian case from Section 10.1.4. In the bootstrap method, the equations (13.36)-(13.37) are solved sequentially for $n=1, \ldots, N-1$, with the two equations on step $n$ used to determine two unknown model parameters $\left(\lambda_{n}, b_{n}\right)$. For example, the following algorithm could be used.

1. Set $\left(\lambda_{n}, b_{n}\right), n=1, \ldots, N-1$, to some reasonable starting values, e.g. set $\lambda_{n}$ 's to (properly scaled) volatilities obtained by calibrating a pure Gaussian model as in Section 10.1.4, and $b_{n}=\widehat{b}_{S_{n}}$.
2. Set $n=1$.
3. For given $n,\left(\lambda_{\imath}^{*}, b_{i}^{*}\right)$ are known for $i=1, \ldots, n-1$. Note we use a star to denote calibrated values of the model parameters.
4. Calculate $\bar{x}(t)$ (Lemma 13.1.6) and $\bar{y}(t)$ (Proposition 13.1.4) for $t \in \left[0, T_{n}\right]$ using $\left(\lambda_{i}^{*}, b_{i}^{*}\right), i=1, \ldots, n-1$, and the initial guess for $\left(\lambda_{n}, b_{n}\right)$ from Step 1. Note that $\bar{x}(t), \bar{y}(t)$ implicitly depend on $n$ as their definition depends on the swap rate/annuity measure used.
5. Calculate $\lambda_{S_{n}}(t), b_{S_{\prime \prime}}(t)$ for $t \in\left[0, T_{n-1}\right]$ from $\left(\lambda_{i}^{*}, b_{i}^{*}\right), i=1, \ldots, n-1$, using (13.34)-(13.35).
6. Make another guess for $\left(\lambda_{n}, b_{n}\right)$.
7. Update $\lambda_{S_{n}}(t), b_{S_{n}}(t)$ for $t \in\left(T_{n-1}, T_{n}\right]$ from $\left(\lambda_{i}^{*}, b_{i}^{*}\right), i=1, \ldots, n-1$, using (13.34)-(13.35).
8. Calculate $\bar{\lambda}_{S^{\prime \prime}}, \bar{b}_{S_{n}}$ using Proposition 13.1.10.
9. Compare ( $\bar{\lambda}_{S_{n}}, \bar{b}_{S_{n}}$ ) to ( $\hat{\lambda}_{S_{n}}, \hat{b}_{S_{n}}$ ). If not equal within given tolerance, go to Step 6. Otherwise, proceed to Step 10.
10. As we have reached acceptable convergence between $\left(\bar{\lambda}_{S_{n}}, \bar{b}_{S_{n}}\right)$ and ( $\hat{\lambda}_{S_{n}}, \hat{b}_{S_{n}}$ ), set the calibrated model parameter values to the latest trial values, $\left(\lambda_{n}^{*}, b_{n}^{*}\right)=\left(\lambda_{n}, b_{n}\right)$.
11. Update $n \rightarrow n+1$. If $n \leq N-1$ go to Step 3. Otherwise, conclude.

It may appear more accurate to make Step 4 a part of the calibration loop (for each $n$ ), with $\bar{x}(t), \bar{y}(t), t \in\left(T_{n-1}, T_{n}\right]$, and dependent quantities such as $\partial S_{n}(t, \bar{x}(t), \bar{y}(t)) / \partial x$, etc. computed using the current guess for $\left(\lambda_{n}, b_{n}\right)$ in each loop iteration (together with already-calibrated values $\left(\lambda_{i}^{*}, b_{i}^{*}\right)$, $i=1, \ldots, n-1)$. While such extensions are relatively straightforward, our experience shows that the quality of the calibration is rarely improved enough to justify the additional complexity.

It is possible to add regularity terms to Step 9 of the algorithm above, to ensure that the resulting model parameters do not behave irregularly as functions of $t$. See Section 13.1.8.2 below for an example of this.

### 13.1.8 Mean Reversion Calibration

With the volatility calibration out of the way, let us discuss what to do with the remaining model parameter, the mean reversion function $\varkappa(t)$ in (13.10). We start with a short review of the effects of mean reversion.

### 13.1.8.1 Effects of Mean Reversion

Let us first consider a few simple examples as a way of building intuition about the effects of mean reversion on market values of various securities. For simplicity, we use continuously compounded rates as convenient proxies for Libor and swap rates, and consider a pure Gaussian model with constant volatility and mean reversion,

$$
\sigma_{r}(t) \equiv \sigma_{r}, \quad \varkappa(t) \equiv \varkappa
$$

A continuously compounded forward yield over a period $[T, M]$, observed at time $t$, is given by ${ }^{4}$

$$
F(t, T, M)=-\frac{1}{M-T} \ln \frac{P(t, M)}{P(t, T)}
$$

[^42]According to Proposition 10.1.7, the forward yield can be expressed as a function of the state variables and parameters of the model,

$$
\begin{aligned}
F(t, T, M) & =\frac{G(t, M)-G(t, T)}{M-T} x(t)+\frac{1}{2} \frac{G(t, M)^{2}-G(t, T)^{2}}{M-T} y(t), \\
G(t, u) & =\frac{1-e^{-\varkappa(u-t)}}{\varkappa} .
\end{aligned}
$$

Recall that $y(t)$ is deterministic in the Gaussian case, so the standard deviation of $F(T, T, M)$ is equal to

$$
\begin{aligned}
\operatorname{Stdev}(F(T, T, M)) & =\frac{G(T, M)-G(T, T)}{M-T}(\operatorname{Var}(x(T)))^{1 / 2} \\
& =\frac{1-e^{-\varkappa(M-T)}}{\varkappa(M-T)}(\operatorname{Var}(x(T)))^{1 / 2}
\end{aligned}
$$

For two maturities $M_{1}, M_{2}, T \leq M_{1} \leq M_{2}$, we therefore observe that

$$
\frac{\operatorname{Stdev}\left(F\left(T, T, M_{2}\right)\right)}{\operatorname{Stdev}\left(F\left(T, T, M_{1}\right)\right)}=\left(\frac{1-e^{-\varkappa\left(M_{2}-T\right)}}{\varkappa\left(M_{2}-T\right)}\right) /\left(\frac{1-e^{-\varkappa\left(M_{1}-T\right)}}{\varkappa\left(M_{1}-T\right)}\right),
$$

i.e the ratio of standard deviations of two forward yields with the same expiry $T$ but different tenors is independent of the volatility parameter $\sigma_{r}$, and solely determined by the mean reversion parameter $\varkappa$. Since standard deviations of forward yields can loosely be thought of as proxies for implied swaption volatilities, we observe that the mean reversion parameter changes the relative levels of implied volatilities of swaptions with the same expiry but different underlying swap tenors. Specifically, for a fixed level of volatility $\sigma_{r}$, an increase in mean reversion makes the volatility of a longer-tenor swaption decrease relative to the volatility of a shorter-tenor swaption, assuming both have the same expiry date.

With the mean reversion effect above in mind, consider (say) a caplet and a swaption with the same expiry, and imagine an experiment in which the mean reversion is changed, but the volatility $\sigma_{r}$ is adjusted to keep the implied volatility of the caplet unchanged. It should be clear from the discussion above that as mean reversion increases, the swaption volatility will decrease. If, instead, the market volatility of the swaption is kept constant by adjusting $\sigma_{r}$ for each level of mean reversion, then the caplet volatility will increase when the mean reversion increases. Such complementarity of effects of the mean reversion $\varkappa$ and volatility $\sigma_{r}$ allows us, in principle, to set both in such a way that we match the market-implied volatilities of both the caplet and the swaption.

For an alternative look at the effect of mean reversion, let us consider two forward rates with different fixing dates, $F\left(T_{1}, T_{1}, M_{1}\right)$ and $F\left(T_{2}, T_{2}, M_{2}\right)$ with $T_{1} \leq T_{2}$. Observing that

$$
x(t)=\sigma_{r} \int_{0}^{t} e^{-\varkappa(t-u)} d W(u)
$$

it is easy to establish that

$$
\begin{aligned}
& \operatorname{Corr}\left(F\left(T_{1}, T_{1}, M_{1}\right), F\left(T_{2}, T_{2}, M_{2}\right)\right) \\
& \quad=\operatorname{Corr}\left(x\left(T_{1}\right), x\left(T_{2}\right)\right)=e^{-\varkappa\left(T_{2}-T_{1}\right)}\left(\frac{1-e^{-2 \varkappa T_{1}}}{1-e^{-2 \varkappa T_{2}}}\right)^{1 / 2}
\end{aligned}
$$

This correlation, which we can call inter-temporal correlation as the forward yields are observed at different times, depends only on $\varkappa$, i.e. on the autocorrelation properties of the process $x(t)$. At a level of $\varkappa=0$, the intertemporal correlation is $\left(T_{1} / T_{2}\right)^{1 / 2}$ and decreases to 0 as $\varkappa \rightarrow \infty$.

The dependence of inter-temporal correlation on mean reversion is potentially useful for calibration purposes. To see this, consider a hypothetical security, a basket option on a set of rates with different expiry dates, with a payoff

$$
\max \left\{F\left(T_{n}, T_{n}, M_{n}\right), n=1, \ldots, N-1\right\}
$$

It is well-known (and intuitively obvious) that prices of basket options are decreasing functions of correlation, hence we expect the price of the contract above to be increasing in the mean reversion $\varkappa$, ceteris paribus. If such a basket option were traded in the market, the mean reversion could in principle be implied from its price.

While basket options on forward yields are, of course, not traded outright, the example above is not as far-fetched as it might seem since a Bermudan swaption (see Section 5.12) gives the holder a right to exercise into one of several different swaps that, critically, are observed on different exercise dates. The Bermudan swaption is therefore conceptually similar to a basket option on a strip of swap rates, with each rate fixing on its own fixing date. The implications of this analogy, and the effect of mean reversion on inter-temporal correlations, will be exploited later in developing the local projection method for Bermudan swaptions (see Chapter 19).

### 13.1.8.2 Calibrating Mean Reversion to Volatility Ratios

We consider the qG model with the volatility function (13.39). In Section 13.1.7 we developed the volatility calibration algorithm for this model under the assumption that the mean reversion function $\varkappa(t)$ was already available. Since it is often beneficial to calibrate different model parameters separately, the calibration of mean reversion should then ideally not require the knowledge of the volatility $\sigma_{r}(t, x, y)$ of the model.

As presented earlier, the qG model volatility calibration involves a strip of swap rates $\left\{S_{n}(\cdot)\right\}_{n=1}^{N-1}$, with the rate $S_{n}(\cdot)$ fixing on $T_{n}$ and having $\mu(n)$ periods; these swap rates define the volatility function $\sigma_{r}(t, x, y)$ in (13.39)
and serve as volatility calibration targets. As indicated in Section 13.1.8.1, the ratio of volatilities of two swaptions with the same expiry date is more-or-less independent of volatility, suggesting the usage of a second strip ${ }^{5}$ of swap rates to define targets for mean reversion calibration as the ratios between these rates and the original ones.

We assume that a second strip of swap rates is given, with the $n$-th rate fixing on $T_{n}$ and having $\nu(n)$ periods, where $\nu(n) \neq \mu(n)$. We use extended notations to distinguish the two strips, with

$$
\left\{S_{n, \mu(n)}(\cdot)\right\}_{n=1}^{N-1}
$$

used for the original strip and

$$
\left\{S_{n, \nu(n)}(\cdot)\right\}_{n=1}^{N-1}
$$

used for the additional one. The mean reversion is then calibrated to the pairwise ratios of implied volatilities of $S_{n, \nu(n)}\left(T_{n}\right)$ and $S_{n, \mu(n)}\left(T_{n}\right)$, $n=1, \ldots, N-1$. The following result forms the basis of mean reversion calibration.

Proposition 13.1.11. In the quasi-Gaussian model with local short rate volatility function (13.9), the ratio of variances of two swap rates fixing on $T_{n}$ with $m_{1}$ and $m_{2}$ periods, respectively, is approximately given by either

$$
\begin{equation*}
\frac{\operatorname{Var}\left(S_{n, m_{1}}\left(T_{n}\right)\right)}{\operatorname{Var}\left(S_{n, m_{2}}\left(T_{n}\right)\right)} \approx \frac{\int_{0}^{T_{n}}\left(\sigma_{r}^{0}(t) \frac{\partial S_{n, m_{1}}}{\partial x}(t, 0,0)\right)^{2} d t}{\int_{0}^{T_{n}}\left(\sigma_{r}^{0}(t) \frac{\partial S_{n, m_{2}}}{\partial x}(t, 0,0)\right)^{2} d t} \tag{13.40}
\end{equation*}
$$
^eq13-40

where $\sigma_{r}^{0}(t)$ is defined in (13.12), or

$$
\begin{equation*}
\frac{\operatorname{Var}\left(S_{n, m_{1}}\left(T_{n}\right)\right)}{\operatorname{Var}\left(S_{n, m_{2}}\left(T_{n}\right)\right)} \approx \frac{\int_{0}^{T_{n}}\left(\frac{\partial S_{n, m_{1}}}{\partial x}(t, 0,0)\right)^{2} d t}{\int_{0}^{T_{n}}\left(\frac{\partial S_{n, m_{2}}}{\partial x}(t, 0,0)\right)^{2} d t} \tag{13.41}
\end{equation*}
$$
^eq13-41

Proof. The first formula is obtained by using $x(t)=0, y(t)=0$ in (13.16); see also Remark 13.1.7. The second one follows from the first under the approximation of the time-dependent volatility $\sigma_{r}^{0}(t)$ by a constant,

$$
\begin{equation*}
\sigma_{r}^{0}(t) \approx \sigma_{r}^{0}(0) \tag{13.42}
\end{equation*}
$$
^eq13-42

[^43]Remark 13.1.12. In practice, using the simpler approximation (13.41) instead of (13.40) does not reduce the accuracy of mean reversion calibration much. However, if the more accurate formula (13.40) is preferred, one can use an estimate of $\sigma_{r}^{0}(t)$ obtained from, for example, a pre-calibration of a pure Gaussian model.

Proposition 13.1.11 suggests an algorithm for mean reversion calibration. Recalling formulas (13.17), (13.3) we notice that the ratios

$$
\begin{equation*}
\left\{\frac{\operatorname{Var}\left(S_{n, \nu(n)}\left(T_{n}\right)\right)}{\operatorname{Var}\left(S_{n, \mu(n)}\left(T_{n}\right)\right)}\right\}_{n=1}^{N-1} \tag{13.43}
\end{equation*}
$$
^eq13-43

as computed in (13.41) depend on the mean reversion parameter $\varkappa(t)$ only, allowing us to set up an optimization problem in which $\varkappa(t)$ is chosen to match the ratios (13.43) to their market-implied values.

For a possible calibration algorithm, suppose that we discretize $\varkappa(t)$ on the time grid,

$$
\varkappa(t)=\sum_{n=1}^{N-1} \varkappa_{n} \times 1_{\left\{t \in\left(T_{n-1}, T_{n}\right]\right\}}+\varkappa_{N} \times 1_{\left\{t \in\left(T_{N-1}, \infty\right)\right\}}
$$

It is generally not advisable to find these mean reversions only by bestfitting to variance ratios, as the function $\varkappa(t)$ would likely end up being quite irregular. To prevent this, we suggest the inclusion of a regularization term penalizing non-stationary behavior. While there are multiple ways of doing this, one simple approach would be to set optimal mean reversion levels $\left\{\varkappa_{n}^{*}\right\}_{n=1}^{N}$ as the solution to the following optimization problem,

$$
\begin{align*}
\left\{\varkappa_{n}^{*}\right\}=\operatorname{argmin} & \left\{\sum _ { n = 1 } ^ { N - 1 } \left(\frac{\operatorname{Var}\left(S_{n, \nu(n)}\left(T_{n}\right)\right)}{\operatorname{Var}\left(S_{n, \mu(n)}\left(T_{n}\right)\right)}\left(\left\{\varkappa_{n}\right\}\right)\right.\right. \\
& \left.\left.-\frac{\widehat{\operatorname{Var}}\left(S_{n, \nu(n)}\left(T_{n}\right)\right)}{\widehat{\operatorname{Var}}\left(S_{n, \mu(n)}\left(T_{n}\right)\right)}\right)^{2}+w \sum_{n=1}^{N-1}\left(\varkappa_{n+1}-\varkappa_{n}\right)^{2}\right\} \tag{13.44}
\end{align*}
$$
^eq13-44

where

$$
\frac{\operatorname{Var}\left(S_{n, \nu(n)}\left(T_{n}\right)\right)}{\operatorname{Var}\left(S_{n, \mu(n)}\left(T_{n}\right)\right)}\left(\left\{\varkappa_{n}\right\}\right)=\frac{\int_{0}^{T_{\|}}\left(\frac{\partial S_{\|, \nu(\|)}}{\partial x}(t, 0,0)\right)^{2} d t}{\int_{0}^{T_{n}}\left(\frac{\partial S_{n, \mu(m)}}{\partial x}(t, 0,0)\right)^{2} d t}
$$

Here $w>0$ is a user-specified regularization weight, and the $\widehat{\operatorname{Var}}\left(S_{n, m}\right)$ 's are market-implied variances of swap rates. This minimization problem is easily solved by standard non-linear optimization methods, to be discussed in some detail later (in particular in Section 14.5).

At this point the reader might, of course, wonder whether target variances of swap rates, as required for the mean reversion calibration, could indeed be observed in the market. The answer is yes: in general, a market-implied variance of a swap rate can be calculated directly from values of options on the swap rate (i.e. swaptions) across a collection of strikes, as discussed in detail in Chapter 16. In the most commonly used linear local short rate volatility ( 13.29 ) case, we may simplify matters further. Indeed, in this model the market volatility parameters of swap rates are the displaced log-normal volatilities $\widehat{\lambda}_{S_{n, m}}$ (see (13.33) and Section 13.1.7), in which case the market-implied variance of a swap rate can be approximated by

$$
\begin{equation*}
\widehat{\operatorname{Var}}\left(S_{n, m}\left(T_{n}\right)\right) \approx\left(S_{n, m}(0) \hat{\lambda}_{S_{n, m}}\right)^{2} T_{n} \tag{13.45}
\end{equation*}
$$
^eq13-45

Hence, in this case the mean reversion calibration targets could be specified directly in terms of market-implied displaced log-normal volatilities. In any case, when paired up with the volatility calibration procedure from Section 13.1.7, numerical solution of the optimization problem (13.44) ultimately allows us to calibrate the model to market-implied values of two separate swaption strips.

### 13.1.8.3 Calibrating Mean Reversion to Inter-Temporal Correlations

In the previous section, we took advantage of the observation from Section 13.1.8.1 that a ratio of standard deviations of two rates with the same expiry in a one-factor constant-volatility Gaussian model is independent of volatility. In this section, we develop a different calibration method for mean reversion, based on our earlier observation that inter-temporal correlations between forward yields are also nearly independent of volatility.

We focus on correlations (for the rest of this section, we omit the "intertemporal" qualifier for brevity) of the original strip of swap rates $\left\{S_{n}(\cdot)\right\}_{n=1}^{N-1}$. While a different swap rate strip could be used, the practical importance of such a generalization is limited.

Proposition 13.1.13. Let $\sigma_{r}^{0}(t)$ be given by (13.12). The correlation between two rates $S_{n_{1}}\left(T_{n_{1}}\right)$ and $S_{n_{2}}\left(T_{n_{2}}\right), n_{1} \leq n_{2}$, in the quasi-Gaussian model with the general local volatility function (13.9) can be approximated by either

$$
\begin{align*}
\operatorname{Corr} & \left(S_{n_{1}}\left(T_{n_{1}}\right), S_{n_{2}}\left(T_{n_{2}}\right)\right) \\
\approx & \int_{0}^{T_{n_{1}}}\left(\frac{\partial S_{n_{1}}}{\partial x}(t, 0,0)\right)\left(\frac{\partial S_{n_{2}}}{\partial x}(t, 0,0)\right) \sigma_{r}^{0}(t)^{2} d t \\
& \times\left(\int_{0}^{T_{n_{1}}}\left(\frac{\partial S_{n_{1}}}{\partial x}(t, 0,0)\right)^{2} \sigma_{r}^{0}(t)^{2} d t\right)^{-1 / 2} \\
& \times\left(\int_{0}^{T_{n_{2}}}\left(\frac{\partial S_{n_{2}}}{\partial x}(t, 0,0)\right)^{2} \sigma_{r}^{0}(t)^{2} d t\right)^{-1 / 2} \tag{13.46}
\end{align*}
$$
^eq13-46

or

$$
\begin{align*}
& \operatorname{Corr}\left(S_{n_{1}}\left(T_{n_{1}}\right), S_{n_{2}}\left(T_{n_{2}}\right)\right) \approx \int_{0}^{T_{n_{1}}}\left(\frac{\partial S_{n_{1}}}{\partial x}(t, 0,0)\right)\left(\frac{\partial S_{n_{2}}}{\partial x}(t, 0,0)\right) d t \\
& \times\left(\int_{0}^{T_{n_{1}}}\left(\frac{\partial S_{n_{1}}}{\partial x}(t, 0,0)\right)^{2} d t\right)^{-1 / 2}\left(\int_{0}^{T_{n_{2}}}\left(\frac{\partial S_{n_{2}}}{\partial x}(t, 0,0)\right)^{2} d t\right)^{-1 / 2} \tag{13.47}
\end{align*}
$$
^eq13-47

Proof. By Proposition 13.1.2,

$$
d S_{n_{i}}(t)=\frac{\partial S_{n_{i}}}{\partial x}(t, x(t), y(t)) \sigma_{r}(t, x(t), y(t)) d W^{A_{n_{i}}}(t)
$$

in $\mathrm{Q}^{A_{n_{i}}}, i=1,2$. In the risk-neutral measure, the SDE for $S_{n_{i}}(t)$ will have a stochastic drift; however, for the purposes of calculating the correlations in question, we ignore drift contributions. Thus, in the risk-neutral measure $Q$ we have approximately that

$$
\begin{aligned}
& \mathrm{E}\left(S_{n_{1}}\left(T_{n_{1}}\right) S_{n_{2}}\left(T_{n_{2}}\right)\right)-S_{n_{1}}(0) S_{n_{2}}(0) \\
& \quad=\mathrm{E}\left(S_{n_{1}}\left(T_{n_{1}}\right) S_{n_{2}}\left(T_{n_{1}}\right)\right)-S_{n_{1}}(0) S_{n_{2}}(0) \\
& \quad=\mathrm{E} \int_{0}^{T_{n_{1}}} \frac{\partial S_{n_{1}}}{\partial x}(t, x(t), y(t)) \frac{\partial S_{n_{2}}}{\partial x}(t, x(t), y(t)) \sigma_{r}(t, x(t), y(t))^{2} d t
\end{aligned}
$$

Using (13.12) and approximating all functions of state variables with their values at $x=y=0$, we obtain

$$
\begin{aligned}
\mathrm{E}\left(S_{n_{1}}\left(T_{n_{1}}\right) S_{n_{2}}\left(T_{n_{2}}\right)\right. & )-S_{n_{1}}(0) S_{n_{2}}(0) \\
& \approx \int_{0}^{T_{n_{1}}}\left(\frac{\partial S_{n_{1}}}{\partial x}(t, 0,0)\right)\left(\frac{\partial S_{n_{2}}}{\partial x}(t, 0,0)\right) \sigma_{r}^{0}(t)^{2} d t
\end{aligned}
$$

Similar formulas hold for $\mathrm{E}\left(S_{n_{i}}\left(T_{n_{i}}\right)^{2}\right), i=1,2$, and (13.46) follows. The result (13.47) follows from (13.46) by approximating the deterministic volatility function $\sigma_{r}^{0}(t)$ with a constant,

$$
\sigma_{r}^{0}(t) \approx \sigma_{r}^{0}(0)
$$

口
If a "market-implied" correlation matrix

$$
\widehat{\chi}_{n_{1}, n_{2}}=\operatorname{Corr}\left(S_{n_{1}}\left(T_{n_{1}}\right), S_{n_{2}}\left(T_{n_{2}}\right)\right), \quad 1 \leq n_{1} \leq n_{2}<N-1,
$$

is somehow known - extracted, for example, from Bermudan swaption prices that, per discussion in Section 13.1.8.1, depend strongly on such correlations - then Proposition 13.1.13 can be used to calibrate the mean reversion function $\varkappa(t)$ to this matrix. The formula (13.47) is independent of the volatility term $\sigma_{r}(t, x, y)$, and the mean reversion calibration can precede volatility calibration. Alternatively, the more accurate formula (13.46) could be used with the volatility obtained from a pre-calibration of a pure Gaussian model.

### 13.1.8.4 Final Comments on Mean Reversion Calibration

The reader has undoubtedly noticed that the formulas for mean reversion calibration were derived using rather crude approximations. Refinements are certainly possible, but achieving a high level of accuracy for mean reversion calibration was never our objective. Indeed, while it is possible to execute a global calibration to vanilla option in which the mean reversion and the volatility function are calibrated together using numerical valuation methods (such as the PDE method, see later in the chapter), we believe that in setting mean reversion, market information should be a rough guide rather than a "hard" calibration target. For this reason, we recommend using mean reversion calibration to match ratios of volatilities or inter-temporal correlations in an approximate sense only; subsequently, we can apply much more precise volatility calibration to recover our main targets, namely the implied volatilities of the primary swaption strip.

Calibrating a model with relatively limited set of parameters to market inevitably leads to time-dependent parameters and, subsequently, questions about the stationarity of the resulting volatility structure. In most applications, a completely time-stationary model would yield a clearly unacceptable fit to market data, and a certain degree of non-stationarity is unavoidable. However, as far as mean reversion calibration is concerned, we often advocate using a constant mean reversion function $\varkappa(t) \equiv \varkappa$ in the calibration routines developed previously. For example, we can use a constant mean reversion $\varkappa$ to roughly match the volatility ratios of caplets and swaptions, and then calibrate a time-dependent volatility function to match the volatilities of swaptions in the swaption strip exactly. To encourage even more time stationarity in the model, it is also possible to set (through an optimizer) the mean reversion $\varkappa$ in such a way that the volatility parameters of the calibrated model (i.e. the $\lambda_{n}$ 's) are as close to constant as possible.

Let us note that there is one instance of the application of the model where we recommend time-dependent mean reversion, namely when the
model is used as part of the local projection method.. The method is developed in detail later in the book (see e.g. Sections 18.4, 19.2, 20.1.3, 20.2.1), but it can loosely be described as using a "small" model, such as the qG model, as a local or instrument-specific proxy for a "big" model, such as a Libor market model (see Chapter 14). In this case the dynamics of the volatility structure are defined by the "big" and, hopefully, realistic model, and the local model is effectively just a mechanism to reduce numerical complexity of valuation.

Finally, we note that it is, of course, also possible to use the qG model without explicitly calibrating the mean reversion parameter. For example, it could be exposed as an "exotic risk" parameter and set exogenously by a trader to reflect his estimation of the market prices of Bermudan swaptions or other exotic securities. Such practice is, we believe, quite common.

### 13.1.9 Numerical Methods

### 13.1.9.1 Direct Integration

We start the discussion of numerical methods for pricing derivatives in the quasi-Gaussian model by deriving an approximation to the density of the state variables $x(T)$ and $y(T), T>0$. Our approximation is constructed to be suitable for small $T$; we find that it has good accuracy for $T$ around 1-2 years or less, depending on the level of volatility. While usable for valuing European-style derivatives by direct integration - a method generally preferable to PDE or Monte Carlo methods when available - the real utility of having a probability density comes in improving the accuracy of the PDE method, as described in Sections 2.8.2 and 2.8.3.

Consider a contract with a payoff $V(x(T), y(T))$ at time $T$. Its value at time 0 is given by

$$
V_{0}=P(0, T) \mathrm{E}^{T}(V(x(T), y(T))),
$$

where $\mathrm{E}^{T}$ (as always) denotes expectation in the $T$-forward measure $\mathrm{Q}^{T}$. Accordingly, we seek the density of $x(T), y(T)$ in the $T$-forward measure. As we focus on short times to maturity, it suffices to replace $y(T)$ with a deterministic approximation $\bar{y}(T)$ from Proposition 13.1.4. The local volatility (13.9) and mean reversion $\varkappa(t)$ are, as a rule, piecewise constant in time; thus, in the small- $T$ regime, they can be assumed to be independent of $t$. It is also safe to ignore the dependence of the volatility function on $y$. With this in mind, we define

$$
\begin{align*}
\varkappa & =\varkappa(0)  \tag{13.48}\\
v(x) & =\sigma_{r}(0, x, 0) / \sigma_{r}^{0}(0)
\end{align*}
$$
^eq13-48

The following result holds.

Theorem 13.1.14. Let us define $\pi(x)$ by the ODE

$$
v(x)\left(\pi^{\prime}(x)\right)^{2}+2 \varkappa G_{2}(T) \pi(x)(x / \pi(x))^{\prime}-1=0, \quad x \in \mathbb{R}, \quad \pi(0)=0
$$

where $G_{2}(T)$ is defined by (see (13.3) for the definition of $h(t)$ )

$$
\begin{equation*}
G_{2}(T)=\int_{0}^{T} h(s)^{2} d s \tag{13.49}
\end{equation*}
$$
^eq13-49

and the prime denotes differentiation with respect to $x$. Set

$$
\varpi(x)=x / \pi(x)
$$

and let $\Psi(T, x)$ be the $C D F$ of $x(T)$ in the $T$-forward measure. Then

$$
\begin{aligned}
\Psi(T, x) \approx \Phi\left(\frac{x}{\sigma_{r}^{0}(0) \varpi(x) \sqrt{G_{2}(T)}}\right) & \\
& +\sigma_{r}^{0}(0) \sqrt{G_{2}(T)} \varpi^{\prime}(x) \phi\left(-\frac{x}{\sigma_{r}^{0}(0) \varpi(x) \sqrt{G_{2}(T)}}\right)
\end{aligned}
$$

where $\Phi(z), \phi(z)$ are the standard Gaussian CDF and PDF, respectively.
Proof. While lengthy and somewhat technical, the proof is instructive as it shows a general approach to deriving short-time densities for local volatility term structure models. Full details are shown in Appendix 13.A of this chapter.
Remark 13.1.15. If $\sigma_{r}(0, x, 0)=\mathrm{const}$ (pure Gaussian case), then $v(x)=1$, and $\pi(x)=x$ is a solution to the ODE. Therefore, $\varpi(x)=1$, and we recover the Gaussian CDF of $x(T)$ as expected. The function $\varpi(x)=x / \pi(x)$ measures the deviation of the model from the Gaussian case; $\sigma_{r}^{0}(0) \varpi(x)$ can be thought of as an "effective term volatility" at $x$.

Equipped with Theorem 13.1.14, we can recover an approximation for the density

$$
\psi(T, x)=\frac{\partial \Psi}{\partial x}(T, x)
$$

which allows us to value (short-dated) derivative contracts by numerical integration. Specifically, the value of a contract with a payoff $V(x(T), y(T))$ at time $T$ in the quasi-Gaussian local volatility model is approximately equal to

$$
V_{0} \approx P(0, T) \int_{-\infty}^{\infty} V(x, \bar{y}(T)) \psi(T, x) d x
$$

where $\bar{y}(T)$ is given in Proposition 13.1.4. Again, the primary use of the results in Theorem 13.1.14 is to improve on the finite difference method, using the results of Sections 2.8.2 and 2.8.3.

### 13.1.9.2 Finite Difference Methods

We consider the model in its general local volatility form (13.10),

$$
\begin{aligned}
& d x(t)=(y(t)-\varkappa(t) x(t)) d t+\sigma_{r}(t, x(t), y(t)) d W(t) \\
& d y(t)=\left(\sigma_{r}(t, x(t), y(t))^{2}-2 \varkappa(t) y(t)\right) d t
\end{aligned}
$$

Using methods from Chapter 2, a PDE for the value of a security as a function of $x, y$ can be easily derived from these SDEs. We find it more convenient, however, to transform the variables first to replace $y(t)$ with a locally deterministic process $u(t)$ that is drift-free on average.

Recall the definition (13.12), and define the deterministic function $\bar{y}(t)$ by (as in Proposition 13.1.4)

$$
d \bar{y}(t)=\left(\sigma_{r}^{0}(t)^{2}-2 \varkappa(t) \bar{y}(t)\right) d t, \quad \bar{y}(0)=0
$$

or

$$
\bar{y}(t)=h(t)^{2} \int_{0}^{t} \sigma_{r}^{0}(s)^{2} h(s)^{-2} d s
$$

We define a new, normalized auxiliary variable $u(t)$ by

$$
\begin{equation*}
u(t)=y(t)-\bar{y}(t) \tag{13.50}
\end{equation*}
$$
^eq13-50

so that the state process $(x(t), u(t))$ satisfies

$$
\begin{align*}
& d x(t)=(u(t)+\bar{y}(t)-\varkappa(t) x(t)) d t+\sigma_{r}(t, x(t), u(t)+\bar{y}(t)) d W(t)  \tag{13.51}\\
& d u(t)=\left(\left(\sigma_{r}(t, x(t), u(t)+\bar{y}(t))^{2}-\sigma_{r}^{0}(t)^{2}\right)-2 \varkappa(t) u(t)\right) d t \tag{13.52}
\end{align*}
$$
^eq13-51
^eq13-52

subject to $x(0)=u(0)=0$. Values of zero-coupon discount bonds can easily be re-expressed in terms of the new state variables,
$P(t, T, x, u)=\frac{P(0, T)}{P(0, t)} \exp \left(-G(t, T) x-\frac{1}{2} G(t, T)^{2} u-\frac{1}{2} G(t, T)^{2} \bar{y}(t)\right)$.
The new parameterization of the qG model reduces nicely to the (Gaussian) case of deterministic volatility, in the following sense. If $\sigma_{r}(t, x, y)$ is independent of $x$ and $y$, then

$$
\begin{equation*}
\sigma_{r}(t, x, y) \equiv \sigma_{r}^{0}(t) \tag{13.53}
\end{equation*}
$$
^eq13-53

and the SDE for the state variable $u(t)$ becomes

$$
d u(t)=-2 \varkappa(t) u(t) d t, \quad u(0)=0
$$

with the unique solution

$$
\begin{equation*}
u(t) \equiv 0 \tag{13.54}
\end{equation*}
$$
^eq13-54

Thus, in the pure Gaussian case, the system of SDEs (13.51) reduces to a single SDE, in line with the way a one-factor Gaussian model was developed in Section 10.1.2.2.

Aesthetic reasons aside, the change of variable from $y(t)$ to $u(t)$ improves the numerical properties of a discretization of the PDE. Specifically, the variable $y$ (or $u$ ) does not have a diffusion term so, at least in the $y$ direction, the PDE is convection-dominated, in the sense described in Section 2.6. Removing most of the drift outside of the time-stepping scheme alleviates some of the numerical issues associated with such $\mathrm{PDEs}^{6}$.

The PDE associated with the dynamics (13.51) is derived in the standard way. Let $V(t, x, u)$ be the value, at time $t$, of a derivative with a payoff $V(x(T), u(T))$ at time $T$, given that $x(t)=x, u(t)=u$,

$$
V(t, x, u)=\mathrm{E}\left(e^{-\int_{t}^{T} r(s) d s} V(x(T), u(T)) \mid x(t)=x, u(t)=u\right)
$$

Then the function $V(t, x, u)$ satisfies the PDE

$$
\begin{align*}
\frac{\partial V}{\partial t}(t, x, u)+(\mathcal{L} V)(t, x, u) & =(f(0, t)+x) V(t, x, u), \quad 0 \leq t<T  \tag{13.55}\\
V(T, x, u) & =V(x, u)
\end{align*}
$$
^eq13-55

where

$$
\mathcal{L}=\mathcal{L}_{x}+\mathcal{L}_{u},
$$

and

$$
\begin{aligned}
\mathcal{L}_{x} & =(u+\bar{y}(t)-\varkappa(t) x) \frac{\partial}{\partial x}+\frac{1}{2} \sigma_{r}(t, x, u+\bar{y}(t))^{2} \frac{\partial^{2}}{\partial x^{2}} \\
\mathcal{L}_{u} & =\left(\left(\sigma_{r}(t, x, u+\bar{y}(t))^{2}-\sigma_{r}^{0}(t)^{2}\right)-2 \varkappa(t) u\right) \frac{\partial}{\partial u}
\end{aligned}
$$

The PDE (13.55) has two space dimensions and no mixed derivatives, and can be solved numerically by the Douglas-Rachford ADI method outlined in Section 2.10. The fact that $u(t)$ has no diffusion term does not complicate the situation much; even without upwinding (which we nevertheless recommend), we have observed no noticeable deterioration in the stability or accuracy of the scheme. There is some evidence that a 5-point discretization in the $u$ direction may improve precision slightly, as may semi-Lagrangian schemes (see e.g. Chen and Forsyth [2007]); as standard methods produce adequate results, we consider such improvements optional and leave them to the reader to explore.

[^44]The simple nature of the process for $u(t)$ typically allows one to use a rather coarse discretization in the $u$ direction. For instance, a typical setting might involve $n_{t}=100$ time steps, $n_{x}=150$ steps in $x$ direction, and $n_{u}=10$ steps in the $u$ direction. With these settings, most instruments are priced to basis point precision. It should be noted that $n_{u}$ can be chosen to reflect the degree to which the model deviates from the pure Gaussian case. To elaborate, consider a volatility term of the form (13.29). If $b=0$ the model is Gaussian, in which case only one discretization point ( $u_{0}=0$ ) is required, as is clear from (13.54). As $b$ is increased, the model becomes increasingly non-Gaussian, and an ever larger number of points in $u$ direction are required to maintain adequate precision. In other words, a practical scheme would set $u$ as an increasing function of the skew parameter $b$.

The choice of the domain for $(x(t), u(t))$ follows standard prescriptions from Chapter 2. Boundaries in the $x$ dimension are most easily obtained under the Gaussian approximation to the short rate state dynamics

$$
\begin{equation*}
d x(t) \approx(\bar{y}(t)-\varkappa(t) x(t)) d t+\sigma_{r}^{0}(t) d W(t), \quad x(0)=0 \tag{13.56}
\end{equation*}
$$
^eq13-56

using the formula (10.31) for the size of the grid while calculating $\mathrm{E}(x(T))$, $\operatorname{Var}(x(T))$ from (13.56) (see also footnote 7 in Chapter 10). For a slightly more refined approach, we note that with the linear local volatility (13.29) (or under a linear approximation in $x$ to the general volatility function $\left.\sigma_{r}(t, x, y)\right)$ the distribution of $x(T)$ is closer to a displaced log-normal than to a Gaussian; we can then set the boundaries by moment-matching a displaced log-normal variable to the distribution of $x(T)$ and using appropriate quantiles. We leave it to the reader to fill in missing details.

A slightly more interesting question is how to dimension the grid in the $u$ direction. As $u(t)$ does not have its own diffusion term, the randomness in $u(T)$ comes from the stochasticity of the drift in (13.52) which, to first order, is driven by $x(t)$. To extract this dependence, we apply a linear approximation to the volatility function $\sigma_{r}(t, x, y)$,

$$
\sigma_{r}(t, x, u+\bar{y}(t))^{2}-\sigma_{r}^{0}(t)^{2} \approx 2 \sigma_{r}^{0}(t)\left(\partial \sigma_{r}(t, 0,0) / \partial x\right) x
$$

in the SDE (13.52), which allows us to integrate (13.52) in an approximate sense,

$$
u(T) \approx 2 h(T)^{2} \int_{0}^{T} x(t) \sigma_{r}^{0}(t)\left(\partial \sigma_{r}(t, 0,0) / \partial x\right) h(t)^{-2} d t
$$

Then, under the Gaussian approximation (13.56) to the dynamics of $x(t)$, we see that $u(T)$ is also Gaussian with an easily calculated variance, which allows us to set the boundaries in the $u$ direction using an analog to the formula (10.31).

Once the grid boundaries are determined, appropriate spatial boundary conditions need to be specified. This is straightforward, with the ideas from

Section 10.1.5 easily transferable to the quasi-Gaussian model. For instance, a reasonable boundary specification is to assume linearity in $V$ as a function of $u$ at the $u$-boundaries (as in Section 2.2.2), while using the PDE itself (as in Section 10.1.5.2) to establish the boundary conditions at the $x$-boundaries.

### 13.1.9.3 Monte Carlo Simulation

Application of the Monte Carlo method to the quasi-Gaussian model is straightforward and can follow general guidelines from Chapter 3. As with the PDE method, there may be some accuracy gains associated with using the state variable $u(t)$ in (13.50) instead of $y(t)$.

Consider the problem of computing the value of a derivative with the payoff $V(x(T), u(T))$ at time $T$. Let

$$
0=t_{0}<t_{1}<\ldots<t_{N}=T, \quad \Delta_{n}=t_{n}-t_{n-1}
$$

be the discretization of the time domain. By applying a standard Euler discretization (see Section 3.2.3) to (13.51), the following stepping scheme is obtained,

$$
\begin{aligned}
\widehat{x}_{n} & =\widehat{x}_{n-1}+\left(\widehat{u}_{n-1}+\bar{y}\left(t_{n-1}\right)-\varkappa\left(t_{n-1}\right) \widehat{x}_{n-1}\right) \Delta_{n} \\
& +\sigma_{r}\left(t_{n-1}, \widehat{x}_{n-1}, \widehat{u}_{n-1}+\bar{y}\left(t_{n-1}\right)\right) Z_{n} \sqrt{\Delta_{n}} \\
\widehat{u}_{n} & =\widehat{u}_{n-1}+\Delta_{n} \\
& \times\left(\sigma_{r}\left(t_{n-1}, \widehat{x}_{n-1}, \widehat{u}_{n-1}+\bar{y}\left(t_{n-1}\right)\right)^{2}-\sigma_{r}^{0}\left(t_{n-1}\right)^{2}-2 \varkappa\left(t_{n-1}\right) \widehat{u}_{n-1}\right)
\end{aligned}
$$

where $\widehat{x}_{0}=\widehat{u}_{0}=0,\left\{Z_{n}\right\}_{n=1}^{N}$ is a collection of i.i.d. standard Gaussian random variables, and $\left\{\widehat{x}_{n}, \widehat{u}_{n}\right\}_{n=0}^{N}$ is an approximation to $\left\{x\left(t_{n}\right), u\left(t_{n}\right)\right\}_{n=0}^{N}$. Of course, more advanced discretization schemes are possible, as explained in Chapter 3. Some of the ideas of Section 10.1.6 are also applicable here, including the observation that the model can be simulated under either the terminal measure or the spot measure to avoid the bias involved in timediscretizing the continuously compounded money market account $e^{\int_{0}^{r} r(s) d s}$ in the risk-neutral measure.

### 13.1.9.4 Single-State Approximations

The extra state variable, and the resultant requirement of a two-dimensional PDE scheme (see (13.55)) for an essentially one-factor model, is the price one has to pay for the flexibility of volatility specification in the qG model. This price is relatively modest in practice, but does make the model slightly slower than a classical one-factor short rate model, and also makes it somewhat more challenging to use as a building block for more complicated models, such as equity or FX-linked interest rate hybrid models. In this section we
briefly outline a few ideas for reducing the dimensionality of the model to one state only.

A very simple idea that can be traced to Hagan and Woodward [1999a] (in the multi-dimensional setting) is to force the $y(t)$ variable to be deterministic in the SDE for $x(t)$ :

$$
d x(t)=(y(t)-\varkappa(t) x(t)) d t+\sigma_{r}(t, x(t)) d W(t), \quad x(0)=0
$$

where now $y(t)$ is deterministic. Then, using generic machinery of Section 11.3, we can fit $y(t)$ to the initial yield curve via forward induction. Of course, tractability of bond reconstruction formulas is then largely lost, although much of the intuition behind the model is retained and we still maintain separate control over the at-the-money volatility structure (via mean reversion) and the volatility smile (via the local volatility function). Also, the bond reconstruction formulas from Proposition 13.1.1 can be considered as approximations in the new model, and potentially could be used to speed up volatility calibration.

A straight deterministic approximation is a rather blunt tool and would most likely not deliver the level of accuracy we require. A more refined approach for replacing the stochastic variable $y(t)$ in the qG model involves using its projection on the variable $x(t)$, as proposed in Kramin [2008]. To develop this idea in a bit more detail, consider a qG model with local volatility (13.10). Focusing first on calculating the following one-dimensional risk-neutral expectation

$$
\mathrm{E}(V(x(T)))
$$

with the payoff $V(x)$ is a function of $x(T)$ only, the ideas behind Markovian projection (Gyöngy's theorem, see Theorem A.1.1) give us the following exact result.

Proposition 13.1.16. The undiscounted expected value of a payoff $V(x(T))$ in the model (13.10) is equal to

$$
\begin{equation*}
\mathrm{E}(V(x(T)))=\mathrm{E}(V(\widetilde{x}(T))), \tag{13.57}
\end{equation*}
$$
^eq13-57

where the process $\widetilde{x}(t)$ satisfies

$$
\begin{equation*}
d \widetilde{x}(t)=(\widetilde{y}(t, \widetilde{x}(t))-\varkappa(t) \widetilde{x}(t)) d t+\widetilde{\sigma}_{r}(t, \widetilde{x}(t)) d W(t) \tag{13.58}
\end{equation*}
$$
^eq13-58

with

$$
\begin{equation*}
\widetilde{y}(t, x)=\mathrm{E}(y(t) \mid x(t)=x), \quad \tilde{\sigma}_{r}(t, x)^{2}=\mathrm{E}\left(\sigma_{r}(t, x(t), y(t))^{2} \mid x(t)=x\right) \tag{13.59}
\end{equation*}
$$
^eq13-59

The equality (13.57) does not hold for the more realistic, and useful, case of calculating discounted expected values of a payoffs that depend on both $x$ and $y$. Nor is there much theoretical justification for using this
projection when calculating expected values of payoffs that depend on values of the state variables at multiple times. Nevertheless, there is some empirical evidence that the approximations work reasonably well in practice. With that in mind, let us define a generic one-state approximate quasi-Gaussian local volatility model by

$$
\begin{align*}
d x(t) & =(\widetilde{y}(t, x(t))-\varkappa(t) x(t)) d t+\widetilde{\sigma}_{r}(t, x(t)) d W(t)  \tag{13.60}\\
P(t, T) & =\frac{P(0, T)}{P(0, t)} \exp \left(-G(t, T) x(t)-\frac{1}{2} G(t, T)^{2} \widetilde{y}(t, x(t))\right) \\
r(t) & =f(0, t)+x(t)
\end{align*}
$$
^eq13-60

with $\widetilde{y}(t, x), \widetilde{\sigma}_{r}(t, x)$ given by (13.59). Then, the value $V(t, x)$ of a given security at time $t$ in state $x$ satisfies the following one-dimensional PDE

$$
\begin{equation*}
\frac{\partial V}{\partial t}+(\widetilde{y}(t, x)-\varkappa(t) x) \frac{\partial V}{\partial x}+\frac{1}{2} \widetilde{\sigma}_{r}(t, x)^{2} \frac{\partial^{2} V}{\partial x^{2}}=(f(0, t)+x) V \tag{13.61}
\end{equation*}
$$
^eq13-61

Assuming that we can evaluate all terms efficiently, solving the onedimensional PDE (13.61) is typically quicker than solving the PDE (13.55) for the real model. There is little, if any, benefit in applying the approximation to the Monte Carlo method.

Needless to say, the PDE (13.61) should only be considered for problems inside the domain of applicability of the approximation (13.60). In general, we would expect the approximation to work reasonably well for low to moderate volatilities and maturities (up to, say, 20 years), and deteriorate for longer maturities and/or for large volatilities. Kramin [2008] reports good results across a wide maturity spectrum.

To effectively use (13.60), we need to compute/approximate $\widetilde{y}(t, x)$ and $\tilde{\sigma}_{r}(t, x)$. Typically the volatility term $\sigma_{r}(t, x, y)$ either does not depend on $y$ at all (see e.g. (13.29)), or depends on $y$ in a close-to-linear fashion, due to the low variance of $y(t)$ compared to $x(t)$. In both cases, the following simple approximation

$$
\tilde{\sigma}_{r}(t, x)=\sigma_{r}(t, x, \tilde{y}(t, x))
$$

appears to be justified.
To calculate $\widetilde{y}(t, x)$, we recall the definition of $y(t)$,

$$
y(t)=h(t)^{2} \int_{0}^{t} \sigma_{r}(s, x(s), y(s))^{2} h(s)^{-2} d s
$$

Conditioning on $x(t)$ and replacing $y(s)$ with $\bar{y}(s)$ in the argument of $\sigma_{r}$, where $\bar{y}(\cdot)$ is defined in (13.23), we obtain

$$
\widetilde{y}(t, x) \approx h(t)^{2} \int_{0}^{t} \mathrm{E}\left(\sigma_{r}(s, x(s), \bar{y}(s))^{2} \mid x(t)=x\right) h(s)^{-2} d s
$$

Invoking approximate linearity of $\sigma_{r}(s, x, y)^{2}$ in $x$ we obtain

$$
\widetilde{y}(t, x) \approx h(t)^{2} \int_{0}^{t} \sigma_{r}(s, \mathrm{E}(x(s) \mid x(t)=x), \bar{y}(s))^{2} h(s)^{-2} d s
$$

Under the Gaussian approximation

$$
\mathrm{E}(x(s) \mid x(t)=x) \approx \frac{\operatorname{Var}(x(s))}{\operatorname{Var}(x(t))} x \approx \frac{h(s)^{2} \int_{0}^{s} \sigma_{r}^{0}(u)^{2} h(u)^{-2} d u}{h(t)^{2} \int_{0}^{t} \sigma_{r}^{0}(u)^{2} h(u)^{-2} d u} x=\frac{\bar{y}(s)}{\bar{y}(t)} x
$$

so we obtain

$$
\begin{equation*}
\widetilde{y}(t, x) \approx h(t)^{2} \int_{0}^{t} \sigma_{r}(s,(\bar{y}(s) / \bar{y}(t)) x, \bar{y}(s))^{2} h(s)^{-2} d s \tag{13.62}
\end{equation*}
$$
^eq13-62

A direct application of (13.62) is rather costly: with $n_{t}$ discretized points in $t$ direction and $n_{x}$ points in $x$ direction, the cost of computing $\widetilde{y}(\cdot, \cdot)$ for all $t, x$ on the grid is $O\left(n_{t}^{2} n_{x}\right)$, i.e. higher than for solving the PDE (13.61) itself. One remedy for this issue is to approximate $\sigma_{r}(s, x, \bar{y}(s))^{2}$ by a firstor second-order polynomial in $x$ for each $s$, to obtain a polynomial approximation to $\widetilde{y}(s, x)$ with the coefficients computed at $O\left(n_{t}\right)$ cost. Alternatively, we can derive a recursive update equation for $\widetilde{y}(t, x)$ with the additional advantage that it does not rely on approximate linearity of $\sigma_{r}(s, x, y)^{2}$ in $x$ (unlike (13.62)), a condition that, although generally desirable as we pointed out before, is not necessarily satisfied in all applications. We recall the equation (13.11) satisfied by $y(t)$, and discretize it for the time step $\left[t_{n}, t_{n+1}\right]$. Using short-hand notations $x_{n}=x\left(t_{n}\right)$, etc., we obtain

$$
y_{n+1}=y_{n}+\left(\sigma_{r}\left(t_{n+1}, x_{n+1}, \bar{y}_{n+1}\right)^{2}-2 \varkappa(t) y_{n}\right) \Delta_{n}, \quad, \Delta_{n}=t_{n+1}-t_{n},
$$

where $\sigma_{r}(t, x, \bar{y})^{2}$ is evaluated at the right point of the interval for reasons that will be clear momentarily. Next, conditioning on $x_{n}, x_{n+1}$, we obtain

$$
\begin{aligned}
& \mathrm{E}\left(y_{n+1} \mid x_{n}, x_{n+1}\right)=\mathrm{E}\left(y_{n} \mid x_{n}, x_{n+1}\right) \\
& +\left(\sigma_{r}\left(t_{n+1}, x_{n+1}, \bar{y}_{n+1}\right)^{2}-2 \varkappa(t) \mathrm{E}\left(y_{n} \mid x_{n}, x_{n+1}\right)\right) \Delta_{n}
\end{aligned}
$$

By the Markov property

$$
\mathrm{E}\left(y_{n} \mid x_{n}, x_{n+1}\right)=\mathrm{E}\left(y_{n} \mid x_{n}\right),
$$

so

$$
\begin{aligned}
\mathrm{E}\left(y_{n+1} \mid x_{n}, x_{n+1}\right)= & \mathrm{E}\left(y_{n} \mid x_{n}\right) \\
& +\left(\sigma_{r}\left(t_{n+1}, x_{n+1}, \bar{y}_{n+1}\right)^{2}-2 \varkappa(t) \mathrm{E}\left(y_{n} \mid x_{n}\right)\right) \Delta_{n},
\end{aligned}
$$

which gives us a way to obtain $\mathrm{E}\left(y_{n+1} \mid x_{n}, x_{n+1}\right)$ from $\mathrm{E}\left(y_{n} \mid x_{n}\right)$. To get at the quantity that we want, namely $\mathrm{E}\left(y_{n+1} \mid x_{n+1}\right)$, we average over $x_{n}$,

$$
\begin{aligned}
& \mathrm{E}\left(y_{n+1} \mid x_{n+1}\right) \\
& \quad=\int \mathrm{E}\left(y_{n+1} \mid x_{n}=x, x_{n+1}\right) \mathrm{Q}\left(x_{n} \in d x \mid x_{n+1}\right) \\
& \quad=\int\left(\mathrm{E}\left(y_{n} \mid x_{n}=x\right)\right. \\
& \left.\quad+\left(\sigma_{r}\left(t_{n+1}, x_{n+1}, \bar{y}_{n+1}\right)^{2}-2 \varkappa(t) \mathrm{E}\left(y_{n} \mid x_{n}=x\right)\right) \Delta_{n}\right) \\
& \quad \times \mathrm{Q}\left(x_{n} \in d x \mid x_{n+1}\right)
\end{aligned}
$$

and, after rearranging some terms, we obtain a recursive formula for $\mathrm{E}\left(y_{n+1} \mid x_{n+1}\right)$,

$$
\begin{aligned}
\mathrm{E}\left(y_{n+1} \mid x_{n+1}\right)=\left(1-2 \varkappa(t) \Delta_{n}\right) \int \mathrm{E}\left(y_{n} \mid\right. & \left.x_{n}=x\right) \mathrm{Q}\left(x_{n} \in d x \mid x_{n+1}\right) \\
& +\sigma_{r}\left(t_{n+1}, x_{n+1}, \bar{y}_{n+1}\right)^{2} \Delta_{n}
\end{aligned}
$$

For small $\Delta_{n}$, the density

$$
\mathrm{Q}\left(x_{n} \in d x \mid x_{n+1}\right)
$$

is approximately Gaussian, and the required integral can be quickly computed numerically with just a few terms, giving us an algorithm of numerical complexity $O\left(n_{t} n_{x}\right)$.

Finally, we point out that the model (13.60) could be made exactly arbitrage-free by introducing a time-dependent deterministic component in its drift that is fit numerically to the initial yield curve, in line with the discussion at the beginning of this section - but of course at the cost of losing analytical tractability.

### 13.2 One-Factor Quasi-Gaussian Model with Stochastic Volatility

The most general one-factor quasi-Gaussian model specification allows for the short rate volatility to be a stochastic process, see (13.2). While so far we only considered the case of deterministic dependence of the volatility on state variables of the model, we now proceed to generalize the setup to include stochastic volatility.

### 13.2.1 Definition

Introduction of a stochastic variance process (see Chapter 8) in the specification of $g(\cdot)$ in (13.2) leads to a stochastic volatility quasi-Gaussian model. In particular, defining $z(t)$ to be the familiar CIR process,

$$
d z(t)=\theta\left(z_{0}-z(t)\right) d t+\eta(t) \sqrt{z(t)} d Z(t), \quad\langle d Z(t), d W(t)\rangle=0
$$

we obtain a stochastic volatility qG model by specifying the volatility structure of the form

$$
\begin{equation*}
g(t, \omega)=\sqrt{z(t)} g(t, x(t), y(t)) \tag{13.63}
\end{equation*}
$$
^eq13-63

where $g(t, x, y)$ is a function of $t, x, y$ only. With the standard definition

$$
\sigma_{r}(t, x, y)=g(t, x, y) h(t)
$$

the model is defined by the collection of SDEs

$$
\begin{align*}
& d x(t)=(y(t)-\varkappa(t) x(t)) d t+\sqrt{z(t)} \sigma_{r}(t, x(t), y(t)) d W(t)  \tag{13.64}\\
& d y(t)=\left(z(t) \sigma_{r}(t, x(t), y(t))^{2}-2 \varkappa(t) y(t)\right) d t \\
& d z(t)=\theta\left(z_{0}-z(t)\right) d t+\eta(t) \sqrt{z(t)} d Z(t)
\end{align*}
$$
^eq13-64

subject to

$$
x(0)=y(0)=0, \quad z(0)=z_{0}=1, \quad\langle d Z(t), d W(t)\rangle=0
$$

When specifying the local volatility function, it was natural to use piecewise constant functions for various parameters (see (13.39)), and we do the same with the volatility of variance,

$$
\begin{equation*}
\eta(t)=\sum_{n=1}^{N-1} \eta_{n} 1_{\left\{t \in\left(T_{n-1}, T_{n}\right]\right\}} . \tag{13.65}
\end{equation*}
$$
^eq13-65

The bond reconstitution formulas in the model (13.64) are the same as for the local volatility case; as follows from Proposition 13.1.1, they are the same for any quasi-Gaussian model. In particular, the zero-coupon discount bond formulas do not depend on the stochastic volatility process $z(t)$, and thus the model is a "true" stochastic volatility model, i.e. its stochastic volatility is unspanned and cannot be hedged by discount bonds. We remind the reader of the discussion of this topic in Section 11.2.3 and note that the model (13.64) has the lowest possible number of state variables - three - for an unspanned stochastic volatility term structure model, see Collin-Dufresne and Goldstein [2002b].

We note in passing that the assumption of zero correlation $\langle d Z(t), d W(t)\rangle=0$ is a technical restriction helpful for developing efficient calibration formulas. It does not, however, restrict the range of available volatility smiles in the model, as the local volatility term can be used to control the slope of the smile. See also our discussion in Section 13.2.5.

### 13.2.2 Swap Rate Dynamics

Many results obtained in the local volatility case extend naturally to incorporate stochastic volatility, including Proposition 13.1.2, Lemma 13.1.3 and Proposition 13.1.8. The following analog to Corollary 13.1.9 is particularly useful for calibration.

Proposition 13.2.1. Under the assumption of linear local short rate volatility (13.29), the dynamics of a swap rate $S(t)$ in the stochastic volatility quasi-Gaussian model (13.64) are given approximately by

$$
\begin{aligned}
d S(t) & =\sqrt{z(t)} \lambda_{S}(t)\left(b_{S}(t) S(t)+\left(1-b_{S}(t)\right) S(0)\right) d W^{A}(t) \\
d z(t) & =\theta\left(z_{0}-z(t)\right) d t+\eta(t) \sqrt{z(t)} d Z(t) \\
z(0) & =z_{0}=1, \quad\left\langle d Z(t), d W^{A}(t)\right\rangle=0
\end{aligned}
$$

where $\lambda_{S}(t)$ and $b_{S}(t)$ are given by (13.34)-(13.35).
The dynamics in Proposition 13.2.1 are easily recognized to be those of a time-dependent SV model, see Chapter 9. As time averaging methods are available for stochastic volatility models (see Section 9.3), the following proposition, an analog to Proposition 13.1.10, should not come as a surprise. See also Theorem 9.3.1, Corollary 9.3.5, and Theorem 9.3.6.

Proposition 13.2.2. In the setting of the stochastic volatility quasiGaussian model (13.64) with the linear local volatility (13.29), consider a $T$-maturity swaption on a swap rate $S(T)$. For the purpose of European option pricing, the dynamics of $S(t)$ in its annuity measure can be approximated by the following time-homogeneous stochastic volatility model,

$$
\begin{aligned}
d S(t) & =\sqrt{z(t)} \bar{\lambda}_{S}\left(\bar{b}_{S} S(t)+\left(1-\bar{b}_{S}(t)\right) S(0)\right) d W^{A}(t) \\
d z(t) & =\theta\left(z_{0}-z(t)\right) d t+\bar{\eta}_{S} \sqrt{z(t)} d Z(t)
\end{aligned}
$$

where

- The effective volatility of variance $\bar{\eta}_{S}$ is given by

$$
\begin{equation*}
\bar{\eta}_{S}^{2}=\frac{\int_{0}^{T} \eta(t)^{2} \rho_{S}(t) d t}{\int_{0}^{T} \rho_{S}(t) d t} \tag{13.66}
\end{equation*}
$$
^eq13-66

with the weight function $\rho_{S}(t)$ given by

$$
\rho_{S}(r)=\int_{r}^{T} \int_{s}^{T} \lambda_{S}(t)^{2} \lambda_{S}(s)^{2} e^{-\theta(t-s)} e^{-2 \theta(s-r)} d t d s
$$

- The effective skew $\bar{b}_{S}$ is given by

$$
\begin{equation*}
\bar{b}_{S}=\int_{0}^{T} b_{S}(t) w_{S}(t) d t \tag{13.67}
\end{equation*}
$$
^eq13-67

with the weight function $w_{S}(t)$ given by

$$
\begin{aligned}
& w_{S}(t)=\frac{v_{S}(t)^{2} \lambda_{S}(t)^{2}}{\int_{0}^{T} v_{S}(u)^{2} \lambda_{S}(u)^{2} d u} \\
& v_{S}(t)^{2}=z_{0}^{2} \int_{0}^{t} \lambda_{S}(s)^{2} d s+z_{0} e^{-\theta t} \int_{0}^{t} \lambda_{S}(s)^{2} e^{-\theta s} \int_{0}^{s} \eta(u)^{2} e^{2 \theta u} d u d s
\end{aligned}
$$

- The effective volatility $\bar{\lambda}_{S}$ is given by the solution to the equation

$$
\begin{equation*}
\Psi_{\bar{z}}\left(\frac{h^{\prime \prime}\left(\zeta_{S}\right)}{h^{\prime}\left(\zeta_{S}\right)} \bar{\lambda}_{S}, 0 ; T\right)=\Psi_{\overline{z \lambda^{2}}}\left(\frac{h^{\prime \prime}\left(\zeta_{S}\right)}{h^{\prime}\left(\zeta_{S}\right)}, 0 ; T\right) \tag{13.68}
\end{equation*}
$$
^eq13-68

where (see Theorem 9.3.1)

$$
\begin{aligned}
\zeta_{S} & =z_{0} \int_{0}^{T} \lambda_{S}(t)^{2} d t \\
\Psi_{\overline{z \lambda^{2}}}(v, 0 ; T) & =\mathrm{E}\left(\exp \left(v \int_{0}^{T} \lambda_{S}(t)^{2} z(t) d t\right)\right) \\
\Psi_{\bar{z}}(v, 0 ; T) & =\mathrm{E}\left(\exp \left(v \int_{0}^{T} z(t) d t\right)\right) \\
h(x) & =\frac{S(0)}{\bar{b}_{S}}\left(2 \Phi\left(\bar{b}_{S} \sqrt{x} / 2\right)-1\right)
\end{aligned}
$$

The functions $\lambda_{S}(t), b_{S}(t)$ are given by (13.34)-(13.35).

### 13.2.3 Volatility Calibration

With a swaption strip $\left\{S_{n}(\cdot)\right\}_{n=1}^{N-1}$ given as in Sections 13.1.6 and 13.1.7, the volatility calibration algorithm can proceed along the same principles as in Section 13.1.7, where swap rate distributions in the stochastic volatility qG model can be found from the constant-parameter displaced SV SDEs in Proposition 13.2.2. As before, we assume that a collection of market parameters ( $\widehat{\lambda}_{S_{n}}, \widehat{b}_{S_{n}}, \widehat{\eta}_{S_{n}}$ ), $n=1, \ldots, N-1$, is given; in practice, these parameters may be obtained by fitting a vanilla SV model to swaptions of a given expiry/tenor across strikes, a procedure described in more detail in Section 16.1.4.

While it is easy to modify the algorithm of Section 13.1.7 to introduce one more variable to solve for $\left(\eta_{n}\right)$ in Step 6 for each $n$, the calibration algorithm
is typically more stable if we first solve for the volatility of variance function $\eta(t)$, i.e. find $\left\{\eta_{n}^{*}\right\}$ for all $n$, and then follow the algorithm from Section 13.1.7 to solve recursively in $n$ for ( $\lambda_{n}^{*}, b_{n}^{*}$ ), using slightly modified formulas from Proposition 13.1.10 as given in Proposition 13.2.2. For completeness, we repeat the algorithm with the necessary modifications.

1. Set $\left(\lambda_{n}, b_{n}\right), n=1, \ldots, N-1$, to some reasonable starting values, e.g. set $\lambda_{n}$ 's to (properly scaled) volatilities obtained by calibrating a pure Gaussian model as in Section 10.1.4, and $b_{n}=\widehat{b}_{S_{n}}$.
2. Solve for $\eta_{n}^{*}$ for $n=1, \ldots, N-1$, using (13.66). For weights $\rho_{S_{n}}(t)$, use $\lambda_{S_{n}}(t)$ as computed from the first guess for $\lambda_{n}$ 's (using (13.34)) obtained on Step 1.
3. Set $n=1$.
4. For given $n,\left(\lambda_{i}^{*}, b_{i}^{*}\right)$ are known for $i=1, \ldots, n-1$. Note we use a star to denote calibrated values of the model parameters.
5. Calculate $\bar{x}(t)$ (Lemma 13.1.6) and $\bar{y}(t)$ (Proposition 13.1.4) for $t \in \left[0, T_{n}\right]$ using $\left(\lambda_{i}^{*}, b_{i}^{*}\right), i=1, \ldots, n-1$, and the initial guess for $\left(\lambda_{n}, b_{n}\right)$ from Step 1. Note that $\bar{x}(t), \bar{y}(t)$ implicitly depend on $n$ as their definition depends on the swap rate/annuity measure used.
6. Calculate $\lambda_{S_{n}}(t), b_{S_{n}}(t)$ for $t \in\left[0, T_{n-1}\right]$ from $\left(\lambda_{i}^{*}, b_{i}^{*}\right), i=1, \ldots, n-1$, using (13.34)-(13.35).
7. Make another guess for $\left(\lambda_{n}, b_{n}\right)$.
8. Update $\lambda_{S_{n}}(t), b_{S_{n}}(t)$ for $t \in\left(T_{n-1}, T_{n}\right]$ from $\left(\lambda_{i}^{*}, b_{i}^{*}\right), i=1, \ldots, n-1$, using (13.34)-(13.35).
9. Calculate $\bar{\lambda}_{S_{n}}, \bar{b}_{S_{n}}$ using Proposition 13.2.2.
10. Compare ( $\bar{\lambda}_{S_{n}}, \bar{b}_{S_{n}}$ ) to ( $\widehat{\lambda}_{S_{n}}, \widehat{b}_{S_{n}}$ ). If not equal within given tolerance, go to Step 7. Otherwise, proceed to Step 11.
11. As we have reached acceptable convergence between $\left(\bar{\lambda}_{S_{n}}, \bar{b}_{S_{n}}\right)$ and $\left(\hat{\lambda}_{S_{n}}, \hat{b}_{S_{n}}\right)$, set the calibrated model parameter values to the latest trial values, $\left(\lambda_{n}^{*}, b_{n}^{*}\right)=\left(\lambda_{n}, b_{n}\right)$.
12. Update $n \rightarrow n+1$. If $n \leq N-1$ go to Step 4. Otherwise, conclude.

As in Section 13.1.7, Step 5 can be performed inside of the calibration loop with some positive impact on the quality of calibration at a (moderate) cost of extra complexity.

### 13.2.4 Mean Reversion Calibration

The stochastic variance scaling $z(t)$ has some impact on the inter-temporal correlations or volatility ratios of swaptions of different maturities, an effect we discuss further in Section 20.2.4. Still, as mean reversion calibration is not meant to be overly precise, we can continue to use formulas developed in Section 13.1.8 that, in the context of a stochastic volatility qG model, imply that we roundly ignore any such effects.

While on the topic of mean reversions, let us not forget the parameter $\theta$, the mean reversion of the variance process. The parameter $\theta$ fundamentally determines the decay of volatility smile curvature as a function of option maturity $T$; a good setting of this parameter will help keep the parameter $\eta(t)$ from depending excessively on calendar time $t$. In general, we inherit the same $\theta$ as used for the vanilla SV model calibration to European swaptions. a subject we discuss in depth in Section 16.1.4.

### 13.2.5 Non-Zero Correlation

A question sometimes arises whether it is too restrictive to assume (as we do) that the correlation between the Brownian motions driving the curve factor and the stochastic volatility is zero. Empirical evidence from all major fixed income markets generally suggests that correlations between interest rates and their (short-dated) volatilities are small; see e.g. the analysis in Chen and Scott [2001]. Moreover, the assumption of zero correlation has little impact on the range of volatility smiles that the model can produce, as the skew term in the local volatility can produce the necessary tilting of the smile. Nor does it affect hedging implications of the model as long as minimum variance hedging is employed; see our discussion in Section 8.9. In our view, the zero-correlation constraint has little consequence in practice, but brings substantial technical benefits, particularly the ability to shift pricing numeraires without affecting the form of the stochastic variance process. If desired, non-zero correlation can still be accommodated, as most numerical schemes - including averaging formulas - are easy to adapt, as we briefly explain in Remark 9.3.7. One should be mindful, however, of the fact that under non-zero correlation, measure changes introduce a ratedependent term in the drift of the stochastic variance (see Proposition 8.3.9), which requires additional considerations when deriving approximations for European swaptions, say. If interested, the reader can attack the problem along the same lines as in Chapter 15, where the case of non-zero correlation is considered in a context of a different model (the Libor market model).

### 13.2.6 PDE and Monte Carlo Methods

The PDE method of Section 13.1.9.2 can be extended to cover the stochastic volatility case, using the techniques of Section 9.4 for the stochastic volatility part. The resulting PDE will involve three spatial dimensions, which can be handled by the Craig-Sneyd scheme. The same is true for the Monte Carlo method, where a combination of ideas from Sections 13.1.9.3 and 9.5 cover most practical issues.

### 13.3 Multi-Factor Quasi-Gaussian Model

### 13.3.1 General Multi-Factor Model

Multi-factor quasi-Gaussian models combine the flexibility of volatility specification of multi-factor models (see Chapter 12) with the ability to generate volatility skews and smiles. Practical multi-factor quasi-Gaussian models are relatively new (see Andreasen [2005]), but could provide a compelling alternative to the Libor market model in Chapter $14^{7}$, the current de-facto market standard for multi-factor models.

Following the steps that lead to the one-factor quasi-Gaussian model, the multi-factor qG model is obtained by imposing a separability condition on the volatility structure of a multi-factor HJM model. Specifically, let us consider the forward rate process

$$
\begin{equation*}
d f(t, T)=\sigma_{f}(t, T, \omega)^{\top}\left(\left(\int_{t}^{T} \sigma_{f}(t, u, \omega) d u\right) d t+d W(t)\right) \tag{13.69}
\end{equation*}
$$
^eq13-69

where $\sigma_{f}(t, T, \omega)$ is a $d$-dimensional stochastic process, and $W(t)$ a $d$ dimensional Brownian motion in the risk-neutral measure. Let us assume that $\sigma_{f}(t, T, \omega)$ is separable, in the sense that it can be written as

$$
\begin{equation*}
\sigma_{f}(t, T, \omega)=g(t, \omega) h(T) \tag{13.70}
\end{equation*}
$$
^eq13-70

where $g(t, \omega)$ is a $d \times d$ stochastic matrix-valued process, and $h(t)$ is a $d$-dimensional deterministic vector-valued function of time. Then we define

$$
H(t)=\operatorname{diag}(h(t))=\left(\begin{array}{cccc}
h_{1}(t) & 0 & \ddots & 0 \\
0 & h_{2}(t) & \ddots & \ddots \\
\ddots & \ddots & \ddots & 0 \\
0 & \ddots & 0 & h_{d}(t)
\end{array}\right)
$$

Let us assume further that $h_{i}(t) \neq 0, i=1, \ldots, d$, for all $t$, whereby $H(t)$ is then invertible, and so we can define a diagonal $d \times d$ matrix $\varkappa(t)$ by

$$
\begin{equation*}
\varkappa(t)=-\frac{d H(t)}{d t} H(t)^{-1} \tag{13.71}
\end{equation*}
$$
^eq13-71

(this is the same as in (12.7)-(12.8)). Moreover, let us define

$$
G(t, T)=\int_{t}^{T} H(u) H(t)^{-1} 1 d u, \quad \sigma_{r}(t, \omega)=g(t, \omega) H(t)
$$

where we use the notation $\mathbf{1}=(1,1, \ldots, 1)^{\top}$ from Section 12.1.1.1.

[^45]Proposition 13.3.1. Consider a general multi-factor quasi-Gaussian model, i.e. an HJM model (13.69) with the separable volatility condition (13.70). Define stochastic processes $x(t), y(t)$ by

$$
\begin{align*}
d x(t) & =(y(t) \mathbf{1}-\varkappa(t) x(t)) d t+\sigma_{r}(t, \omega)^{\top} d W(t)  \tag{13.72}\\
d y(t) & =\left(\sigma_{r}(t, \omega)^{\top} \sigma_{r}(t, \omega)-\varkappa(t) y(t)-y(t) \varkappa(t)\right) d t
\end{align*}
$$
^eq13-72

where $x(t) \in \mathbb{R}^{d}, y(t) \in \mathbb{R}^{d \times d}$, and $x(0)=0, y(0)=0$. Then, zero-coupon discount bonds are given by

$$
P(t, T)=P(t, T, x(t), y(t)),
$$

with

$$
P(t, T, x, y)=\frac{P(0, T)}{P(0, t)} \exp \left(-G(t, T)^{\top} x-\frac{1}{2} G(t, T)^{\top} y G(t, T)\right) .
$$

In addition, the instantaneous forward rates are given by

$$
\begin{equation*}
f(t, T)=f(0, T)+1^{\top} H(T) H(t)^{-1}(x(t)+y(t) G(t, T)) \tag{13.73}
\end{equation*}
$$
^eq13-73

with the short rate

$$
r(t)=f(0, t)+1^{\top} x(t) .
$$

Proof. Follows closely that of Proposition 12.1.2.

### 13.3.2 Local and Stochastic Volatility Parameterization

While a pure local volatility specification of the multi-factor qG model is certainly possible, for brevity let us proceed directly to a more general setting where we have both local and stochastic volatility. Following Section 13.2, we start by specifying a one-dimensional process $z(t)$ by

$$
\begin{equation*}
d z(t)=\theta\left(z_{0}-z(t)\right) d t+\eta(t) \sqrt{z(t)} d Z(t), \quad z(0)=z_{0}=1 \tag{13.74}
\end{equation*}
$$
^eq13-74

with $\langle d Z(t), d W(t)\rangle=0$. Inspired by the one-dimensional case, we would like to specify a model with the volatility structure of the type

$$
\begin{equation*}
\sigma_{r}(t, \omega)^{\top}=\sqrt{z(t)} \sigma_{x}(t, x(t), y(t))^{\top}, \tag{13.75}
\end{equation*}
$$
^eq13-75

where $\sigma_{x}(t, x, y)$ is a multi-dimensional local volatility function responsible for inducing the skews in volatility smiles of swaptions. However, it is not entirely obvious how to parametrize $\sigma_{x}(t, x, y)$ sensibly, as the volatility function is not only responsible for skews but also for the general volatility structure of the model, including volatilities and correlations of all the rates.

Fortunately, the ideas of Section 12.1.7 could be fruitfully recycled and extended here, as suggested by Andreasen [2005] (see also Cheyette [1991]).

Recalling the definition of benchmark rates from Section 12.1.7, we take $d$ benchmark tenors $\delta_{1}<\ldots<\delta_{d}$, and define $d$ "rolling" benchmark rates $f_{i}(t)=f\left(t, t+\delta_{i}\right), i=1, \ldots, d$. Ideally, it would be convenient if the qG model specification was such that the dynamics of the benchmark rates $f_{i}(t)$, $i=1, \ldots, d$, were of the familiar form

$$
\begin{equation*}
d f_{i}(t)=\sqrt{z(t)} \lambda_{i}^{f}(t)\left(\alpha_{i}^{f}(t)+b_{i}^{f}(t) f_{i}(t)\right) d U_{i}(t)+O(d t), \quad i=1, \ldots, d \tag{13.76}
\end{equation*}
$$
^eq13-76

where $\left\{U_{i}(t)\right\}_{i=1}^{d}$ is a $d$-dimensional vector of Brownian motions with the correlation matrix $X^{f}(t)=\left\{\chi_{i, j}(t)\right\}$. The following proposition shows how model parameters need to be set for these dynamics to hold.

Proposition 13.3.2. Let us define the $d \times d$ matrix-valued process $H^{f}$ by

$$
H^{f}(t)=\left(\begin{array}{c}
h\left(t+\delta_{1}\right)^{\top} \\
\vdots \\
h\left(t+\delta_{d}\right)^{\top}
\end{array}\right)
$$

and $\sigma^{f}$ by

$$
\begin{aligned}
& \sigma^{f}(t, f(t)) \\
& \quad=\operatorname{diag}\left(\left(\lambda_{1}^{f}(t)\left(\alpha_{1}^{f}(t)+b_{1}^{f}(t) f_{1}(t)\right), \ldots, \lambda_{d}^{f}(t)\left(\alpha_{d}^{f}(t)+b_{d}^{f}(t) f_{d}(t)\right)\right)^{\top}\right)
\end{aligned}
$$

where $f(t)=\left(f_{1}(t), \ldots, f_{d}(t)\right)^{\top}$. Also, let $D^{f}(t)$ be specified by $X^{f}(t)= D^{f}(t)^{\top} D^{f}(t)$. In (13.72), let us set

$$
\begin{align*}
\sigma_{r}(t, \omega)^{\top} & =\sqrt{z(t)} \sigma_{x}(t, x(t), y(t))^{\top}  \tag{13.77}\\
\sigma_{x}(t, x(t), y(t))^{\top} & =H(t) H^{f}(t)^{-1} \sigma^{f}(t, f(t)) D^{f}(t)^{\top}
\end{align*}
$$
^eq13-77

where $\sigma_{x}$ is a function of $x(t), y(t)$ because $f_{i}$ 's are, see (13.73). Then the $q G$ model in Proposition 13.3 . 1 is consistent with the benchmark rate dynamics of equation (13.76).

Proof. From (13.73) and using the definition of the vector $f(t)$, we obtain

$$
f(t)=f(0)+H^{f}(t) H(t)^{-1}(x(t)+y(t) G(t, T))
$$

Then, with the help of (13.72),

$$
\begin{aligned}
d f(t) & =O(d t)+H^{f}(t) H(t)^{-1} d x(t) \\
& =O(d t)+H^{f}(t) H(t)^{-1} \sigma_{r}(t, \omega)^{\top} d W(t)
\end{aligned}
$$

Using (13.77), we obtain

$$
\begin{aligned}
d f(t) & =O(d t)+\sqrt{z(t)} H^{f}(t) H(t)^{-1} H(t) H^{f}(t)^{-1} \sigma^{f}(t, f(t)) D^{f}(t)^{\top} d W(t) \\
& =O(d t)+\sqrt{z(t)} \sigma^{f}(t, f(t)) D^{f}(t)^{\top} d W(t)
\end{aligned}
$$

The statement of the proposition follows once we set $d U(t)=D^{f}(t)^{\top} d W(t)$. $\square$

With the parameterization outlined in Proposition 13.3.2, the benchmark rates follow "SV-like" dynamics (13.76), and we can reasonably expect Libor and swap rates to follow similar dynamics, at least approximately. Besides reducing the generic qG model to a familiar class of dynamics, the parameterization in Proposition 13.3.2 also achieves a clear distinction between the effects of the various model parameters: we now have volatility parameters $\left\{\lambda_{i}^{f}(t)\right\}$, rate correlation parameters $\left\{\chi_{i, j}(t)\right\}$, skew parameters $\left\{b_{i}^{f}(t)\right\}$, and the volatility of variance $\eta(t)$.

As was the case for the model in Section 12.1.7, the qG model above has enough flexibility to calibrate to $d$ swaption strips, if we assume that the mean reversions $\varkappa(t)$ and the correlation matrix $X^{f}(t)$ are specified prior to volatility calibration. If the swaption strips are of constant-tenor type - a sensible choice for $d>1$ - then it is natural to set the tenors of the swaptions we want to use in calibration equal to the benchmark tenors $\delta_{i}$, $i=1, \ldots, d$. With the stochastic volatility parameterization (13.74), (13.77), we can calibrate

1. At-the-money volatilities for $d$ swaption strips.
2. Slopes of the volatility smiles for $d$ swaption strips.
3. Curvature of the smile for one swaption strip.

The last point is not as restrictive as it may appear, since the curvatures of the smiles for swaptions of the same expiry but different tenors tend to be fairly similar.

To parametrize the model in a suitable way, let us choose a tenor structure $0=T_{0}<\ldots<T_{N}$ and denote the swap rates/annuities for the $i$-th swaption strip, $i=1, \ldots, d$, in the calibration set by

$$
S_{n, \mu_{i}(n)}(t, x(t), y(t)), \quad A_{n, \mu_{i}(n)}(t, x(t), y(t)), \quad n=1, \ldots, N-1 .
$$

Extending (13.39), it is natural to define, for $i=1, \ldots, d$,

$$
\begin{aligned}
\lambda_{i}^{f}(t) & =\sum_{n=1}^{N-1} \lambda_{i, n} 1_{\left\{t \in\left(T_{n-1}, T_{n}\right]\right\}}, \quad \alpha_{i}^{f}(t)=\sum_{n=1}^{N-1} S_{n, \mu_{i}(n)}(0) 1_{\left\{t \in\left(T_{n-1}, T_{n}\right]\right\}}, \\
b_{i}^{f}(t) & =\sum_{n=1}^{N-1} b_{i, n} D_{i, n} 1_{\left\{t \in\left(T_{n-1}, T_{n}\right]\right\}}, \quad \eta(t)=\sum_{n=1}^{N-1} \eta_{n} 1_{\left\{t \in\left(T_{n-1}, T_{n}\right]\right\}},
\end{aligned}
$$

where the skew scalings $D_{i, n}$ are given by an approximate (as we ignore $\partial y$ terms) derivative of $S_{n, \mu_{i}(n)}$ "in the direction" of $f_{i}$,

$$
D_{i, n}=1^{\top} H\left(t+\delta_{i}\right) H(t)^{-1}\left(\nabla S_{n, \mu_{i}(n)}\right) \text {, }
$$

with

$$
\nabla S \triangleq\left(\partial S / \partial x_{1}, \ldots, \partial S / \partial x_{d}\right)
$$

In summary, the volatility smile parameters in the model are the $(2 d+1) \times (N-1)$ parameters $\left\{\lambda_{i, n}\right\},\left\{b_{i, n}\right\}$ and $\left\{\eta_{n}\right\}$.

### 13.3.3 Swap Rate Dynamics and Approximations

The swap rate dynamics in the multi-factor qG model can be derived and simplified by techniques similar to those we applied earlier in the one-factor case. First, we establish a multi-dimensional counterpart to Proposition 13.1.2, i.e. the exact dynamics of a given swap rate $S(t)$ in the annuity measure corresponding to its annuity $A(t)$, see (13.13)-(13.14).

Proposition 13.3.3. In a multi-factor stochastic volatility quasi-Gaussian model with volatility parameterization (13.75), the dynamics of a swap rate $S(t)$ defined by (13.13) are given by

$$
\begin{equation*}
d S(t)=\sqrt{z(t)}\left((\nabla S) \sigma_{x}^{\top} \sigma_{x}(\nabla S)^{\top}\right)^{1 / 2} d U^{A}(t) \tag{13.78}
\end{equation*}
$$
^eq13-78

where all functions are understood to be evaluated at $(t, x(t), y(t))$, and $U^{A}(t)$ is a one-dimensional Brownian motion in the annuity measure $\mathrm{Q}^{A}$.

Proof. By standard arguments of using Ito's lemma on $S(t)$ and a martingale property of $S(t)$.

Using the Markovian projection method (see Appendix A), the dynamics of (13.78) can be approximated, for the purposes of pricing European options, with

$$
\begin{equation*}
d S(t)=\sqrt{z(t)} \varphi(t, S(t)) d U^{A}(t) \tag{13.79}
\end{equation*}
$$
^eq13-79

where

$$
\varphi(t, s)^{2}=\mathrm{E}^{A}\left((\nabla S) c_{x}(\nabla S)^{\top} \mid S(t)=s\right)
$$

and we have denoted

$$
\begin{equation*}
c_{x}=c_{x}(t, x, y)=\sigma_{x}(t, x, y)^{\top} \sigma_{x}(t, x, y) \tag{13.80}
\end{equation*}
$$
^eq13-80

We expect the local volatility term in (13.79) to mostly control the slope of the volatility smile, hence we look for a linear approximation to $\varphi(t, s)$ in $s$. First, we need to choose the point for expansion. One can use $x(t)=0, y(t)=0$ as a decent choice; however, using $\mathrm{E}^{A}(x(t)), \mathrm{E}^{A}(y(t))$ or approximations thereof is, as always, preferable.

Proposition 13.3.4. Let

$$
\begin{equation*}
\sigma_{x}^{0}(t)=\sigma_{x}(t, 0,0) \tag{13.81}
\end{equation*}
$$
^eq13-81

Then

$$
\mathrm{E}^{A}(y(t)) \approx \bar{y}(t) \triangleq H(t)\left(\int_{0}^{t} H(s)^{-1} \sigma_{x}^{0}(s)^{\top} \sigma_{x}^{0}(s) H(s)^{-1} d s\right) H(t)
$$

An approximation $\bar{x}(t)$ to $\mathrm{E}^{A}(x(t))$ is given by

$$
\bar{x}(t)=H(t)\left(\int_{0}^{t} H(s)^{-1}\left(\bar{y}(s) 1-\sigma_{x}^{0}(s)^{\top} \sigma_{x}^{0}(s) G_{A}(s)\right) d s\right), \quad t \in\left[0, T_{1}\right]
$$

where (recall that $S(t)=S_{0, N}(t), A(t)=A_{0, N}(t)$ )

$$
G_{A}(s)=\frac{1}{A(0)} \sum_{n=0}^{N-1} \tau_{n} P\left(0, T_{n+1}\right) G\left(s, T_{n+1}\right) .
$$

Proof. The result for $\mathrm{E}^{A}(y(t))$ follows after approximating $\sigma_{r}(t, \omega)$ in (13.72) by $\sigma_{x}^{0}(t)$ and then proceeding as in the proof of Proposition 13.1.4 (see also Section 12.1.1.1).

Using the same approximation and replacing $y(t)$ with $\bar{y}(t)$, we obtain from (13.72) the following SDE for $x_{g}(t)$, a Gaussian approximation to $x(t)$,

$$
d x_{g}(t)=\left(\bar{y}(t) 1-\varkappa(t) x_{g}(t)\right) d t+\sigma_{x}^{0}(t)^{\top} d W(t)
$$

For a given $T>0$, in the $T$-forward measure,

$$
d x_{g}(t)=\left(\bar{y}(t) 1-\sigma_{x}^{0}(t)^{\top} \sigma_{x}^{0}(t) G(t, T)-\varkappa(t) x_{g}(t)\right) d t+\sigma_{x}^{0}(t)^{\top} d W^{T}(t)
$$

where $d W^{T}$ is a driftless Brownian motion; hence

$$
\mathrm{E}^{T}\left(x_{g}(t)\right)=H(t)\left(\int_{0}^{t} H(s)^{-1}\left(\bar{y}(s) 1-\sigma_{x}^{0}(s)^{\top} \sigma_{x}^{0}(s) G(s, T)\right) d s\right)
$$

Then,

$$
\begin{aligned}
\mathrm{E}^{A}(x(t)) \approx & \mathrm{E}^{A}\left(x_{g}(t)\right) \\
= & \frac{1}{A(0)} \mathrm{E}\left(\beta(t)^{-1} A(t) x_{g}(t)\right) \\
= & \frac{1}{A(0)} \sum_{n=0}^{N-1} \tau_{n} \mathrm{E}\left(x_{g}(t) \beta(t)^{-1} P\left(t, T_{n+1}\right)\right) \\
= & \sum_{n=0}^{N-1} \frac{\tau_{n} P\left(0, T_{n+1}\right)}{A(0)} \mathrm{E}^{T_{n+1}}\left(x_{g}(t)\right) \\
= & \sum_{n=0}^{N-1} \frac{\tau_{n} P\left(0, T_{n+1}\right)}{A(0)} H(t) \\
& \times \int_{0}^{t} H(s)^{-1}\left(\bar{y}(s) 1-\sigma_{x}^{0}(s)^{\top} \sigma_{x}^{0}(s) G\left(s, T_{n+1}\right)\right) d s
\end{aligned}
$$

and the result follows.
In preparation for our next result, let us define $S_{g}(t)$ to be the Gaussian approximation to $S(t)$, i.e. a process with the dynamics given by (13.78) where all functions are evaluated at ( $t, 0,0$ ) and there is no stochastic volatility,

$$
d S_{g}(t)=\left.\left((\nabla S)\left(\sigma_{x}^{0}\right)^{\top} \sigma_{x}^{0}(\nabla S)^{\top}\right)^{1 / 2}\right|_{(t, 0,0)} d U^{A}(t)
$$

The dynamics of $S_{g}(t)$ can alternatively be represented using a multidimensional Brownian motion,

$$
\begin{equation*}
d S_{g}(t)=\nabla S(t, 0,0)\left(\sigma_{x}^{0}\right)^{\top} d W^{A}(t) \tag{13.82}
\end{equation*}
$$
^eq13-82

Theorem 13.3.5. Let $\bar{x}(t)$ and $\bar{y}(t)$ be as given in Proposition 13.3.4, and let $c_{x}(t)$ be defined as in (13.80). For pricing European swaptions, the dynamics of the swap rate $S(t)$ (defined by (13.13)) in the multi-factor quasi-Gaussian model with the volatility structure given by (13.77) can be approximated by the following time-dependent Heston dynamics,

$$
d S(t) \approx \sqrt{z(t)} \lambda_{S}(t)\left(b_{S}(t) S(t)+\left(1-b_{S}(t)\right) S(0)\right) d U^{A}(t)
$$

where

$$
\begin{aligned}
\lambda_{S}(t)= & \frac{1}{S(0)}\left((\nabla S) c_{x}(\nabla S)^{\top}\right)^{1 / 2}(t, \bar{x}(t), \bar{y}(t)) \\
b_{S}(t)= & S(0)\left(\frac{1}{2} \frac{(\nabla S) d_{x}(\nabla S)^{\top}}{(\nabla S)\left(c_{x}(\nabla S)^{\top}(\nabla S) c_{x}\right)(\nabla S)^{\top}}\right. \\
& \left.+\frac{(\nabla S)\left(c_{x}\left(\nabla^{2} S\right) c_{x}\right)(\nabla S)^{\top}}{(\nabla S)\left(c_{x}(\nabla S)^{\top}(\nabla S) c_{x}\right)(\nabla S)^{\top}}\right)\left.\right|_{(t, \bar{x}(t), \bar{y}(t))}
\end{aligned}
$$

Here we have denoted

$$
\begin{gather*}
d_{x}=\sum_{l=1}^{d}\left(c_{x}(\nabla S)^{\top}\right)_{l}\left(\frac{\partial \sigma_{x}^{\top}}{\partial x_{l}} \sigma_{x}+\sigma_{x}^{\top} \frac{\partial \sigma_{x}}{\partial x_{l}}\right) \\
\frac{\partial \sigma_{x}^{\top}}{\partial x_{l}}=H(t) H^{f}(t)^{-1} \operatorname{diag}\left(\lambda_{1}^{f}(t) b_{1}^{f}(t) \frac{h_{l}\left(t+\delta_{1}\right)}{h_{l}(t)}, \ldots\right. \\
\left.\ldots, \lambda_{d}^{f}(t) b_{d}^{f}(t) \frac{h_{l}\left(t+\delta_{d}\right)}{h_{l}(t)}\right) D^{f}(t)^{\top} \tag{13.83}
\end{gather*}
$$
^eq13-83

and

$$
\nabla^{2} S \triangleq\left\{\partial^{2} S / \partial x_{i} \partial x_{j}\right\}_{i, j=1}^{d}
$$

Proof. In (13.79), consider the conditional expected value

$$
\mathrm{E}^{A}\left((\nabla S) c_{x}(\nabla S)^{\top} \mid S(t)=s\right)
$$

Expanding the integrand $(\nabla S) c_{x}(\nabla S)^{\top}$ around $(t, \bar{x}(t), \bar{y}(t))$ to first order in $x$, we obtain

$$
\begin{aligned}
& \left((\nabla S) c_{x}(\nabla S)^{\top}\right)(t, x(t), y(t)) \approx\left((\nabla S) c_{x}(\nabla S)^{\top}\right)(t, \bar{x}(t), \bar{y}(t)) \\
& \quad+\left.\nabla\left((\nabla S) c_{x}(\nabla S)^{\top}\right)\right|_{(t, \bar{x}(t), \bar{y}(t))}(x(t)-\bar{x}(t))
\end{aligned}
$$

Then

$$
\begin{aligned}
& \mathrm{E}^{A}\left((\nabla S) c_{x}(\nabla S)^{\top} \mid S(t)=s\right) \approx\left((\nabla S) c_{x}(\nabla S)^{\top}\right)(t, \bar{x}(t), \bar{y}(t)) \\
& +\left.\nabla\left((\nabla S) c_{x}(\nabla S)^{\top}\right)\right|_{(t, \bar{x}(t), \bar{y}(t))} \mathrm{E}^{A}(x(t)-\bar{x}(t) \mid S(t)=s) .
\end{aligned}
$$

Using a Gaussian approximation for the conditional expected value,

$$
\begin{aligned}
\mathrm{E}^{A}(x(t)-\bar{x}(t) \mid S(t)=s) & \approx \mathrm{E}^{A}\left(x_{g}(t)-\bar{x}(t) \mid S_{g}(t)=s\right) \\
& =\frac{\operatorname{Cov}\left(S_{g}(t), x_{g}(t)\right)}{\operatorname{Var}\left(S_{g}(t)\right)}(s-S(0)) \\
& =\frac{\left(\sigma_{x}^{0}\right)^{\top} \sigma_{x}^{0}(\nabla S(t, 0,0))^{\top}}{(\nabla S(t, 0,0))\left(\sigma_{x}^{0}\right)^{\top} \sigma_{x}^{0}(\nabla S(t, 0,0))^{\top}}(s-S(0)) \\
& \left.\approx \frac{c_{x}(\nabla S)^{\top}}{(\nabla S) c_{x}(\nabla S)^{\top}}\right|_{(t, \bar{x}(t), \bar{y}(t))}(s-S(0))
\end{aligned}
$$

For any $l, l=1, \ldots, d$, we have

$$
\begin{aligned}
& \frac{\partial}{\partial x_{l}}\left((\nabla S) c_{x}(\nabla S)^{\top}\right) \\
& \quad=\left(\left(\nabla^{2} S\right)_{l}\right)^{\top} c_{x}(\nabla S)^{\top}+(\nabla S)\left(\frac{\partial c_{x}}{\partial x_{l}}\right)(\nabla S)^{\top}+(\nabla S) c_{x}\left(\nabla^{2} S\right)_{l} \\
& \quad=(\nabla S)\left(\frac{\partial \sigma_{x}^{\top}}{\partial x_{l}} \sigma_{x}+\sigma_{x}^{\top} \frac{\partial \sigma_{x}}{\partial x_{l}}\right)(\nabla S)^{\top}+2(\nabla S) c_{x}\left(\nabla^{2} S\right)_{l}
\end{aligned}
$$

where $\left(\nabla^{2} S\right)_{l}$ is the $l$-th column of the matrix $\nabla^{2} S$. Thus

$$
\begin{aligned}
\nabla\left((\nabla S) c_{x}(\nabla S)^{\top}\right) & c_{x}(\nabla S)^{\top} \\
& =(\nabla S) d_{x}(\nabla S)^{\top}+2(\nabla S)\left(c_{x}\left(\nabla^{2} S\right) c_{x}\right)(\nabla S)^{\top}
\end{aligned}
$$

where the matrix $d_{x}$ is given in the statement of the theorem, and

$$
\begin{aligned}
& \left.\nabla\left((\nabla S) c_{x}(\nabla S)^{\top}\right)\right|_{(t, \bar{x}(t), \bar{y}(t))} \mathrm{E}^{A}(x(t)-\bar{x}(t) \mid S(t)=s) \\
& \quad \approx\left(\frac{(\nabla S) d_{x}(\nabla S)^{\top}}{(\nabla S) c_{x}(\nabla S)^{\top}}\right. \\
& \left.\quad+2 \frac{(\nabla S)\left(c_{x}\left(\nabla^{2} S\right) c_{x}\right)(\nabla S)^{\top}}{(\nabla S) c_{x}(\nabla S)^{\top}}\right)\left.\right|_{(t, \bar{x}(t), \bar{y}(t))}(s-S(0))
\end{aligned}
$$

Thus, we obtain (all terms on the right-hand side evaluated at $(t, \bar{x}(t), \bar{y}(t))$ )

$$
\begin{aligned}
& \varphi(t, s)^{2} \approx(\nabla S) c_{x}(\nabla S)^{\top} \\
& +\left(\frac{(\nabla S) d_{x}(\nabla S)^{\top}}{(\nabla S) c_{x}(\nabla S)^{\top}}+2 \frac{(\nabla S)\left(c_{x}\left(\nabla^{2} S\right) c_{x}\right)(\nabla S)^{\top}}{(\nabla S) c_{x}(\nabla S)^{\top}}\right)(s-S(0))
\end{aligned}
$$

so that, linearizing around $s \approx S(0)$,

$$
\begin{aligned}
\varphi(t, s) & \approx\left((\nabla S) c_{x}(\nabla S)^{\top}\right)^{1 / 2}+\frac{1}{2\left((\nabla S) c_{x}(\nabla S)^{\top}\right)^{1 / 2}} \\
& \times\left(\frac{(\nabla S) d_{x}(\nabla S)^{\top}}{(\nabla S) c_{x}(\nabla S)^{\top}}+2 \frac{\left.(\nabla S)\left(c_{x}\left(\nabla^{2} S\right) c_{x}\right)^{\top} \nabla S\right)^{\top}}{(\nabla S) c_{x}(\nabla S)^{\top}}\right)(s-S(0))
\end{aligned}
$$

Setting

$$
\lambda_{S}(t)=\frac{\varphi(t, S(0))}{S(0)}, \quad b_{S}(t)=S(0) \frac{\frac{\partial \varphi}{\partial s}(t, S(0))}{\varphi(t, S(0))}
$$

the main statement of the theorem follows. Finally, from the definition of $\sigma_{x}(t, x, y)$ (see (13.77)),
$\frac{\partial \sigma_{x}^{\top}}{\partial x_{l}}=H(t) H^{f}(t)^{-1} \operatorname{diag}\left(\lambda_{1}^{f}(t) b_{1}^{f}(t) \frac{\partial f_{1}(t)}{\partial x_{l}}, \ldots, \lambda_{d}^{f}(t) b_{d}^{f}(t) \frac{\partial f_{d}(t)}{\partial x_{l}}\right) D^{f}(t)^{\top}$
and the expression (13.83) follows from

$$
\frac{\partial f_{i}(t)}{\partial x_{l}}=h_{l}\left(t+\delta_{i}\right) / h_{l}(t) .
$$

Remark 13.3.6. Using averaging techniques, the time-dependent Heston dynamics could be easily translated into time-independent ones. The derivation and the result essentially mimic those for the one-dimensional case, see Proposition 13.2.2.

### 13.3.4 Volatility Calibration

As explained in Section 13.3.2, the model (13.72), (13.74), (13.77) has enough degrees of freedom to calibrate to the smiles of $d$ swaption strips if the vector function $h(t)$, as well as the time-dependent correlation matrix of benchmark rates $X^{f}(t)$, are specified exogenously. For multi-factor quasiGaussian models, the strips are usually taken to be constant-tenor strips with swap tenors matching benchmark rate tenors. With $d=4$ or 5 factors being a typical choice of dimensionality, a calibration to 4 or 5 swaption strips essentially recovers the whole universe of swaption volatilities, so there is little need to choose calibration targets in a product-specific way.

As with the one-factor stochastic volatility qG model, we favor splitting the calibration into two main steps. First, we calibrate the volatility of volatility curve $\eta(t)$ to the market-implied curvatures of the smiles or, better yet, to average curvatures of volatility smiles across swap tenors, as we only have the flexibility of making $\eta(t)$ time-specific, not tenor-specific. After that, the main calibration is performed, matching the overall levels and slopes of volatility smiles of $d$ swaption strips. We omit the details as they follow closely the algorithm of Section 13.2.3, with the only difference being that on each time step, the calibration problem involves $d$ swaptions and not 1 . Since all formulas are closed-form, the calibration is essentially instantaneous.

### 13.3.5 Mean Reversions, Correlations, and Numerical Schemes

In the multi-factor context, the time-dependent "loadings" vector $h(t)$ essentially defines the interpolation rule, i.e. how the volatilities and correlations of non-benchmark rates are obtained from those of the benchmark rates. We advocate choosing $d$ fixed values of mean reversions and using them for all cases - note that they should all be different, since the inverse of matrix $H_{f}(t)$ is required to exist. For example, a reasonable choice is to span the interval $[0,1]$ with mean reversions while always including the point 0 , i.e. set

$$
\varkappa(t)=\operatorname{diag}\left((0.015,0.15,0.30,1.20)^{\top}\right)
$$

for a 4 -dimensional model, corresponding to benchmark tenors

$$
\left\{\delta_{1}, \ldots, \delta_{4}\right\}=\{6 m, 2 y, 10 y, 30 y\}
$$

However, in principle at least, the mean reversions can be calibrated as well, giving us additional $d$ strips to calibrate to. Formulas for mean reversion calibration could be derived in the same way as for the one-dimensional case, see Section 13.1.8.

Additionally, the correlation matrix between benchmark rates, $X^{f}(t)$, could in principle be used in calibration, particularly when valuing products with strong correlation sensitivity. In this case, to capture market-implied correlation information in the model, one sometimes chooses to best-fit
market-observed prices of CMS spread options by tweaking benchmark rate correlations. We discuss this in more details in the context of LM models in Section 14.5.9. For now we just note that spread option values exhibit a correlation smile, i.e. the dependence of implied correlation on the strike of the spread option (see Section 17.4.2), so the choice of the strike to calibrate to should be carefully considered.

Finally, some brief words on numerical implementation. For $d>1$, PDE methods quickly become impractical - even for the simple case of $d=2$, there are 3 auxiliary ( $y$ ) variables to take care of, pushing the dimension of the PDE to 5 which is prohibitively expensive in virtually all applications. However, by using tricks such as "freezing" or projecting some of the auxiliary variables, as in Section 13.1.9.4, a PDE scheme for $d=2$ or $d=3$ could possibly be made viable.

With Monte Carlo methods, the usual considerations apply, and no special tricks beyond those for a one-dimensional stochastic volatility model are required.

## 13. A Appendix: Density Approximation

We prove Theorem 13.1.14 in a number of steps. Denoting for brevity

$$
\sigma=\sigma^{0}(0),
$$

and using the notations of Section 13.1.9.1, we can write down the approximate risk-neutral model dynamics as

$$
\begin{align*}
& d x(t)=(\bar{y}(t)-\varkappa x(t)) d t+\sigma v(x(t)) d W(t)  \tag{13.84}\\
& d \bar{y}(t)=\left(\sigma^{2}-2 \varkappa \bar{y}(t)\right) d t \tag{13.85}
\end{align*}
$$
^eq13-84
^eq13-85

## 13.A.1 Simplified Forward Measure Dynamics

As a first step we simplify the dynamics of the state process.
Proposition 13. A.1. For small time $T$, the distribution of $x(T)$ in the $T$-forward measure can be approximated by the distribution of $\widetilde{x}(T)$, with the dynamics of $\widetilde{x}(t)$ given by

$$
\begin{equation*}
d \widetilde{x}(t)=-\varkappa \widetilde{x}(t) d t+\sigma v(\widetilde{x}(t)) d W^{T}(t) \tag{13.86}
\end{equation*}
$$
^eq13-86

where $W^{T}(t)$ is a Brownian motion in the $T$-forward measure, and $v(x)$ is the same as in (13.84) (and defined by (13.48)).

Remark 13.A.2. Note that the statement is only about approximating the marginal distribution of $x(T)$ with $\widetilde{x}(T)$, not the dynamics of $x(\cdot)$ with $\widetilde{x}(\cdot)$ in the $T$-forward measure.

Proof. The process

$$
d W^{T}(t)=d W(t)+\sigma v(x(t)) G(t, T) d t
$$

is a driftless Brownian motion in the $T$-forward measure, hence

$$
d x(t)=\left(\bar{y}(t)-\sigma^{2} v(x(t))^{2} G(t, T)-\varkappa x(t)\right) d t+\sigma v(x(t)) d W^{T}(t)
$$

Also,

$$
\bar{y}(t)=\sigma^{2} \int_{0}^{t} e^{-2 \varkappa(t-s)} d s=\sigma^{2} G_{2}(t)
$$

using notation introduced in (13.3). Thus

$$
\begin{aligned}
& x(T)=\sigma^{2} \int_{0}^{T} e^{-\varkappa(T-t)}\left(G_{2}(t)-v(x(t))^{2} G(t, T)\right) d t \\
& \quad+\sigma \int_{0}^{T} e^{-\varkappa(T-t)} v(x(t)) d W^{T}(t)
\end{aligned}
$$

On the other hand, the instantaneous forward rate $f(t, T)$ is a martingale in the $T$-forward measure,

$$
\mathrm{E}^{T}(f(T, T))=f(0, T)
$$

which implies that

$$
\mathrm{E}^{T}(x(T))=0
$$

Therefore

$$
\mathrm{E}^{T}\left(\int_{0}^{T} e^{-\varkappa(T-t)}\left(G_{2}(t)-v(x(t))^{2} G(t, T)\right) d t\right)=0
$$

where the equality is only approximate as we replaced $y(t)$ with $\bar{y}(t)$, but is as accurate as the approximation of $y(t)$ with $\bar{y}(t)$. We replace the equality with a somewhat stronger condition

$$
\int_{0}^{T} e^{-\varkappa(T-t)}\left(G_{2}(t)-v(x(t))^{2} G(t, T)\right) d t=0
$$

leading to

$$
x(T)=\sigma \int_{0}^{T} e^{-\varkappa(T-t)} v(x(t)) d W^{T}(t)
$$

which is equivalent to (13.86).

## 13.A. 2 Effective Volatility

From now on we consider the distribution of $x(T)$ to be given by (13.86), and we drop the tilde to simplify the notations. The (undiscounted) value of a call option on $x(t)$ with strike $k$ is denoted by

$$
\begin{equation*}
c(t, k)=\mathrm{E}^{T}(x(t)-k)^{+} \tag{13.87}
\end{equation*}
$$
^eq13-87

By the Bachelier formula (Remark 7.2.9), this function is known explicitly for $v(x) \equiv 1$ (since (13.86) is a Gaussian SDE then), and we denote it $c_{0}$ :

$$
c_{0}(t, k, \sigma)=\left(x_{0}-k\right) \Phi\left(\frac{x_{0}-k}{\sigma \sqrt{G_{2}(t)}}\right)+\sigma \sqrt{G_{2}(t)} \phi\left(\frac{x_{0}-k}{\sigma \sqrt{G_{2}(t)}}\right)
$$

where $\phi(z)$ is the standard Gaussian PDF, $\Phi(z)$ is the standard Gaussian CDF,$x_{0}=x(0)=0$, and $G_{2}(t)$ is defined in (13.49). Using this expression as the base case, we look for the approximate value of $c(t, k)$, for $x(t)$ governed by (13.86), of the form (compare to the methods of Section 7.5)

$$
\begin{equation*}
c(t, k)=c_{0}(t, k, \sigma \varpi(k)) \tag{13.88}
\end{equation*}
$$
^eq13-88

Here $\varpi(k)$ has the meaning of the effective term volatility. For notational convenience, we also define

$$
\begin{equation*}
\zeta(t, k)=\frac{x_{0}-k}{\sigma \varpi(k) \sqrt{G_{2}(t)}} \tag{13.89}
\end{equation*}
$$
^eq13-89

Then

$$
\begin{equation*}
c(t, k)=\left(x_{0}-k\right) \Phi(\zeta(t, k))+\sigma \varpi(k) \sqrt{G_{2}(t)} \phi(\zeta(t, k)) \tag{13.90}
\end{equation*}
$$
^eq13-90

In the next few sections, we obtain an expression for $\zeta(t, k)$ in the small- $\sigma$ limit. To do so, we, firstly, derive a PDE for $c(t, k)$ defined by (13.87). Then, we substitute the expression (13.90) into the PDE to derive an equation on $\zeta(t, k)$. We drop the terms of order $O\left(\sigma^{2}\right)$ and smaller, and solve the simplified equation for $\zeta(t, k)$. Finally, we obtain a CDF and a PDF of $x(t)$ by differentiating $c(t, k)$ in strike.

## 13.A. 3 The Forward Equation for Call Options

In this section we derive a PDE for $c(t, k)$ in the variables $t$ (time to expiry) and $k$ (strike), just as we did in Proposition 7.4.2 for vanilla local volatility models. Let $\psi(t, k)$ be the density of $x(t)$ and $\Psi(t, k)$ its CDF. We use subscripts to denote partial derivatives, and primes to denote derivatives of functions of a single variable.

Proposition 13.A.3. The function $c(t, k)$ defined by (13.87) for $x(t)$ following (13.86) satisfies

$$
\begin{equation*}
c_{t}(t, k)=-\varkappa\left[c(t, k)-k c_{k}(t, k)\right]+\frac{\sigma^{2}}{2} v(k)^{2} c_{k k}(t, k) \tag{13.91}
\end{equation*}
$$
^eq13-91

with the initial condition

$$
c(0, k)=\delta(k)
$$

where $\delta(k)$ is the Dirac delta function at 0 .
Proof. Follows that of Proposition 7.4.2 closely, with the use of the following identity

$$
\begin{aligned}
\mathrm{E}\left(x(t) 1_{\{x(t)>k\}}\right)=\mathrm{E}\left((x(t)-k) 1_{\{x(t)>k\}}\right)+k \mathrm{E}( & \left.1_{\{x(t)>k\}}\right) \\
& =c(t, k)-k c_{k}(t, k)
\end{aligned}
$$

## 13.A. 4 Asymptotic Expansion

Lemma 13.A.4. The following holds for $c(t, k)$ as defined by (13.90),

$$
\begin{align*}
\frac{c_{t}}{\phi(\zeta)} & =\sigma \varpi(k)\left(\sqrt{G_{2}(t)}\right)^{\prime}  \tag{13.92}\\
\frac{c-k c_{k}}{\phi(\zeta)} & =\left(\varpi(k)-k \varpi^{\prime}(k)\right) \sigma \sqrt{G_{2}(t)}  \tag{13.93}\\
\frac{c_{k k}}{\phi(\zeta)} & =-\zeta_{k}+\sigma \sqrt{G_{2}(t)}\left(\varpi^{\prime \prime}(k)-\varpi^{\prime}(k) \zeta_{k} \zeta\right) \tag{13.94}
\end{align*}
$$
^eq13-92
^eq13-93
^eq13-94

where $\varpi, \zeta$ are defined by (13.88), (13.89).
Proof. Follows by straightforward differentiation of $c(t, k)$ defined by (13.90) with the help of the identities

$$
\begin{equation*}
\phi^{\prime}(z)=-z \phi(z), \quad \phi^{\prime \prime}(z)=\left(z^{2}-1\right) \phi(z) \tag{13.95}
\end{equation*}
$$
^eq13-95

Proposition 13.A.5. If the function $\varpi(k)$ is such that

$$
\pi(k) \triangleq \frac{k}{\varpi(k)}
$$

satisfies the ODE

$$
v(k)^{2}\left(\pi^{\prime}(k)\right)^{2}+2 \varkappa G_{2}(t) \pi(k)\left(\frac{k}{\pi(k)}\right)^{\prime}-1=0, \quad k \in \mathbb{R}
$$

with the boundary condition

$$
\pi(0)=0
$$

then

$$
c(t, k)=c_{0}(t, k, \sigma \varpi(k))+O\left(\sigma^{2}\right)
$$

i.e. $c_{0}(t, k, \sigma \varpi(k))$ is an approximation to $c(t, k)$ (given in (13.90)) to the first order in $\sigma, \sigma \rightarrow 0$.

Proof. Substituting (13.92)-(13.94) into the PDE (13.91) and keeping only the terms of order $O\left(\sigma^{2}\right)$ we obtain,

$$
\begin{aligned}
\sigma \varpi(k) \frac{d}{d t} \sqrt{G_{2}(t)}=- & \varkappa\left(\varpi(k)-k \varpi^{\prime}(k)\right) \sigma \sqrt{G_{2}(t)} \\
& +\frac{\sigma^{2}}{2} v(k)^{2}\left(-\zeta^{\prime}(k)-\sigma \sqrt{G_{2}(t)} \varpi^{\prime}(k) \zeta^{\prime}(k) \zeta(k)\right) .
\end{aligned}
$$

Dividing by $\sigma$ and using the fact that

$$
\begin{equation*}
\sigma \varpi(k) \sqrt{G_{2}(t)} \zeta(k)=-k, \tag{13.96}
\end{equation*}
$$
^eq13-96

we obtain

$$
\begin{align*}
\varpi(k) \frac{d}{d t} \sqrt{G_{2}(t)}=-\varkappa(\varpi(k) & \left.-k \varpi^{\prime}(k)\right) \sqrt{G_{2}(t)} \\
& +\frac{\sigma}{2} v(k)^{2}\left(\left(\frac{k \varpi^{\prime}(k)}{\varpi(k)}-1\right) \zeta^{\prime}(k)\right) . \tag{13.97}
\end{align*}
$$
^eq13-97

By definition of $\pi(k)$, with the help of (13.96), and using the definition of $G_{2}(t)$,

$$
\begin{gathered}
\zeta^{\prime}(k)=-\frac{1}{\sigma \sqrt{G_{2}(t)}} \pi^{\prime}(k), \quad \frac{k \varpi^{\prime}(k)}{\varpi(k)}-1=-\pi^{\prime}(k) \varpi(k) \\
\frac{d}{d t} \sqrt{G_{2}(t)}=\frac{d G_{2}(t) / d t}{2 \sqrt{G_{2}(t)}}, \quad \frac{d}{d t} G_{2}(t)=e^{-2 \varkappa t}=1-2 \varkappa G_{2}(t)
\end{gathered}
$$

which, substituted into (13.97), gives us, after some simplifications,

$$
2 \varkappa G_{2}(t) \frac{k \varpi^{\prime}(k)}{\varpi(k)}+v(k)^{2}\left(\pi^{\prime}(k)\right)^{2}=1
$$

In addition,

$$
\frac{k \varpi^{\prime}(k)}{\varpi(k)}=\pi(k)\left(\frac{k}{\pi(k)}\right)^{\prime}
$$

so, finally,

$$
v(k)^{2}\left(\pi^{\prime}(k)\right)^{2}+2 \varkappa G_{2}(t) \pi(k)\left(\frac{k}{\pi(k)}\right)^{\prime}-1=0
$$

To obtain the boundary conditions on $\pi(k)$, we recall that $\varpi(k)=k / \pi(k)$ and, as $\varpi(k)$ has to be bounded at $k=0$, we must have $\pi(k)=0$.

## 13.A. 5 Proof of Theorem 13.1.14

The statement of the theorem follows by using Proposition 13.A. 3 to simplify the model dynamics to (13.86), and then differentiating $c_{0}(T, x, \sigma \varpi(x))$ with $\varpi(x)$ from Proposition 13.A. 5 once with respect to $x$ to obtain the approximate CDF of $x(T)$,

$$
\Psi(T, x)=\frac{\partial c_{0}(T, x, \sigma \varpi(x))}{\partial x}+1 .
$$

We omit tedious but straightforward details.

