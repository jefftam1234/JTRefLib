---
tags: [quant/knowledge, quant/books, interest-rates]
---

## Fundamentals of Interest Rate Modeling

The purpose of this brief chapter is twofold. First, we introduce notations to characterize prices and yields of basic fixed income market securities. In addition to providing the foundation for a more expansive discussion of fixed income markets (which we shall undertake in Chapter 5), this part of the chapter serves to identify and characterize a number of probability measures that are of fundamental importance in models for the term structure of interest rates. A brief discussion of measures used in a two-currency setting is also provided.

In the second part of the chapter, we discuss general characteristics of models with dynamics driven by vector-valued Brownian motions. This analysis leads to the fundamental class of Heath-Jarrow-Morton (HJM) (see Heath et al. [1992]) models of continuously compounded forward rates. Among other special cases, we discuss in some detail tractable HJM models with Gaussian volatility structure, and provide some results for the case where such models are Markovian. These discussions continue in Chapters 10 through 12 where we consider one- and multi-factor short rate models, and in the Chapter 13 where we introduce the important class of quasi-Gaussian HJM models with local and stochastic volatility.

### 4.1 Fixed Income Notations

### 4.1.1 Bonds and Forward Rates

As in earlier chapters, let $P(t, T)$ denote the time $t$ price of a zero-coupon bond (also known as a discount bond) delivering for certain $\$ 1$ at maturity $T \geq t$. Suppose we are interested in purchasing at some future time $T$ a zero-coupon bond maturing at $T+\tau, \tau>0$. At time $t<T$, the price of such a bond can be locked in by i) purchasing at time $t$ one ( $T+\tau$ )-maturity zero-coupon bond; and ii) selling short ("shorting") $P(t, T+\tau) / P(t, T)$
$T$-maturity zero-coupon bonds. The time $t$ cost of executing this strategy is zero,

$$
-1 \cdot P(t, T+\tau)+P(t, T+\tau) / P(t, T) \cdot P(t, T)=0
$$

but a flow of

$$
-P(t, T+\tau) / P(t, T)
$$

will take place at time $T$ as the $T$-maturity short position matures. This is compensated by an inflow of $\$ 1$ at time $T+\tau$. In other words, our trading strategy effectively fixes the time $T$ purchase price of the $(T+\tau)$-maturity bond at

$$
P(t, T, T+\tau) \triangleq P(t, T+\tau) / P(t, T), \quad \tau>0
$$

a quantity known as the time $t$ forward price for the zero-coupon bond spanning $[T, T+\tau]$.

It is often convenient to characterize a forward bond price by a discount rate. One such rate is the continuously compounded forward yield $y(t, T, T+ \tau)$, defined by

$$
\begin{equation*}
e^{-y(t, T, T+\tau) \tau}=P(t, T, T+\tau) \tag{4.1}
\end{equation*}
$$
^eq4-1

The time between the maturity of the forward bond and the expiry of the forward contract, i.e. $\tau$, is often called the tenor of the forward bond or the forward yield. In the definition of the continuously compounded yield lies an implicit, and idealized, assumption of continuous reinvestment of investment proceeds. Most actual market quotes, however, are based on discrete-time compounding of proceeds. Accordingly, we define a simple forward rate $L(t, T, T+\tau)$ as

$$
\begin{equation*}
1+\tau L(t, T, T+\tau)=1 / P(t, T, T+\tau) . \tag{4.2}
\end{equation*}
$$
^eq4-2

Again, $\tau$ is the tenor of the forward rate. For an arbitrary set of dates $T=T_{0}<T_{1}<T_{2}<\ldots<T_{n}$, notice that forward bond prices can be recovered from forward rates by simple compounding,

$$
P\left(t, T_{n}\right) / P(t, T)=\prod_{i=1}^{n} \frac{1}{1+\left(T_{i}-T_{i-1}\right) L\left(t, T_{i-1}, T_{i}\right)}
$$

Unless we state otherwise, throughout this book we shall typically make the assumption that spot rates $L(T, T, T+\tau)$ are the Libor (London Interbank Offered Rate) rates quoted in the interbank market. Libor rates are quoted on values of $\tau$ ranging from one week $(\tau=1 / 52)^{1}$ to 12 months $(\tau=1)$, and form the basis for a number of floating-rate derivative contracts, such as interest rate swaps and Eurodollar futures. We shall examine these securities in more detail in Chapter 5.

[^59]In the limit $\tau \downarrow 0$,

$$
L(t, T, T+\tau) \rightarrow f(t, T)
$$

where the quantity $f(t, T)$ is the time $t$ instantaneous forward rate to time $T$. We think of $f(t, T)$ as the forward rate spanning $[T, T+d T]$, observed at time $t$. The relation between instantaneous forward rates and bond prices is given by the continuous compounding formula

$$
\begin{equation*}
P(t, T, T+\tau)=\exp \left(-\int_{T}^{T+\tau} f(t, u) d u\right) \tag{4.3}
\end{equation*}
$$
^eq4-3

such that

$$
\begin{equation*}
f(t, T)=-\frac{\partial \ln P(t, T)}{\partial T} \tag{4.4}
\end{equation*}
$$
^eq4-4

and, from (4.1),

$$
y(t, T, T+\tau)=\tau^{-1} \int_{T}^{T+\tau} f(t, u) d u, \quad f(t, T)=\lim _{\tau \downarrow 0} y(t, T, T+\tau)
$$

We also notice the relationship

$$
f(t, T)=\frac{\partial(y(t, t, T)(T-t))}{\partial T}=y(t, t, T)+(T-t) \frac{\partial y(t, t, T)}{\partial T} .
$$

The quantity

$$
\begin{equation*}
r(t) \triangleq f(t, t) \tag{4.5}
\end{equation*}
$$
^eq4-5

is an $\mathcal{F}_{t}$-measurable random variable known as the short rate or sometimes the spot rate. Loosely speaking, we can think of $r(t)$ as the overnight rate in effect at time $t$.

Finally, let us note that interest rates of various flavors and related quantities are typically quoted in percentage points or sometimes in basis points, where 1 basis point $=1 / 100$ of one percent.

### 4.1.2 Futures Rates

Through the market for Eurodollar futures (see Chapter 5), investors can enter into securities that will pay at time $T$ an amount of

$$
\begin{equation*}
1-L(T, T, T+\tau) \tag{4.6}
\end{equation*}
$$
^eq4-6

At time 0, a Eurodollar futures contract can be entered into at no upfront cost, but with an implicit obligation of the holder to pay at time $T$ per unit of notional

$$
1-F(0, T, T+\tau)
$$

in return for the payout (4.6). Here, $F(t, T, T+\tau)$ is the time $t$ simple futures rate for the period $[T, T+\tau]$. Importantly, the futures rate is marked to market (or resettled) each day, with the day's change in the futures rate immediately credited to or debited from the contract holder's account with the futures exchange. Specifically, after holding the contract for a period of $\Delta=1$ day, the futures contract holder would thus experience a cash flow of

$$
\begin{aligned}
(1-F(\Delta, T, T+\tau))-(1-F & (0, T, T+\tau)) \\
& =-(F(\Delta, T, T+\tau)-F(0, T, T+\tau))
\end{aligned}
$$

Continuing the mark-to-market process to maturity shows that the total amount of cash flow received by the holder on $[0, T]$ is

$$
\begin{equation*}
-(F(T, T, T+\tau)-F(0, T, T+\tau))=-(L(T, T, T+\tau)-F(0, T, T+\tau)) \tag{4.7}
\end{equation*}
$$
^eq4-7

where we have used the fact that $F(T, T, T+\tau)$ must equal $L(T, T, T+\tau)$ to avoid a delivery arbitrage.

The fact that the net cash flow payment (4.7) on a Eurodollar futures contract has been made incrementally on a daily basis has important valuation consequences, and causes the futures rate to differ from the forward rate defined earlier. For instance, under a scenario of rising interest rates, the holder of a Eurodollar futures contract must make payments to the futures exchange. As rates are rising, the contract holder will be faced with a high-rate - and thus unfavorable - borrowing environment for funding these payments. Conversely, when interest rates fall, the reinvestment of received funds will take place at increasingly low rates. Due to the adverse behavior of funding costs and reinvestment gains, we would expect the purchaser of a Eurodollar futures contract to pay less for these instruments than for a comparable instrument without daily mark-to-market. Consequently, we would expect the futures rate to be above the corresponding forward rate. We shall quantify this effect in Section 4.5.1 and, with more advanced models, in Section 16.8.

We notice that we can define instantaneous futures rates $q(t, T)$ in the same fashion as we defined instantaneous forward rates:

$$
q(t, T)=\lim _{\tau \downarrow 0} F(t, T, T+\tau)
$$

### 4.1.3 Annuity Factors and Par Rates

Most fixed income securities involve multiple cash flows taking place on a pre-set schedule of dates, often referred to as a tenor structure,

$$
0 \leq T_{0}<T_{1}<\ldots<T_{N} .
$$

Given a tenor structure, for any two integers $k, m$ satisfying $0 \leq k<N$, $m>0$, and $k+m \leq N$, we can define an annuity factor $A_{k, m}$ by

$$
\begin{equation*}
A_{k, m}(t)=\sum_{n=k}^{k+m-1} P\left(t, T_{n+1}\right) \tau_{n}, \quad \tau_{n}=T_{n+1}-T_{n} \tag{4.8}
\end{equation*}
$$
^eq4-8

Annuity factors provide for compact notation when pricing coupon-bearing securities. For instance, a security making $m$ coupon payments of $c \tau_{n}$ at all $T_{n+1}, n=k, \ldots, k+m-1$, is easily seen to have time $t$ value of

$$
c A_{k, m}(t), \quad t \leq T_{k}
$$

If the security also involves a back-end return of notional at time $T_{k+m}$ (as is the case for a regular coupon-bearing bond), the pricing expression is

$$
\begin{equation*}
c A_{k, m}(t)+P\left(t, T_{k+m}\right), \tag{4.9}
\end{equation*}
$$
^eq4-9

where we assume that the bond has been normalized to have a unit notional. The time $t$ forward price to $T_{k}$ of the security (4.9) is

$$
c A_{k, m}(t) / P\left(t, T_{k}\right)+P\left(t, T_{k+m}\right) / P\left(t, T_{k}\right)
$$

the value of the coupon $c$ for which this expression equals 1 is known as the forward par rate or, when used in the context of swap pricing, as the forward swap rate. With $S_{k, m}(t)$ denoting the time $t$ swap rate, we apparently have

$$
\begin{equation*}
S_{k, m}(t)=\frac{P\left(t, T_{k}\right)-P\left(t, T_{k+m}\right)}{A_{k, m}(t)}, \quad t \leq T_{k} . \tag{4.10}
\end{equation*}
$$
^eq4-10

From the definition of $L\left(t, T_{n}, T_{n+1}\right)$ in (4.2), a little thought shows that the numerator of the expression for $S_{k, m}(t)$ can be expanded into a weighted sum of forward rates, leading to the alternative expression

$$
\begin{equation*}
S_{k, m}(t)=\frac{\sum_{n=k}^{k+m-1} \tau_{n} P\left(t, T_{n+1}\right) L_{n}(t)}{A_{k, m}(t)}, \quad t \leq T_{k} \tag{4.11}
\end{equation*}
$$
^eq4-11

where we have introduced the useful shorthand

$$
L_{n}(t) \triangleq L\left(t, T_{n}, T_{n+1}\right)
$$

It follows that the forward swap rate can be loosely interpreted as a weighted average of simple forward rates on the specified tenor structure. We note for the future that the time $T_{k}$ is sometimes referred to as the fixing date, or expiry, of the swap rate $S_{k, m}$, while the length of the corresponding swap, $T_{k+m}-T_{k}$, is sometimes called the tenor of the swap rate.

### 4.2 Fixed Income Probability Measures

As discussed in Chapter 1, selecting an equivalent martingale measure is largely a matter of choosing a numeraire, an asset price process used to
re-normalize the prices of other traded securities. For later reference, this section lists and names a number of important numeraires and measures used in fixed income pricing. Throughout the section, we assume that the market is complete, and we use $V(t)$ to denote the time $t$ price of a derivative security making an $\mathcal{F}_{T}$-measurable payment of $V(T)$.

### 4.2.1 Risk Neutral Measure

The numeraire defining the risk-neutral measure Q is the continuously compounded money market account $\beta(t)$, satisfying the locally deterministic SDE

$$
\begin{equation*}
d \beta(t)=r(t) \beta(t) d t, \quad \beta(0)=1 \tag{4.12}
\end{equation*}
$$
^eq4-12

where $r(t)$ is the short rate, $r(t)=f(t, t)$. Solving this equation yields

$$
\beta(t)=e^{\int_{0}^{t} r(u) d u}
$$

From the results of Chapter 1, in the absence of arbitrage the numerairedeflated process $V(t) / \beta(t)$ must be a martingale, implying the derivative security valuation formula

$$
\begin{equation*}
V(t) / \beta(t)=\mathrm{E}_{t}^{\mathrm{Q}}(V(T) / \beta(T)), \quad t \leq T, \tag{4.13}
\end{equation*}
$$
^eq4-13

or equivalently

$$
\begin{equation*}
V(t)=\mathrm{E}_{t}^{\mathrm{Q}}\left(e^{-\int_{t}^{\dot{T}} r(u) d u} V(T)\right) \tag{4.14}
\end{equation*}
$$
^eq4-14

If we apply (4.14) to the special case of $V(T)=1$, we obtain a fundamental bond pricing formula. We highlight the importance of this result by listing it in a lemma.

Lemma 4.2.1. In the absence of arbitrage, the time $t$ price $P(t, T)$ of a $T$-maturity zero-coupon bond is

$$
\begin{equation*}
P(t, T)=\mathrm{E}_{t}^{\mathrm{Q}}\left(e^{-\int_{t}^{T} r(u) d u}\right) \tag{4.15}
\end{equation*}
$$
^eq4-15

It follows from Lemma 4.2.1 that specification of the dynamics of $r(t)$ under $Q$ suffices to determine the prices of discount bonds at all times and maturities. Models that are based on such a direct specification of $r(t)$ dynamics are known as short rate models and are the subject of Chapters 10 through 12. Notice the resemblance between expressions (4.3) and (4.15); if $r(t)$ is deterministic, the two expressions will agree as $r(u)=f(t, u)$, $u \geq t$. If $r(t)$ is random, one may wonder whether this result will hold in expectation. The answer to this is negative, i.e.

$$
\begin{equation*}
f(t, u) \neq \mathrm{E}_{t}^{\mathrm{Q}}(r(u)) \tag{4.16}
\end{equation*}
$$
^eq4-16

provided $r$ is random. We prove this in the section below. Under certain idealized conditions, however, equality holds in (4.16) provided $f(t, u)$ is replaced by the futures rate $q(t, u)$. The exact result is as follows.

Lemma 4.2.2. Assume that mark-to-market takes place continuously. Under regularity conditions on the short rate $r(\cdot)$ - it suffices that $r(\cdot)$ is positive and bounded - the futures rate $F(\cdot, T, T+\tau)$ is a Q -martingale, and

$$
\begin{equation*}
F(t, T, T+\tau)=\mathrm{E}_{t}^{\mathrm{Q}}(L(T, T, T+\tau)) \tag{4.17}
\end{equation*}
$$
^eq4-17

Proof. Over a small interval $[t, t+d t]$, we have earlier shown that the cash proceeds from a futures contract are proportional to

$$
d F(t, T, T+\tau)=F(t+d t, T, T+\tau)-F(t, T, T+\tau)
$$

Suppose that we hold the futures contract up to some arbitrary horizon $t<\mathcal{T} \leq T$ at which point we exit (e.g., by selling the futures contract). Deflating all cash proceeds from this strategy with the numeraire $\beta(t)$ and integrating provides us with the time $t$ value of the futures contract as

$$
V_{\mathrm{fut}}(t)=\beta(t) \mathrm{E}_{t}^{\mathrm{Q}}\left(\int_{t}^{\mathcal{T}} \beta(s)^{-1} d F(s, T, T+\tau)+\beta(\mathcal{T})^{-1} V_{\mathrm{fut}}(\mathcal{T})\right)
$$

As it is always costless to enter into a futures contract, $V_{\text {fut }}(t)=V_{\text {fut }}(\mathcal{T})=0$ by definition, so for arbitrary $t<\mathcal{T} \leq T$ we must have (since $\beta(t)$ is positive)

$$
\begin{equation*}
\mathrm{E}_{t}^{\mathrm{Q}}\left(\int_{t}^{\mathcal{T}} \beta(s)^{-1} d F(s, T, T+\tau)\right)=0 \tag{4.18}
\end{equation*}
$$
^eq4-18

Provided that $\beta(s)^{-1}$ is almost surely positive (which is the case if $r$ is bounded), the fact that (4.18) holds for arbitrary horizons $t<\mathcal{T} \leq T$ shows that

$$
\mathrm{E}_{t}^{\mathrm{Q}}(d F(s, T, T+\tau))=0, \quad t \leq s \leq T
$$

which demonstrates that $F$ is a Q-martingale. The result (4.17) then immediately follows. $\square$

Equation (4.17) states that the futures rate is the Q-expectation of the time $T$ spot rate $L(T, T, T+\tau)$. A similar relation must then hold for the instantaneous futures rate, i.e.

$$
\begin{equation*}
q(t, u)=\mathrm{E}_{t}^{\mathrm{Q}}(r(u)) \tag{4.19}
\end{equation*}
$$
^eq4-19

as stated earlier².
Lemma 4.2.2 was first proven by non-probabilistic methods in Cox et al. [1981], who employed a direct, and quite instructive, hedging argument to

[^60]show the result. The assumption of continuous resettlement in the lemma may appear idealized, but the difference between daily and continuous settlement is quite small, as shall be demonstrated in Chapter 16. Explicit modeling of discrete resettlement is nevertheless quite straightforward, and basically involves shifting measure, from the risk-neutral measure to the so-called spot measure, defined below. We return to this issue in Chapter 16.

### 4.2.2 $T$-Forward Measure

The $T$-forward measure $Q^{T}$ was introduced in Jamshidian [1991b] (see also Geman et al. [1995]), and uses a $T$-maturity zero-coupon bond as the numeraire asset. As is customary, we let $\mathrm{E}^{T}(\cdot)$ denote expectations in measure $\mathrm{Q}^{T}$, such that

$$
V(t) / P(t, T)=\mathrm{E}_{t}^{T}(V(T) / P(T, T)), \quad t \leq T .
$$

As obviously $P(T, T)=1$, this expression simplifies to the convenient form

$$
\begin{equation*}
V(t)=P(t, T) \mathrm{E}_{t}^{T}(V(T)) \tag{4.20}
\end{equation*}
$$
^eq4-20

Comparison of (4.20) and (4.14) shows that shifting to the $T$-forward measure in a sense decouples the expectation of the terminal payout $V(T)$ from that of the numeraire. As we shall see, this is often very convenient when we attempt to construct analytical formulas for prices of certain simple interest rate derivatives. From the results of Section 1.3, we note that the explicit connection between the risk-neutral and $T$-forward measures is given by the density

$$
\begin{equation*}
\mathrm{E}_{t}^{\mathrm{Q}}\left(\frac{d \mathrm{Q}^{T}}{d \mathrm{Q}}\right)=\frac{P(t, T) / P(0, T)}{\beta(t)} \tag{4.21}
\end{equation*}
$$
^eq4-21

As $P(t, T+\tau)$ is the price of a traded asset, from the definition of the $T$-forward measure it follows that forward bond prices

$$
P(t, T, T+\tau)=P(t, T+\tau) / P(t, T)
$$

are martingales in the $T$-forward measure. We highlight a related result for forward rates below.

Lemma 4.2.3. In the absence of arbitrage the forward Libor rate $L(t, T, T+ \tau$ ) is a martingale under $\mathrm{Q}^{T+\tau}$, such that

$$
\begin{equation*}
L(t, T, T+\tau)=\mathrm{E}_{t}^{T+\tau}(L(T, T, T+\tau)), \quad t \leq T . \tag{4.22}
\end{equation*}
$$
^eq4-22

Proof. By definition (see (4.2))

$$
L(t, T, T+\tau)=\tau^{-1}(P(t, T) / P(t, T+\tau)-1)
$$

As $P(t, T) / P(t, T+\tau)$ is a martingale under $\mathrm{Q}^{T+\tau}$, so is $L(t, T, T+\tau)$. The result follows.

Taking the limit $\tau \downarrow 0$ and setting $T=u$ yields

$$
\begin{equation*}
f(t, u)=\mathrm{E}_{t}^{u}(f(u, u))=\mathrm{E}_{t}^{u}(r(u)), \tag{4.23}
\end{equation*}
$$
^eq4-23

which should be compared to the result (4.19).

### 4.2.3 Spot Measure

When working with a multitude of forward rates on a tenor structure $0=T_{0}<T_{1}<\ldots<T_{N}$, it is often convenient to introduce a numeraire that can be extended to arbitrary horizons by compounding. While the continuously compounded money market account $\beta$ would accomplish this, working with a continuously compounded numeraire is inherently awkward in a setting with a discrete tenor structure. As an alternative, we can introduce a discrete-time equivalent of the continuously compounded money market account to be the value of the following trading strategy. At time $0, \$ 1$ is invested in $1 / P\left(0, T_{1}\right) T_{1}$-maturity discount bonds, returning the amount

$$
1 / P\left(0, T_{1}\right)=1+\tau_{0} L\left(0,0, \tau_{0}\right)
$$

at time $T_{1}$. This amount is then reinvested ("rolled") at time $T_{1}$ in $T_{2}$ maturity bonds, returning

$$
1 / P\left(0, T_{1}\right) \cdot 1 / P\left(T_{1}, T_{2}\right)=\left(1+\tau_{0} L\left(0,0, T_{1}\right)\right)\left(1+\tau_{1} L\left(T_{1}, T_{1}, T_{2}\right)\right)
$$

at time $T_{2}$. Repeating this re-investment strategy at each date in the tenor structure gives rise to an asset price process $B(t)$, where $B(0)=1$ and

$$
\begin{equation*}
B(t)=\prod_{n=0}^{i}\left(1+\tau_{n} L_{n}\left(T_{n}\right)\right) P\left(t, T_{i+1}\right), \quad T_{i}<t \leq T_{i+1} \tag{4.24}
\end{equation*}
$$
^eq4-24

where we used the already introduced short-hand notation $L_{n}(t)= L\left(t, T_{n}, T_{n+1}\right)$. The process $B(t)$ is effectively a rolling certificate of deposit, and can be interpreted as a discrete-time equivalent of $\beta(t) . B(t)$ will approach $\beta(t)$ as the time spacing of the tenor structure is made increasingly fine.

The measure induced by $B(t)$ is known as the spot measure (or sometimes spot Libor measure), denoted $\mathrm{Q}^{B}$. With $\mathrm{E}^{B}(\cdot)$ denoting expectations in this measure we have

$$
V(t)=\mathrm{E}_{t}^{B}\left(V(T) \frac{B(t)}{B(T)}\right)
$$

where

$$
\begin{gathered}
\frac{B(t)}{B(T)}=\prod_{n=i+1}^{j}\left(1+\tau_{n} L_{n}\left(T_{n}\right)\right)^{-1} \frac{P\left(t, T_{i+1}\right)}{P\left(T, T_{j+1}\right)} \\
T_{i}<t \leq T_{i+1}, \quad T_{j}<T \leq T_{j+1}
\end{gathered}
$$

The similarity between the discrete and continuous money market accounts makes the spot Libor measure resemble the risk-neutral measure in many ways. For example, as we recall from Lemma 4.2.2, the risk-neutral measure is characterized by the fact that a continuously resettled futures rate is a martingale. In close parallel, the futures rate that is marked to market (resettled) discretely on dates $T_{0}, \ldots, T_{N}$ turns out to be a martingale in the spot Libor measure. We show this in Section 16.8.

### 4.2.4 Terminal and Hybrid Measures

One advantage of the spot measure over an arbitrary $T$-forward measure is the fact that the numeraire asset $B(t)$ will remain alive throughout the span of the tenor structure $\left\{T_{n}\right\}_{n=0}^{N}$. This property of $B(t)$ is necessary for the valuation of securities which may mature (randomly) at any date in the tenor structure. Securities of this type include for instance barrier options (such as range accruals) and options with early exercise rights (such as Bermudan swaptions). On the other hand, if we pick the $T$-forward measure corresponding to the last maturity in the tenor structure, $T=T_{N}$, this also yields a numeraire asset - the $T_{N}$-maturity zero-coupon bond - that is certain to remain alive at all dates in the tenor structure. The measure induced by $P\left(t, T_{N}\right)\left(\mathrm{Q}^{T_{N}}\right)$ is often referred to as the terminal measure. For a security $V$ maturing at a date $T \leq T_{N}$ we get, from the usual martingale pricing formula,

$$
\begin{equation*}
V(t)=P\left(t, T_{N}\right) \mathrm{E}_{t}^{T_{N}}\left(V(T) / P\left(T, T_{N}\right)\right), \quad t \leq T \leq T_{N} . \tag{4.25}
\end{equation*}
$$
^eq4-25

In (4.25) it is useful to notice that $V(T) / P\left(T, T_{N}\right)$ is the time $T_{N}$ proceeds of rolling at time $T$ the security payout $V(T)$ into a zero-coupon bond maturing at time $T_{N}$, effectively aligning the maturity of the numeraire and the cash flow date of the underlying asset. As an alternative, $V(T)$ could be rolled into the spot numeraire asset $B(T)$, leading to a $T_{N}$ payout of $V(T) / B(T) \cdot B\left(T_{N}\right)$. This gives rise to the equivalent formula

$$
\begin{equation*}
V(t)=P\left(t, T_{N}\right) \mathrm{E}_{t}^{T_{N}}\left(V(T) B\left(T_{N}\right) / B(T)\right), \quad t \leq T \leq T_{N} . \tag{4.26}
\end{equation*}
$$
^eq4-26

We note that this formula can also be derived from the basic relationship between the measures $\mathrm{Q}^{B}$ and $\mathrm{Q}^{T_{N}}$ by simply noting that

$$
P\left(T, T_{N}\right) \mathrm{E}_{T}^{T_{N}}\left(B\left(T_{N}\right) / B(T)\right)=B(T) \mathrm{E}_{T}^{B}\left(B\left(T_{N}\right) / B(T) / B\left(T_{N}\right)\right)=1,
$$

such that, by iterated conditional expectations ${ }^{3}$,[^61]

$$
\begin{aligned}
V(t) & =P\left(t, T_{N}\right) \mathrm{E}_{t}^{T_{N}}\left(V(T) / P\left(T, T_{N}\right)\right) \\
& =P\left(t, T_{N}\right) \mathrm{E}_{t}^{T_{N}}\left(V(T) \mathrm{E}_{T}^{T_{N}}\left(B\left(T_{N}\right) / B(T)\right)\right) \\
& =P\left(t, T_{N}\right) \mathrm{E}_{t}^{T_{N}}\left(V(T) B\left(T_{N}\right) / B(T)\right)
\end{aligned}
$$
as before.
As mentioned, equations (4.25) and (4.26) effectively involve reinvestment of the proceeds $V(T)$ to align cash payment with the numeraire $P\left(t, T_{N}\right)$. If the numeraire expires before the derivative security, we can apply the same reinvestment idea, but this time to the numeraire asset. Consider for instance a derivative security maturing at time $T_{N}$ (paying $V\left(T_{N}\right)$ ), and suppose we wish to extend the $T$-forward measure to price this option. For instance, we can define a numeraire asset as follows:

$$
\widetilde{P}(t, T)= \begin{cases}P(t, T), & t \leq T \\ B(t) / B(T), & t>T\end{cases}
$$

This asset corresponds to an investment strategy where we i) at time 0 purchase the $T$-maturity zero-coupon bond; and ii) at time $T$ invest the proceeds from the zero-coupon bond ( $\$ 1$ ) in the spot measure numeraire asset (4.24). Letting $\widetilde{\mathrm{Q}}^{T}$ denote the measure induced by $\widetilde{P}(t, T)$, we can write

$$
\begin{aligned}
V(t) & =\widetilde{P}(t, T) \widetilde{\mathrm{E}}_{t}^{T}\left(V\left(T_{N}\right) / \widetilde{P}\left(T_{N}, T\right)\right) \\
& =\widetilde{P}(t, T) \widetilde{\mathrm{E}}_{t}^{T}\left(V\left(T_{N}\right) B(T) / B\left(T_{N}\right)\right), \quad T<T_{N}, t<T_{N},
\end{aligned}
$$

where $\widetilde{\mathrm{E}}^{T}$ is the expectations operator for the measure $\widetilde{\mathrm{Q}}^{T}$. If also $t \leq T$, this expression becomes (compare to (4.26))

$$
V(t)=P(t, T) \widetilde{\mathrm{E}}_{t}^{T}\left(V\left(T_{N}\right) B(T) / B\left(T_{N}\right)\right), \quad t \leq T<T_{N},
$$

which, in effect, uses $B(T) / B\left(T_{N}\right)$ to discount $V\left(T_{N}\right)$ back to time $T$. The equivalent result in the $T$-forward measure is

$$
V(t)=P(t, T) \mathrm{E}_{t}^{T}\left(B(T) \mathrm{E}_{T}^{B}\left(V\left(T_{N}\right) / B\left(T_{N}\right)\right)\right)
$$

Notice that if $V$ matures at time $T$, rather than at $T_{N}$, we simply have

$$
V(t)=P(t, T) \widetilde{\mathrm{E}}_{t}^{T}(V(T))=P(t, T) \mathrm{E}_{t}^{T}(V(T))
$$

which is obvious from the definition of the numeraire $\widetilde{P}(t, T)$.
The measure $\widetilde{\mathrm{Q}}^{T}$ is by construction a hybrid between the spot measure and the $T$-forward measure. Obviously, many other such measures exist, corresponding to different reinvestment strategies of expiring numeraire assets.

### 4.2.5 Swap Measures

Being a linear combination of zero-coupon bonds (see (4.8)), an annuity factor $A_{k, m}(t)$ on a tenor structure qualifies as a numeraire asset. The measure $\mathrm{Q}^{k, m}$ induced by this numeraire is known as a swap measure or an annuity measure. In the absence of arbitrage we have

$$
V(t)=A_{k, m}(t) \mathrm{E}_{t}^{k, m}\left(V(T) / A_{k, m}(T)\right)
$$

where $\mathrm{E}^{k, m}(\cdot)$ denotes expectation under $\mathrm{Q}^{k, m}$.
Lemma 4.2.4. In the absence of arbitrage, the forward swap rate $S_{k, m}(t)$ is a martingale in measure $\mathrm{Q}^{k, m}$.

Proof. By definition

$$
S_{k, m}(t)=\frac{P\left(t, T_{k}\right)-P\left(t, T_{k+m}\right)}{A_{k, m}(t)}
$$

As the numeraire deflated assets $P\left(t, T_{k}\right) / A_{k, m}(t)$ and $P\left(t, T_{k+m}\right) / A_{k, m}(t)$ must both be martingales, so must be their difference.

As we shall see later, swap measures are very useful for analytical manipulations of price formulas for options on swaps.

### 4.3 Multi-Currency Markets

While this book is primarily dedicated to the study of single-currency interest rate derivatives, occasionally it will be necessary to consider certain effects associated with trading in a multi-currency economy. For instance, in Chapter 6 we touch upon issues of yield curve constructions in non-domestic currencies, and in Chapter 16 we discuss the important practical case where a derivative pays out in a foreign currency, but has a payout function that depends on one or more domestic interest rate variables. This brief section provides background material and notation required for these and other cross-currency applications.

### 4.3.1 Notations and FX Forwards

We consider two economies, a "domestic" economy and a "foreign" economy. Let $P_{d}(t, T)$ and $P_{f}(t, T)$ denote time $t$ zero-coupon bond prices in the domestic and foreign economies, respectively. So, $P_{f}(t, T)$ (say) is the time $t$ price, in foreign currency, of one unit of foreign currency delivered for certain at time $T$. Translation of values in foreign currency to domestic currency takes place at a foreign exchange (FX) rate of $X(t)$, measured in units of domestic currency per unit of foreign currency. In other words, the value $\widetilde{P}_{d}$ to a domestic investor of one foreign zero-coupon bond is

$$
\widetilde{P}_{d}(t, T)=X(t) P_{f}(t, T)
$$

The quantity

$$
X_{T}(t)=\frac{\widetilde{P}_{d}(t, T)}{P_{d}(t, T)}=X(t) \frac{P_{f}(t, T)}{P_{d}(t, T)}
$$

is known as the forward $F X$ rate to time $T$. Then name is motivated by the following arbitrage strategy:

- Buy one foreign zero-coupon bond, at a cost of $\widetilde{P}_{d}(t, T)$ in domestic currency.
- Finance the purchase by selling short domestic zero-coupon bonds on a notional of $\widetilde{P}_{d}(t, T) / P_{d}(t, T)$.
With no outlay at time $t$, the strategy will generate a net cash flow at time $T$ of one unit of foreign currency and $-X_{T}(t)$ units of domestic currency, such that the trading strategy in effect has locked in a time $t$ a future time $T$ exchange rate of $X_{T}(t)$.


### 4.3.2 Risk Neutral Measures

Let $\beta_{d}(t)$ and $\beta_{f}(t)$ be the continuously compounded money market accounts in the domestic and foreign economies, respectively. $\beta_{d}(t)$ and $\beta_{f}(t)$ induce two separate risk-neutral measures, denoted $\mathrm{Q}^{d}$ and $\mathrm{Q}^{f}$; let us investigate how these measures are related. If $g(T)$ is a random payout at time $T$ made in foreign currency, in a complete market the value (in units of foreign currency) of this payout to a foreign investor is, from standard principles,

$$
\begin{equation*}
V_{f}(t)=\beta_{f}(t) \mathrm{E}_{t}^{f}\left(g(T) \beta_{f}(T)^{-1}\right) \tag{4.27}
\end{equation*}
$$
^eq4-27

where $\mathrm{E}_{t}^{f}$ denotes expectations in the foreign risk-neutral measure $\mathrm{Q}^{f}$. For a domestic investor, the payout of $g(T)$ must be translated to domestic currency units at a rate of $X(T)$, making the effective domestic payout function $g(T) X(T)$. Thereby,

$$
\begin{equation*}
V_{d}(t)=\beta_{d}(t) \mathrm{E}_{t}^{d}\left(g(T) X(T) \beta_{d}(T)^{-1}\right) \tag{4.28}
\end{equation*}
$$
^eq4-28

where $\mathrm{E}_{t}^{d}$ denotes expectations in measure $\mathrm{Q}^{d}$. Importantly, the expressions in (4.27) and (4.28) are linked by the spot exchange rate, as the absence of a cross-currency arbitrage dictates that

$$
V_{d}(t)=X(t) V_{f}(t)
$$

or

$$
\begin{equation*}
\beta_{d}(t) \mathrm{E}_{t}^{d}\left(g(T) X(T) \beta_{d}(T)^{-1}\right)=X(t) \beta_{f}(t) \mathrm{E}_{t}^{f}\left(g(T) \beta_{f}(T)^{-1}\right) \tag{4.29}
\end{equation*}
$$
^eq4-29

We use this result to establish the following lemma:

Lemma 4.3.1. The domestic and foreign risk-neutral probability measures $\mathrm{Q}^{d}$ and $\mathrm{Q}^{f}$ are related by the density process

$$
\mathrm{E}_{t}^{d}\left(\frac{d \mathrm{Q}^{f}}{d \mathrm{Q}^{d}}\right)=\frac{\beta_{f}(t) X(t)}{\beta_{d}(t) X(0)}, \quad t \geq 0 .
$$

Proof. For an $\mathcal{F}_{T}$-measurable variable $Y(T)=g(T) X(T) \beta_{d}(T)^{-1}$ satisfying regularity conditions, a rearrangement of the basic relation (4.29) yields

$$
\mathrm{E}_{t}^{d}(Y(T))=X(t) \frac{\beta_{f}(t)}{\beta_{d}(t)} \mathrm{E}_{t}^{f}\left(\frac{Y(T)}{X(T)} \frac{\beta_{d}(T)}{\beta_{f}(T)}\right),
$$

From the results of Section 1.3, the density relating measures $\mathrm{Q}^{d}$ and $\mathrm{Q}^{f}$ is then as given in the lemma.

With $\beta_{f}(t) X(t) / \beta_{d}(t)$ being a martingale in the domestic risk-neutral measure, we note that if $X(t)$ is an Ito process, it must take the form

$$
d X(t)=X(t)\left(r_{d}(t)-r_{f}(t)\right) d t+\sigma_{X}(t)^{\top} d W(t)
$$

where $r_{d}(t)-r_{f}(t)$ is the spread between domestic and foreign short rates, $W(t)$ is a (vector-valued) $\mathrm{Q}^{d}$-Brownian motion, and $\sigma_{X}(t)$ is some adapted stochastic process satisfying regularity conditions.

### 4.3.3 Other Measures

Having established the Radon-Nikodym derivative relating the domestic and foreign risk-neutral measures, relations between various other domestic and foreign probability measures are easy to establish. For instance, the following result is easily proven the same way as Lemma 4.3.1.

Lemma 4.3.2. Let $\mathrm{E}_{t}^{T, d}$ denote expectations in the domestic $T$-forward probability measure. The domestic and foreign $T$-forward probability measures $\mathrm{Q}^{T, d}$ and $\mathrm{Q}^{T, f}$ are related by the density process

$$
\mathrm{E}_{t}^{T, d}\left(\frac{d \mathrm{Q}^{T, f}}{d \mathrm{Q}^{T, d}}\right)=\frac{P_{f}(t, T) P_{d}(0, T) X(t)}{P_{d}(t, T) P_{f}(0, T) X(0)}=\frac{X_{T}(t)}{X_{T}(0)}, \quad t \geq 0 .
$$

We highlight the fact that the forward FX rate $X_{T}(t)$ is a $\mathrm{Q}^{T, d}$-martingale satisfying $X_{T}(T)=X(T)$. For an $\mathcal{F}_{T}$-measurable variable $Y(T)$, we thereby have the convenient expression

$$
\mathrm{E}_{t}^{T, d}(Y(T))=X_{T}(t) \mathrm{E}_{t}^{T, f}\left(\frac{Y(T)}{X(T)}\right)
$$

where $\mathrm{E}_{t}^{T, f}$ denotes expectations in the foreign $T$-forward probability measure.

### 4.4 The HJM Analysis

Having defined notations and established basic arbitrage relationships, let us turn to assigning dynamics to the many quantities we have introduced so far. We shall here follow the Heath et al. [1992] (Heath-Jarrow-Morton, or HJM) approach, where all information in the economy is assumed to originate with a finite number of Brownian motions. The resulting class of models is quite broad, and in much of the rest of this book we shall deal with ways to reduce the general HJM model to specific, and tractable, special cases. For now, however, we concentrate on a general analysis, although we keep our treatment fairly informal.

### 4.4.1 Bond Price Dynamics

In the HJM framework, we concern ourselves with the modeling of how an entire continuum of $T$-indexed bond prices $P(\cdot, T)$ jointly evolves our time, starting from a known condition $P(0, T)$. We consider models of finite horizon, i.e. with $T \in[0, \mathcal{T}], \mathcal{T}<\infty$, and specialize to a filtration generated by a $d$-dimensional Brownian motion. We assume that a risk-neutral measure Q exists and is unique. Let $W(t)$ be an adapted $d$-dimensional Q-Brownian motion, and define deflated bond values as $P_{\beta}(t, T)=P(t, T) / \beta(t)$, where $\beta(t)$ as always is the continuously rolled money market account. In the absence of arbitrage, $P_{\beta}(t, T)$ is a martingale in the risk-neutral measure, and the martingale representation theorem then implies that

$$
\begin{equation*}
d P_{\beta}(t, T)=-P_{\beta}(t, T) \sigma_{P}(t, T)^{\top} d W(t), \quad t \leq T \tag{4.30}
\end{equation*}
$$
^eq4-30

where $\sigma_{P}(t, T)=\sigma_{P}(t, T, \omega)$ is a $d$-dimensional stochastic process adapted to the filtration generated by $W$. We assume that $\sigma_{P}(t, T)$ is regular enough for $P_{\beta}(t, T)$ to be a square-integrable martingale. Also, as the bond $P(t, T)$ must equal $\$ 1$ at $t=T$ ("pull to par"), we impose the consistency condition

$$
\sigma_{P}(T, T)=0
$$

Using (4.12) and Ito's lemma, it follows from (4.30) that

$$
\begin{equation*}
d P(t, T) / P(t, T)=r(t) d t-\sigma_{P}(t, T)^{\top} d W(t) \tag{4.31}
\end{equation*}
$$
^eq4-31

where $r(t)$ is the short rate process. Equation (4.31) defines the class of $d$-dimensional HJM models.

Another application of Ito's lemma shows that forward bond prices $P(t, T, T+\tau)=P(t, T+\tau) / P(t, T)$ must satisfy

$$
\begin{align*}
d P(t, T, T+\tau) / P(t, T, T+\tau) & =-\left[\sigma_{P}(t, T+\tau)-\sigma_{P}(t, T)\right]^{\top} \sigma_{P}(t, T) d t \\
& -\left[\sigma_{P}(t, T+\tau)-\sigma_{P}(t, T)\right]^{\top} d W(t) \tag{4.32}
\end{align*}
$$
^eq4-32

In the $T$-forward measure $\mathrm{Q}^{T}, P(t, T, T+\tau)$ is a martingale (see Section 4.2.2), and

$$
\begin{equation*}
d P(t, T, T+\tau) / P(t, T, T+\tau)=-\left[\sigma_{P}(t, T+\tau)-\sigma_{P}(t, T)\right]^{\top} d W^{T}(t) \tag{4.33}
\end{equation*}
$$
^eq4-33

where $W^{T}(t)$ is a $\mathrm{Q}^{T}$-Brownian motion. Comparison of (4.32) and (4.33) shows that

$$
\begin{equation*}
d W^{T}(t)=d W(t)+\sigma_{P}(t, T) d t \tag{4.34}
\end{equation*}
$$
^eq4-34

which by Girsanov's theorem identifies the density process for the measure shift between $\mathrm{Q}^{T}$ and Q in the HJM setting:

$$
\begin{equation*}
\varsigma(t)=\mathrm{E}_{t}^{\mathrm{Q}}\left(\frac{d \mathrm{Q}^{T}}{d \mathrm{Q}}\right)=\mathcal{E}\left(-\int_{0}^{t} \sigma_{P}(u, T)^{\top} d W(u)\right) \tag{4.35}
\end{equation*}
$$
^eq4-35

or

$$
d \varsigma(t) / \varsigma(t)=-\sigma_{P}(t, T)^{\top} d W(t)
$$

This result could, of course, have been established from the first principles as well - see equation (4.21).

### 4.4.2 Forward Rate Dynamics

Traditionally, HJM models are stated in terms of instantaneous forward rates, rather than bond prices. Besides eliminating the need to consider the short rate $r$, this also reveals a number of fundamental properties of the class of HJM models. By Ito's lemma, in measure Q,

$$
d \ln P(t, T)=O(d t)-\sigma_{P}(t, T)^{\top} d W(t)
$$

where for convenience we have omitted writing out the drift term. Differentiating the right- and left-hand sides of this equation with respect to $T$, we get from equation (4.4),

$$
d f(t, T)=\mu_{f}(t, T) d t+\sigma_{f}(t, T)^{\top} d W(t)
$$

where

$$
\begin{equation*}
\sigma_{f}(t, T)=\frac{\partial}{\partial T} \sigma_{P}(t, T) \tag{4.36}
\end{equation*}
$$
^eq4-36

and $\mu_{f}(t, T)$ is listed below.
Lemma 4.4.1. The process for $f(t, T)$ in the $T$-forward measure is

$$
\begin{equation*}
d f(t, T)=\sigma_{f}(t, T)^{\top} d W^{T}(t) \tag{4.37}
\end{equation*}
$$
^eq4-37

In the risk-neutral measure, the process is

$$
\begin{align*}
d f(t, T) & =\sigma_{f}(t, T)^{\top} \sigma_{P}(t, T) d t+\sigma_{f}(t, T)^{\top} d W(t) \\
& =\sigma_{f}(t, T)^{\top} \int_{t}^{T} \sigma_{f}(t, u) d u d t+\sigma_{f}(t, T)^{\top} d W(t) \tag{4.38}
\end{align*}
$$
^eq4-38

Proof. The SDE (4.37) follows directly from the martingale relation (4.23). The risk-neutral process (4.38) then can be derived from the relations (4.34) and (4.35), with the second equality following from (4.36). $\square$

The equation (4.38) is often considered to be the main result of Heath et al. [1992]. It demonstrates that an HJM model is fully specified once the forward rate diffusion coefficients $\sigma_{f}(t, T)$ have been specified for all $t$ and $T$. Note that HJM models take initial forward rates $f(0, T)$ as exogenous inputs, ensuring that these models are automatically consistent with discount bond prices at time 0 . This is true irrespective of the choice of $\sigma_{f}(t, T)$, which can be set freely (subject to regularity conditions) from either empirical analysis, or from a calibration to market prices of fixed income derivatives.

While it is convenient that HJM models are automatically calibrated to initial bond prices, a number of other features of the general HJM model are less attractive. Particularly problematic is the sheer dimensionality of the model: to describe the time $t$ state of a discount bond curve spanning $[t, T]$, we need to keep track of a continuum of forward rates $\{f(t, u), t \leq u \leq T\}$. By Lemma 4.4.1 the forward rate curve follows an infinite-dimensional diffusion process, leaving us with an infinite number of state variables to diffuse. In practice, the implementation of an HJM model will require either making special assumptions about the $\sigma_{f}$ process that permit a finite-dimensional Markovian representation of the forward rate curve; or moving from infinitesimal forward rates to continuously compounded forward rates that span time-buckets of finite length. Chapters 10 through 13 and Section 4.5.2 below give examples of the former idea, and Andersen [1995] discusses the latter approach in a Monte Carlo setting. An idea closely related to the discussion in Andersen [1995] is to built a model around a finite set of simple (Libor) forward rates on a fixed tenor structure. This approach has a number of computational and theoretical advantages, and is the subject of Chapter 14. For now, we note that any arbitrage-free interest rate model set in a filtration generated exclusively by Brownian motions must be a special case of an HJM model. In particular, any such model must correspond to a particular choice of $\sigma_{f}(t, T)$.

### 4.4.3 Short Rate Process

As discussed earlier, specification of a short rate process is, in principle, sufficient to completely specify a full yield curve model. In the HJM framework, it follows from (4.38) that the short rate $r(t)$ in measure Q is

$$
r(t)=f(t, t)=f(0, t)+\int_{0}^{t} \sigma_{f}(u, t)^{\top} \int_{u}^{t} \sigma_{f}(u, s) d s d u+\int_{0}^{t} \sigma_{f}(u, t)^{\top} d W(u)
$$

The process for $r(t)$ is generally not Markovian, as can be seen by focusing on the path-dependent term

$$
D(t)=\int_{0}^{t} \sigma_{f}(u, t)^{\top} d W(u)
$$

for which we must have

$$
\begin{align*}
D(T)=D(t)+ & \int_{t}^{T} \sigma_{f}(u, T)^{\top} d W(u) \\
& +\left\{\int_{0}^{t} \sigma_{f}(u, T)^{\top} d W(u)-\int_{0}^{t} \sigma_{f}(u, t)^{\top} d W(u)\right\} \tag{4.39}
\end{align*}
$$
^eq4-39

Thereby

$$
\mathrm{E}^{Q}(D(T) \mid D(t)) \neq \mathrm{E}_{t}^{Q}(D(T))
$$

unless the bracketed term in (4.39) is either non-random, or a deterministic function of $D(t)$ (which is generally not the case).

An interesting area of investigation concerns the conditions under which either $r(t)$ is outright Markov ${ }^{4}$ or, less restrictively, can be written as

$$
r(t)=h(t, x(t))
$$

for a deterministic function $h$ and a finite-dimensional Markovian vector of state variables $x(t)$. Definitive results are given in Björk [2001], building on earlier (and considerably less abstract) work by Jamshidian [1991b], Cheyette [1991], and Ritchken and Sankarasubramanian [1995]. Section 4.5.2 and Chapter 13 list some of the results of these papers.

### 4.5 Examples of HJM Models

### 4.5.1 The Gaussian Model

In the HJM bond price dynamics (4.31), we now assume that $\sigma_{P}(t, T)$ is a bounded ( $d$-dimensional) deterministic function of $t$ and $T$. It follows from (4.32) and (4.33) that forward bond prices are then log-normally distributed in both Q and $\mathrm{Q}^{T}$. The forward rate process in Q is

$$
\begin{equation*}
d f(t, T)=\sigma_{f}(t, T)^{\top} \sigma_{P}(t, T) d t+\sigma_{f}(t, T)^{\top} d W(t), \quad \sigma_{f}(t, T)=\frac{\partial}{\partial T} \sigma_{P}(t, T) \tag{4.40}
\end{equation*}
$$
^eq4-40

which implies that $r(T)=f(T, T)$ is Gaussian with Q-moments

$$
\begin{aligned}
\mathrm{E}_{t}^{\mathrm{Q}}(f(T, T)) & =\int_{t}^{T} \sigma_{f}(u, T)^{\top} \sigma_{P}(u, T) d u \\
\operatorname{Var}_{t}^{\mathrm{Q}}(f(T, T)) & =\int_{t}^{T} \sigma_{f}(u, T)^{\top} \sigma_{f}(u, T) d u
\end{aligned}
$$

[^62]The simple form of the Gaussian HJM model makes it quite tractable, permitting analytical price formulas for a number of European options and futures contracts ${ }^{5}$. While the Gaussian HJM model suffers from the drawback of allowing negative forward and spot rates, analytical results derived in the model are often very useful in gaining a deeper understanding of a given contract, even if ultimately a more realistic model will be required for serious pricing purposes. Indeed, results derived for the Gaussian HJM model can often be used as a starting point for development of closed-form approximations in other models; we shall see many examples of this later in the book.

For illustration, we list a few select analytical results below. More formulas can be found in numerous sources, including Chapter II in Andersen [1996], Jamshidian [1991b], and Jamshidian [1993], to name a few.

Proposition 4.5.1 (Option on Zero-Coupon Bond). Consider a European call option paying at maturity $T$ the amount

$$
V(T)=\left(P\left(T, T^{*}\right)-K\right)^{+}, \quad T^{*}>T
$$

In the Gaussian HJM model (4.40), we have

$$
\begin{equation*}
V(t)=P\left(t, T^{*}\right) \Phi\left(d_{+}\right)-P(t, T) K \Phi\left(d_{-}\right) \tag{4.41}
\end{equation*}
$$
^eq4-41

where

$$
\begin{aligned}
d_{ \pm} & =\frac{\ln \left(P\left(t, T^{*}\right) /(K P(t, T))\right) \pm v / 2}{\sqrt{v}} \\
v & =\int_{t}^{T}\left|\sigma_{P}\left(u, T^{*}\right)-\sigma_{P}(u, T)\right|^{2} d u
\end{aligned}
$$

Proof. In the $T$-forward measure $\mathrm{Q}^{T}$ we have, from (4.20),

$$
\begin{aligned}
V(t) & =P(t, T) \mathrm{E}_{t}^{T}\left(\left(P\left(T, T^{*}\right)-K\right)^{+}\right) \\
& =P(t, T) \mathrm{E}_{t}^{T}\left(\left(P\left(T, T, T^{*}\right)-K\right)^{+}\right)
\end{aligned}
$$

From the discussion in Section 4.4.1 we know that $P\left(t, T, T^{*}\right)$ is a $\mathrm{Q}^{T}$ martingale characterized by the SDE

$$
d P\left(t, T, T^{*}\right) / P\left(t, T, T^{*}\right)=-\left[\sigma_{P}\left(t, T^{*}\right)-\sigma_{P}(t, T)\right]^{\top} d W^{T}(t)
$$

where $W^{T}$ is a $d$-dimensional $Q^{T}$-Brownian motion. As this is just a GBM process with time-dependent coefficients, the Black-Scholes-Merton results in Section 1.9 apply and lead to (4.41).

[^63]Proposition 4.5.2 (Caplet). Consider a European call option paying at $T+\tau$ the amount (a caplet)

$$
V(T+\dot{\tau})=\tau(L(T, T, T+\tau)-K)^{+}, \quad \tau>0 .
$$

In the Gaussian HJM model (4.40), we have

$$
\begin{gathered}
V(t)=\tau P(t, T+\tau)\left(L(t, T, T+\tau) \Phi\left(d_{+}\right)-K \Phi\left(d_{-}\right)\right) \\
d_{ \pm}=\frac{\ln (L(t, T, T+\tau) / K) \pm v / 2}{\sqrt{v}}, v=\int_{t}^{T}\left|\sigma_{P}(u, T+\tau)-\sigma_{P}(u, T)\right|^{2} d u
\end{gathered}
$$

Proof. From Lemma 4.2.3 we know that $L(t, T, T+\tau)$ is a martingale in the $(T+\tau)$-forward measure. An application of Ito's lemma to the definition (4.2) reveals that

$$
d L(t, T, T+\tau) / L(t, T, T+\tau)=\left[\sigma_{P}(t, T+\tau)-\sigma_{P}(t, T)\right]^{\top} d W^{T}(t)
$$

and the result follows immediately along the same lines as in the proof of Proposition 4.5.1. $\square$

Proposition 4.5.3 (Futures Rate). In the Gaussian HJM model (4.40), futures rates are given by

$$
\begin{equation*}
F(t, T, T+\tau)=\tau^{-1}\left((1 / P(t, T, T+\tau)) e^{\Omega(t, T)}-1\right) \tag{4.42}
\end{equation*}
$$
^eq4-42

where

$$
\Omega(t, T)=\int_{t}^{T}\left[\sigma_{P}(u, T+\tau)-\sigma_{P}(u, T)\right]^{\top} \sigma_{P}(u, T+\tau) d u
$$

Proof. From Lemma 4.2.2,

$$
\begin{align*}
F(t, T, T+\tau) & =\mathrm{E}_{t}^{\mathrm{Q}}(L(T, T, T+\tau)) \\
& =\tau^{-1} \mathrm{E}_{t}^{\mathrm{Q}}(1 / P(T, T+\tau)-1) \\
& =\tau^{-1} \mathrm{E}_{t}^{\mathrm{Q}}(G(T)-1) \tag{4.43}
\end{align*}
$$
^eq4-43

where we have introduced an auxiliary variable

$$
G(t) \triangleq P(t, T) / P(t, T+\tau)=1 / P(t, T, T+\tau) .
$$

Ito's lemma shows that (see also (4.32)) in measure Q

$$
\begin{aligned}
d G(t) / G(t) & =\left[\sigma_{P}(t, T+\tau)-\sigma_{P}(t, T)\right]^{\top} \sigma_{P}(t, T+\tau) d t \\
& +\left[\sigma_{P}(t, T+\tau)-\sigma_{P}(t, T)\right]^{\top} d W(t)
\end{aligned}
$$

such that

$$
\mathrm{E}_{t}^{\mathrm{Q}}(G(T))=G(t) e^{\Omega(t, T)}=(1 / P(t, T, T+\tau)) e^{\Omega(t, T)},
$$

where $\Omega$ is as given above. The result of Proposition 4.5.3 then follows directly from (4.43).

In any rational model $\Omega(t, T) \geq 0$, such that $F(t, T, T+\tau) \geq L(t, T, T+ \tau)$, consistent with the qualitative discussion in Section 4.1.2. As shown in Chapter II of Andersen [1996], the spread (also known as futures convexity) between futures and forward rates can be decomposed into two components: i) a term originating from the mark-to-market mechanism of a futures contract; and ii) a term originating from the fact that a futures contract unlike a regular forward rate agreement - pays out the rate at the date it settles (at time $T$ ) rather than one period ahead (at time $T+\tau$ ). Andersen [1996], Chapter II, additionally contains a number of numerical examples examining typical futures-forward spreads, and also investigates the pricing of options on futures rates.

Section 16.8 looks in detail into pricing interest rate futures under more advanced models.

### 4.5.2 Gaussian HJM Models with Markovian Short Rate

Although quite tractable, the Gaussian HJM model generally does not allow for a finite-dimensional Markovian representation, and typically does not imply Markov-diffusive behavior of the short rate. As shown in Carverhill [1994], the short rate can be made Markovian, however, by imposing certain conditions on the deterministic forward rate volatility function $\sigma_{f}(t, T)$. To explore this, first recall from Section 4.4.3 the relation

$$
r(t)=f(0, t)+\int_{0}^{t} \sigma_{f}(u, t)^{\top} \int_{u}^{t} \sigma_{f}(u, s) d s d u+\int_{0}^{t} \sigma_{f}(u, t)^{\top} d W(u)
$$

where now $\sigma_{f}$ is deterministic. Consider imposing the special choice

$$
\begin{equation*}
\sigma_{f}(t, T)=g(t) h(T) \tag{4.44}
\end{equation*}
$$
^eq4-44

where $h$ is a positive real function and $g: \mathbb{R} \rightarrow \mathbb{R}^{d \times 1}$ can take any sign. For this case we have

$$
\sigma_{P}(t, T)=\int_{t}^{T} \sigma_{f}(t, u) d u=g(t) \int_{t}^{T} h(u) d u
$$

and

$$
\begin{align*}
r(t) & =f(0, t)+h(t) \int_{0}^{t} g(u)^{\top} g(u)\left(\int_{u}^{t} h(s) d s\right) d u+h(t) \int_{0}^{t} g(u)^{\top} d W(u)  \tag{4.45}\\
& \triangleq f(0, t)+h(t) \int_{0}^{t} m_{f}(t, u) d u+h(t) \int_{0}^{t} g(u)^{\top} d W(u)
\end{align*}
$$
^eq4-45

Importantly, the term

$$
D(t)=\int_{0}^{t} \sigma_{f}(u, t)^{\top} d W(u)=h(t) \int_{0}^{t} g(u)^{\top} d W(u)
$$

is now Markov, since

$$
D(T)=h(T) \int_{0}^{T} g(u)^{\top} d W(u)=\frac{h(T)}{h(t)} D(t)+h(T) \int_{t}^{T} g(u)^{\top} d W(u)
$$

which should be compared to the general (non-Markov) expression (4.39).
To show that the short rate is Markovian, we differentiate (4.45) with respect to $t$, yielding

$$
\begin{align*}
d r(t)= & \frac{\partial f(0, t)}{\partial t} d t+h^{\prime}(t)\left(\int_{0}^{t} m_{f}(t, u) d u+\frac{D(t)}{h(t)}\right) d t \\
& +h(t) \frac{\partial}{\partial t} \int_{0}^{t} m_{f}(t, u) d u d t+h(t) g(t)^{\top} d W(t) \\
= & \frac{\partial f(0, t)}{\partial t} d t+h^{\prime}(t)\left(\frac{r(t)-f(0, t)}{h(t)}\right) d t \\
& +h(t)^{2} \int_{0}^{t} g(u)^{\top} g(u) d u d t+h(t) g(t)^{\top} d W(t) \\
= & \left(\frac{\partial f(0, t)}{\partial t}-\frac{h^{\prime}(t)}{h(t)} f(0, t)+h(t)^{2} \int_{0}^{t} g(u)^{\top} g(u) d u+\frac{h^{\prime}(t)}{h(t)} r(t)\right) d t \\
& +h(t) g(t)^{\top} d W(t) \tag{4.46}
\end{align*}
$$
^eq4-46

where the second equality follows from rearrangement of (4.45). This leads to the following result.

Proposition 4.5.4. In the $d$-dimensional Gaussian HJM model, when (4.44) holds the short rate satisfies an SDE of the type

$$
d r(t)=(a(t)-\varkappa(t) r(t)) d t+\sigma_{r}(t)^{\top} d W(t)
$$

where $\varkappa: \mathbb{R} \rightarrow \mathbb{R}$ and $\sigma_{r}: \mathbb{R} \rightarrow \mathbb{R}^{d \times 1}$ are deterministic functions of time, and

$$
\begin{aligned}
a(t) & =\frac{\partial f(0, t)}{\partial t}+\varkappa(t) f(0, t)+\int_{0}^{t} e^{-2 \int_{a}^{t} \varkappa(s) d s} \sigma_{r}(u)^{\top} \sigma_{r}(u) d u \\
& =\frac{\partial f(0, t)}{\partial t}+\varkappa(t) f(0, t)+\int_{0}^{t} \sigma_{f}(u, t)^{\top} \sigma_{f}(u, t) d u
\end{aligned}
$$

Proof. First, by way of defining $\varkappa$ and $\sigma_{r}$, we set

$$
h(T)=e^{-\int_{0}^{T} \varkappa(s) d s} ; \quad g(t)=e^{\int_{0}^{t} \varkappa(s) d s} \sigma_{r}(t)
$$

such that $h^{\prime}(t) / h(t)=-\varkappa(t)$ and

$$
\sigma_{f}(t, T)=e^{-\int_{t}^{T} \varkappa(s) d s} \sigma_{r}(t)
$$

The result of Proposition 4.5.4 then follows directly by insertion into (4.46).

### 4.5.3 Log-Normal HJM Models

To avoid the negative forward rates inherent in Gaussian HJM models, it is tempting to consider forward rate specifications of the type

$$
\begin{equation*}
\sigma_{f}(t, T)=f(t, T) \sigma(t, T) \tag{4.47}
\end{equation*}
$$
^eq4-47

where $\sigma(t, T)$ is deterministic and bounded. In the $T$-forward measure

$$
d f(t, T)=f(t, T) \sigma(t, T)^{\top} d W^{T}(t)
$$

such that $f(t, T)$ is log-normally distributed. While avoiding negative rates, the specification (4.47) has severe technical problems: in Q , forward rates will explode to infinity with non-zero probability. Attempts to apply the valuation formula (4.15) will thereby result in all zero-coupon bond prices being zero, implying obvious arbitrage opportunities. To suggest a rationale for the exploding rates, consider the Q-dynamics
$d f(t, T)=\left(f(t, T) \sigma(t, T)^{\top} \int_{t}^{T} f(t, u) \sigma(t, u) d u\right) d t+f(t, T) \sigma(t, T)^{\top} d W(t)$.
Loosely speaking, the drift-term is proportional to forward rates squared, which, in the light of the linear growth condition in Theorem 1.6.1, may cause us to suspect problems with the existence of a non-exploding solution. Morton [1988] confirms this rigorously.

One solution to the explosion problem involves enforcing a strict upper bound on $\sigma_{f}(t, T)$, as in

$$
\sigma_{f}(t, T)=\min (f(t, T), M) \sigma(t, T),
$$

where $M$ is a large positive constant. For the one-factor case ( $d=1$ ) Heath et al. [1992] demonstrate that this specification will ensure non-negative forward rates ${ }^{6}$ and will prevent rate explosions. Nevertheless, the model is clearly awkward in its dependence on the arbitrary constant $M$. A more satisfying solution is discussed in Chapter 14, where we show that the explosion problem can be circumvented by working with simply - rather than continuously - compounded forward rates. A related issue in short rate models is also discussed in Chapter 11.

[^64]
