---
tags: [quant/knowledge, quant/books]
---
# Front and Back Matter — Lorenzo Bergomi, *Stochastic Volatility Modeling*

[Book Index](index.md)

This file preserves the Preface, Epilogue, and Bibliography of the book — content that lived only in the original `WholeBook.md` monolith and `epilogue, footnote.md`, and is not reproduced in the split `Ch01.md`–`Ch12.md` files or in `index.md`. Extracted 2026-07-06 prior to deleting those two files, once the 12-chapter rewrite of `index.md` was complete and fully validated.

Footnotes [^0] and [^1] below belong to the Preface specifically (all other footnotes, [^2]–[^174], are already inlined within their respective chapter files and are not repeated here).

---

## Preface

> Pour soulever un poids si lourd Sisyphe, il faudrait ton courage! Bien qu'on ait du cœur à l'ouvrage, L'Art est long et le Temps est court.

Baudelaire, Les Fleurs du mal

Tu quid? Quousque sub alio moveris? Impera et dic quod memoriae tradatur, aliquid et de tuo profer.

Seneca, Letters to Lucilius, XXXIII

> C'est icy un livre de bonne foy, lecteur.

Montaigne, Essais

- This, Reader, is an honest book. It warns you at the outset that it is not a treatise on stochastic volatility.

Nor is it a mathematical finance textbook - there are treatises and textbooks galore on the shelves of bookshops and university libraries.

Rather, my intention has been to explain how stochastic volatility - and which kind of stochastic volatility - can be used to address practical issues arising in the modeling of derivatives.

Modeling in finance is an engineering field: while our task as engineers is to frame problems in mathematical terms and solve them using sophisticated machinery whenever necessary, the problems themselves originate in the form of embarrassingly practical trading questions.

I have been fortunate to have spent my career as a quant in an institution - Société Générale - and in a field - equity derivatives - that the derivatives community has come to associate. I have tried to convey the experience thus gained. My main objective has been to clearly state the motivation for each question I address and to represent the thought process one follows when designing or using a model.

A model's (ir)relevance is measured by its (in)ability to account for all nonlinearities of a derivative payoff with respect to hedging instruments, and to adequately reflect the former in the prices it produces.

Thus have I often lamented the propensity of "pricing quants" or "expectation calculators" to envision their task as that of producing (real) numbers auspiciously called prices, and to favor analytical tractability or computational speed at the expense of model relevance. What is the point in ultrafast mispricing? [^0]

This book is intended to be read in sequence. It assumes familiarity with the basic concepts and models of quantitative finance. The motivation for stochastic volatility is the subject of the first chapter; this is followed by a chapter on local volatility, both a special breed of stochastic volatility and a market model, which I survey as such.

I urge you, Reader, to read it fully, as I introduce notions and discuss issues that are referred to repeatedly throughout the book. After a warm-up with forward-start options we embark on stochastic volatility. Jump and Lévy models are briefly dealt with at the end of Chapter 10.

This work was written mostly at night, during hours normally devoted to rest and sleep. Thus, dear Reader, I beg your forgiveness. You will do me a great service by reporting typos, inaccuracies, downright errors, to lb.svbook@gmail.com.

For everything I have learned so far I have to be thankful to more people than my memory can remember. I wish to thank the practitioners and academics that I have met regularly at conferences and seminars, and also fellow workers at Société Générale: quants, traders, structurers.

I am especially indebted to my coworkers in the Quantitative Research team, past and present. They will find here a reflection of the very many discussions we have had over the years.

The generous help of colleagues - quants and traders - in proofreading the manuscript is also gratefully acknowledged, chief among them Florent Bersani and Pierre Henry-Labordère whose eagle eyes have helped clear the text of many imprecisions.

Finally, I would like to express my gratitude to Rena, Elisa and Chiara for their patience.

Thank you Rena for occasionally filling the loneliness of nightly writing sessions with music. This has resulted for the author in some unlikely pairings, e.g. local volatility and the piano part of the Kreutzer sonata, the Heston model and the Chopin Barcarolle, Rachmaninoff's Third and options on variance, etc.

Paris, Summer 2015

---

## Epilogue

> Bornons ici cette carrière. Les longs ouvrages me font peur. Loin d'épuiser une matière, On n'en doit prendre que la fleur.

La Fontaine, Fables, VI
Et pour ceux qui joignent le bon sens avec l'étude, lesquels seuls je souhaite pour mes juges [...]

Descartes, Discours de la méthode
Quant au mouvement en lui-même, je vous le déclare avec humilité, nous sommes impuissants à le définir.

Balzac, La Peau de chagrin

- This excursion, Reader, ends here. I hope the journey has been pleasant and instructive.

We attempt to build market models for hedging instruments - the underlying along with vanilla options, or a subset of them, or variance swaps, ot yet other convex payoffs - such that (a) the P&L of a hedged option position is of the typical gamma/theta form, with payoff-independent break-even levels, (b) these break-even levels are - at least partially - in our control.

We could consider higher-order contributions but, in practice, managing the risk of all second-order greeks at a book level is a formidable task already.

Vanilla options are not quite independent instruments, however; they are related to one another and to the underlying itself through their terminal condition. At order one in volatility of volatility, for example, the SSR of homogeneous stochastic volatility models is related to the decay of the AMTF skew.

Characterizing the restrictions that the initial configuration of vanilla implied volatilities places on their future evolution is an unsettled issue, that may be more easily addressed by modeling other convex payoffs, for example the power payoffs of Chapter 4.

Even achieving objective (a) is not as straightforward, as highlighted by our analysis of mixed models in Chapter 12. Expressing the instantaneous volatility as the product of local and stochastic volatility components, a seemingly reasonable and innocuous ansatz, generally leads to non-functional models, for lack of a proper breakdown of the carry P&L.

This serves as a reminder that derivatives modeling does not start with the assumption of a process for $S_{t}$, but with clearly articulated modeling objectives relating to observable quantities, regardless of how these objectives are achieved mathematically.

A commonplace statement one has been hearing at recent quant conferences is that the age of modeling is now over, that quantitative finance has become but a tedious form of accounting, the preserve of a new order of adjusters - xxA-quants: CVA-quants, FVA-quants, KVA-quants, etc.

Focusing on derivatives' risks, on their proper modeling and the meaningfulness of model-generated prices - the things being adjusted - is the symptom of an old-fashioned, outdated, mindset.

I respectfully dissent. Product risks and modeling choices are still begging for a proper understanding; forty years after Black-Scholes, work on the next generation of models has just started.

What about calibration, a deceptive notion we should strive to abolish? No one says that they "calibrate" the spot value in the Black-Scholes model: a model should naturally take as inputs the market values of instruments used as hedges. Then, parameters should be chosen so as to generate the desired break-even levels for the carry P&L.

Should they be calibrated to market prices? This is meaningful only if the difference between their calibrated and realized values can be materialized as the P&L of an actual trading strategy, a rare occurrence.

Therefore, Reader, you will do well to resist the compulsion of calibration and the addictive psychological reward that comes with it.

With these last words of encouragement, Reader, I bid you farewell.

---

## Bibliography

[1] Andreasen, J., Huge B.: Random grids, Risk Magazine, July, pp. 66-71, 2011.
[2] Andersen, L., Andreasen, J.: Jump-Diffusion Processes: Volatility Smile Fitting and Numerical Methods for Pricing, available at SSRN: http://ssrn.com/abstract=171438.
[3] Avellaneda, M., Levy, A., Paras, A.: Pricing and hedging derivative securities in markets with uncertain volatilities, Applied Mathematical Finance 2(2), pp. 73-88, 1995.
[4] Avellaneda, M., Paras, A.:Pricing and hedging derivative securities in markets with uncertain volatilities: the Lagrangian uncertain volatility model, Applied Mathematical Finance, 3, pp. 21-52, 1996.
[5] Backus, D., Foresi, S., Wu, L.: Accounting for biases in Black-Scholes, available at SSRN: http://ssrn.com/abstract=585623, 1997.
[6] Berestycki, H., Busca, J., Florent, I.: Asymptotics and calibration of local volatility models, Quantitative Finance, 2(1), pp. 61-69, 2002.
[7] Berestycki, H., Busca, J., Florent, I.: Computing the implied volatility in stochastic volatility models, Communications on Pure and Applied Mathematics, Vol. LVII, pp. 1352-1373, 2004.
[8] Bergomi, L.: Smile dynamics, Risk Magazine, September, pp. 117-123, 2004. Also available at SSRN: http://ssrn.com/abstract=1493294.
[9] Bergomi, L.: Smile dynamics II, Risk Magazine, October, pp. 67-73, 2005. Also available at SSRN: http://ssrn.com/abstract=1493302.
[10] Bergomi, L.: Smile dynamics III, Risk Magazine, October, pp. 90-96, 2008. Also available at SSRN: http://ssrn.com/abstract=1493308.
[11] Bergomi, L.: Smile dynamics IV, Risk Magazine, December, pp. 94-100, 2009. Also available at SSRN: http://ssrn.com/abstract=1520443.
[12] Bergomi, L.: Correlations in asynchronous markets, Risk Magazine, November, pp. 76-82, 2010. Also available at SSRN: http://ssrn.com/abstract=1635866.
[13] Bergomi, L., Guyon, J.: Stochastic volatility's orderly smiles, Risk Magazine, May, pp. 60-66, 2012. Also available at SSRN: http://ssrn.com/abstract=1967470.
[14] Bick, A.: Quadratic-variation-based dynamic strategies, Management Science, 41(4), pp. 722-732, 1995.
[15] Björk, T., Blix, M., Landén, C.: On finite dimensional realizations for the term structure of futures prices, International Journal of Theoretical and Applied Finance, 9(03), pp. 281-314, 2006.
[16] Bos, M., Vandermark, S.: Finessing fixed dividends, Risk Magazine, September, pp. 157-158, 2002.
[17] Bouchaud, J. Ph., Potters, M., Sestovic, D.: Hedge your Monte Carlo, Risk Magazine, March, pp. 133-136, 2001. Also available at SSRN: http://ssrn.com/abstract=238868.
[18] Breeden, D., Litzenberger, R.: Prices of state contingent claims implicit in option prices, Journal of Business, 51, pp. 621-651, 1978.
[19] Buehler, H.: Volatility and dividends - volatility modeling with cash dividends and simple credit risk, available at SSRN: http://ssrn.com/abstract=1141877.
[20] Buehler, H.: Consistent variance curve models, Finance and Stochastics, 10(2), pp. 178-203, 2006. Also available at SSRN: http://ssrn.com/abstract=687258.
[21] Carmona, R., Nadtochiy, S.: Local volatility dynamic models, Finance and Stochastics (13), pp. 1-48, 2009.
[22] Carr, P., Chou, A.: Breaking barriers, Risk Magazine, September, pp. 139-145, 1997.
[23] Carr, P., Madan, D.: Determining volatility surfaces and option values from an implied volatility smile, Quantitative Analysis in Financial Markets, Vol II, M. Avellaneda, ed, pp. 163-191, 1998.
[24] Carr, P., Lewis, K.: Corridor variance swaps, Risk Magazine, February, pp. 67-72, 2004.
[25] Carr, P., Madan, D.: Towards a theory of volatility trading, Volatility, Risk Publications, Robert Jarrow, ed., pp. 417-427, 2002.
[26] Carr, P., Geman, H., Madan, D., Yor, M.: Stochastic volatility for Lévy processes, Mathematical Finance 13(3), pp. 345-382, 2003.
[27] Carr, P., Lee, R.: Robust replication of volatility derivatives, available at: http://www.math.uchicago.edu/~rl/rrvd.pdf, 2009.
[28] Carr, P., Lee, R. Hedging variance options on continuous semimartingales, Finance and Stochastics, 14 (2), pp. 179-207, 2010.
[29] Castagna, A., Mercurio, F.: The vanna-volga method for implied volatilities, Risk Magazine, January, pp. 106-111, 2007.
[30] Cheyette, O.: Markov representation of the Heath-farrow-Morton model, available at SSRN: http://ssrn.com/abstract=6073.
[31] Cherny, A., Dupire, B.: On certain distributions associated with the range of martingales, Optimality and Risk - Modern Trends in Mathematical Finance, pp. 29-38, 2010.
[32] Chriss, N., Morokoff, W.: Market risk of variance swaps, Risk Magazine, October, pp. 55-59, 1999.
[33] Cont, R., Tankov, P.: Financial modeling with jump processes, Chapman & Hall / CRC Press, Financial Mathematics Series, 2003.
[34] Cox, A., Wang, J.: Optimal robust bounds for variance options, http://arxiv.org/abs/1308.4363, 2013.
[35] Davis, M. H., Panas, V. G., Zariphopoulou, T.: European option pricing with transaction costs, SIAM Journal on Control and Optimization, 31(2), pp. 470-493, 1993.
[36] De Marco, S., Henry-Labordère, P.: Linking vanillas and VIX options: a constrained martingale optimal transport problem, available at SSRN: http://ssrn.com/abstract=2354898.
[37] Derman, E., Kani, I.: Riding on a smile, Risk Magazine, February, pp. 32-39, 1994.
[38] Duanmu, Z.: Rational pricing of options on realized volatility - the Black Scholes way, Global Derivatives conference, Madrid, 2004.
[39] Derman, E., Kani, I.: Stochastic implied trees: arbitrage pricing with stochastic term and strike structure of volatility, International Journal of Theoretical and Applied Finance, 1(01), pp. 61-110, 1998.
[40] Dupire, B.: Pricing with a smile, Risk Magazine, January, pp. 18-20, 1994.
[41] Dupire, B.: Functional Ito calculus and volatility hedge, Global Derivatives conference, Paris, 2009.
[42] Dupire, B.: Arbitrage bounds for volatility derivatives as free boundary problem, available at http://www.math.kth.se/pde_finance/presentations/Bruno.pdf, 2005.
[43] Durrleman, V.: From implied to spot volatilities, Finance and Stochastics, 14 (4), pp. 157-177, 2010.
[44] Durrleman, V., El Karoui, N.: Coupling smiles, Quantitative Finance, 8(6), pp. 573-590, 2008. Also available at SSRN: http://ssrn.com/abstract=1005332.
[45] Fisher, T., Tataru, G.: Non-parametric stochastic local volatility modeling, Global Derivatives conference, Paris, 2010.
[46] Fukasawa, M.: The normalizing transformation of the implied volatility smile, Mathematical Finance, 22(4), pp. 753-762, 2012. Preprint version available at http://arxiv.org/abs/1008.5055, 2010.
[47] Galichon, A., Henry-Labordère, P., Touzi, N.: A stochastic control approach to no-arbitrage bounds given marginals, with an application to lookback options, Annals of Applied Probability, 24(1), pp. 312-336, 2014.
[48] Gatheral, J.: The volatility surface: a practitioner's guide, Wiley Finance, 2006.
[49] Gatheral, J., Jacquier, A.: Arbitrage-free SVI volatility surfaces, Quantitative Finance, 14(1), pp. 59-71, 2014. Also available at SSRN: http://ssrn.com/abstract=2033323.
[50] Green, R.C., Jarrow, R.A.: Spanning and completeness in markets with contingent claims, Journal of Economic Theory, 41, pp. 202-210, 1987.
[51] Guyon, J., Henry-Labordère, P.: From spot volatilities to implied volatilities, Risk Magazine, June, pp. 79-84, 2011. Also available at SSRN: http://ssrn.com/abstract=1663878.
[52] Guyon, J., Henry-Labordère, P.: Being particular about calibration, Risk Magazine, January, pp. 92-97, 2012.
[53] Guyon, J., Henry-Labordère, P.: Nonlinear option pricing, Chapman & Hall/CRC Press, Financial Mathematics Series, 2013.
[54] Gyöngi, I.: Mimicking the one-dimensional marginal distributions of processes having an Ito differential, Probability Theory and Related Fields, 71, pp. 501516, 1986.
[55] Hagan, P.S., Kumar, D., Lesniewski, A.S., Woodward, D.E.: Managing smile risk, Wilmott Magazine, September, pp. 84-108, 2002.
[56] Henry-Labordère, P.: Analysis, geometry, and modeling in finance: advanced methods in option pricing, Chapman & Hall/CRC Press, Financial Mathematics Series, 2008.
[57] Henry-Labordère, P.: Automated option pricing: numerical methods, available at SSRN: http://ssrn.com/abstract=1968344, 2011.
[58] Henry-Labordère, P.: Calibration of local stochastic volatility models to market smiles: a Monte-Carlo approach, Risk Magazine, September, pp. 113-117, 2009. Also available at SSRN: http://ssrn.com/abstract=1493306.
[59] Henry-Labordère, P.: Vega decomposition of exotics on vanillas: a Monte-Carlo approach, available at SSRN: http://ssrn.com/abstract=2229990.
[60] Heston, S.: A closed-form solution for options with stochastic volatility with applications to bond and currency options, Review of Financial Studies, 6(2), pp. 327-343, 1993.
[61] Higham, N. J.: Computing the nearest correlation matrix - a problem from finance, IMA Journal of Numerical Analysis, 22(3), pp. 329-343, 2002.
[62] Hobson, D.: The Skorokhod embedding problem and model-independent bounds for option prices, Paris-Princeton Lectures on Mathematical Finance, pp. 267318, 2010.
[63] Hobson, D., Klimmek, M.: Model-independent hedging strategies for variance swaps, Finance and Stochastics, 16(4), pp. 611-649, 2012.
[64] Hunt, P., Kennedy, J., Pelsser, A.: Markov-functional interest rate models, Finance and Stochastics, 4(4), pp. 391-408, 2000.
[65] Lee, R.: Weighted variance swap, Encyclopedia of Quantitative Finance, 2010. Also available at http://math.uchicago.edu/~rl/EQF_weightedvarianceswap.pdf.
[66] Lee, R.: Implied and local volatilities under stochastic volatility, International Journal of Theoretical and Applied Finance, 4(01), pp. 45-89, 2001.
[67] Lee, R.: The moment formula for implied volatility at extreme strikes, Mathematical Finance, 14(3), pp. 469-480, 2004.
[68] Leland, H. E.: Option pricing and replication with transaction costs, The Journal of Finance, 40(5), pp. 1283-1301, 1985.
[69] Lewis, A.: Option valuation under stochastic volatility, Finance Press, 2000.
[70] Lipton, A.: The vol smile problem, Risk Magazine, February, pp. 61-65, 2002.
[71] Lyons, T.: Uncertain volatility and the risk-free synthesis of derivatives, Applied Mathematical Finance, 2(2), pp. 117-133, 1995.
[72] Marcinkiewicz, J.: Sur une propriété de la loi de Gauss, Mathematische Zeitschrift, 44, pp. 612-618, 1939.
[73] Matytsin, A.: Perturbative analysis of volatility smiles, Columbia Practitioners Conference on the Mathematics of Finance, New York, 2000.
[74] Nachman, D.: Spanning and completeness with options, Review of Financial Studies, 3, 31, pp. 311-328, 1988.
[75] Neuberger, A.: The log contract, Journal of Portfolio Management 20, 1994.
[76] Nicolay, D.: Asymptotic Chaos Expansions in Finance: Theory and Practice, Springer Finance, Springer Finance Lecture Notes, 2014.
[77] Papanicolaou, A.: Extreme-strike comparisons and structural bounds for SPX and VIX Options, available at SSRN: http://ssrn.com/abstract=2532020.
[78] Piterbarg, V.: Time to smile, Risk Magazine, May, pp. 71-75, 2005.
[79] Romano, M., Touzi, N.: Contingent claims and market completeness in a stochastic volatility model, Mathematical Finance, 7(4), pp. 399-412, 1997.
[80] Schönbucher, P. J.: A market model for stochastic implied volatility, Phil. Trans. Roy. Soc., Ser. A 357, pp. 2071-2092, 1999.
[81] Schweizer, M., Wissel, J.: Term structures of implied volatilities: absence of arbitrage and existence results, Mathematical Finance, 18(1), pp. 77-114, 2008.
[82] Sepp, A.: Efficient Numerical PDE Methods to Solve Calibration and Pricing Problems in Local Stochastic Volatility Models, available at http://kodu.ut.ee/~spartak/.
[83] Strikwerda, J.: Finite Difference Schemes and Partial Differential Equations, Society for Industrial and Applied Mathematics, 2nd edition, 2007.
[84] Whalley, A. E., Wilmott P.: An asymptotic analysis of an optimal hedging model for option pricing with transaction costs. Mathematical Finance 7(3), pp. 307-324, 1997.
[85] Willard, G.A.: Calculating prices and sensitivities for path-independent derivative securities in multifactor models. Journal of Derivatives, 5(1), pp. 45-61, 1997.
[86] Wissel, J.: Arbitrage-free market models for option prices, Working Paper Series, Working Paper No. 248, FINRISK, 2007.

[^0]: Kind souls will allege it is still preferable to ultraslow mispricing.
[^1]: Cover art: in Albrecht Dürer's famous Melencolia I print (detail), a disgruntled quant sits amid the tools of his trade some unfamiliar to the author. The scale and hourglass are symbols, respectively, for no-arbitrage and the convex order. Rows and columns of the magic square all sum up to 34, as do numbers in the four quadrants, as well as in the center quadrant. The middle cells of the bottom row spell out the year the engraving was made: 1514.
