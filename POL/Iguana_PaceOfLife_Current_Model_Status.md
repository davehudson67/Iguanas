# Iguana actuarial senescence modelling: current model status, rationale and next steps

## Purpose of this note

This document summarises where the iguana survival/senescence analysis currently stands, why the model has taken its present form, what we have learned from the model-development process, and how the current preliminary results should be interpreted.

The aim of the analysis is not simply to fit a Gompertz survival curve, but to separate two biologically distinct components of adult mortality:

1. **mortality level at a representative adult age**, and  
2. **the rate at which mortality increases with age (actuarial senescence)**.

The main biological predictors of interest are **sex**, **cohort density**, and **supplementary feeding**. The current model allows each of these focal predictors to influence both components of the mortality trajectory.

The results below are **preliminary 100,000-iteration results**, intended to assess model structure, parameterisation and convergence before the final long MCMC run. They should therefore not yet be treated as final inferential results.

---

# 1. Data and latent life-history structure

The current analysis contains **2,937 individuals** across eight island populations.

The survival model accommodates uncertainty in both birth and death time. This is important because many individuals were first encountered after birth and most individuals are right-censored rather than observed dying.

## 1.1 Birth-time uncertainty

Birth times are handled using three sources of information.

### Known-age individuals

For **597 individuals**, age information provides a relatively precise birth interval.

### SVL-informed birth intervals

For individuals without known age, first recorded snout–vent length (SVL) is used to constrain plausible age and therefore birth year.

The current dataset contains:

- **2,312 individuals** with standard SVL-informed birth intervals;
- **3 individuals** whose SVL-derived birth interval initially implied a birth later than an earlier recorded capture.

For those three animals, first capture was treated as hard evidence that the animal already existed, so the upper birth bound was constrained by first capture.

### Individuals without SVL information

A further **25 individuals** lack both known age and usable SVL information. These retain much broader birth uncertainty.

Thus, **28 individuals** in total are flagged as having particularly weak or adjusted birth histories:

- 25 with no SVL information;
- 3 with SVL/first-capture conflict.

These individuals will later be removed in a targeted sensitivity analysis.

---

# 2. Correction to the previous birth/death handling

An important data-preparation issue was identified during the model audit.

The previous workflow applied a blanket correction when the upper birth interval overlapped the upper death/censoring interval. This altered the birth interval for **23 individuals**.

Closer examination showed that **all 23 were right-censored individuals**.

For a right-censored animal, the censoring date is evidence that the animal was still alive at that time; it is not a death date. Therefore, modifying the birth interval merely because it overlapped the right-censor boundary was unnecessary.

The current setup reconstructs the intended birth intervals directly and does **not** apply this blanket correction.

As a result:

- 23 birth intervals were restored by one year;
- all 23 changes involved right-censored animals;
- no impossible interval-censored birth/death histories remained.

This revised treatment is now the basis of the current model.

---

# 3. Why the Gompertz model was reparameterised

The original models used the conventional Gompertz hazard:

\[
h(t)=a_0\exp(bt),
\]

where:

- \(a_0\) is the mortality hazard extrapolated to age zero;
- \(b\) controls how rapidly mortality increases with age.

Although mathematically valid, this parameterisation created an interpretational and statistical problem.

Most of the data concern animals observed during adult life. Consequently, the model was being asked to estimate an age-zero mortality intercept that is largely outside the most informative region of the data. Changes in \(a_0\) could therefore compensate strongly for changes in \(b\).

This made it harder to distinguish:

- baseline mortality level;
- the rate of senescence.

## 3.1 Reference-age parameterisation

The model has therefore been reparameterised as:

\[
h(t)=a_{10}\exp[b(t-10)].
\]

Here:

- \(a_{10}\) is mortality hazard at **age 10 years**;
- \(b\) remains the Gompertz actuarial-senescence parameter.

This is exactly the same Gompertz family. It is a reparameterisation rather than a different biological model.

The conventional age-zero parameter used internally by the custom Gompertz distribution is recovered as:

\[
a_0=a_{10}\exp(-10b).
\]

Age 10 was chosen because it lies in a strongly data-supported part of the life course and is biologically interpretable as an adult reference age.

## 3.2 What the reparameterisation achieved

This change substantially improved separation of mortality level and senescence.

In the current time-pooled model:

\[
\mathrm{cor}(\log a_{10,\mathrm{fed}},\log b_{\mathrm{fed}})
=-0.041.
\]

This is effectively negligible and is a major improvement over the strong intercept/slope compensation present in earlier age-zero formulations.

---

# 4. Current biological model structure

The current primary model treats **sex, cohort density and feeding** as focal biological predictors.

Rather than deciding in advance that a predictor may affect only mortality level or only senescence, each focal variable is allowed to influence **both** components where biologically plausible.

The model is therefore:

\[
\log(a_{10,i}) =
\log(a_{10,\mathrm{fed}})
+ \beta_{\mathrm{sex},a}X_{\mathrm{sex},i}
+ \beta_{\mathrm{dens},a}X_{\mathrm{dens},i}
+ \beta_{\mathrm{feed},a}X_{\mathrm{feed},i}
+ u_{\mathrm{time},i}
+ u_{\mathrm{island},i},
\]

and

\[
\log(b_i) =
\log(b_{\mathrm{fed}})
+ \beta_{\mathrm{sex},b}X_{\mathrm{sex},i}
+ \beta_{\mathrm{dens},b}X_{\mathrm{dens},i}
+ \beta_{\mathrm{feed},b}X_{\mathrm{feed},i}.
\]

Thus the focal hypotheses are:

- does sex alter mortality at age 10?
- does sex alter the rate of actuarial senescence?
- does early-life/cohort density alter mortality at age 10?
- does cohort density alter actuarial senescence?
- is supplementary feeding associated with mortality level?
- is supplementary feeding associated with senescence rate?

---

# 5. Sex effect

Sex is coded female = 0, male = 1 and centred for model fitting.

Sex is allowed to affect both:

- \(a_{10}\), mortality at age 10;
- \(b\), the rate of actuarial senescence.

This is biologically preferable to forcing sex to influence only one part of the mortality curve.

An effect on \(b\) can also be viewed as an age-dependent sex difference, because it changes how quickly the male/female mortality difference develops across age.

---

# 6. Cohort density

Cohort density is intended to capture the environment experienced by an individual's birth cohort.

Density is estimated from the minimum-number-alive information for each island and year, divided by usable island area. It is then log-transformed and standardised.

The resulting density covariate therefore represents approximately a **one-standard-deviation shift in log cohort density**.

For years in which direct density information is unavailable, an island-level mean is used.

Importantly, density is now allowed to affect both:

- mortality level at age 10;
- the rate of senescence.

This was a deliberate revision from earlier models in which density had been allowed to affect mortality level only.

That earlier restriction turned out to matter because the clearest density signal in the current model occurs on the senescence parameter rather than on \(a_{10}\).

---

# 7. Supplementary feeding and the island-level design

Feeding is applied at island level rather than individually.

Five islands are in the fed category and three are unfed.

Therefore, individual animals are **not independent treatment replicates** for feeding. The effective replication for feeding is at the island/block level.

This is critical for interpretation.

Any feeding result should therefore be described as an **island-level feeding association**, rather than a straightforward individual-level causal treatment effect.

## 7.1 Feeding reparameterisation

An earlier centred coding of feeding produced a very strong posterior dependency between the overall \(b\) intercept and the feeding effect on \(b\):

\[
\mathrm{cor}(\log b,\beta_{\mathrm{feed},b})\approx -0.96.
\]

A diagnostic showed why.

The senescence rate among fed islands was estimated well, whereas the corresponding rate among the three unfed islands was much less precise.

The model was therefore reparameterised so that the fed group is an explicit observed reference:

- fed = 0;
- unfed = -1.

Under this coding:

\[
\beta_{\mathrm{feed},b}
=
\log(b_{\mathrm{fed}})
-
\log(b_{\mathrm{unfed}}).
\]

The same coding is used for \(a_{10}\).

This is biologically easier to interpret and greatly improves parameter separation.

In the current model:

\[
\mathrm{cor}(\log b_{\mathrm{fed}},\beta_{\mathrm{feed},b})=0.130,
\]

so the previous near-perfect posterior ridge has disappeared.

---

# 8. Residual island heterogeneity

A residual island random effect is included on \(a_{10}\).

This allows baseline adult mortality to differ among islands for reasons not explained by:

- feeding;
- cohort density;
- sex.

The island random effects are constrained to have an exact zero mean so that the model intercept retains a meaningful population-level interpretation.

There is currently **no island random effect on \(b\)** in the primary model.

This is deliberate.

With only eight islands, and feeding itself varying only among islands, a separate unrestricted island effect on senescence would compete directly with the feeding contrast for limited between-island information.

An `island_b` formulation remains a useful sensitivity model, but is not currently part of the primary specification.

---

# 9. Temporal cohort effects

The model includes residual temporal variation in \(a_{10}\).

This represents unmeasured differences among birth cohorts that are not captured by measured density.

The initial model used independent five-year cohort blocks.

A support audit showed very sparse representation of the earliest historical periods:

| Original cohort period | Birth-interval midpoint count |
|---|---:|
| 1946–1950 | 0 |
| 1951–1955 | 0 |
| 1956–1960 | 9 |
| 1961–1965 | 29 |
| 1966–1970 | 54 |
| 1971–1975 | 79 |
| 1976–1980 | 106 |
| 1981–1985 | 285 |
| 1986–1990 | 342 |
| 1991–1995 | 473 |
| 1996–2000 | 552 |
| 2001–2005 | 466 |
| 2006–2010 | 334 |
| 2011 onward | 208 |

Because it was not realistic to estimate independent five-year effects for periods represented by zero or very few midpoint birth assignments, the temporal structure was simplified.

The current temporal blocks are:

1. 1946–1965  
2. 1966–1970  
3. 1971–1975  
4. 1976–1980  
5. 1981–1985  
6. 1986–1990  
7. 1991–1995  
8. 1996–2000  
9. 2001–2005  
10. 2006–2010  
11. 2011–2019/20

The first four raw historical blocks are therefore pooled, as are the latest sparse raw periods.

This pooling was chosen on the basis of data support, not because it improved a particular biological result.

---

# 10. Effect of temporal pooling

The pooled temporal model produces almost identical focal biological estimates to the preceding model.

For example:

| Parameter | Previous fed-reference model | Time-pooled model |
|---|---:|---:|
| \(\beta_{\mathrm{sex},a}\) | 0.439 | 0.444 |
| \(\beta_{\mathrm{sex},b}\) | 0.128 | 0.131 |
| \(\beta_{\mathrm{dens},a}\) | 0.086 | 0.092 |
| \(\beta_{\mathrm{dens},b}\) | -0.163 | -0.163 |
| \(b_{\mathrm{fed}}\) | 0.206 | 0.204 |

This stability is reassuring.

However, pooling the sparse early blocks did **not** solve the slow mixing of the temporal variance parameter:

- previous `sigma_time`: Rhat ≈ 1.13, ESS ≈ 90;
- pooled `sigma_time`: Rhat ≈ 1.12, ESS ≈ 87.

This indicates that the slow mixing of the temporal variance is not simply a consequence of the sparse early blocks.

The pooled structure is nevertheless retained because it is a more defensible description of the actual temporal resolution in the data.

---

# 11. Preliminary results from the current 100k model

The current results come from:

- 100,000 MCMC iterations;
- 30,000 burn-in;
- 3 chains;
- thinning interval of 10.

These are **model-development results**, not yet the final inferential run.

## 11.1 Sex and mortality at age 10

\[
\beta_{\mathrm{sex},a}=0.444
\]

95% credible interval:

\[
0.157,\ 0.739.
\]

Exponentiating gives:

\[
\exp(0.444)=1.56.
\]

Thus males currently have an estimated mortality hazard at age 10 approximately **1.56 times that of females**, holding the other model terms constant.

The 95% credible interval on this hazard ratio is approximately:

\[
1.17-2.09.
\]

Posterior probability that the effect is positive:

\[
P(\beta_{\mathrm{sex},a}>0)\approx0.999.
\]

This is currently one of the clearest results in the model.

---

# 12. Sex and actuarial senescence

For the senescence parameter:

\[
\beta_{\mathrm{sex},b}=0.132.
\]

95% credible interval:

\[
-0.022,\ 0.279.
\]

The corresponding multiplicative effect on \(b\) is approximately:

\[
\exp(0.132)=1.14.
\]

Thus the current posterior suggests that male Gompertz senescence may be around **14% faster** than female senescence.

However, the 95% credible interval still just overlaps no effect:

\[
0.978-1.32.
\]

The posterior probability of a positive effect is approximately:

\[
P(\beta_{\mathrm{sex},b}>0)=0.954.
\]

The current interpretation is therefore:

> There is strong evidence that males have higher adult mortality at age 10, together with substantial but less conclusive evidence that male mortality also accelerates more rapidly with age.

---

# 13. Cohort density and mortality at age 10

For mortality level:

\[
\beta_{\mathrm{dens},a}=0.093
\]

with 95% credible interval:

\[
-0.144,\ 0.328.
\]

On the multiplicative scale:

\[
\exp(0.093)\approx1.10.
\]

The credible interval is approximately:

\[
0.87-1.39.
\]

There is therefore currently **little clear evidence that cohort density simply shifts mortality hazard at age 10**.

---

# 14. Cohort density and actuarial senescence

The clearest new result from the revised model is the effect of density on \(b\):

\[
\beta_{\mathrm{dens},b}=-0.164.
\]

95% credible interval:

\[
-0.238,\ -0.082.
\]

The multiplicative effect is:

\[
\exp(-0.164)=0.849.
\]

Thus a one-standard-deviation increase in cohort density is associated with an estimated senescence-rate multiplier of approximately **0.85**.

In other words, higher-density cohorts are currently estimated to have a Gompertz senescence rate about **15% lower per +1 SD cohort density**.

The 95% credible interval for this multiplier is approximately:

\[
0.788-0.922.
\]

The posterior probability that the effect is negative is effectively 1.0.

This parameter also mixes reasonably well in the 100k test:

- Rhat ≈ 1.01;
- effective sample size ≈ 961.

This result is particularly important because the earlier model structure had constrained the density effect on \(b\) to zero. The current biologically symmetric model shows that the strongest density association may concern **the rate at which mortality increases with age**, rather than mortality level at a single adult reference age.

This finding will need to survive the final long run and sensitivity analyses before being treated as robust.

---

# 15. Feeding and mortality at age 10

The feeding effect on \(a_{10}\) is:

\[
\beta_{\mathrm{feed},a}\approx-0.112
\]

at the posterior median.

The corresponding fed/unfed mortality ratio is approximately:

\[
0.89.
\]

However, uncertainty is extremely large:

\[
95\%\,CrI \approx 0.20-4.69.
\]

There is therefore currently **no clear evidence for an island-level feeding association with mortality at age 10**.

---

# 16. Feeding and actuarial senescence

For senescence:

\[
\beta_{\mathrm{feed},b}\approx0.172
\]

at the posterior median.

This corresponds to a fed/unfed senescence-rate ratio of approximately:

\[
1.19.
\]

The 95% credible interval is broad:

\[
0.53-2.43.
\]

The posterior probability of a positive feeding effect is only approximately 0.61.

Thus there is currently **no compelling evidence that feeding increases or decreases actuarial senescence**.

The group-specific estimates reinforce this point.

For fed islands:

\[
b_{\mathrm{fed}}\approx0.204,
\]

95% credible interval:

\[
0.183-0.226.
\]

This parameter is estimated well:

- Rhat ≈ 1.00;
- ESS ≈ 990.

For unfed islands:

\[
b_{\mathrm{unfed}}\approx0.172
\]

at the posterior median, with a much broader distribution.

Its uncertainty reflects the limited replication among only three unfed islands.

The current feeding result should therefore be interpreted primarily as:

> The available data estimate the senescence trajectory of fed islands quite precisely, but provide much weaker information about the contrast with unfed islands. Consequently, the feeding association remains highly uncertain.

---

# 17. Temporal pattern in mortality level

The temporal random effects show a substantial long-term pattern rather than purely small fluctuations around zero.

Approximate posterior means are:

| Birth-cohort period | Temporal effect on log \(a_{10}\) |
|---|---:|
| 1946–1965 | +2.88 |
| 1966–1970 | +1.63 |
| 1971–1975 | +0.75 |
| 1976–1980 | +0.30 |
| 1981–1985 | +0.30 |
| 1986–1990 | +0.72 |
| 1991–1995 | +0.87 |
| 1996–2000 | -0.07 |
| 2001–2005 | -2.09 |
| 2006–2010 | -2.57 |
| 2011 onward | -2.73 |

This suggests a pronounced secular decline in baseline adult mortality across birth cohorts.

At present this is represented as exchangeable temporal cohort heterogeneity rather than a fitted linear temporal trend.

The temporal variance parameter remains the slowest-mixing part of the current 100k model:

\[
\sigma_{\mathrm{time}}\approx2.15,
\]

with:

- Rhat ≈ 1.12;
- ESS ≈ 87.

This is not yet satisfactory for final inference and is one reason the full longer MCMC run is required.

The temporal effect pattern will also motivate a later sensitivity model in which a long-term cohort trend is represented explicitly, potentially together with residual temporal variation.

---

# 18. Current convergence assessment

Several focal biological parameters already mix well in the 100k run:

| Parameter | Rhat | ESS |
|---|---:|---:|
| \(b_{\mathrm{fed}}\) | 1.00 | 990 |
| \(\beta_{\mathrm{sex},a}\) | 1.01 | 1311 |
| \(\beta_{\mathrm{sex},b}\) | 1.01 | 1024 |
| \(\beta_{\mathrm{dens},a}\) | 1.01 | 1389 |
| \(\beta_{\mathrm{dens},b}\) | 1.01 | 961 |
| detection \(p\) | 1.00 | 10224 |

The main parameters still requiring substantially more sampling are:

- `sigma_time`;
- the feeding effects;
- derived unfed-group parameters;
- some mortality-level intercept/random-effect quantities.

This pattern is consistent with limited block-level replication for feeding and with a difficult hierarchical temporal variance parameter.

---

# 19. Why no additional interactions are currently in the primary model

The model already contains the key biologically relevant age-dependent effects through \(b\).

For example:

- sex effect on \(b\) represents sex-dependent ageing;
- density effect on \(b\) represents density-dependent ageing;
- feeding effect on \(b\) represents feeding-dependent ageing.

Thus the core biological interactions with age are already built into the model.

Potential additional interactions include:

- feeding × density;
- feeding × sex;
- density × sex.

These are not currently included because they would add substantial complexity and, particularly for feeding × density, could be difficult to distinguish from island-level contrasts.

A feeding × density interaction is biologically interesting because supplementation could potentially buffer density-dependent resource limitation, but it is better treated as a targeted secondary hypothesis than automatically added to the primary model.

---

# 20. Why species is not currently included

The two iguana forms in the dataset are biologically similar and earlier model development provided little evidence for a strong independent species effect.

Species is also strongly structured by island, meaning that it competes with island-level effects for the same limited among-island information.

For those reasons, species is not currently included in the primary model.

It can be revisited if there is a specific biological hypothesis that requires it.

---

# 21. Why the current model is preferable to the earlier model-development versions

The current model is preferable for several reasons.

## 21.1 Its parameters are biologically interpretable

Mortality level is defined at age 10 rather than as an extrapolated age-zero hazard.

## 21.2 Focal predictors are treated symmetrically

Sex, density and feeding can each influence both mortality level and senescence.

We are therefore no longer deciding whether a biological predictor belongs on \(a\) or \(b\) based primarily on whether a particular model happened to mix.

## 21.3 Feeding is parameterised at a real observed group

The fed-reference parameterisation removed the near-perfect posterior dependency between the senescence intercept and feeding contrast.

## 21.4 Sparse historical time periods are not over-resolved

The model no longer tries to estimate independent five-year temporal effects for periods containing effectively no birth-cohort information.

## 21.5 The focal biological results are stable across structural changes

The sex and density estimates have remained very similar when:

- the Gompertz intercept was moved to age 10;
- feeding was reparameterised;
- early temporal blocks were pooled.

In particular, the negative density effect on \(b\) has been approximately:

\[
-0.16
\]

throughout these revised formulations.

That stability is encouraging, although it is not a substitute for final convergence and sensitivity analysis.

---

# 22. Current biological interpretation

At this stage, the model suggests a potentially coherent pace-of-life pattern.

### Sex

Males appear to:

- have clearly higher adult mortality at age 10;
- probably experience somewhat faster actuarial senescence.

### Cohort density

Higher cohort density currently shows:

- little evidence for a strong shift in mortality at age 10;
- a clear negative association with the Gompertz senescence rate.

This implies that cohort environment may influence the **shape of the mortality trajectory** rather than simply increasing or decreasing mortality uniformly across adult life.

Potential mechanisms could eventually include developmental effects, selective disappearance, cohort quality, or broader life-history responses to density.

Those mechanisms should not yet be distinguished on the basis of the present model alone.

### Supplementary feeding

There is currently no clear feeding association with either:

- mortality level at age 10;
- senescence rate.

However, the feeding contrast is intrinsically imprecise because feeding is replicated at island level, with only five fed and three unfed islands.

The appropriate conclusion at present is therefore uncertainty, rather than evidence that feeding has no biological effect.

---

# 23. Important caution about interpreting \(a_{10}\) and \(b\) independently

A predictor that affects both \(a_{10}\) and \(b\) changes the entire mortality trajectory.

For example, a positive effect on mortality at age 10 combined with a negative effect on \(b\) could generate convergence or crossing of hazard curves at different ages.

Therefore, the final biological interpretation should not rely only on tables of \(\beta_a\) and \(\beta_b\).

Once the final posterior is available, the analysis should generate age-specific posterior hazard and survival curves for:

- males versus females;
- low versus high cohort density;
- fed versus unfed islands.

These plots will show how the estimated effects combine over the life course.

---

# 24. What happens next

The immediate next step is to run the current time-pooled, fed-reference, age-10 model for substantially longer:

- **500,000 MCMC iterations**;
- **100,000 burn-in**;
- **3 chains**;
- **thin = 10**.

No further terms will be added to the primary model before this run.

The purpose of the long run is to determine whether:

1. `sigma_time` converges adequately;
2. the feeding parameters mix sufficiently for their uncertainty to be interpreted;
3. the strong negative density effect on senescence remains stable;
4. the sex effects remain stable;
5. the low posterior correlations produced by the new parameterisation are maintained.

If the long run is satisfactory, this will become the **primary model**.

After that, a deliberately limited set of sensitivity analyses will be run:

- exclude the 28 weak/adjusted birth histories;
- leave one island out at a time when assessing feeding;
- fit a temporal-trend sensitivity to test whether the pronounced cohort pattern is better represented by an explicit secular trend plus residual temporal variation;
- consider one model allowing residual island variation in \(b\), to check whether the senescence results depend on excluding island-specific ageing rates.

Potential higher-order interactions such as feeding × density will be considered only as secondary hypothesis-driven analyses rather than added automatically to the primary model.

---

# Brief summary for discussion

The model-development process has moved from a conventional age-zero Gompertz model toward a more interpretable **reference-age pace-of-life model**.

The current formulation separates:

- mortality hazard at age 10;
- actuarial senescence rate.

Sex, cohort density and feeding are each allowed to influence both components.

Preliminary results suggest:

- **higher male adult mortality**;
- **probable faster male senescence**;
- **little clear density effect on mortality level at age 10**;
- **a strong negative association between cohort density and actuarial senescence**;
- **no clear feeding effect, with considerable uncertainty due to limited island-level replication**.

These patterns have remained strikingly stable through the recent model reparameterisations.

The next step is therefore not to add further complexity, but to run the current model for **500,000 iterations**, assess final convergence, and then subject the resulting biological conclusions to a small set of targeted sensitivity analyses.
