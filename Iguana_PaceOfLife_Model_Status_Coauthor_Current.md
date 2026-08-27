# Iguana pace-of-life survival modelling: current status, model development, preliminary results and final checks

## Purpose of this note

This document summarises the current state of the iguana survival/senescence analysis, including:

- the biological questions being addressed;
- the sequence of modelling decisions that led to the present model;
- the treatment of uncertain birth and death times;
- the reference-age Gompertz reparameterisation;
- the treatment of sex, supplementary feeding, cohort density, island effects and temporal effects;
- the preliminary results obtained so far;
- the recent audit of the cohort-density covariate;
- the small extension now being tested in the latest model;
- the remaining checks required before treating the analysis as final.

The overall aim is to separate two biologically distinct components of adult mortality:

1. **mortality level at a representative adult age**, and  
2. **the rate at which mortality increases with age (actuarial senescence)**.

The main biological predictors of interest are:

- **sex**;
- **cohort density**;
- **supplementary feeding**.

The latest model is currently running as a 100,000-iteration development fit after revising the cohort-density construction to reduce bias associated with the early years of monitoring.

The numerical posterior results presented below refer to the immediately preceding 100,000-iteration model, before that density correction. They are therefore **preliminary**, but they provide the benchmark against which the current run will be compared.

---

# 1. Data and latent life-history structure

The current analysis contains **2,937 individuals** from eight island populations.

The model explicitly accounts for uncertainty in both birth and death time.

This is essential because:

- many animals entered the monitoring programme after birth;
- birth year is often not known exactly;
- many individuals were last seen alive rather than observed dying;
- study entry occurred at different times among islands.

The model therefore treats birth time and, where appropriate, death time as latent variables rather than fixing them to a single guessed year.

---

# 2. Birth-time reconstruction

Birth information comes from three levels of evidence.

## 2.1 Known-age individuals

For **597 individuals**, age information gives a relatively narrow birth interval.

## 2.2 SVL-informed birth intervals

For animals without known age, first recorded snout–vent length (SVL) is used to constrain plausible age.

Quantile regressions relating age to SVL are used to define lower and upper plausible ages, which are then translated into birth-year intervals.

This provides standard SVL-informed birth intervals for **2,312 individuals**.

Three additional individuals produced an SVL-derived birth interval that conflicted with an earlier known capture.

For those animals, first capture is treated as hard evidence that the animal was already alive, so the inferred birth window is constrained accordingly.

## 2.3 Individuals without usable SVL information

A further **25 individuals** lack both known age and usable first-SVL information.

Their birth intervals are therefore much broader.

In total, **28 individuals** are flagged as having particularly weak or adjusted birth histories:

- 25 with no usable SVL information;
- 3 with an SVL/first-capture conflict.

These 28 individuals are retained in the primary analysis but will be removed in a targeted sensitivity analysis.

---

# 3. Correction to an earlier birth/death rule

During the model audit, an important issue was identified in the earlier data-preparation workflow.

A previous rule shortened the upper birth interval whenever it overlapped the upper death/censoring interval.

This affected **23 individuals**.

All 23 were subsequently shown to be **right-censored individuals**.

For a right-censored individual, the censoring date means only that the individual was known to be alive up to that point. It is not a death date.

There was therefore no justification for altering the birth interval merely because it overlapped the censoring boundary.

The current workflow reconstructs birth intervals directly and does not apply that blanket adjustment.

After the correction:

- 23 birth intervals were restored by one year;
- all 23 were right-censored animals;
- no genuinely impossible interval-censored birth/death histories remained.

This revised treatment is retained in all current models.

---

# 4. Survival model

The survival process is modelled using a Gompertz hazard.

The conventional form is

\[
h(t)=a_0\exp(bt),
\]

where:

- \(a_0\) is the extrapolated mortality hazard at age zero;
- \(b\) controls the increase in mortality with age.

The parameter \(b\) is therefore the key parameter describing **actuarial senescence**.

---

# 5. Why the Gompertz model was reparameterised

The original analyses used the age-zero Gompertz intercept \(a_0\).

That is mathematically valid, but it is not ideal for these data because most information comes from adult animals.

Consequently, \(a_0\) is estimated by extrapolating away from the best-supported part of the life course.

This generated strong compensation between:

- mortality level;
- senescence rate.

To improve biological interpretation and posterior geometry, the model was rewritten around age 10:

\[
h(t)=a_{10}\exp[b(t-10)].
\]

Here:

- \(a_{10}\) is mortality hazard at age 10;
- \(b\) is still the Gompertz senescence parameter.

The conventional age-zero parameter is recovered internally as

\[
a_0=a_{10}\exp(-10b).
\]

Thus this is not a different survival model; it is a more interpretable parameterisation of the same Gompertz hazard.

Age 10 was chosen because it lies within a strongly represented part of the observed life course.

The reparameterisation worked well.

In the latest completed development model,

\[
\mathrm{cor}(\log a_{10,\mathrm{fed}},\log b_{\mathrm{fed}})
=-0.041,
\]

showing that adult mortality level and senescence are now almost completely separated statistically.

---

# 6. Current biological model

The focal predictors are allowed to affect both Gompertz components.

Conceptually,

\[
\log(a_{10,i}) =
\alpha_a
+ \beta_{\mathrm{sex},a}X_{\mathrm{sex},i}
+ \beta_{\mathrm{dens},a}D_{j,t}
+ \beta_{\mathrm{feed},a}X_{\mathrm{feed},j}
+ u_{\mathrm{time},t}
+ u_{\mathrm{island},j},
\]

and

\[
\log(b_i) =
\alpha_b
+ \beta_{\mathrm{sex},b}X_{\mathrm{sex},i}
+ \beta_{\mathrm{dens},b}D_{j,t}
+ \beta_{\mathrm{feed},b}X_{\mathrm{feed},j}.
\]

The current coding is deliberately simple:

- **female = 0, male = 1**;
- **unfed = 0, fed = 1**;
- cohort density is a standardised continuous variable.

The latest implementation retains the well-informed fed group as the reference group for the Gompertz intercepts, while the feeding contrast itself is still the ordinary fed-versus-unfed categorical comparison.

---

# 7. Why the focal variables are included on both Gompertz components

An important part of the model-development process was moving away from assigning biological predictors to only one mortality component.

There is no strong biological reason to assume in advance that sex, density or feeding can influence mortality level but not senescence, or vice versa.

The primary model therefore asks two separate questions for each focal predictor:

1. does it alter mortality at age 10?
2. does it alter the rate at which mortality increases with age?

This is central to the pace-of-life interpretation.

A predictor acting on \(b\) is inherently an **age-dependent effect**, because it changes the slope of the log-hazard with age.

---

# 8. Sex

Sex is entered directly as:

- female = 0;
- male = 1.

It is allowed to influence both:

- mortality at age 10;
- actuarial senescence.

This means that the model can distinguish a simple adult mortality difference from a difference in the ageing trajectory itself.

---

# 9. Supplementary feeding

Feeding is treated as a binary island-level predictor:

- unfed = 0;
- fed = 1.

The important design point is that feeding is applied at the **island/block level** rather than individually.

Five islands are fed and three are unfed.

Therefore, thousands of individual iguanas do not constitute thousands of independent feeding replicates.

The effective replication for the feeding comparison is the number of islands.

This limits precision and means that feeding results should be described as **island-level associations**, rather than as straightforward individual-level causal treatment effects.

The fed islands are much better represented for estimating the Gompertz senescence trajectory than the three unfed islands, which explains much of the uncertainty in the feeding contrast.

---

# 10. Island and temporal effects

The model includes two additional sources of heterogeneity on mortality at age 10.

## 10.1 Residual island effect

Each island has a residual effect on \(a_{10}\).

This represents persistent island-level mortality differences not explained by:

- feeding;
- density;
- sex.

The island effects are constrained to average zero.

## 10.2 Shared temporal cohort effect

The model also contains a temporal birth-cohort effect on \(a_{10}\).

This effect is **shared across all islands**.

It represents broad environmental or historical conditions affecting animals born in the same general period.

This is conceptually different from density.

The structure is therefore:

- **island effect:** varies among islands, constant through time;
- **time effect:** varies through time, common to all islands;
- **density:** varies jointly by island and cohort year.

That distinction has become particularly important during the recent density audit.

---

# 11. Temporal pooling

The original temporal model used separate five-year birth-cohort blocks.

An audit showed very little support for several early historical blocks.

The approximate birth-window midpoint support was:

| Cohort period | Number of midpoint assignments |
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
| 2011–2019 | 208 |

The temporal structure was therefore simplified to:

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
11. 2011–2019  

This pooling was based on data support, not on the direction of any biological coefficient.

---

# 12. Effect of temporal pooling

The key biological estimates were almost unchanged after pooling.

For example:

| Parameter | Previous model | Time-pooled model |
|---|---:|---:|
| sex effect on \(a_{10}\) | 0.439 | 0.444 |
| sex effect on \(b\) | 0.128 | 0.131 |
| density effect on \(a_{10}\) | 0.086 | 0.092 |
| density effect on \(b\) | -0.163 | -0.163 |
| \(b\) in the well-informed fed group | 0.206 | 0.204 |

This stability is encouraging.

Pooling did not materially improve the temporal variance parameter:

- before pooling: `sigma_time` Rhat ≈ 1.13, ESS ≈ 90;
- after pooling: `sigma_time` Rhat ≈ 1.12, ESS ≈ 87.

This indicates that the slow mixing of the temporal variance is not simply caused by sparse early blocks.

The pooled structure is nevertheless retained because it better reflects the temporal resolution actually supported by the data.

---

# 13. Preliminary results from the latest completed 100k model

The latest completed development model used:

- 100,000 iterations;
- 30,000 burn-in;
- 3 chains;
- thin = 10.

These results are the benchmark against which the current corrected-density model will be assessed.

---

# 14. Sex effect on mortality at age 10

The estimated male-versus-female effect was:

\[
\beta_{\mathrm{sex},a}=0.444
\]

with 95% credible interval:

\[
0.157,\ 0.739.
\]

On the multiplicative scale:

\[
\exp(0.444)=1.56.
\]

Thus males had an estimated mortality hazard at age 10 approximately **56% higher than females**, holding the remaining model terms constant.

The 95% credible interval for the hazard ratio was approximately:

\[
1.17-2.09.
\]

The posterior probability of a positive effect was approximately:

\[
P(\beta_{\mathrm{sex},a}>0)=0.999.
\]

This is currently one of the strongest and most stable results.

---

# 15. Sex effect on actuarial senescence

The estimated sex effect on \(b\) was:

\[
\beta_{\mathrm{sex},b}=0.132.
\]

The 95% credible interval was:

\[
-0.022,\ 0.279.
\]

On the multiplicative scale:

\[
\exp(0.132)\approx1.14.
\]

Thus the posterior suggests that male senescence may be approximately **14% faster** than female senescence.

The 95% interval on this ratio was approximately:

\[
0.98-1.32.
\]

The posterior probability that the effect is positive was approximately **0.954**.

The current interpretation is therefore:

> Males show clear evidence of higher adult mortality, together with substantial but not yet conclusive evidence of faster actuarial senescence.

---

# 16. Preliminary cohort-density result before the new audit correction

The previous density covariate was based on island-year **minimum number alive (MNA)** divided by usable island area.

Density was log-transformed and standardised before entering the model.

The previous model produced little evidence that density affected mortality at age 10:

\[
\beta_{\mathrm{dens},a}=0.093
\]

with 95% credible interval approximately:

\[
-0.144,\ 0.328.
\]

However, the density effect on senescence was much stronger:

\[
\beta_{\mathrm{dens},b}=-0.164.
\]

The 95% credible interval was:

\[
-0.238,\ -0.082.
\]

On the multiplicative scale:

\[
\exp(-0.164)=0.849.
\]

This means that a one-standard-deviation increase in the previous cohort-density covariate was associated with an estimated **15% lower senescence rate**.

The corresponding 95% interval for the senescence-rate multiplier was approximately:

\[
0.788-0.922.
\]

The posterior probability of a negative density effect was essentially 1.

The parameter also mixed well:

- Rhat ≈ 1.01;
- ESS ≈ 961.

This was therefore the strongest new result to emerge from the revised model structure.

Because it was potentially important, the construction of the density covariate was audited before moving to the final long run.

---

# 17. MNA does not assume complete detection

An important clarification is that MNA does **not** assume every individual on an island was detected during every survey.

An animal is counted as known alive between its first and last recorded year, even if it was missed during intermediate surveys.

MNA is therefore a conservative index of the minimum number known to have been alive.

However, MNA is not automatically equivalent to true abundance.

In particular, it can be biased low during the early years of a monitoring programme because animals that were already alive but had not yet been encountered do not enter the marked/known population until their first capture.

This creates the possibility of an **ascertainment ramp**.

---

# 18. Density audit: evidence of a monitoring start-up effect

The audit showed strong evidence that several islands exhibit exactly this pattern.

For the first ten survey years, the correlation between MNA and time since monitoring began was:

| Island | Correlation between MNA and monitoring time |
|---|---:|
| Bitter | 0.902 |
| FFRC | 0.948 |
| Gaulin | 0.982 |
| Leaf | 0.917 |
| NAdder | -0.042 |
| Noddy | 0.672 |
| UCay | 0.981 |
| WhiteB | 0.322 |

At the same time, the fraction of MNA represented by newly encountered animals generally declined strongly through the first monitoring years.

Examples of mean entrant fraction in survey years 1–3 versus year 6 onward were:

| Island | Years 1–3 | Year 6+ |
|---|---:|---:|
| Bitter | 0.851 | 0.449 |
| FFRC | 0.728 | 0.133 |
| Gaulin | 0.572 | 0.235 |
| Leaf | 0.700 | 0.169 |
| NAdder | 0.678 | 0.233 |
| Noddy | 0.872 | 0.698 |
| UCay | 0.645 | 0.114 |
| WhiteB | 0.462 | 0.295 |

Likewise, median MNA in the first three survey years versus survey year 6 onward was:

| Island | Median MNA years 1–3 | Median MNA year 6+ |
|---|---:|---:|
| Bitter | 7 | 54 |
| FFRC | 5 | 89 |
| Gaulin | 19 | 113 |
| Leaf | 2 | 211 |
| NAdder | 53 | 54 |
| Noddy | 33 | 90 |
| UCay | 52 | 162 |
| WhiteB | 51 | 80.5 |

NAdder is an important exception, with little evidence of a large early MNA ramp.

Taken together, the pattern on several other islands is consistent with the marked population accumulating through the early monitoring period rather than MNA immediately representing a stable proxy for true abundance.

---

# 19. Why the shared temporal effect cannot fully solve this problem

The model already contains a temporal cohort effect, but this is common across all islands.

The MNA start-up problem is different.

Monitoring began at different times on different islands.

Therefore, an artificial density ramp associated with “years since monitoring began” occurs at different calendar dates on different islands.

A single shared temporal effect cannot simply absorb this structure.

This is exactly why density and time need to remain conceptually separate:

- **time:** common environmental/historical variation across all islands;
- **density:** local island × cohort-year demographic environment;
- **monitoring start-up:** island-specific ascertainment process.

The density covariate therefore needed to be revised directly.

---

# 20. Audit of observed versus imputed density

The audit also showed that many cohort years occur before direct annual MNA data are available.

At the midpoint of each individual's birth interval, the proportion receiving directly observed survey-year MNA was approximately:

| Island | Observed MNA | Island-mean imputation |
|---|---:|---:|
| Bitter | 4% | 96% |
| FFRC | 33% | 67% |
| Gaulin | 52% | 48% |
| Leaf | 68% | 32% |
| NAdder | 6% | 94% |
| Noddy | 11% | 89% |
| UCay | 46% | 54% |
| WhiteB | 10% | 90% |

Across complete possible birth intervals, the mean proportion of candidate cohort years with directly observed density was approximately **0.42**.

Thus the density covariate contains two distinct types of information:

1. observed annual island-specific density during monitored years;
2. an island-level typical value when annual density is unavailable.

This distinction is now made explicit in the revised construction.

---

# 21. Current extension: corrected cohort-density model

The model currently running retains the same biological and survival structure but changes how cohort density is constructed.

The aim is to preserve genuine island × year density variation while preventing the early accumulation of marked animals from being interpreted as true low density.

## 21.1 Definition of the start-up period

A single reproducible rule is used across islands:

> **The first five actual survey years on each island are treated as the monitoring/ascertainment phase.**

This avoids choosing a different cut-off by inspection for each island.

The first-five-year rule is intentionally conservative and will itself be checked in sensitivity analysis.

## 21.2 Established monitoring period

From **survey year 6 onward**, annual MNA-derived density is treated as sufficiently established to contribute annual variation.

Thus:

- established surveyed year → use that island-year's observed MNA density;
- first five survey years → do not use the raw early MNA value;
- pre-monitoring year → no raw annual MNA is assumed;
- unsurveyed gap → no raw annual MNA is assumed.

## 21.3 Established island baseline

For each island, a typical established density is calculated from survey years 6 onward.

The mean is calculated on the `log1p(density)` scale and transformed back where necessary.

For years without trusted annual density, the model uses this **established island mean**.

This includes:

- years before monitoring began;
- the first five monitoring years;
- unsurveyed gaps.

The rationale is not that density was literally constant in those years.

Rather, where the data cannot support annual density variation, it is preferable to represent the island's typical density than to insert an artificial start-up trajectory or zero density.

## 21.4 Standardisation

The mean and standard deviation used to standardise density are now calculated only from **trusted established annual observations**.

This is an important improvement.

A one-standard-deviation density effect now reflects the scale of real observed annual density variation during established monitoring, rather than a distribution dominated by repeated imputed values and artificially low start-up MNA.

---

# 22. Structure of the model currently running

The current candidate model can therefore be summarised as follows.

### Mortality at age 10

\[
\log(a_{10,i}) =
\alpha_a
+ \beta_{\mathrm{sex},a}X_{\mathrm{sex},i}
+ \beta_{\mathrm{dens},a}D^*_{j,t}
+ \beta_{\mathrm{feed},a}X_{\mathrm{feed},j}
+ u_{\mathrm{time},t}
+ u_{\mathrm{island},j}.
\]

### Actuarial senescence

\[
\log(b_i) =
\alpha_b
+ \beta_{\mathrm{sex},b}X_{\mathrm{sex},i}
+ \beta_{\mathrm{dens},b}D^*_{j,t}
+ \beta_{\mathrm{feed},b}X_{\mathrm{feed},j}.
\]

where \(D^*_{j,t}\) is now:

- observed annual density for island \(j\), year \(t\), when monitoring is established;
- otherwise the established mean density of island \(j\).

The temporal effect remains shared across islands.

The island random effect remains specific to island.

The model still contains latent birth time, left truncation, interval/right censoring, and the detection component.

---

# 23. Preliminary feeding results

The previous time-pooled model showed no clear feeding association.

For mortality at age 10, the posterior median fed/unfed ratio was approximately:

\[
0.89,
\]

with a very broad 95% credible interval:

\[
0.20-4.69.
\]

For senescence, the posterior median fed/unfed ratio was approximately:

\[
1.19,
\]

with 95% credible interval:

\[
0.53-2.43.
\]

The posterior probability of a positive feeding effect on senescence was only about **0.61**.

The appropriate interpretation is therefore not “feeding has no effect”, but rather:

> the available island-level replication is insufficient to estimate the feeding contrast precisely.

This is why leave-one-island-out analysis will be especially important for feeding.

---

# 24. Preliminary temporal pattern

The pooled temporal effects showed a substantial long-term pattern in mortality level.

Approximate posterior means were:

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
| 2011–2019 | -2.73 |

This looks more like a pronounced secular change than purely exchangeable noise around zero.

The current primary model still treats these as exchangeable period effects.

A later sensitivity model will test whether a model containing an explicit long-term temporal trend plus residual temporal deviations provides a more natural representation.

---

# 25. Current convergence status

In the latest completed 100k model, several focal biological parameters were already mixing well:

| Parameter | Rhat | ESS |
|---|---:|---:|
| \(b\) in the well-informed fed group | 1.00 | 990 |
| sex effect on \(a_{10}\) | 1.01 | 1311 |
| sex effect on \(b\) | 1.01 | 1024 |
| density effect on \(a_{10}\) | 1.01 | 1389 |
| density effect on \(b\) | 1.01 | 961 |
| detection probability | 1.00 | 10224 |

The principal convergence problem was the temporal variance:

- `sigma_time` ≈ 2.15;
- Rhat ≈ 1.12;
- ESS ≈ 87.

Feeding-related quantities were also less well sampled because of weak island-level replication.

The new density model is therefore being run for 100k first so that a meaningful change in the density coefficient is not hidden inside a very long computation.

---

# 26. What we are looking for in the current run

The most important comparison is the new estimate of:

\[
\beta_{\mathrm{dens},b}.
\]

The previous estimate was approximately:

\[
-0.164.
\]

Three broad outcomes are possible.

## 26.1 The effect remains strongly negative and similar in magnitude

This would be highly encouraging.

It would indicate that the density–senescence association is not primarily driven by the early monitoring ascertainment ramp.

The result would then deserve serious biological interpretation.

## 26.2 The effect remains negative but becomes weaker or less precise

This would suggest that part of the original signal was contributed by the density construction, but that an underlying biological association may remain.

The effect size and uncertainty would then need to be interpreted conservatively.

## 26.3 The effect disappears or changes substantially

This would indicate that the previous density–senescence result was strongly dependent on the MNA start-up structure.

That would be an important result in itself and would prevent over-interpretation of an artefactual signal.

The current run is therefore a genuine robustness test, not simply an attempt to preserve the previous result.

---

# 27. Why we are not adding more complexity yet

The current model already captures the central pace-of-life questions.

Effects on \(b\) are intrinsically age-dependent effects.

For example:

- sex on \(b\) asks whether males and females age at different rates;
- density on \(b\) asks whether cohort environment changes the ageing rate;
- feeding on \(b\) asks whether fed and unfed islands differ in ageing rate.

Possible interactions such as:

- feeding × density;
- feeding × sex;
- density × sex

may be biologically interesting, but are not currently being added to the primary model.

This avoids creating a large exploratory model in a dataset where feeding replication is limited to eight islands.

Such interactions can be investigated later as targeted secondary hypotheses if justified.

---

# 28. Why species is not currently included

The two iguana forms represented in the data are biologically similar.

Earlier model development provided little evidence for a substantial independent species effect.

Species is also strongly confounded with island.

Including species alongside island-level effects would therefore consume limited between-island information without a strong prior biological reason.

Species is consequently not part of the current primary model, although it could be revisited if a specific species-level hypothesis becomes important.

---

# 29. Why the current model is now close to a final primary model

The model-development process has progressively addressed the main statistical and biological issues.

### Birth uncertainty

Handled explicitly using latent birth intervals informed by known age, SVL and first-capture evidence.

### Death uncertainty

Handled using interval censoring and right censoring rather than converting uncertain death times into fixed events.

### Delayed study entry

Handled through left truncation.

### Mortality/senescence confounding

Greatly reduced by moving from age-zero mortality to an age-10 reference hazard.

### Sex

Represented as a direct female-versus-male categorical effect on both mortality level and senescence.

### Feeding

Represented as a block-level fed-versus-unfed contrast, with inference appropriately limited by island-level replication.

### Density

Represented as an island × cohort-year process, now being corrected so that early marked-population accumulation is not mistaken for real low density.

### Island heterogeneity

Represented as persistent residual island variation in mortality level.

### Temporal heterogeneity

Represented as a shared cohort-period effect across islands.

### Sparse historical periods

Handled by pooling poorly supported temporal blocks.

The current model is therefore biologically coherent and substantially more interpretable than the earlier versions.

---

# 30. Remaining checks before final inference

The current corrected-density 100k run is the immediate next checkpoint.

If it behaves sensibly, the remaining analyses should be deliberately limited and hypothesis-driven.

## 30.1 Confirm the corrected density result

Compare the new density coefficients with the previous values, particularly:

\[
\beta_{\mathrm{dens},b}\approx-0.164.
\]

## 30.2 Check sensitivity to the start-up rule

The five-survey-year rule is defensible and reproducible, but it is still an analytical choice.

A targeted sensitivity should therefore assess whether the density conclusion changes materially under a nearby rule, for example:

- three start-up survey years;
- five start-up survey years;
- or an established-monitoring-only density analysis.

The aim is not to search for the cut-off producing the strongest effect, but to show that the inference is not driven by a particular reasonable definition.

## 30.3 Remove the 28 weak/adjusted birth histories

Refit after excluding:

- 25 animals without usable SVL information;
- 3 animals with SVL/first-capture conflicts.

## 30.4 Leave-one-island-out feeding analysis

Because feeding is an island-level treatment, rerun the feeding analysis leaving out one island at a time.

This will show whether the apparent feeding contrast is disproportionately driven by a single island.

## 30.5 Residual island variation in senescence

Fit one sensitivity model containing an island-level random effect on \(b\).

This asks whether the density and feeding senescence effects persist after allowing islands to have persistent differences in ageing rate.

This is a sensitivity rather than the primary model because there are only eight islands.

## 30.6 Alternative temporal representation

Fit a model containing an explicit secular cohort trend plus residual temporal variation.

This will test whether the strong long-term temporal pattern is better described as a trend rather than purely exchangeable block variation.

## 30.7 Final age-specific predictions

Once the final posterior is accepted, generate mortality and survival trajectories across age for:

- males versus females;
- lower versus higher cohort density;
- fed versus unfed islands.

This is important because a predictor may affect both \(a_{10}\) and \(b\), meaning that its effect on mortality can change with age.

---

# 31. Final MCMC run

The 500,000-iteration final run has deliberately not yet been started for the corrected-density model.

The current 100k fit is first being used to verify that the revised density construction behaves sensibly.

Once the density definition and model structure are accepted, the candidate primary model will be run for:

- **500,000 iterations**;
- **100,000 burn-in**;
- **3 chains**;
- **thin = 10**.

The final assessment will focus on:

- Rhat;
- effective sample size;
- posterior stability across chains;
- stability of the focal biological coefficients;
- the temporal variance;
- feeding uncertainty;
- sensitivity-analysis consistency.

---

# 32. Current biological picture

The emerging biological picture is currently as follows.

## Sex

There is strong evidence that males have higher adult mortality than females.

There is also substantial evidence that males may experience faster actuarial senescence, although the 95% interval still marginally overlaps no difference.

## Cohort density

Before correcting the early MNA ascertainment issue, higher cohort density was strongly associated with lower actuarial senescence.

The effect was approximately a 15% reduction in \(b\) per one-standard-deviation increase in density.

This is potentially the most interesting new result, but it is currently being subjected to a stricter density reconstruction before being accepted.

## Supplementary feeding

No clear feeding effect has yet been detected on either mortality level or senescence.

The uncertainty is substantial and consistent with the fact that feeding is replicated across only eight islands.

## Time

There is a strong long-term cohort pattern in baseline adult mortality that is common across islands in the current model.

Its form may ultimately be better described by an explicit temporal trend.

---

# Brief co-author summary

The analysis has developed into a **reference-age Gompertz pace-of-life model** that separates mortality level at age 10 from actuarial senescence.

Sex, cohort density and supplementary feeding are each allowed to influence both mortality components.

The latest completed model suggests:

- substantially higher male adult mortality;
- probable faster male senescence;
- little clear cohort-density effect on mortality at age 10;
- a strong negative density–senescence association;
- no clear feeding effect, with substantial uncertainty due to island-level replication.

The main new issue identified during model checking concerns the cohort-density covariate.

MNA itself does not assume complete detection, but several islands show a pronounced rise in MNA during the first monitoring years while the fraction of newly encountered animals falls. This strongly suggests that early MNA partly reflects accumulation into the marked population.

Because monitoring began at different dates on different islands, the shared temporal effect cannot fully absorb this artefact.

The model currently running therefore uses annual MNA-derived density only once monitoring is established. The first five survey years, pre-monitoring years and unsurveyed gaps receive each island's established typical density instead. Density is standardised using only the trusted established annual observations.

The immediate question is whether the previous density–senescence coefficient of approximately **−0.16** survives this correction.

If it does, the result becomes substantially more convincing.

If it weakens or disappears, that will show that part of the previous association was generated by the structure of the monitoring record.

Once this 100k corrected-density run is assessed, a small number of pre-specified sensitivity analyses will be conducted, followed by the final 500k MCMC fit and age-specific posterior mortality/survival predictions.
