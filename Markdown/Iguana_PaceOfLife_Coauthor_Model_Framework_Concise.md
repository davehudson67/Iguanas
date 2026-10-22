# Iguana pace-of-life survival model
## Current modelling framework, rationale and biological interpretation

**Status:** Candidate primary model structure following model and data audits. Results below are from the current 100,000-iteration development fit; the intention is to agree the biological/model structure with co-authors before running the definitive longer MCMC.

---

## 1. Aim

The central question is whether individual and cohort-level conditions are associated with variation in **pace of life**, separating two components of mortality:

1. **adult mortality level**, represented by mortality hazard at age 10 (\(a_{10}\)); and
2. **actuarial senescence**, represented by the Gompertz ageing parameter (\(b\)), which describes how rapidly mortality increases with age.

The current model is therefore designed to ask:

- Do males and females differ in adult mortality and/or the rate of senescence?
- Do cohorts born in unusually high- or low-density years for their own island differ in subsequent adult mortality and senescence?
- Is there a long-term change in adult mortality among birth cohorts?
- Is there residual variation among islands and among cohort periods after these measured effects are accounted for?

The emphasis is on estimating **associations with mortality trajectories**, rather than assigning causal effects to density or other environmental variables.

---

## 2. Core survival framework

Birth and death are incompletely observed for many individuals, so the model retains uncertainty in the life history rather than assigning a single estimated age or lifespan.

Birth time is treated as latent within an individual-specific interval informed by known age where available and by the SVL-age relationship otherwise. Death is handled using interval- or right-censoring, and the survival likelihood is left-truncated so that animals already alive when monitoring began are correctly conditioned on having survived to the observable period.

A detection component also accounts for opportunities to observe an animal while it was alive. Thus uncertainty in birth date, lifespan and observation history is propagated into the mortality parameters and cohort covariates.

The current analysis includes **2,937 individuals from eight islands**.

### Why use mortality at age 10?

The standard Gompertz model can be written as

\[
h(t)=a_0\exp(bt),
\]

but \(a_0\) is mortality extrapolated to age zero and can be strongly correlated with \(b\). We instead use the mathematically equivalent reference-age form

\[
h(t)=a_{10}\exp[b(t-10)].
\]

Age 10 was chosen **before interpreting the current biological effects**, because it lies within a well-supported part of the age distribution:

- median estimated age at first observation: **10.5 years**;
- median age at last observation: **15.5 years**;
- median midpoint of observed age ranges: **13.5 years**;
- among known-birth animals, the median midpoint of the observed age range was **9.5 years**.

Thus \(a_{10}\) is a meaningful adult mortality parameter rather than an extrapolated juvenile intercept. The reparameterisation also successfully separated the two Gompertz components: in the no-feeding reference-age model the posterior correlation between \(\log(a_{10})\) and \(\log(b)\) was only **0.011**.

---

## 3. Why the current covariate structure was chosen

### Sex

Sex is included on both components of the mortality trajectory:

\[
\log(a_{10,i}) \supset \beta_{\mathrm{sex},a}\,\mathrm{sex}_i
\]

and

\[
\log(b_i) \supset \beta_{\mathrm{sex},b}\,\mathrm{sex}_i.
\]

Females are the reference category (0) and males are coded 1. This allows sex to affect both adult mortality level and the subsequent rate of actuarial ageing.

### Cohort density: within-island rather than between-island

Density is derived from minimum-number-alive (MNA) by island and year, divided by usable habitat area and analysed on the \(\log(1+\mathrm{density})\) scale.

Auditing the raw MNA series showed a pronounced **monitoring start-up/ascertainment ramp** on several islands: the earliest survey years often contained very low MNA and high proportions of newly encountered animals. We therefore do not treat the first five survey years on an island as reliable annual density estimates. Annual density deviations are used from survey year 6 onward.

A second audit showed that, among trusted established observations, **88.7% of the variation in density was between islands and only 11.3% was within islands through time**. An overall density covariate would therefore be dominated by comparisons such as “dense island versus sparse island”, potentially conflating density with persistent island differences.

The current model instead uses an island-centred density variable:

\[
D^{W}_{j,t}
=
\log(1+D_{j,t})
-
\overline{\log(1+D)}_j.
\]

This asks a much cleaner question:

> **For a given island, do animals from cohorts born in a denser- or sparser-than-usual established year show different subsequent mortality trajectories?**

The deviation is standardised using the variation in trusted established annual observations. A +1 unit change therefore represents a **one-SD increase in density relative to that island's own typical density**.

For pre-monitoring years, the five start-up survey years and unsurveyed gaps, the within-island deviation is set to zero. This does **not** assert that density was exactly average; it means that we do not invent evidence for an annual departure from the island mean where reliable annual information is unavailable. This is deliberately conservative and may attenuate the estimated density effect.

### Long-term temporal trend

The earlier model allowed cohort periods to differ through exchangeable random effects on \(a_{10}\). Those effects revealed a striking long-term pattern: older cohorts had much higher estimated mortality levels and recent cohorts much lower mortality levels. A descriptive line through those period estimates explained about **88%** of their variation.

The current model therefore separates:

1. an explicit linear **birth-cohort trend**, expressed per decade, and
2. residual cohort-period effects around that trend.

Conceptually,

\[
\log(a_{10,i}) =
\alpha_a
+\beta_{\mathrm{sex},a}X_i
+\beta_{\mathrm{dens},a}D^W_{j,t}
+\beta_{\mathrm{time}}T_i
+u_{\mathrm{period},t}
+u_{\mathrm{island},j}.
\]

The temporal trend is applied to mortality level at age 10. Residual period effects retain flexibility for particular historical periods that depart from the long-term trend.

This formulation substantially reduced the large monotonic pattern previously carried by the period random effects. After adding the trend, most residual temporal effects were approximately between **-0.65 and +0.63** on the log-hazard scale rather than spanning roughly **+3 to -2.8** from old to recent cohorts. There remains some structured period variation, particularly around the 1970s and late 1980s-1990s, but the broad secular pattern is now represented explicitly rather than being forced into exchangeable random effects.

### Residual island variation

A random island effect remains on \(a_{10}\) to account for persistent differences in adult mortality among the eight island populations that are not explained by sex, within-island cohort density or temporal effects.

The primary model does not currently include an unrestricted island random effect on \(b\). With only eight islands, such an effect is considerably harder to estimate and would represent an additional hypothesis about persistent population differences in senescence. It remains a possible publication-stage sensitivity rather than part of the current primary structure.

---

## 4. Why feeding is not in the candidate primary model

Supplementary feeding was initially considered because islands can be classified broadly as fed or unfed. However, the exposure is applied at the island/block level rather than independently to animals, giving only five nominally fed and three unfed island blocks.

More importantly, the binary classification may be a poor proxy for actual exposure. Islands differ in size, access, tourist visitation and iguana habituation, so not all animals on a nominally fed island are likely to experience equivalent supplementation.

This became both a biological and statistical concern. When feeding was included, its posterior was broad and the senescence contrast mixed poorly. We therefore fitted the same model without feeding. Removing it produced **almost no change** in the sex or density estimates. For example:

| Effect | Feeding included | Feeding removed |
|---|---:|---:|
| Sex -> \(a_{10}\) | 0.366 | 0.365 |
| Sex -> \(b\) | 0.175 | 0.175 |
| Density -> \(a_{10}\) | -1.15 | -1.15 |
| Density -> \(b\) | -0.090 | -0.092 |

Thus feeding was not acting as an important adjustment variable for the main relationships, while removing it made the model cleaner and more defensible. Feeding can still be reported as an exploratory analysis if there is strong biological interest from the co-author group.

---

## 5. Current candidate model

The current candidate primary model is:

\[
\begin{aligned}
\log(a_{10,i}) =\;&
\alpha_a
+\beta_{\mathrm{sex},a}\mathrm{sex}_i
+\beta_{\mathrm{dens},a}D^W_{j,t}\\
&+\beta_{\mathrm{time}}T_i
+u_{\mathrm{period},t}
+u_{\mathrm{island},j},
\end{aligned}
\]

and

\[
\log(b_i)=
\alpha_b
+\beta_{\mathrm{sex},b}\mathrm{sex}_i
+\beta_{\mathrm{dens},b}D^W_{j,t}.
\]

The individual mortality trajectory is then

\[
h_i(t)=a_{10,i}\exp[b_i(t-10)].
\]

In words, the model separates differences in **adult mortality level** from differences in the **rate of ageing**, while allowing both sex and within-population cohort density to influence each component. Long-term cohort change, shorter-term period departures and persistent island heterogeneity are modelled separately on adult mortality.

---

## 6. What the current model is telling us

The current 100k development fit gives a coherent biological pattern.

| Effect | Median \(\beta\) | 95% credible interval | Multiplicative effect | Interpretation |
|---|---:|---:|---:|---|
| Male -> age-10 mortality | 0.315 | 0.034 to 0.581 | 1.37 | ~37% higher age-10 mortality in males |
| Male -> senescence | 0.216 | 0.092 to 0.349 | 1.24 | ~24% faster actuarial senescence in males |
| +1 SD within-island density -> age-10 mortality | -0.447 | -0.632 to -0.237 | 0.64 | ~36% lower age-10 mortality |
| +1 SD within-island density -> senescence | -0.131 | -0.219 to -0.030 | 0.88 | ~12% lower actuarial senescence |

Posterior support for these directions is strong: approximately **98.6%** for higher male age-10 mortality, **>99.9%** for faster male senescence, **>99.9%** for lower age-10 mortality at higher within-island cohort density, and **99.3%** for slower senescence at higher within-island cohort density.

### Sex

The sex result has been notably stable through the model-development checks. The current interpretation is that males appear to experience both:

- a higher adult mortality level; and
- a steeper increase in mortality with age.

This is therefore a difference in the overall mortality trajectory rather than simply a shift in baseline risk.

### Density

The density result is particularly interesting because it survives the move from an overall island-density measure to a **within-island cohort comparison**.

The current result is not simply that naturally dense islands contain longer-lived animals. Instead, within the established monitoring period, cohorts associated with **higher-than-usual density for their own island** subsequently show:

- lower mortality at age 10; and
- a slower increase in mortality with age.

This is consistent with a cohort or environmental association with subsequent pace of life. It should not yet be interpreted causally. MNA is a population-history-derived proxy rather than an independent experimental measure of density, and the amount of established annual information differs among islands.

Possible biological mechanisms — including cohort quality, selective processes, developmental environment or population dynamics — should therefore be discussed only after the statistical pattern and the biology of the populations are considered together.

### Time

There is also a substantial long-term cohort signal in adult mortality. Later-born cohorts appear to have lower age-10 mortality than earlier-born cohorts, with additional residual deviations among periods.

The important point at this stage is not the precise linear coefficient, but that explicitly modelling the secular trend produces a much more natural temporal structure: the residual period effects no longer have to reproduce the entire old-to-recent mortality gradient.

The causes of this long-term change are currently unresolved. Potential explanations could include real environmental change, management change, population processes or changes in observation/monitoring, and should be discussed with the co-author group before being given a biological interpretation.

---

## 7. Checks that materially informed the current framework

Four checks materially changed or strengthened the final modelling story:

1. **Birth-history audit and reference-age parameterisation.**  
   Birth intervals and right-censoring were checked carefully, and the Gompertz model was expressed at age 10 rather than age zero. This made the adult mortality parameter biologically interpretable and greatly reduced dependence between mortality level and senescence.

2. **Density start-up audit.**  
   Early MNA values showed clear monitoring/ascertainment ramps on several islands. The first five actual survey years are therefore not used as annual density deviations. A 3-versus-5-year start-up sensitivity retained the same broad direction of density effects.

3. **Feeding removal.**  
   Feeding is weakly replicated and actual individual exposure is uncertain. Removing it left the principal sex and density estimates essentially unchanged, supporting the cleaner no-feeding framework.

4. **Between- versus within-island density and temporal audit.**  
   Most overall density variation was between islands, while the temporal random effects showed a strong secular pattern. The current model therefore asks a within-island density question and explicitly represents the long-term cohort trend. The density association remains after making those changes.

These checks are important because they show that the current structure was driven by the design and diagnostics rather than by selecting whichever model produced the strongest statistical result.

---

## 8. Current interpretation and remaining caution

At this stage, the clearest emerging picture is:

- **males have a faster pace of mortality**, with both higher adult mortality and faster actuarial senescence;
- **higher-than-usual cohort density within an island is associated with a slower mortality trajectory**, affecting both adult mortality level and senescence;
- **adult mortality has changed markedly among birth cohorts through time**, beyond the measured density effect;
- meaningful **residual island and period differences** remain.

The density finding is potentially the most novel but should be worded carefully. The model supports a within-island association, but the density index is derived from long-term capture histories and reliable annual density information is concentrated in established monitoring years. For unobserved/start-up years the model assumes no known annual deviation, which is conservative but reduces information.

Likewise, the temporal trend should initially be treated as an important empirical pattern rather than assigned to a particular ecological mechanism.

---

## 9. Proposed next step

I suggest treating this as the **candidate primary model**.

But I think we need to consider these questions:

- Does within-island MNA-derived density appropriately represent the cohort environment we want to discuss?
- What historical or management processes could plausibly explain the strong long-term cohort trend?
- Is there a compelling biological reason to retain feeding, species, or additional island-level senescence structure in the primary analysis rather than as secondary sensitivities?

Once the structure is agreed we can run with longer chains for final posterior estimation, followed by a small number of pre-specified publication-stage checks if needed. The final presentation should focus on posterior mortality trajectories and survival/hazard curves rather than on coefficients alone.

---

### Working summary

The current framework is a latent-life-history Gompertz model parameterised at age 10, with sex and **within-island cohort-density deviations** acting on both adult mortality and actuarial senescence, an explicit birth-cohort trend and residual period/island heterogeneity acting on adult mortality; no feeding term in the candidate primary model.

The main current result is that **male mortality is both higher and more strongly age-dependent, while animals associated with higher-than-usual density in their birth cohort show lower adult mortality and slower actuarial senescence. A strong long-term decline in adult mortality across birth cohorts is also evident and warrants biological interpretation.**
