# Final iguana survival model: technical and plain-language explanation

## Purpose of this document

This document describes the **final primary model** fitted to the nine-site Allen Cays rock iguana dataset through 2019. It explains both the statistical formulation and the biological meaning of each component.

The model was designed to answer four linked questions:

1. How does mortality change with **age**?
2. Do males and females differ in **adult mortality level** and/or the **rate of senescence**?
3. Is survival associated with **density around the time of birth** and **birth cohort**?
4. Does mortality change during calendar periods of **high/regular provisioning**, while allowing for site-specific background temporal change and the major disturbance at FFRC associated with Hurricane Matthew?

At the same time, the model explicitly accounts for imperfect and highly variable **detection**, uncertainty in **birth year**, uncertain **sex** for a minority of animals, **delayed entry** into the study, and interval/right censoring of death.

This document describes the primary model only. The provisioning, Matthew and density sensitivity models and posterior predictive checks are separate validation analyses.

---

# 1. Dataset entering the model

The final cleaned dataset contained:

- **4,152 individuals**
- **9 sites**
- data through **2019**
- **74 confirmed deaths**
- **4,078 right-censored individuals**
- **1,879 observed females**
- **1,983 observed males**
- **290 individuals with latent sex**
  - 39 with conflicting recorded F/M assignments
  - 251 never assigned a sex
- **2,117 individuals** whose entire possible birth interval was covered by the density history and therefore contributed to estimation of the primary birth-density effect
- **170 reconstructed field-visit detection occasions** used after first capture
- **64,195 post-entry individual × visit detection opportunities**

The nine modeled populations were:

- Bitter
- FFRC
- Gaulin
- Leaf North
- Leaf South
- North Adderly
- Noddy
- U Cay
- White Bay

Leaf North and Leaf South were treated as separate biological populations rather than pooled as a single "Leaf" site.

---

# 2. The model in one sentence

For each iguana, we estimated an uncertain birth date and age at death, modeled age-specific mortality with a **left-truncated Gompertz curve**, allowed sex and birth density to alter both the level and rate of senescence, allowed provisioning and Hurricane Matthew to modify mortality in the calendar years when they occurred, and simultaneously modeled whether a living animal was actually detected at each field visit.

---

# 3. Plain-language overview

An iguana does not enter the study at birth. It might first be caught at age 3, 10, 20 or older. We therefore cannot treat first capture as the start of its biological lifespan.

Instead, the model says:

> "This animal was already some uncertain age when we first caught it. We know it survived up to that first capture, so start its likelihood there, but keep its mortality risk appropriate to its true age."

The underlying ageing process is Gompertz: mortality increases approximately exponentially with age. The model then asks whether that age trajectory differs according to sex and density around birth.

Provisioning is handled differently. It is not assigned permanently to an individual. An animal can live through years with no high-intensity provisioning and later live through years with high/regular provisioning. Its underlying ageing curve continues normally, while its instantaneous mortality hazard is multiplied by the provisioning effect only during the calendar years when high provisioning is present at its site.

Finally, not seeing an iguana does not mean it is dead. A separate detection model estimates the probability of catching a living animal at each actual field visit.

---

# 4. Time and notation

For individual \(i\):

- \(B_i\) = latent birth time
- \(A_i\) = age at death
- \(D_i = B_i + A_i\) = calendar/model time of death
- \(L_i\) = age at first live capture
- \(j(i)\) = site occupied by individual \(i\)
- \(s_i\) = sex indicator, with female \(=0\), male \(=1\)
- \(z_i\) = standardized within-site density covariate around the inferred birth year
- \(c_i\) = centered birth-cohort covariate

Calendar time is represented internally on an adjusted scale beginning in 1943, but all biological interpretations are in ordinary calendar years and ages.

---

# 5. Latent birth time

Birth date/year is not known exactly for most individuals.

The model treats birth time as a latent variable:

\[
B_i \sim \mathrm{Uniform}(B_{i,\min}, B_{i,\max}).
\]

The interval \([B_{i,\min}, B_{i,\max}]\) comes from the data-preparation stage.

For animals with reliable known-age information, these intervals are tight. For other animals, the interval is inferred principally from the SVL-age calibration and capture history.

This is important because age is not inserted into the survival model as a single guessed value. Uncertainty in birth date propagates through the survival analysis.

### Plain-language interpretation

Rather than saying:

> "This iguana was definitely born in 1997",

the model can say:

> "Given the available age and size information, its birth date is somewhere in this defensible interval, and we will average the survival inference over that uncertainty."

---

# 6. Latent sex

For animals whose recorded sex was unambiguous, sex is treated as observed.

The cleaning rule was:

- F records plus unknown records only \(\rightarrow\) female
- M records plus unknown records only \(\rightarrow\) male
- both F and M recorded \(\rightarrow\) sex remains uncertain
- all records unknown \(\rightarrow\) sex remains uncertain

For uncertain individuals:

\[
s_i \sim \mathrm{Bernoulli}(0.5).
\]

The \(0.5\) prior is intentionally weak. We did not use the observed capture sex ratio as the prior because detection itself differs between sexes.

### Plain-language interpretation

The model does not majority-vote contradictory sex records. If an animal has genuinely conflicting male and female records, the uncertainty is carried into the analysis.

---

# 7. Baseline Gompertz ageing process

The underlying age-specific mortality hazard is parameterized around **age 10**:

\[
h_{0i}(a)
=
a_{10,i}
\exp\left[b_i(a-10)\right].
\]

Here:

- \(a\) is true age
- \(a_{10,i}\) is the individual's mortality hazard at age 10 before calendar-period multipliers are applied
- \(b_i\) controls how rapidly mortality accelerates with age

The usual Gompertz model is often written as

\[
h(a)=a_0e^{ba},
\]

where \(a_0\) is mortality at age zero. We instead parameterize it at age 10 because mortality at age 10 is much more biologically interpretable for these long-lived iguanas.

The two forms are algebraically equivalent:

\[
a_0
=
a_{10}e^{-10b}.
\]

No new ageing distribution is created by this reparameterization; it is simply a more useful reference point on the same Gompertz curve.

---

# 8. Mortality level at age 10

The individual-specific age-10 hazard is:

\[
\log(a_{10,i})
=
\log(a_{10})
+
\beta_{\mathrm{sex},a}s_i
+
\beta_{\mathrm{dens},a}z_i
+
\beta_{\mathrm{cohort}}c_i
+
u^{(a)}_{j(i)}.
\]

The terms mean:

- \(\log(a_{10})\): population reference mortality at age 10
- \(\beta_{\mathrm{sex},a}\): male-female difference in mortality level at age 10
- \(\beta_{\mathrm{dens},a}\): effect of birth-period density on age-10 mortality
- \(\beta_{\mathrm{cohort}}\): residual effect of birth cohort
- \(u^{(a)}_j\): residual site effect on age-10 mortality

The birth-cohort variable is centered on 1990 and scaled in decades:

\[
c_i
=
\frac{\mathrm{birth\ year}_i-1990}{10}.
\]

Therefore \(\beta_{\mathrm{cohort}}\) is the log-hazard change associated with being born one decade later, after accounting for the other modeled effects.

---

# 9. Rate of senescence

The Gompertz senescence parameter is modeled separately:

\[
\log(b_i)
=
\log(b)
+
\beta_{\mathrm{sex},b}s_i
+
\beta_{\mathrm{dens},b}z_i.
\]

Thus sex and birth density can affect not only mortality at age 10 but also **how quickly mortality accelerates thereafter**.

This distinction is central to the biological interpretation.

A covariate can therefore produce:

- higher mortality in early adulthood but slower subsequent senescence,
- little difference at age 10 but much faster senescence,
- or effects in the same direction on both components.

### Why there is not one constant sex hazard ratio

Because sex changes \(b_i\), the male:female hazard ratio can change with age.

For otherwise equivalent individuals:

\[
\frac{h_M(a)}{h_F(a)}
=
\exp(\beta_{\mathrm{sex},a})
\,
\exp\left[(b_M-b_F)(a-10)\right].
\]

So a near-zero sex effect at age 10 can coexist with a strong sex difference at older ages.

---

# 10. Birth-density covariate

Density was constructed from minimum-number-alive histories within each site.

For each site-year:

1. MNA was transformed as \(\log(1+\mathrm{MNA})\).
2. The site's long-term mean was removed.
3. The resulting within-site deviations were standardized.

The density covariate therefore represents **relative density within a site**, rather than comparing the absolute abundance of differently sized islands.

Conceptually:

\[
z_{j,y}
=
\mathrm{standardize}
\left[
\log(1+\mathrm{MNA}_{j,y})
-
\overline{\log(1+\mathrm{MNA}_{j,\cdot})}
\right].
\]

Density was interpolated between observed survey years.

The individual's density value is evaluated at its latent birth year:

\[
z_i = z_{j(i),B_i}.
\]

Because \(B_i\) is uncertain, the relevant density year can change during MCMC.

## Conservative primary support rule

The primary model allows the density coefficient to be informed only by individuals for whom **the entire possible birth interval** lies within years supported by the density history.

There were 2,117 such individuals.

For the other 2,035 animals:

\[
z_i=0
\]

for purposes of the density term.

Here, zero means the **within-site average on the centered density scale**; it does not mean zero animals.

Those individuals still contribute fully to survival, detection, sex, site and other model components. They simply do not provide information about the density coefficient in the primary analysis.

---

# 11. Residual site effect on mortality

Sites may differ in persistent ways not represented by the measured covariates.

We therefore include:

\[
u^{(a)}_j
=
\sigma_{\mathrm{island},a}
\left(
u^{*(a)}_j-\bar u^{*(a)}
\right),
\]

with

\[
u^{*(a)}_j\sim N(0,1).
\]

Centering the site effects makes their average exactly zero.

### Plain-language interpretation

The model does not force all islands to have identical baseline mortality. It estimates persistent residual differences between sites while keeping the overall intercept interpretable.

---

# 12. Calendar-varying mortality

The full hazard is:

\[
h_i(a,y)
=
a_{10,i}
\exp\left[b_i(a-10)\right]
\exp\left[
\beta_F F_{j(i),y}
+
\beta_M M_{j(i),y}
+
\gamma_{j(i)}P_y
\right].
\]

The first part,

\[
a_{10,i}\exp[b_i(a-10)],
\]

is the age-driven Gompertz curve.

The second part,

\[
\exp\left[
\beta_F F_{j,y}
+
\beta_M M_{j,y}
+
\gamma_jP_y
\right],
\]

is a calendar-year multiplier.

This allows an individual's mortality to change when its environmental exposure changes, without resetting its biological ageing clock.

---

# 13. High/regular provisioning

The primary provisioning covariate is:

\[
F_{j,y}
=
\begin{cases}
1, & \text{high/regular provisioning at site }j\text{ in year }y\\
0, & \text{otherwise.}
\end{cases}
\]

Thus:

\[
\exp(\beta_F)
\]

is the proportional change in instantaneous mortality hazard during high/regular provisioning periods, conditional on all the other components of the model.

The primary history was:

| Site | Primary high/regular provisioning history |
|---|---|
| Bitter | 0 through 2005; 1 from 2006–2019 |
| FFRC | 0 |
| Gaulin | 0 in the primary contrast; lower/intermittent provisioning is not classified as high |
| Leaf North | 1 throughout the modeled monitoring period |
| Leaf South | 0 through 2011; 1 from 2012–2019 |
| North Adderly | 0 |
| Noddy | 0 |
| U Cay | 1 throughout the modeled monitoring period |
| White Bay | 0 in the primary contrast; lower/intermittent provisioning is not classified as high |

This is deliberately a **site × calendar-year exposure**.

It is **not** an individual "fed versus unfed" label.

### Example

An iguana on Bitter that survives across 2005–2006 is modeled as experiencing:

\[
F_{\mathrm{Bitter},2005}=0
\]

and

\[
F_{\mathrm{Bitter},2006}=1.
\]

Its Gompertz age curve continues smoothly across those years. Only the calendar-period multiplier changes.

---

# 14. Hurricane Matthew term

For the primary model:

\[
M_{j,y}
=
1
\]

only for FFRC from 2017 onward, and is zero otherwise.

The corresponding multiplier is:

\[
\exp(\beta_M).
\]

This term is intended to stop a known major FFRC disturbance from being silently absorbed into other temporal terms.

It should be interpreted as an **association with the post-2016 FFRC period**, not as a clean experimental estimate of the causal effect of Hurricane Matthew.

An acute 2016–2018 coding is being examined separately as a sensitivity analysis.

---

# 15. Site-specific temporal slopes

Provisioning histories are partly confounded with calendar time. For example, provisioning begins at Bitter during the study, so a simple before/after comparison could mistakenly attribute any unrelated temporal change at Bitter to feeding.

To make the provisioning estimate more conservative, each site is allowed its own residual temporal slope:

\[
\gamma_j
=
\sigma_{\mathrm{period}}
\left(
\gamma^*_j-\bar\gamma^*
\right),
\]

where:

\[
\gamma^*_j\sim N(0,1).
\]

Calendar time is centered on 2010:

\[
P_y
=
\frac{y-2010}{10}.
\]

Therefore:

\[
\gamma_jP_y
\]

represents a site-specific residual change in log mortality per decade.

The slopes are centered so that:

\[
\sum_j \gamma_j = 0.
\]

## Why there is no unrestricted global calendar-period trend

Age, birth cohort and calendar period are structurally related:

\[
\mathrm{period}
=
\mathrm{birth\ cohort}
+
\mathrm{age}.
\]

A free global age effect, birth-cohort effect and calendar-period effect would therefore create an age-period-cohort identifiability problem.

We retain:

- Gompertz age
- birth cohort
- centered site-specific temporal deviations

but do not add a separate unrestricted global calendar-period slope.

---

# 16. Exact piecewise cumulative hazard

Because provisioning, Matthew status and site-period effects vary through calendar time, survival cannot be evaluated by simply assigning one fixed exposure value to an animal.

Within a given calendar year \(y\), define:

\[
m_{j,y}
=
\exp\left[
\beta_FF_{j,y}
+
\beta_MM_{j,y}
+
\gamma_jP_y
\right].
\]

For an individual with Gompertz parameters \(a_{10,i}\) and \(b_i\), the cumulative hazard between ages \(a_1\) and \(a_2\), when that whole interval lies within the same calendar exposure segment, is:

\[
H_i(a_1,a_2 \mid y)
=
m_{j,y}
\frac{a_{10,i}}{b_i}
\left[
e^{b_i(a_2-10)}
-
e^{b_i(a_1-10)}
\right].
\]

If an animal's post-entry life spans several calendar years, the model divides it at calendar-year boundaries and sums these pieces:

\[
H_i(L_i,A_i)
=
\sum_r
H_i(a_{r,\mathrm{start}},a_{r,\mathrm{end}}\mid y_r).
\]

This is calculated exactly rather than approximated by assuming a single average exposure over the animal's life.

---

# 17. Correct delayed entry / left truncation

An individual only enters the survival study when it is first caught alive.

If its age at first capture is \(L_i\), we already know:

\[
A_i>L_i.
\]

The appropriate death-age density is therefore conditional on survival to entry:

\[
f(A_i \mid A_i>L_i)
=
h_i(A_i)
\exp\left[
-\int_{L_i}^{A_i} h_i(u)\,du
\right].
\]

Equivalently:

\[
f(A_i \mid A_i>L_i)
=
\frac{f(A_i)}{S(L_i)}.
\]

The division by \(S(L_i)\) renormalizes the remaining probability distribution so that:

\[
\int_{L_i}^{\infty}
f(a\mid A_i>L_i)\,da
=
1.
\]

This is the key correction relative to the earlier shifted-Gompertz implementation.

## What the earlier implementation effectively did

The earlier implementation shifted age:

\[
a^*=a-L_i,
\]

which effectively restarted the Gompertz age clock at first capture.

Thus an animal first captured at 20 could initially receive the same Gompertz hazard as an animal first captured at 5.

## What the final model does

The Gompertz curve remains anchored to true biological age.

An animal entering at age 20 starts further up the same ageing curve than an animal entering at age 5, while the likelihood is conditioned on the fact that both survived to their respective entry ages.

---

# 18. An important consequence of conditioning on entry

The survival likelihood is conditional on an individual having survived to its first capture.

Therefore the model estimates **post-entry mortality**.

Calendar-varying provisioning or disturbance histories before that individual's first observed capture are not used to reconstruct the probability that it survived long enough to enter the sample.

This is deliberate: conditioning on first capture is what removes the left-truncation bias.

Birth-period density and birth cohort remain fixed individual covariates because they characterize the individual's early-life/background history, but the survival likelihood begins at observed entry.

---

# 19. Death time and censoring

Calendar/model death time is:

\[
D_i=B_i+A_i.
\]

Observed death information is linked to this latent death time through interval censoring.

For individuals with a confirmed death:

\[
D_i \in [D_{i,\min},D_{i,\max}].
\]

For individuals last observed alive:

\[
D_i>D_{i,\mathrm{last\ alive}}.
\]

Thus animals do not need an exactly known death date to contribute to the survival likelihood.

Animals permanently translocated out of the study population are treated as right-censored rather than dead.

---

# 20. Reconstruction of field visits

Detection opportunities were not defined simply as one occasion per site-year.

The capture-date audit showed that some site-years contained clearly separated field bouts, sometimes months apart.

A new field visit was therefore defined when the next capture day at the same site and in the same year occurred **more than 14 days after the preceding capture day**.

The 14-day choice was based on the empirical audit:

- >14 days produced 179 candidate field visits
- >21 days produced 178 candidate field visits

so the reconstructed visit count had essentially stabilized.

After excluding first-capture occasions where appropriate, **170 field visits** contributed to the detection model.

---

# 21. Which detection opportunities are included?

For each animal:

- there are **no artificial non-detections before first capture**
- the field visit on which the animal is first captured is excluded from the detection likelihood
- subsequent visits are genuine recapture opportunities
- a later visit in the same calendar year can be used
- no opportunities are created after known permanent removal
- no opportunities are created after an observed death

This produced **64,195 post-entry individual × field-visit opportunities**.

---

# 22. Detection model

For individual \(i\) at field visit \(e\), define age:

\[
a_{ie}
=
\max(0,t_e-B_i).
\]

The age effect on detection is:

\[
g(a_{ie})
=
\min\left(1,\frac{a_{ie}}{15}\right).
\]

Thus detection probability increases linearly with age up to 15 years and then plateaus.

Conditional on being alive:

\[
\mathrm{logit}(p_{ie})
=
\alpha_p
+
\beta_{p,\mathrm{sex}}s_i
+
\beta_{p,\mathrm{age}}g(a_{ie})
+
u^{(p)}_{j(i)}
+
v^{(p)}_e.
\]

Here:

- \(\alpha_p\) = reference detection intercept
- \(\beta_{p,\mathrm{sex}}\) = sex effect on detection
- \(\beta_{p,\mathrm{age}}\) = juvenile-to-adult detection effect
- \(u^{(p)}_j\) = site detection effect
- \(v^{(p)}_e\) = field-visit detection effect

Both site and field-visit effects are centered random effects.

---

# 23. Linking survival and detection

The model explicitly checks whether an animal is still alive at a visit:

\[
I_{ie}
=
I(D_i\ge t_e).
\]

The observed detection probability is:

\[
p^{\mathrm{obs}}_{ie}
=
I_{ie}p_{ie}.
\]

Then:

\[
Y_{ie}
\sim
\mathrm{Bernoulli}
\left(
p^{\mathrm{obs}}_{ie}
\right).
\]

If the latent death time occurs before a survey:

\[
I_{ie}=0
\]

and therefore:

\[
p^{\mathrm{obs}}_{ie}=0.
\]

A later observed recapture conversely forces the latent death time to lie after that recapture.

### Plain-language interpretation

Failure to catch an animal is not automatically treated as death.

The model asks two separate questions:

1. Was the animal still alive?
2. If it was alive, what was the probability that we caught it?

This distinction is essential because detection varied substantially among sites and field visits.

---

# 24. Detection random effects

Persistent site-level detection heterogeneity is modeled as:

\[
u^{(p)}_j
=
\sigma_{\mathrm{island},p}
\left(
u^{*(p)}_j-\bar u^{*(p)}
\right),
\]

with:

\[
u^{*(p)}_j\sim N(0,1).
\]

Field-visit heterogeneity is:

\[
v^{(p)}_e
=
\sigma_{\mathrm{survey},p}
\left(
v^*_e-\bar v^*
\right),
\]

with:

\[
v^*_e\sim N(0,1).
\]

The field-visit random effect absorbs variation associated with the actual survey occasion, including unmeasured differences in effort, team composition, weather, animal aggregation and other visit-specific conditions.

We deliberately did **not** use the number caught during a visit as an effort covariate, because capture count itself combines abundance and detectability.

---

# 25. Full model in compact form

The survival component can be summarized as:

\[
B_i
\sim
U(B_{i,\min},B_{i,\max}),
\]

\[
s_i
\sim
\mathrm{Bernoulli}(0.5)
\quad
\text{for latent-sex animals},
\]

\[
\log(a_{10,i})
=
\log(a_{10})
+
\beta_{\mathrm{sex},a}s_i
+
\beta_{\mathrm{dens},a}z_i
+
\beta_{\mathrm{cohort}}c_i
+
u^{(a)}_{j(i)},
\]

\[
\log(b_i)
=
\log(b)
+
\beta_{\mathrm{sex},b}s_i
+
\beta_{\mathrm{dens},b}z_i,
\]

\[
h_i(a,y)
=
a_{10,i}
e^{b_i(a-10)}
e^{
\beta_FF_{j(i),y}
+
\beta_MM_{j(i),y}
+
\gamma_{j(i)}P_y
},
\]

and:

\[
A_i
\sim
\text{left-truncated time-varying Gompertz}
\left(
h_i,\,
A_i>L_i
\right).
\]

The detection component is:

\[
\operatorname{logit}(p_{ie})
=
\alpha_p
+
\beta_{p,\mathrm{sex}}s_i
+
\beta_{p,\mathrm{age}}
\min\left(1,\frac{a_{ie}}{15}\right)
+
u^{(p)}_{j(i)}
+
v^{(p)}_e,
\]

\[
Y_{ie}
\sim
\mathrm{Bernoulli}
\left[
I(D_i\ge t_e)p_{ie}
\right].
\]

---

# 26. Priors

The main priors were intentionally weak relative to the scale of the data.

## Gompertz reference parameters

\[
\log(a_{10})
\sim
N\left[\log(0.03),\,0.75^2\right],
\]

\[
\log(b)
\sim
N\left[\log(0.08),\,0.50^2\right].
\]

## Survival regression coefficients

\[
\beta_{\mathrm{sex},a}
\sim
N(0,0.75^2),
\]

\[
\beta_{\mathrm{sex},b}
\sim
N(0,0.50^2),
\]

\[
\beta_{\mathrm{dens},a}
\sim
N(0,0.75^2),
\]

\[
\beta_{\mathrm{dens},b}
\sim
N(0,0.50^2),
\]

\[
\beta_{\mathrm{cohort}}
\sim
N(0,0.50^2),
\]

\[
\beta_F
\sim
N(0,0.75^2),
\]

\[
\beta_M
\sim
N(0,0.75^2).
\]

## Random-effect scales

\[
\sigma_{\mathrm{period}}
\sim
U(0,0.75),
\]

\[
\sigma_{\mathrm{island},a}
\sim
U(0,2),
\]

\[
\sigma_{\mathrm{island},p}
\sim
U(0,3),
\]

\[
\sigma_{\mathrm{survey},p}
\sim
U(0,3).
\]

## Detection priors

\[
\alpha_p
\sim
N\left[\operatorname{logit}(0.20),1.5^2\right],
\]

\[
\beta_{p,\mathrm{sex}}
\sim
N(0,1^2),
\]

\[
\beta_{p,\mathrm{age}}
\sim
N(0,1^2).
\]

---

# 27. What the key parameters mean

| Parameter | Interpretation |
|---|---|
| \(a_{10}\) | Reference mortality hazard at age 10 before calendar multipliers |
| \(b\) | Reference Gompertz rate of actuarial senescence |
| \(\beta_{\mathrm{sex},a}\) | Male-female contrast in mortality at age 10 |
| \(\beta_{\mathrm{sex},b}\) | Sex effect on the Gompertz senescence parameter |
| \(\beta_{\mathrm{dens},a}\) | Birth-density effect on age-10 mortality |
| \(\beta_{\mathrm{dens},b}\) | Birth-density effect on senescence |
| \(\beta_{\mathrm{cohort}}\) | Residual effect of birth decade |
| \(\beta_F\) | Log hazard ratio during high/regular provisioning |
| \(\beta_M\) | Log hazard ratio during the FFRC post-2016 period |
| \(\sigma_{\mathrm{period}}\) | Magnitude of residual site-specific temporal trends |
| \(\sigma_{\mathrm{island},a}\) | Persistent residual heterogeneity in baseline mortality among sites |
| \(\alpha_p\) | Detection intercept |
| \(\beta_{p,\mathrm{sex}}\) | Sex effect on detection odds |
| \(\beta_{p,\mathrm{age}}\) | Age effect on detection up to age 15 |
| \(\sigma_{\mathrm{island},p}\) | Persistent site heterogeneity in detection |
| \(\sigma_{\mathrm{survey},p}\) | Field-visit heterogeneity in detection |

---

# 28. Reference individual

The reference Gompertz parameters correspond to an individual with:

- female sex (\(s=0\))
- birth density \(z=0\), i.e. the centered within-site mean
- birth cohort centered at 1990
- zero residual site mortality effect

Calendar multipliers are applied separately.

In a non-high-provisioning, non-Matthew year at the 2010 period reference, the reference hazard at age 10 is therefore \(a_{10}\), before any non-zero site-period or site mortality deviation.

---

# 29. MCMC fitting

The final primary model was fitted in NIMBLE using three independent chains.

Each chain used:

- **150,000 iterations**
- **50,000 burn-in**
- **thin = 10**
- **10,000 retained samples**

giving **30,000 retained posterior draws** across the three chains.

Birth time and lifespan were sampled jointly with adaptive factor slice samplers to improve mixing in the strongly coupled latent life-history parameters.

The three chains were run independently with separate seeds.

---

# 30. Convergence of the final primary model

The completed long run showed:

- no monitored parameter with \(\hat R>1.05\)
- most \(\hat R\) values at approximately 1.00–1.01
- effective sample sizes ranging from several hundred to many thousands for the main parameters

The model was therefore considered sufficiently converged for substantive interpretation.

---

# 31. Main posterior estimates from the final primary run

The following are posterior medians with 95% credible intervals unless otherwise stated.

| Parameter | Median | 95% CrI |
|---|---:|---:|
| \(a_{10}\) | 0.0383 | 0.0268–0.0510 |
| \(b\) | 0.0450 | 0.0313–0.0673 |
| \(\beta_{\mathrm{sex},a}\) | -0.0699 | -0.258–0.118 |
| \(\beta_{\mathrm{sex},b}\) | 0.685 | 0.451–0.982 |
| \(\beta_{\mathrm{dens},a}\) | 0.371 | 0.236–0.513 |
| \(\beta_{\mathrm{dens},b}\) | -0.322 | -0.471 to -0.173 |
| \(\beta_{\mathrm{cohort}}\) | 0.0742 | -0.0564–0.317 |
| \(\beta_F\) | 0.403 | -0.047–0.875 |
| \(\beta_M\) | -0.383 | -1.70–0.797 |
| \(\sigma_{\mathrm{period}}\) | 0.163 | 0.0124–0.586 |
| \(\sigma_{\mathrm{island},a}\) | 0.267 | 0.0967–0.615 |
| \(\alpha_p\) | -2.75 | -2.98 to -2.53 |
| \(\beta_{p,\mathrm{sex}}\) | 0.322 | 0.254–0.390 |
| \(\beta_{p,\mathrm{age}}\) | 1.39 | 1.23–1.55 |
| \(\sigma_{\mathrm{island},p}\) | 0.683 | 0.282–1.47 |
| \(\sigma_{\mathrm{survey},p}\) | 1.37 | 1.21–1.56 |

---

# 32. Interpretation of the sex effect

The age-10 sex coefficient is close to zero:

\[
\beta_{\mathrm{sex},a}
\approx -0.07.
\]

Thus there is little evidence that males and females differ strongly in mortality around age 10.

However:

\[
\beta_{\mathrm{sex},b}
\approx 0.685,
\]

so the male senescence parameter is approximately:

\[
e^{0.685}
\approx
1.98
\]

times the female value for otherwise equivalent reference animals.

This produces increasingly divergent mortality with age.

Using posterior-median reference values, the approximate male:female hazard ratios are:

| Age | Male:female hazard ratio |
|---:|---:|
| 10 | 0.93 |
| 12 | 1.02 |
| 20 | 1.45 |
| 30 | 2.26 |
| 40 | 3.52 |

### Plain-language interpretation

The primary sex signal is not that males simply have higher mortality throughout adulthood.

Instead, males and females have similar mortality around early adulthood, but **male mortality accelerates more rapidly with age**.

This is evidence for faster actuarial senescence in males.

---

# 33. Interpretation of birth density

For a one-standard-deviation increase in the supported within-site birth-density covariate:

\[
e^{\beta_{\mathrm{dens},a}}
=
e^{0.371}
\approx
1.45.
\]

Thus age-10 mortality is approximately 45% higher at the reference age.

At the same time:

\[
e^{\beta_{\mathrm{dens},b}}
=
e^{-0.322}
\approx
0.72.
\]

Thus the Gompertz senescence parameter is lower.

### Plain-language interpretation

Individuals associated with higher-density birth periods show a combination of:

- higher mortality around age 10
- slower subsequent acceleration of mortality

Because density affects both \(a_{10}\) and \(b\), there is no single density hazard ratio that applies at every age.

This pattern should not yet be assigned a mechanism. The primary density analysis intentionally uses only fully supported birth windows, and a less-conservative density-support sensitivity analysis is being run.

---

# 34. Interpretation of provisioning

The final primary estimate was:

\[
\mathrm{HR}_{\mathrm{high\ provisioning}}
=
e^{\beta_F}
\approx
1.50.
\]

The 95% credible interval was:

\[
0.954\text{–}2.40.
\]

The posterior probability that the hazard ratio exceeds 1 was:

\[
P(\mathrm{HR}>1)=0.964.
\]

### Plain-language interpretation

The model assigns approximately a **96% posterior probability that mortality hazard is higher during periods of high/regular provisioning**, with a posterior median corresponding to about a **50% higher instantaneous hazard**.

The central 95% credible interval still slightly overlaps 1, so the result is better described as strong posterior evidence for a positive association than as a definitive fixed 50% increase.

The inference is also not equivalent to a randomized feeding experiment.

The strongest within-site before/after information comes from Bitter, with additional but weaker transition information from Leaf South, while Leaf North and U Cay contribute high-provisioning populations without an observed within-study transition.

---

# 35. Interpretation of the Matthew term

For FFRC after 2016:

\[
\mathrm{HR}_{\mathrm{Matthew}}
\approx
0.68
\]

with a very wide 95% credible interval:

\[
0.18\text{–}2.22.
\]

The posterior probability of an increase in mortality was only approximately:

\[
P(\mathrm{HR}>1)=0.279.
\]

### Plain-language interpretation

The available data do not provide a precise estimate of a post-Matthew mortality effect.

The coefficient is retained as an adjustment for an important known disturbance, but it should not be interpreted as evidence that Matthew reduced mortality or as a clean causal estimate of hurricane impact.

---

# 36. Interpretation of birth cohort

The posterior for the residual birth-cohort effect overlapped zero:

\[
\beta_{\mathrm{cohort}}
=
0.074
\quad
(95\%\,\mathrm{CrI}\approx -0.056,\ 0.317).
\]

Thus there is no strong residual evidence that individuals born in later decades have systematically higher or lower mortality after accounting for density, sex, site and the other modeled components.

---

# 37. Interpretation of detection

The fitted detection model showed:

\[
\beta_{p,\mathrm{sex}}
\approx
0.322,
\]

so males had higher detection odds than females:

\[
e^{0.322}
\approx
1.38.
\]

The age effect was strongly positive:

\[
\beta_{p,\mathrm{age}}
\approx
1.39.
\]

Detection therefore increases substantially from young animals to age 15, after which the modeled age effect plateaus.

There is also substantial variation among field visits:

\[
\sigma_{\mathrm{survey},p}
\approx
1.37.
\]

### Plain-language interpretation

Capture probability is far from constant.

It varies with sex, age, island and the particular field visit. Explicitly modeling this heterogeneity prevents variation in field conditions or animal catchability from being mistaken for mortality.

---

# 38. What the model does **not** claim

The model does **not** assume that:

- every unseen animal is dead
- first capture is birth
- all individuals were the same age when entering the study
- provisioning is a fixed individual trait
- all sites have identical baseline mortality
- all field visits have equal detectability
- density is known for birth years outside observed support
- high provisioning was experimentally randomized among islands
- Hurricane Matthew is a controlled natural experiment
- targeted removal of large males can be estimated from these data

The large-male-removal effect was deliberately omitted because the available documentary evidence was not strong enough to construct a defensible site × year exposure history.

---

# 39. Strengths of the final formulation

The main advances over the earlier analysis are:

### Correct biological age under delayed entry

The Gompertz ageing clock is no longer restarted at first capture.

### Explicit uncertainty in birth year

Age uncertainty propagates through the survival likelihood rather than being replaced by a single age estimate.

### Separate mortality level and senescence effects

Sex and density can influence \(a_{10}\) and \(b\) independently.

### Calendar-varying provisioning

An animal's exposure can change through its lifetime rather than being assigned permanently as "fed" or "unfed".

### Protection against temporal confounding

Centered site-specific period slopes allow unrelated temporal drift instead of forcing all before/after changes into the provisioning coefficient.

### Explicit detection process

Survival and recapture are modeled jointly.

### Actual field-visit occasions

Distinct field trips separated by long gaps are no longer collapsed into one annual occasion.

### Conservative density support

The primary density effect is estimated only where the possible birth interval is fully supported by the density history.

---

# 40. Remaining assumptions and limitations

## Gompertz form

The model assumes that, conditional on covariates, adult mortality follows an exponential Gompertz age trajectory.

## Proportional calendar effects

Provisioning, Matthew and residual site-period effects act multiplicatively on the hazard. The primary model does not allow the provisioning effect itself to vary with age or sex.

## Site-level provisioning exposure

Provisioning is represented at site × year level. Individual feeding behavior and within-site spatial exposure are not directly modeled.

## Pre-entry survival is conditioned away

The model estimates mortality after first capture. It cannot use animals that died before ever being observed, and it does not estimate how provisioning before first capture affected the probability of entering the sample.

## Density reconstruction

Density is derived from MNA histories rather than a perfect census and is incompletely supported in early historical periods.

## Growth calibration transfer

Reliable known-age calibration is concentrated in the northern *inornata* populations. Applying the resulting age-size information to *figginsi* therefore entails a transferability assumption that should be examined in sensitivity analysis.

## Detection simplification

The age effect on detection is linear to age 15 and flat thereafter. More complex age or seasonal detection structures are possible, but were not included in the primary model.

---

# 41. Validation and sensitivity analyses

The primary model is now frozen.

The following analyses are intended to assess robustness rather than to search for a preferred result:

1. **Any provisioning vs none**  
   Replace the primary high/regular-versus-other exposure with any provisioning versus no provisioning.

2. **Acute Matthew period**  
   Replace FFRC post-2016 with an acute 2016–2018 disturbance window.

3. **Partial density support**  
   Relax the requirement that the entire birth interval has density support.

4. **Growth/birth-window sensitivity**  
   Examine the assumption that the northern known-age SVL relationship can be transferred to the southern *figginsi* populations.

5. **Posterior predictive checks**  
   Test whether replicated datasets generated from the fitted model reproduce important features of the observed recapture histories.

These are validation analyses. The primary model specification is not being altered in response to whether those alternatives produce a more desirable result.

---

# 42. Short lay summary

We followed more than four thousand iguanas whose exact ages and death dates were often uncertain. The model allows each animal to have an uncertain birth date and asks how its risk of death changes as it ages.

Mortality is modeled using a Gompertz curve, meaning that the risk of death can rise increasingly rapidly at older ages. Crucially, this curve is tied to the animal's true age rather than to the number of years since it was first caught.

The model also recognizes that we do not catch every living iguana on every field trip. It separately estimates capture probability using the animal's sex and age, the island, and the particular field visit.

Sex and density around birth can alter both mortality in early adulthood and how quickly mortality accelerates later in life. Provisioning is treated differently: it changes through calendar time. If an animal lived on an island before and after high-intensity feeding began, the model treats it as unexposed before the change and exposed afterwards.

The final analysis suggests that males do not have markedly greater mortality around age 10, but their mortality rises much faster as they age. Higher density around birth is associated with higher mortality around age 10 but slower subsequent senescence. There is also substantial posterior support for higher mortality during high/regular provisioning periods: the median estimated hazard was about 1.5 times that during other periods, although uncertainty remains and the analysis is observational rather than experimental.

---

# 43. Suggested concise Methods description

> We fitted a Bayesian hierarchical survival–recapture model in NIMBLE. Age-specific mortality followed a Gompertz hazard parameterized at age 10, with individual birth times treated as latent within empirically derived birth intervals. The survival distribution was correctly left-truncated at each individual's age at first capture, thereby conditioning the likelihood on survival to study entry without resetting the Gompertz age scale. Sex and within-site density at the inferred birth year were allowed to affect both age-10 mortality and the Gompertz senescence parameter, while birth cohort and a centered site random effect affected age-10 mortality. Calendar-varying high/regular provisioning and the FFRC post-2016 period entered as proportional-hazard multipliers, alongside centered site-specific temporal slopes. Detection was modeled separately on reconstructed field visits, defined by gaps of >14 days between capture bouts, using sex, an age effect increasing to age 15, centered site effects and centered field-visit random effects. First-capture occasions and all pre-entry non-detections were excluded from the detection likelihood. Death times were interval- or right-censored as appropriate, and sex was treated as latent for individuals with conflicting or absent sex assignments. Three independent MCMC chains of 150,000 iterations were run, with 50,000 burn-in iterations and thinning by 10.

---

# 44. Suggested concise provisioning interpretation

> High/regular provisioning was associated with a posterior median mortality hazard ratio of 1.50 relative to other site-years (95% CrI 0.95–2.40), with a 0.964 posterior probability that the hazard ratio exceeded one. Because provisioning was modeled as a site- and calendar-specific exposure, individuals could move from unexposed to exposed periods as provisioning histories changed. The estimate should nevertheless be interpreted as an observational association rather than a randomized treatment effect, particularly because much of the within-site temporal contrast is provided by Bitter and, to a lesser extent, Leaf South.

---

# 45. Conceptual takeaway

The final model separates four processes that were previously easier to conflate:

\[
\boxed{
\text{biological ageing}
+
\text{fixed individual history}
+
\text{calendar-varying environment}
+
\text{observation process}
}
\]

Specifically:

\[
\boxed{
\text{Gompertz age trajectory}
}
\]

captures ageing;

\[
\boxed{
\text{sex + birth density + birth cohort + site}
}
\]

capture fixed or early-life differences;

\[
\boxed{
\text{provisioning + Matthew + site temporal deviation}
}
\]

capture calendar-varying mortality conditions; and

\[
\boxed{
\text{sex + age + site + field visit}
}
\]

capture the probability of observing a living animal.

That separation is the central statistical and biological rationale of the final analysis.
