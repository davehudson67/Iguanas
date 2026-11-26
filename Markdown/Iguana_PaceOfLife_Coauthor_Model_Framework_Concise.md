# Iguana pace-of-life survival model

## Final audited framework and findings — 26 November 2026

This replaces the earlier candidate-model explanation. All planned primary, sensitivity and population-exclusion fits are complete. The companion technical explanation gives the full mathematical specification.

## 1. What are we asking?

Do males and females differ in adult mortality and actuarial senescence? Do conditions around birth, especially relative population density, predict subsequent mortality trajectories? Is mortality associated with tourist provisioning or the period following Hurricane Matthew?

The data comprise 4,152 individuals across nine populations through 2019, including 74 confirmed deaths and 4,078 right-censored individuals. Sex was observed for 1,879 females and 1,983 males; 290 animals had missing or conflicting assignments and their sex was treated as latent.

## 2. How the model works

The Gompertz model separates mortality hazard at age 10 from the rate at which log mortality hazard increases with age. Age 10 is a reference age, not the onset of adulthood or a biological threshold. Sex and within-site birth density can affect both components. Birth cohort and residual site differences affect the reference hazard.

Birth times are uncertain for many animals. Instead of using one estimated age, the model samples birth dates within empirically derived intervals. It conditions survival on reaching first live capture without resetting biological age. Confirmed deaths and right-censoring constrain latent death times.

Failure to recapture an animal is not assumed to mean death. The model jointly estimates detection at reconstructed field visits, accounting for sex, age (the log-odds age contribution plateaus at 15), site, survey occasion and persistent individual catchability. There are 170 reconstructed occasions and 64,195 post-entry individual-by-visit opportunities.

Density is a standardized within-site birth-year index derived from minimum-number-alive histories, rather than a comparison of absolute abundance across islands. Only fully supported birth windows inform the primary density coefficient. Calendar-varying provisioning and the FFRC post-2016 period multiply mortality hazard when the exposure applies. Centered site-specific calendar slopes help accommodate unrelated temporal changes. Provisioning is included in the final primary model.

## 3. Main results

Posterior medians and 95% credible intervals:

| Question | Estimated association | Interpretation |
|---|---|---|
| Male versus female mortality at age 10 | HR 1.03 (0.83–1.27) | Little evidence of a difference at the reference age |
| Male versus female Gompertz rate | Ratio 1.74 (1.41–2.27) | Faster age-related mortality acceleration in males |
| Birth density, +1 SD: age-10 hazard | HR 1.35 (1.17–1.57) | Higher reference mortality |
| Birth density, +1 SD: Gompertz rate | Ratio 0.79 (0.69–0.90) | Slower subsequent acceleration |
| High/regular provisioning vs other periods | HR 1.91 (1.16–3.07) | Higher instantaneous mortality hazard; posterior probability 0.990 |
| FFRC post-Matthew period | HR 0.88 (0.24–3.05) | Highly uncertain |
| Birth cohort | Log HR 0.074 (−0.072–0.361) | No clear residual cohort effect |

Ratios for sex and density are exponentiated coefficients. A Gompertz-rate ratio is not a mortality hazard ratio at every age, and a hazard is not a probability of death. Density therefore does not have one constant effect on survival across ages.

The final results also support substantial individual detection heterogeneity (SD 1.11, 1.04–1.19 on the log-odds scale).

## 4. How robust are these results?

We fitted 14 scenarios: the primary model, alternative provisioning classification, partial birth-density support, an acute post-Matthew period, exclusion of Leaf North and U Cay together, and nine separate leave-one-population-out analyses. Each used three long chains; none had a monitored Rhat above 1.05.

**Sex:** All nine single-population omissions retained a positive male Gompertz-rate effect with credible intervals excluding zero. The result is not driven by Leaf North or U Cay individually. However, removing both together changes the result: the male log rate ratio becomes −0.105 (−0.737–0.479), while the male age-10 log hazard ratio becomes 1.063 (0.443–1.734), equivalent to an age-10 HR of about 2.90. The remaining populations support higher male mortality at age 10 but not a clear difference in senescence rate. We cannot determine whether the change reflects information loss, genuine population differences or undocumented removal.

**Density:** Both effects are stable under partial density support. Their directions persist in every single-population omission, but both intervals include zero without Leaf North; the senescence association is also uncertain after the joint exclusion. We should acknowledge dependence on information from Leaf North.

**Provisioning:** Median HRs range from 1.68 to 2.44 across all analyses, with posterior probabilities of increased hazard of 0.949–0.999. Several credible intervals include one. Classifying any provisioning versus none gives HR 1.85 (0.90–3.31), a different exposure contrast. The primary high/regular association is consistent but not a causal demonstration.

**Matthew:** The acute 2017–2018 contrast gives HR 0.90 (0.22–3.12). The data do not provide a precise estimate. When FFRC is excluded there is no exposed population, so that estimate is not reported.

## 5. Could human removal explain the male result?

The descriptive audit did not identify Leaf North or U Cay as uniquely unusual in male disappearance proportions among animals seen during 2003–2007. Both sites showed a modest male excess, comparable to Gaulin. But disappearance does not establish cause, and we lack a defensible site-year removal history. The targeted and individual exclusions address influence, not whether removal actually occurred. Permanent disappearance cannot reliably be separated from death by these encounter histories alone.

## 6. How well does the model fit?

Posterior predictive checks on 961 known-age, resolved-sex individuals showed reasonable reproduction of several encounter summaries, but underprediction of long recapture spans (at least five and ten years) and U Cay detections. These are genuine limitations. The stability of survival coefficients after removing U Cay is reassuring but does not prove the detection model is adequate in every respect.

## 7. What should the paper say?

The primary analysis supports faster male actuarial senescence, an association between higher birth density and higher reference mortality but slower senescence, and higher mortality hazard during high/regular provisioning. The sex decomposition is sensitive to joint exclusion of two large populations; density inference loses precision without Leaf North. All findings are observational associations. No causal provisioning, density or removal mechanism is established.

The model conditions on survival to entry and is not a direct estimate of juvenile survival. Standardized curves at very young or old ages may extrapolate outside strong data support. Constant-exposure trajectories illustrate fitted contrasts, not real lifelong feeding histories.

## 8. Publication and reproducibility

The final primary posterior is Model_Final_IndivDetect_Audited_3CHAINS_COMPLETE.rds (scripts 37–38). The complete supplementary robustness table and figure come from scripts 51–52; the final audited manuscript figures come from scripts 53–55. Earlier results and candidate-model descriptions are superseded. Final manuscript figure inspection and numerical/text consistency checks remain before submission.
