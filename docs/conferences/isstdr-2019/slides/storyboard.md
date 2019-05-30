# The influence of risk group turnover in STI/HIV epidemics: mechanistic insights from transmission modeling

## From the Abstract

**Results:**

- Across the range of turnover and treatment parameters explored, turnover consistently decreased STI prevalence in the core group.

- In the low-risk group, turnover increased prevalence under low treatment rate, but had the opposite effect under high treatment rate.
// downplay, focus on incidence

- When calibrating to the same STI prevalence, fitted core group partner change rates were higher with turnover than without.
// ratio not number

- Using these fitted parameters, models with turnover then consistently projected a higher tPAF of the core group versus models without.

## Storyboard:

FYI: 12 min + 3 min questions; maximum of 10 slides

### 1. Background

<!-- - (1) Core group theory:
  - defined: risk heterogeneity affects:
    - R0 and initial epidemic speed
    - equilibrium prevalence
    - (other epidemic features)
  - examples of "core groups" (AKA high risk group)
    - ... -->

- (1) {CGT} & Turnover:
  {core group theory here, maybe}
  - defined: movement of individuals b/w risk groups
    - AKA: episodic risk, migration
- (2)
  - [fig]: turnover diagram
  - research questions:
    - influence of turnover on equilibrium prevalence and incidence?
    - influence of turnover on TPAF?

### 2. Methods

- (3) "illustrative" SIR model of STI transmission
  - {developed an approach to turnover, and applied to ...}
  - [fig]: SIR diagram
  - 3 risk groups
  - 1-sex
  - proportional mixing
  - no disease attributable mortality

- (4) Experiments: influence of turnover
  - 1. Fixed parameters:
    - a. vary overall magnitude of turnover -> equilibrium prevalence and incidence
    { b. vary duration of infectiousness }
  - 2. Fitted parameters:
       Fit contact rates to 25% / 5% prevalence in high / low risk groups with & without turnover
    - a. compare fitted contact rates
    - b. estimate TPAF with / out turnover 

### 3. Results

<!-- - (5) [1a] Equilibrium Prevalence in the core group
  - [fig]: high risk prevalence vs turnover
  - declines due to 2 effects:
    1. net movement of infected individuals from high to low; susceptibles from low to high
    2. average contact rate of infected individuals is reduced

- (6) [1a] Equilibrium Prevalence in non-core groups
  - [fig]: low risk prevalence vs turnover
  - at low turnover: prevalence increases (effect #1 dominates)
  - at high turnover: prevalence decreases (effect #2 dominates)
  - [fig]: overall prevalence vs turnover
  - overall prevalence similar to low risk prevalence

- (7) [1b] Influence of treatment
  - [fig]:
      a) high risk prevalence
      b) low risk prevalence
      c) overall prevalence
  - 3 regions:
    - 1. low treatment: same trends as above
    - 2. high treatment: no epidemic (R0 < 1)
    - 3. "boundary": epidemic barely persists
      - turnover decreases incidence and prevalence -->

- (8) [2a] Fitted risk parameters with turnover
  - C_H / C_L:
    - With turnover: 16.9 / 0.28 = 59.6
    - Without turnover: 15.8 / 2.49 = 6.33
  - Turnover "homogenizes risk"
  - Inferred heterogeneity in risk must be higher
    - Contact rates have to "work harder" to achieve same prevalence ratio

- (9) [2b] TPAF of high risk group:
  - before model fitting: higher without turnover than with turnover
    - lower equilibrium prevalence ratio with turnover
  - after model fitting: higher with turnover than without
    - equal equilibrium prevalence ratio (fitted)
    - higher equilibrium ratio of risk heterogeneity (e.g. contact rates)

### 4. Conclusion

- (10) Implications:
  - Modelling turnover is important
    - equilibrium incidence & prevalence
    - fitted parameters:
      - with turnover: inferred risk heterogeneity must be higher
    - Estimated TPAF of high risk group:
      - with turnover: higher
  - How to incorporate turnover:
    - https://github.com/c-uhs/turnover
    - paper preprint
