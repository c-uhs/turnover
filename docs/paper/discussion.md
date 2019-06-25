# Discussion storyboard

refs: Koopman1997, Zhang2012*, Alam2013*, Romero-Severson2013*, Eaton2014, Henry2015* [*

### 1) summary
- new framework for turnover, with data inputs
- influence and mechanisms of turnover on group-specific equilibrium prevalence and incidence
- influence and mechanisms of turnover on inferred heterogeneity and projected TPAF

### 2) New Framework
- differences from prior work:
  - unifying framework for modellers
    - more than 2 risk groups
    - lists possible assumptions and force modeller to choose a combination which best reflects the 'epidemic context' in question (and available data)
      -> makes "definition of turnover explicit" (something like this)
      -> meaning w.r.t. real world is more transparent
      -> assumptions transparent
    - for example: if survey says duration in group is x years ->
    - could be used with other epidemics
    - end with: how to use this system?

### implications (modelling):
- Overall effects
  - turnover affects equilibrium incidence & prevalence
    - overall effect varies with context
    - turnover always decreases "core group" prevalence
    - overall: "risk homogenization"
  - Other works:
    - Stigum1994 ("Migration")
    - Zhang2012, Henry2015 ("Episodic Risk")
- Mechanisms:
  - turnover is "convection" to incidence's "conduction" (continuously mixed compartments)
    - movement of mass between compartments
    - vs: interaction
  - treatment (universal) to achieve epidemic control is reduced by turnover
    - also shown by Henry2015 (maybe)

### implications (epi)
  - TPAF "prior work" (HRG: ~%) , even from other STIs
    - TPAF implies intervention impact among HRG
    - TPAF -> importance for reaching with intervention
      "maxiumum proportion of infections averted with perfect intervention"
  - homogenizing effect of turnover:
    - 2 differnt contexts (with / without turnover):
      - TPAF higher without
    - 2 different models (with / without turnover):
      - risk underestimating heterogeneity -> TPAF (HRG) if we omit turnover
      - may underestimate impact of targeted interventions
        - Henry2015

### 2) limitations
  - all results subject to assumptions & parameters used here
    - proportional mixing
      - but: Henry2015: "major qualitative results still hold" under assort
    - SIR model
      - no "re-S" for STIs like
      - why important (possible direction)
    - no disease history
      - e.g. early infection (Henry2015)
    - no disease attributable mortality
      - challenges with constant group sizes
    - future work: all 3 of these features
  - focus on models at equilibrium
    - see Henry2015 for R0 analysis (i.e. early epidemic dynamics)
    - appendix?
