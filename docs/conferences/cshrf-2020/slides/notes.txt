# Background

- So while of course we may be preoccupied with a certain other pandemic at the moment, sexually transmitted infections continue to represent a major global health burden, including over 1 million new infections per day, with downstream health implications like adverse birth outcomes, increased cancer risk, and increased risk of HIV, which is also mostly sexually transmitted.

- One tool we use to inform STI intervention priorities is mathematical transmission modelling, meaning simulations of populations and their sexual behaviour.

- Based on transmission modelling, we can quantify the proportion of new infections which stem, directly or indirectly, from unmet treatment and prevention needs of various risk groups.

- Such risk groups could be identified by a combination of factors, but many definitions include the number of types of partnerships formed each year, such as married couples, or sex workers.

# Key Modelling Concepts

- two key modelling concepts underpin this work

- the first is Turnover, which is movement of people *between* risk groups. So typically we simulate entry into and exit from sexual activity, but many models do not consider movement of individuals between risk groups.

- for example, individuals may enter into sex work from a lower risk state, and then retire from sex work back into a lower risk state before cessation of sexual activity.

- the second concept is assortative mixing, which is the idea that individuals tend to form partnerships with other people in their risk group. Some models assume that individuals form partnerships randomly, whereas data in many contexts suggest some degree of assortative mixing.

- And in general, we will use dotted arrows to indicate movement of individuals between risk groups, and solid arrows to indicate formation of partnerships between or within risk groups.

# Research Questions

So our research questions are really about how turnover influences two model outputs: 1. equilibrium or steady state STI prevalence, and 2. the tPAF of the highest risk group, that is, the proportion of future infections that could be averted if the highest risk group had perfect treatment and prevention. 

And then we additionally compare the influence of turnover on these two outputs under two different mixing assumptions: random mixing, and assortative mixing.


# Methods

- To answer these questions, we developed a simple susceptible, infectious, recovered (S-I-R) model --- which some of you may be slightly more familiar with now due to COVID-19 --- and we stratified the population by not only these three health states, but three risk groups: a large low risk group, a medium sized medium risk group, and a small high risk group.

- We simulated turnover between the risk groups such that the risk groups did not change size over time, and that the overall amount of turnover could be controlled by a single parameter, which was the duration in the highest risk group.

- Using this framework we then examined the influence of turnover on equilibrium STI prevalence, under both random and assortative mixing conditions.

---

- Next, we calibrated the model to represent a specific epidemic context, meaning STI prevalence in each risk group, by varying the numbers of partnerships formed per year by individuals in each risk group.

- Such model calibration is common when applying epidemic models to answer context-specific research questions, and represents refining uncertain model inputs to reproduce observed data in the outputs, like incidence or prevalence.

- So we calibrated four model variants to the same epidemic context:
  - with versus without turnover
  - and under random versus assortative mixing

- then, finally, for each fitted model, we computed the tPAF, the projected relative importance of the highest risk group.

---
Results
---

# Random mix -- turnover homogenizes STI prevalence

- Our first result is that turnover reduces the ratio of STI prevalence in the highest vs lowest risk groups under random mixing.

- So on the x-axis is turnover, which increases to the right, controlled by duration in the highest risk group. And on the y-axis is the prevalence ratio for the high and low risk groups.

- We can see that increasing turnover reduces the ratio, which can be explained by comparing two simplified models:
  - On the left, without turnover, we have random partnerships between all groups, but infections remain concentrated in the high risk group due to higher total numbers of partnerships in the group.
  - On the right, with turnover, many infections are still acquired in the high risk group, but some of these may move into the low risk group via turnover of infected individuals, which then increases the prevalence in the low risk group, and also decreases prevalence in the high risk group, since individuals entering the high risk group via turnover are more likely uninfected.

- Again, solid arrows denote partnerships, and dotted arrows denote turnover. And please note we've omitted the middle risk group here for simplicity, but it remains in the model.

# Random mix - infer larger risk ratio with turnover → higher tPAF+

- Next, we fit the random mixing model with and without turnover to the same prevalence targets, which included a prevalence ratio of 6.7 between the high and low risk groups.

- In order to reproduce this target ratio, the ratio of partnerships per year between the high and low risk groups needed to be higher with turnover than without, because a greater risk ratio was required to overcome the homogenizing effect of turnover. So after calibration, the ratio of partners per year was only 15.2 without turnover versus 23.9 with turnover.

- Finally, calculating the tPAF of the high risk group over 10 years after model fitting, we find a larger tPAF in the model with turnover than in the model without, because we inferred a greater proportion of transmission associated with the high risk group through model fitting. The difference is about 5%.

- We should note that we represent the risk ratio using partnerships per year, but such heterogeneity could also be associated with other biological, environmental, or individual factors.

- Now let's turn to the assortative mixing case.

# Assort mix - turnover allows infections to "escape" sexual networks

- For our second comparison, we plot the same high vs low risk prevalence ratio, versus increasing turnover, now under assortative mixing. We can see that the decline in prevalence ratio with turnover is larger in the case of assortative mixing versus random mixing.

- Again, we will explain using two simplified models:
  - On the left, without turnover, infections are highly concentrated in the high risk group due a greater number of partnerships, but also more of those partnerships with individuals in the same risk group, making it difficult for infections to escape the concentrated sexual network.
  - On the right, with turnover, the net movement of infectious individuals from high to low risk allows infections to escape the concentrated sexual network. The relative impact of this movement is higher in the context of assortative mixing because the baseline level of risk concentration is higher without turnover.

# Assort mix - higher risk ratio & escaped infections → higher tPAF++

- Repeating the analysis from the random mixing case, we re-calibrated the model with and without turnover to the same STI prevalence targets.

- Again we infer greater risk ratio, or ratio of partnerships per year, with turnover than without.

- However, now the relative difference in estimated tPAF is much larger, on the order of 20% larger with turnover than without, and that's because not only have we inferred a greater risk ratio, but the number of onward infections that are able to escape the high risk sexual network has increased via turnover.

- We should note that the overall tPAF of the high risk group has decreased under assortative mixing, since infections do struggle to escape that sexual network, but turnover of infectious individuals provides a pathway for infections to escape which counteracts this phenomenon, in addition to same risk ratio mechanism relevant to random mixing.

# Implications

- So in sum there are two main implications to this work.

  1. influence of turnover on STI epidemics is larger under assortative vs random mixing, due mainly to the ability of infections to escape concentrated sexual networks via turnover.

  2. If turnover is not adequately captured in calibrated STI epidemic models, the the importance of meeting the needs of higher risk groups for epidemic control may be underestimated.

And these results may be also relevant in other non-STI epidemics.

Finally, in terms of limitations, we should note that all results here are conditional on the modelling assumptions we've made.

- Thank you for listening.


