import sys
import sensitivity
import variants

todo = sys.argv[1]
if todo == 'compare-growth':
  variants.exp_growth(save=True)
if todo == 'compare-hetero':
  variants.exp_hetero(save=True)
if todo == 'compare-turnover':
  variants.exp_turnover(save=True)
if todo == 'fit':
  variants.run_fit(save=True)
if todo == 'compare-tpaf':
  variants.exp_tpaf(save=True)
if todo == 'sensitivity-run':
  sensitivity.run_sims(sys.argv[2:])
if todo == 'sensitivity-plot':
  sensitivity.gen_plots()
