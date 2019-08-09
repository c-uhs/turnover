import sys
import sensitivity
import variants
import flows

todo = sys.argv[1]
if todo == 'compare-turnover':
  variants.exp_turnover()
if todo == 'fit':
  variants.run_fit()
if todo == 'compare-tpaf':
  variants.exp_tpaf()
if todo == 'sensitivity-run':
  sensitivity.run_sims(sys.argv[2])
if todo == 'sensitivity-plot':
  sensitivity.gen_plots()
if todo == 'flows':
  flows.make_figs()
if todo == 'debug': 
  sensitivity.run_sims('debug')
