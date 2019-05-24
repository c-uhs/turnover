import sys
import sensitivity
import batch
# import compare
import variants
# import debug
# import newinf

todo = sys.argv[1]
if todo == 'compare-growth':
  variants.exp_growth(save=True)
if todo == 'compare-hetero':
  variants.exp_hetero(save=True)
if todo == 'compare-turnover':
  variants.exp_turnover(save=True)
if todo == 'compare-tpaf':
  variants.exp_tpaf(save=True)
if todo == 'sensitivity-run':
  sensitivity.run_sims(sys.argv[2:])
if todo == 'sensitivity-plot':
  sensitivity.gen_plots()
# if todo == 'newinf':
#   newinf.main()
# if todo == 'debug':
#   debug.main(sys.argv[2:])
