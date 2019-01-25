import sys
import surface
import batch
import compare

todo = sys.argv[1]
if todo == 'surface-run':
  surface.run_sims(sys.argv[2:])
if todo == 'surface-plot':
  surface.make_plots()
if todo == 'batch':
  batch.run_sims(sys.argv[2:])
if todo == 'compare-growth':
  compare.growth()
if todo == 'compare-hetero':
  compare.hetero()
if todo == 'compare-turnover':
  compare.turnover(*sys.argv[2:])
if todo == 'compare-tpaf':
  compare.tpaf()
print('done')