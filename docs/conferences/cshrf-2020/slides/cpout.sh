turnover=/home/jesse/mishra/proj/turnover
figs=$turnover/outputs/cshrf/asso/figs

for model in {main,asso,both}; do
  cp $figs/sensitivity/1d-${model}-ratio-prevalence-high-low-tau=0.1.pdf fig/prev-ratio-$model.pdf
  cp $figs/compare/asso-tpaf-tpaf-high-all-vs=${model}.pdf fig/tpaf-${model}.pdf
done

for model in {asso,main}; do
  data=$turnover/outputs/cshrf/$model/data/values;
  for fit in {-[fit]-,-}; do
    if [ ${#fit} -gt 1 ]; then fiti=f; else fiti=r; fi
    for tx in {turnover,no-turnover}; do
      txi=$(echo $tx | head -c 1);
      for out in {prevalence,C}; do
        outi=$(echo $out | head -c 1);
        for sub in {ratio-high-low,high,low}; do
          subi=$(echo $sub | head -c 1);
          cp $data/${tx}${fit}${out}-${sub}.txt data/$model/${txi}${fiti}-${outi}${subi};
        done
      done
      cp $data/${tx}${fit}tpaf-high-high.txt data/$model/${txi}${fiti}-tpaf;
    done
  done
done

