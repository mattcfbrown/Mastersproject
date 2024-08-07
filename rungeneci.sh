geneci optimize-ensemble --confidence-list ../GENECI/inferred_networks/test/lists/GRN_CLR.csv \
    --confidence-list ../GENECI/inferred_networks/test/lists/GRN_PIDC.csv \
    --confidence-list ../GENECI/inferred_networks/test/lists/GRN_GENIE3_RF.csv \
    --num-parents 3 --mutation-strength 0.1 \
    --population-size 100 \
    --num-evaluations 20 --cut-off-criteria PercLinksWithBestConf --cut-off-value 0.4 \
    --function Quality --algorithm GA \
    --threads 8 --no-plot-fitness-evolution --no-plot-pareto-front \
    --no-plot-parallel-coordinates --output-dir ../GENECI/geneci_consensus