for (1..15){
    `perl sim_data.pl $ARGV[0]`;
    `perl logit_reg.pl >> sim_rst.best.lambda.$ARGV[0]`;
    `perl naive_enrich.pl >> sim_rst.naive.lambda.$ARGV[0]`;
    `perl em_fmap.pl -msize 1 | tail -1 >> sim_rst.ms1.lambda.$ARGV[0]`;
    `perl em_fmap.pl -it 0.05 | tail -1 >> sim_rst.msa.lambda.$ARGV[0]`;
}
