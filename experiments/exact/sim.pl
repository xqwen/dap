$n = $ARGV[0];

`Rscript sim_pheno.R $n`;

@pops= ("ceu","fin","gbr","tsi","yri");

$g = "gene";
my $con = "";
foreach $p (@pops){
    $file = "pheno_data/gene.$p.pheno";
    open FILE, "$file";
    my @data = <FILE>;
    chomp @data;
    $con .= sprintf "pheno $g $p @data\n";
}

open OUT, ">$g.sbams.dat";
print OUT "$con";
close OUT;
`cat geno_data/$g.meta.geno  >> $g.sbams.dat`;
