$g = $ARGV[0];

@pops= ("ceu","fin","gbr","tsi","yri");


my $con = "";
foreach $p (@pops){
    $file = "pheno_data/$g.$p.pheno";
    open FILE, "$file";
    my @data = <FILE>;
    chomp @data;
    $con .= sprintf "pheno $g $p @data\n";
}

open OUT, ">sbams_data/$g.sbams.dat";
print OUT "$con";
close OUT;
`cat geno_data/$g.meta.geno  >> sbams_data/$g.sbams.dat`;
