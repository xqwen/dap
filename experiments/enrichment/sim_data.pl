
init();
sim_data();
assemble();


sub init {
    
    exit if !defined($ARGV[0]);
    
    @f = <truth/*.truth>;
    unlink @f;

    @f = <truth/*.annot>;
    unlink @f;
    
    @f = <sbams_data/*.dat>;
    unlink @f;
    

    $lambda = $ARGV[0];
}



sub sim_data {

    open OUT, ">batch_sim.cmd";
    open FILE, "gene.list";

    while(<FILE>){
        next if $_ !~ /(ENSG\d+)/;
        print OUT "Rscript sim_pheno.R $1 $lambda\n";
    }

    `openmp_wrapper -d batch_sim.cmd -t 10`;
}


unlink "batch_sim.cmd";;

# assembling
sub assemble {
    open OUT, ">batch_assemble.cmd";
    open FILE, "gene.list";
    while(<FILE>){
	next if $_ !~ /(ENSG\d+)/;
	print OUT "perl assemble.pl $1\n";
    }
    
    `openmp_wrapper -d batch_assemble.cmd -t 10`;
}
    
unlink "batch_assemble.cmd";
