init();
$counter = 1;
$old_logBF = -99999;
$logBF = 0;
while(1){
    printf "  %3d        ",$counter;
    E_step();
    M_step();
    printf "    %9.3f    %9.3f       %9.3f\n", $mu, $lambda,$logBF;
    last if abs($logBF-$old_logBF)<=0.1;
    $old_logBF = $logBF;
    $counter++;

}

get_CI();





sub E_step{

    open FILE, "gene.list";
    open OUT, ">batch_dap.cmd";
    my $count = 0;
    while(<FILE>){
	next if $_ !~ /(ENSG\d+)/;
	if($msize==-1){
	    print OUT "dap -d sbams_data/$1.sbams.dat -g grid -p prior/$1.prior -msize 2 -all  -it $it -st $st  > dap_rst/$1.rst 2>/dev/null\n";
	}else{
	    print OUT "dap -d sbams_data/$1.sbams.dat -g grid -p prior/$1.prior -msize $msize -all -it $it -st $st  > dap_rst/$1.rst 2>/dev/null\n";
	    
	}
	$count++;
    }  
    #print STDERR "First round running ...\n";
    `openmp_wrapper -d batch_dap.cmd -t 10 -all`;
    #print STDERR "Checking insufficient analysis ... \n";
   
    # Adaptive run
    if($msize==-1){

	$out = `grep Warning dap_rst/*.rst`;
	@list = ($out=~/(ENSG\d+)/g);
	$n = scalar(@list);
	#printf STDERR "%d completed, %d need re-analysis\n",$count-$n, $n;
	#print STDERR "Second round re-analysis ... \n";
	foreach $g (@list){
	    #printf STDERR "re-analyzing $g\n";
	    `dap -d sbams_data/$g.sbams.dat -g grid  -p prior/$g.prior  -t 10 -all -it $it -st $st > dap_rst/$g.rst 2>/dev/null`;
	}
    }
    
    #print STDERR "dap completed\n";
	
    $out = `grep LogBF dap_rst/*.rst`;
    # LogBF = 60.34384 
    my @rst = ($out=~/LogBF\s+\=\s+(\S+)\s+/g); 
    $logBF = 0;
 ;
    foreach $v (@rst){
	$logBF += $v;
    }
    #unlink "batch_dap.cmd";
}




sub M_step {

    my %rcd;
    my @list;
    my @files = <truth/*.annot>;
    foreach $f (@files){

        open FILE, $f;
	$f =~ /(ENSG\d+)\./;
	$g = $1;
        while(<FILE>){
            next if $_ !~ /^\s*(\S+)\s+/;
            $snp = "$g\_$1";
            chomp;
            $rcd{$snp}->{header} = sprintf "%15s $_",$g;
            $rcd{$snp}->{inp} = 0;
            push @list, $snp;
        }

    }

    @files = <dap_rst/*.rst>;
    foreach$f (@files){

	open FILE, $f;
        $f =~ /(ENSG\d+)\./;
        $g = $1;
	while(<FILE>){
            next if $_ !~ /\(c\)/;
            my @data = split /\s+/, $_;
            shift @data until $data[0]=~/^\S/;
            $snp = "$g\_$data[1]";
            chomp;
            $rcd{$snp}->{inp} = $data[2];
	}

    }

    open OUT, ">regress.dat";
    foreach $snp (@list){
        print OUT "$rcd{$snp}->{header}  $rcd{$snp}->{inp}\n";
    }
    
    $out = `Rscript enrich_m.R 2>/dev/null`;
    $out =~/^\S+\s+(\S+)\s+(\S+)/;
   
    # update global variables
    $mu = $1;
    $lambda = $2;

    open FILE, "prob.out";
    my %rcd;
    while(<FILE>){

        next if $_ !~ /\d/;
	my @data = split /\s+/, $_;
        shift @data until $data[0]=~/^\S/;

        my $gene = $data[0];
        shift @data;
        $rcd{$gene} .= "@data\n";

    }

    foreach (keys %rcd){
	open OUT, ">prior/$_.prior";
        print OUT $rcd{$_};
        close OUT;
    }
    
    unlink "regress.dat", "prob.out";

   
}




sub init{


    $msize = -1;
    $pi1 = 1e-3;
    $st = 0.01;
    $it = 0.02;

    for $i (0..$#ARGV){
	if($ARGV[$i] eq "-msize" && $ARGV[$i+1]>=1){
	    $msize = $ARGV[++$i];
	    next;
	}

	if($ARGV[$i] eq "-pi1"){
	    $pi1 = $ARGV[++$i];
	    next;
	}
	
	if($ARGV[$i] eq "-st"){
	    $st = $ARGV[++$i];
	    next;
	}

	if($ARGV[$i] eq "-it"){
	    $it = $ARGV[++$i];
	    next;
	}



    }
    
    



    my @pfiles = <prior/*.prior>;
    if($#pfiles>=0){
	unlink @pfiles;
    }
    
    my @files = <truth/*.annot>;
    foreach $f (@files){
	
	open FILE, $f;
        $f =~ /(ENSG\d+)\./;
        open OUT, ">prior/$1.prior";
	
	while(<FILE>){
            next if $_ !~ /^\s*(\S+)\s+/;
	    print OUT "$1  $pi1\n";
        }
	close OUT;
    }

    

}


sub get_CI {
    
    $l1 = $logBF; 
    my $prob = exp($mu)/(1+exp($mu));
    my @files = <truth/*.annot>;
    foreach $f (@files){

        open FILE, $f;
        $f =~ /(ENSG\d+)\./;
        $g = $1;
	open OUT, ">prior/$g.prior";
        while(<FILE>){
	    
	    my @data = split /\s+/, $_;
	    shift @data until $data[0]=~/^\S/;
	    printf OUT "$data[0] %7.3e\n",$prob; 
            
	}
	close OUT;

    }
    E_step();
    
    $l0 = $logBF;
    $se = abs($lambda)/sqrt(2*abs($l1 -$l0));
    #printf "\n%9.3f  %9.3f   %.3f   %.3f  ", $mu, $lambda, $l1,$l0; 
    #printf "\n%9.3f  %9.3f   [%.3f, %.3f]\n", $mu, $lambda, $lambda-1.96*$se,$lambda+1.96*$se; 
    #printf "   [%.3f, %.3f]\n", $lambda-1.96*$se,$lambda+1.96*$se; 
    printf "\n%9.3f  %9.3f   %9.3f   %9.3f\n", $mu, $lambda, $lambda-1.96*$se, $lambda+1.96*$se; 
}


