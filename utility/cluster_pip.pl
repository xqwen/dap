#!/usr/bin/perl

# define the path to cluster_pip.R here
$cluster_pip_R_path = "~/bin/cluster_pip.R";


$cluster_num = -1;
$gene_name = "";

for ($i=0;$i<=$#ARGV;$i++){
    
    if($ARGV[$i] eq "-d"){
	$file_name = $ARGV[++$i];
	next;
    }

    if($ARGV[$i] eq "-c"){
	$cluster_num = $ARGV[++$i];
	next;
    }

    if($ARGV[$i] eq "-g"){
	$gene_name = $ARGV[++$i];
	next;
    }

    if($ARGV[$i] eq "-Rpath"){
	$cluster_pip_R_path = $ARGV[++$i];
	next;
    }

}

if(!defined($cluster_pip_R_path)){

    if(-e "cluster_pip.R"){
	$cluster_pip_R_path = "cluster_pip.R";
    }
    else{
	print STDERR "Error: cannot locate cluster_pip.R\n";
	exit(1);
    }
}






exit if(! -e "$file_name");


open FILE, "$file_name";
$max_model = 0;
while(<FILE>){

    next if $_ !~ /\d/;
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;
    
    if(/\[/){
	next if /NULL/;

	
	my @snps = ($_ =~/\[(\S+)\]/g);

	if(scalar(@snps) > $max_model){
	    $max_model = scalar(@snps);
	}

	foreach $i (0..$#snps){
	    $s = $snps[$i];
	    if(!defined($rcd{$s})){
		$rcd{$s} = {};
		push @list, $s;
	    }
	    $rcd{$s}->{prob} += $data[1];
	    for $j (0..$#snps){
		$p = $snps[$j];
		if ($p ne  $s){
		    $rcd{$s}->{$p} +=  $data[1];
		}
	    }
	    
	}	
    }
    
    if($#data==3){
     	my @data = split /\s+/, $_;
        shift @data until $data[0]=~/^\S/;
    	$snp_bf{$data[1]} = $data[-1];
    }
}

if($#list<1){
    print STDERR "Less than 2 SNPs included in posterior samples ... skip\n";
    exit;
}

$tempfile = ".dist.".time() . "_" . rand();
open OUT, ">$tempfile";
print OUT "prob ";
foreach $s (@list){
    
    print OUT "$s ";
}
print OUT "\n";

foreach $s (@list){
    
    printf OUT "%7.3e ", $rcd{$s}->{prob};
    foreach $p (@list){
	printf OUT "%7.3e ",$rcd{$s}->{$p}
    }
    print OUT "\n";
}
close OUT;


$cluster_num = $max_model if($cluster_num<0);


$rst = `Rscript $cluster_pip_R_path $tempfile $cluster_num`;




#$tempfile2 = "pip_plot.".time() . "_" . rand();

#open OUT, ">$tempfile2";

my @lines = split /\n/, $rst;
my %out;

foreach $l (@lines){
    next if $l !~/\d/;
    my @data = split /\s+/, $l;
    shift @data until $data[0]=~/^\S/;
    if(!defined($output{$data[2]})){
	$output{$data[2]} = {};
    }
    #$output{$data[2]}->{$data[0]} = { con => sprintf("%15s  %8.5f  %d  %8.5f  %7.3f   %8.5f\n",$data[0],$data[1],$data[2],$data[3],$snp_bf{$data[0]},$snp_rbp{$data[0]}),
    $data[0]=~s/^X//;
    $output{$data[2]}->{$data[0]} = { con => sprintf("%15s  %8.5f  %d  %8.5f  %7.3f\n",$data[0],$data[1],$data[2],$data[3],$snp_bf{$data[0]}),
				      pip => $data[1] };
    
}

foreach $c (sort {$a <=> $b} keys %output){
    foreach $s (sort {$output{$c}->{$b}->{pip} <=> $output{$c}->{$a}->{pip}} keys %{$output{$c}}){
	
	$s=~/\.(\d+)/;
	#my @data = split /\s+/, $output{$c}->{$s}->{con};
	#shift @data until $data[0]=~/^\S/;
	#next if $data[0]<0.01;
	
	#print  OUT "$1 $output{$c}->{$s}->{con}";
	print "$output{$c}->{$s}->{con}";
    }
}

#close OUT;

#`Rscript pip_plot.R $tempfile2 $gene_name`;

#unlink $tempfile2;
unlink $tempfile;


