
`openmp_wrapper -d batch_fsingle.cmd -t 10`;

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
	$rcd{$snp}->{gamma} = 0;
	push @list, $snp;
    }
    
}

@files = <fsingle_rst/*.rst>;
foreach$f (@files){
    
    open FILE, $f;
    $f =~ /(ENSG\d+)\./;
    $g = $1;
    while(<FILE>){
	next if $_ !~ /\d/;
	my @data = split /\s+/, $_;
	shift @data until $data[0]=~/^\S/;
	$snp = "$data[0]\_$data[1]";
	chomp;
	$rcd{$snp}->{pval} = $data[-1];
    }


}

open OUT, ">regress_std.dat";
foreach $snp (@list){
    print OUT "$rcd{$snp}->{header}  $rcd{$snp}->{pval}\n";
}

$out = `Rscript enrich_std.R`;
chomp $out;
my @data = split /\s+/, $out;
shift @data until $data[0]=~/^\S/;
shift @data;
printf "%9.3f  %9.3f   %9.3f  %9.3f\n",@data[0..3];

