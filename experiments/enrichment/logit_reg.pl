@files = <truth/*.annot>;
open OUT, ">logit.dat";

foreach $f (@files){

    $f =~ /(ENSG\d+)/;
    $g = $1;
    my %rcd;
    open FILE, "truth/$g.truth";
    while(<FILE>){
	next if $_ !~ /^\s*(\S+)\s+/;
	$rcd{$1} = 1;
    }

    open FILE, "$f";
    while(<FILE>){
	chomp;
	next if $_ !~ /^\s*(\S+)\s+/;
	my $y = 0;
	$y = 1 if(defined($rcd{$1}));
	print OUT "$_ $y\n";
    }
}


$out = `Rscript logitReg.R`;
chomp $out;
my @data = split /\s+/, $out;
shift @data until $data[0]=~/^\S/;
shift @data;
printf "%9.3f  %9.3f   %.3f  %.3f\n",$data[0],$data[1],$data[2],$data[3];
