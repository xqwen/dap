$dir = "./";

for $i (0..$#ARGV){

    if($ARGV[$i] eq "-dir"){
        $dir = $ARGV[++$i];
        next;
    }
}


@files = <$dir/*.vcf>;
foreach $f (@files){
	
	open FILE, "$f";
	while(<FILE>){
		 next if /^\s*\#/;
		 next if $_ !~ /\d/;
        	/^(.*)\t(\S+)$/;
		chomp;
		my $header = $1;
		my $info = $2;
		my @data = split /\s+/, $_;
		$data[0]=~/chr(\S+)/;
		my $chr = $1;
		my $pos = $data[1];
		if(!defined($record{$chr}->{$pos})){
			$record{$chr}->{$pos}->{header} = $header;
			$record{$chr}->{$pos}->{info} = $info;
		}else{
			$record{$chr}->{$pos}->{info} .= "\|".$info;
		}
	}
	print STDERR "Finish processing $f ... \n"; 
}


foreach $chr (1..22){
	foreach $pos (sort {$a <=>$b} keys %{$record{$chr}}){
		my $cinfo = parse_info($record{$chr}->{$pos}->{info});
		print "$record{$chr}->{$pos}->{header}\t$cinfo\n";
	}
}
			



sub parse_info {

	my ($line) = @_;
	my @ldata = split /\|/, $line;
	@index = (0..$#ldata);
	my $con ="";
	my %rcd;
	@rcd{@index}= @ldata;
	my @pip;
	foreach $d (@ldata){
		$d =~/\=(\S+)\[/;
		push @pip, $1;
	}
	my %hash;
	@hash{@index} = @pip;
	$count = 0;
	foreach $k (sort {$hash{$b} <=> $hash{$a}} keys %hash){
		$con .= "$rcd{$k}";
		if($count < $#ldata){
			$con .= "|";
			$count++;
		}
	}
	return($con);
}
