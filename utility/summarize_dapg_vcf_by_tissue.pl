$tissue = "";
$dir = "./";
for $i (0..$#ARGV){

    if($ARGV[$i] eq "-t"){
        $tissue = $ARGV[++$i];
        next;
    }

    if($ARGV[$i] eq "-dir"){

        $dir = $ARGV[++$i];
        next;
    }
}


@files = <$dir/*.rst>;
$count=0;
foreach $f (@files){
    process_fm($f);
}

print STDERR "\nProcessing completes, generating vcf ... \n";  
@chr = (1..22,"X");
foreach $ch (@chr){
    foreach $p (sort {$a<=>$b} keys %{$map{$ch}}){
        $id = $map{$ch}->{$p};
        print "$snp{$id}->{header}\t$snp{$id}->{info}\n";
    }
}


sub process_fm{
    my ($f) = @_;
    $f =~/$dir\/(\S+?)\./;
    my $gene = $1;
    my %cluster;
    open FILE, "grep \\\{ $f | ";
    while(<FILE>){
        s/\{//;
        s/\}//;
        my @data = split /\s+/, $_;
        shift @data until $data[0]=~/^\S/;
        $cluster{$data[0]} = "\[$data[2]:$data[1]\]";
    }

    open FILE, "grep \\\(\\\( $f | ";
    while(<FILE>){
        my @data = split /\s+/, $_;
        shift @data until $data[0]=~/^\S/;
        next if $data[-1] == -1;
        next if $data[2] < 1e-4;
        my $info = "$gene:$data[-1]\@$tissue\=$data[2]".$cluster{$data[-1]};
        if(!defined($snp{$data[1]})){
            $data[1] =~ /chr(\S+)\_(\d+)\_(\S+)\_(\S+)\_b38/;
            $map{$1}->{$2} = $data[1];
            $snp{$data[1]}->{header} = "chr$1\t$2\t$data[1]\t$3\t$4";
            $snp{$data[1]}->{info} = "$info";
        }else{
            $snp{$data[1]}->{info} .="|".$info;
        }


    }

}


