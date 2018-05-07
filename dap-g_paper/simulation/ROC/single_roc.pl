sub read_file {

    my ($id) = @_;
    my %pip_cmp;
    my %pip_truth;
    open FILE, "zcat single_BF/sim.$_.bf.gz |";
    while(<FILE>){
        my @data = split /\s+/, $_;
        shift @data until $data[0]=~/^\S/;
        $data[0]=~/rs(\d+)/;
        $region = int($1/11);
        $index-- if($1%11 ==0);
        #print "$data[0] $region $data[1]\n";
        if(!defined($pip_cmp{$region})){
            $pip_cmp{$region} = $data[2];
        }else{
            if($data[2]>$pip_cmp{$region}){

                $pip_cmp{$region} = $data[2];
            }
        }
    }
    open FILE, "sim_data/sim.$id.truth";
    while(<FILE>){
        next if $_ !~ /\d/;
        my @data = split /\s+/, $_;
        shift @data until $data[0]=~/^\S/;
        $data[1]=~/rs(\d+)/;
        $region = int($1/11);
        $index-- if($1%11 ==0);
        $pip_truth{$region} = 1;
    }

    foreach $r (sort {$a<=>$b} keys %pip_cmp){
        printf "sim%s:r%d\t%7.3e\t%d\n",$id,$r, $pip_cmp{$r}, $pip_truth{$r};
    }

}

for(1..1000){
    read_file($_);
}
