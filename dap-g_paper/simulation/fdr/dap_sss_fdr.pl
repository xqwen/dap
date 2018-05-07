sub read_file {

    my ($id) = @_;

    my %truth;
    open FILE, "sim_data/sim.$id.truth";
    while(<FILE>){
        next if $_ !~ /\d/;
        my @data = split /\s+/, $_;
        shift @data until $data[0]=~/^\S/;
        $truth{$data[1]} = 1;
    }

    my %rcd;
    open FILE, "grep \\\(\\\( dap_out/sim.$id.dap.out |";
    while(<FILE>){
        my @data = split /\s+/, $_;
        shift @data until $data[0]=~/^\S/;
        shift @data;
        if(defined($truth{$data[0]})){
            $rcd{$data[-1]} = 1;
        }
    }
    open FILE,"grep \"\{\" dap_out/sim.$id.dap.out |";
    while(<FILE>){
        my @data = split /\s+/, $_;
        shift @data until $data[0]=~/^\S/;

        $data[0]=~/(\d+)/;
        printf "$id   $1    %7.3e   %d\n",$data[2], $rcd{$1};
    }

}

for(1..1000){
    read_file($_);
}
