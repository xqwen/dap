sub read_file {

    my ($id) = @_;
    my %rcd;




    open FILE, "sim_data/sim.$id.truth";
    while(<FILE>){
        next if $_ !~ /rs/;
        my @data = split /\s+/, $_;
        shift @data until $data[0]=~/^\S/;
        $rcd{$data[1]} = 1;
    }

    open FILE, "grep \\\(\\\( dap_z_out/sim.$id.dap_ss.out |";
    while(<FILE>){
        my @data = split /\s+/, $_;
        shift @data until $data[0]=~/^\S/;
        shift @data;
        printf "$id  $data[0]  $data[1]   %d\n",$rcd{$data[0]};
    }

}

for(1..1000){
    read_file($_);
}
