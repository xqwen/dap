#! /usr/bin/perl

$file = "";
$prob = 0.95;

for $i  (0..$#ARGV){
    if ($ARGV[$i] eq "-d"){
        $file = $ARGV[++$i];
        next;
    }
    if($ARGV[$i] eq "-p"){
        $prob = $ARGV[++$i];
        next;
    }
}


if($file eq ""){
    print STDERR "Error: dap-g output file is not specified.\n";
    exit(1);
}

open FILE, "grep \\\{ $file | ";

while(<FILE>){
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;

    next if $data[2]< $prob;
    $data[0] =~/(\d+)/;
    $sig{$1} = {};
    $cum{$1} = 0;
}



if (scalar(keys %sig)==0){
    printf STDERR "No %d%% credible set can be constructed\n", int($prob*100);
    exit;
}


open FILE, "grep \\\(\\\( $file | ";
while(<FILE>){
    
    my @data = split /\s+/, $_;
    shift @data until $data[0]=~/^\S/;

    if(defined($sig{$data[-1]}) && $cum{$data[-1]}< $prob){
        $sig{$data[-1]}->{$data[1]} = $data[2];
        $cum{$data[-1]} += $data[2];
    }

}


foreach $c (sort {$a <=> $b} keys %sig){
    
    foreach $snp (sort {$sig{$c}->{$b} <=> $sig{$c}->{$a}} keys %{$sig{$c}}){
        printf "%2d\t%10s\t%7.3e\n",$c, $snp, $sig{$c}->{$snp}
    }
}




