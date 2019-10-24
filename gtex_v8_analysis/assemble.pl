# perl assemble.pl master_genotype_vcf query_bed expression_file covariate_file

$geno_vcf = $ARGV[0];
$query_bed = $ARGV[1];


open FILE, "$ARGV[3]";
while(<FILE>){
   next if $_ !~ /\d/;
     
   my @data = split /\s+/, $_;
   shift @data until $data[0]=~/^\S/;

   if (/\s*ID/){
      shift @data;
      @id_list = @data;
      next;
   }

   $var = $data[0];
   shift @data;
   $covar_text .="controlled $var gtex @data\n";
}








$geno_text = "";
open FILE, "bcftools view -R $query_bed $geno_vcf |";
while(<FILE>){

   if(/^\s*\#CHROM/){
      my @data = split /\s+/,$_;

      shift @data until $data[0]=~/^\s*GTEX/;
      @hash{@data} = (0..$#data);
      foreach $id (@id_list){
         push @index, $hash{$id};
      }
      next;
   }

   next if /^\s*\#/;
   next if $_ !~ /\d/;
   my @data = split /\s+/, $_;
   shift @data until $data[0]=~/^\S/;

   $snp = $data[2];
   shift @data for (1..9);
   @sdata = @data[@index];
   my @geno;
   foreach $d (@sdata){
     next if $d !~ /\S\/\S/;
    
     $g = 0 if $d =~ /0\/0/;
     $g = 1 if $d =~ /0\/1/;
     $g = 1 if $d =~ /1\/0/;;
     $g = 2 if $d =~ /1\/1/;;
     $g = "NA" if $d =~ /\D\/\D/;
     push @geno, $g;
   }

   $geno_text .= "geno $snp gtex @geno\n";


}

$expr_text = `cat $ARGV[2]`;


print "$expr_text";
print "$covar_text";
print "$geno_text";
 
