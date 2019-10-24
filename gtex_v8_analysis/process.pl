for ($i=0;$i<=$#ARGV;$i++){
  if($ARGV[$i] eq "-e"){
    $expr_file = $ARGV[++$i];
    next;
   }
   if($ARGV[$i] eq "-g"){
     $geno_vcf = $ARGV[++$i];
     next;
   }
   
   if($ARGV[$i] eq "-c"){
     $covar_file = $ARGV[++$i];
     next;
   }
 
   if($ARGV[$i] eq "-t"){
     $tissue = $ARGV[++$i];
     next;
   }
   
   printf STDERR "Error: unknown option \"$ARGV[$i]\" \n";
   exit;
}  

printf "Expression vcf file: $expr_file\n";
printf "Genotype vcf file: $geno_vcf\n";
printf "Covariate file: $covar_file\n";
printf "Tissue: $tissue\n";

$tissue = "tissue" if !defined($tissue);
`mkdir $tissue` if ! -d $tissue;
open FILE, "zcat $expr_file | ";
open CMD, ">$tissue.assemble.cmd";
while(<FILE>){
  
   next if $_ !~ /ENSG/;
   next if $_ !~ /\d/;
   my @data = split /\s+/, $_;
   shift @data until $data[0]=~/^\S/;
   my $chr = $data[0];
   my $sp = $data[1];
   $pos1 = $sp-1000000;
   $pos1 = 1 if $pos1<1;
   $pos2 = $sp+1000001;
   $gene = $data[3];
   push @gene_list, $gene;
   open OUT , ">$tissue/$gene.cis_snp.bed";
   print OUT "$chr\t$pos1\t$pos2\n";
   @expression = @data[4..$#data];
   close OUT;
   open OUT, ">$tissue/$gene.expr";
   print OUT "pheno $gene gtex @expression\n";
   close OUT;
   print CMD "perl assemble.pl $geno_vcf $tissue/$gene.cis_snp.bed $tissue/$gene.expr $covar_file > $tissue/$gene.sbams.dat\n";
}
   
