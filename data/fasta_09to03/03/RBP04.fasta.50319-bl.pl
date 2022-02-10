#!/usr/bin/env perl
$host = shift;
$instance = shift;
$arg = shift;

#### random sleep, rand() can be a fraction of second
select(undef,undef,undef,rand());

if ($arg) {
  @ids = split(/,/, $arg);
}
else {
  while(1) {
    if (opendir(DDIR, "RBP04.fasta-seq")) { 
      @ids = grep {/^\d+$/} readdir(DDIR);
      last;
    }
    else {
      sleep(1);
    }
  }
}

foreach $id (@ids) {

  next unless (-e "RBP04.fasta-seq/$id");
  next if (-e "RBP04.fasta-seq/$id.lock");
  $cmd = `touch RBP04.fasta-seq/$id.lock`;

  if (50) {
    $cmd = `blastp -outfmt 6 -db ./RBP04.fasta.50319 -seg yes -evalue 0.000001 -max_target_seqs 100000 -num_threads 1 -query RBP04.fasta-seq/$id -out RBP04.fasta-bl/$id`;
    $cmd =                         `/home/aoli/subprogram/cd-hit/cdhit/psi-cd-hit/psi-cd-hit.pl -J parse_blout_multi RBP04.fasta-bl/$id -c 0.3 -ce -1 -aS 0 -aL 0 -G 1 -prog blastp -bs 0 >> RBP04.fasta-blm/$host.$instance`;
  }
  elsif (1) {
    $cmd = `blastp -outfmt 6 -db ./RBP04.fasta.50319 -seg yes -evalue 0.000001 -max_target_seqs 100000 -num_threads 1 -query RBP04.fasta-seq/$id | /home/aoli/subprogram/cd-hit/cdhit/psi-cd-hit/psi-cd-hit.pl -J parse_blout RBP04.fasta-bl/$id -c 0.3 -ce -1 -aS 0 -aL 0 -G 1 -prog blastp -bs 1`;
  }
  else {
    $cmd = `blastp -outfmt 6 -db ./RBP04.fasta.50319 -seg yes -evalue 0.000001 -max_target_seqs 100000 -num_threads 1 -query RBP04.fasta-seq/$id -out RBP04.fasta-bl/$id`;
    $cmd =                         `/home/aoli/subprogram/cd-hit/cdhit/psi-cd-hit/psi-cd-hit.pl -J parse_blout RBP04.fasta-bl/$id -c 0.3 -ce -1 -aS 0 -aL 0 -G 1 -prog blastp -bs 0`;
  }
  $cmd = `rm -f  RBP04.fasta-seq/$id`;
  $cmd = `rm -f  RBP04.fasta-seq/$id.lock`;
}

($tu, $ts, $cu, $cs) = times();
$tt = $tu + $ts + $cu + $cs;
$cmd = `echo $tt >> RBP04.fasta-seq/host.$host.$instance.cpu`;

