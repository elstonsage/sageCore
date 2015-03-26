#! /usr/bin/perl -w
# FIND_blox.pl <marker_file> <LD_file> <D'_cutoff>
# A script to generate a list of LD blocks defined by having
# |D'| exceed a given threshold for all pairs of ADJACENT markers
# within the block
#
# Input files:
# ------------
# marker_file:  Text file with two columns:  marker name and physical
#     map position
# LD_file:  Output from UNIX HaploView, *.LD file
#
# Output file:
# ------------
# LDblox.out:  a list of the LD blocks in physical map order
#
# Rob Igo 3/22/07

# Arguments:  Target directory, LOD score and D' thresholds for high LD
($markfile, $ldfile, $dprimet, $lodt) = @ARGV;
unless (@ARGV > 1) { 
  die "Syntax:  FIND_blox.pl <marker_file> <LD_file> <D_prime_cutoff>\n"; }
unless(defined($dprimet)) { $dprimet = 0.7; }
unless(defined($lodt)) { $lodt = -1; }

###############
# First, load in marker set
open (MARKER, "$markfile") || die "Couldn't open $markfile.\n";
while (<MARKER>) {
  next unless (/\w/);
  s/^\s+//;
  @list = split;
  push (@marklist, $list[0]);
}
close MARKER;

##############
# Next, load |D'| information for adjacent marker pairs
foreach $i (0..($#marklist - 1)) {
  $pair = $marklist[$i] . "_" . $marklist[$i + 1];
  push (@pairlist,$pair);
  $pairhash{$pair}++;
}

open (LD, "$ldfile") || die "Couldn't open $ldfile.\n";
<LD>;  # Skip header line
while (<LD>) {
  next unless (/\w/);
  s/^\s+//;
  @list = split;
  ($m1, $m2, $dprime, $lod) = @list[0..3];
  $pair = $m1 . "_" . $m2;
  next unless (exists($pairhash{$pair}));
  $dprimes{$pair} = $dprime;
  $lods{$pair} = $lod;
}


#############
# Finally, evaluate criteria for each pair of adjacent markers
# and form blocks
foreach $pair (@pairlist) {
  $ldnext = 0;
  if (exists($dprimes{$pair})) {
    if ($dprimes{$pair} >= $dprimet && $lods{$pair} >= $lodt) { $ldnext = 1; }
  }
  else {
    print "Warning:  No |D'| info for $pair.  Will default to no LD.\n"; }
  push(@ldnextlist, $ldnext);
  # print "$dprimes{$pair} $lods{$pair} $ldnext\n";
}

open (OUT, ">LDblox.out") || die "Couldn't write to LDblox.out.\n";
$block = 0;  # Current in-block flag
$blockcount = 1;  # Block designation
@blocksnps = ();  # SNPs in current block
foreach $i (0..$#pairlist) {
  ($m1, $m2) = @marklist[$i, $i + 1];
  if ($ldnextlist[$i]) {
    if ($block) { push(@blocksnps, $m2); }
    else {
      $block = 1;
      @blocksnps = ($m1, $m2);
    }
  }
  else {
    print OUT "Block $blockcount:  ";
    if ($block) {
      print OUT "@blocksnps\n";
      $block = 0;
      @blocksnps = ();
    }
    else { print OUT "$m1\n"; }
    $blockcount++;
  }

  if ($i == $#pairlist) {  # Print out last block or last marker
    print OUT "Block $blockcount:  ";
    if ($block) { print OUT "@blocksnps\n"; }
    else { print OUT "$m2\n"; }
  }
}

close OUT;
