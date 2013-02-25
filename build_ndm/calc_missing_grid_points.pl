# Program to calculate velocities for points at which the grid model hasn't
# been able to calculate real values.
#
# Arguments are the input file, the list of grid points at which new values
# are to be calculated, and the output file.  The grid points to be calculated
# are identified by a column and row number (0 based)
#
# The grid points are calculated by interpolating from adjacent points using
# a simple best fit translation, rotation and scale.  This is only intended to
# extend the model by one or two grid cells at the edge of the valid data 
# where the missing grid data affects points on land.


use strict;
use FindBin;
use lib $FindBin::Bin;
use LeastSquares;

# Maximum range of interpolation - maximum number of cells away from calc point,
#   in latitude or longitude direction
# Minimum number of points required to calculate a new value

my $maxrange = 6;
my $minusept = 4;  


@ARGV==3 || die "Parameters required: input_grid_file missing_points_file output_grid_file\n";


my($dmpin,$pointfile,$dmpout) = @ARGV;

open(DMPIN,"<$dmpin") ||  die "Cannot open $dmpin\n";
open(PNTS,"<$pointfile") || die "Cannot open $pointfile\n";
open(DMPOUT,">$dmpout") ||  die "Cannot open output file $dmpout\n";

my %param;
my $grid=[];
my $fixed = {};

while(<DMPIN>) {
   if ( /^V(\d+)\,(\d+)\:\s+(\S+)\s+(\S+)/ ) {
     $grid->[$1-1]->[$2-1] = [$3,$4] if ($3 ne '*' &&  $4 ne '*' && ! ( $3 eq '0' && $4 eq '0'));
     }
   elsif ( /^(\w+)\:\s+(\S+)/ ) {
     print DMPOUT $_;
     $param{$1} = $2;
     }
   }

my ($lnmin, $lnmax, $nln, $ltmin, $ltmax, $nlt ) 
     = @param{'XMIN','XMAX','NGRDX','YMIN','YMAX','NGRDY'};

my $dln = ($lnmax-$lnmin)/($nln-1);
my $dlt = ($ltmax-$ltmin)/($nlt-1);

# Calculate each point for which information is required ...

while( <PNTS> ) {
   next if ! /^\s*(\d+)\s+(\d+)\s*$/;
   my($r,$c) = ($1,$2);
   if ( ! defined $grid->[$r-1]->[$c-1] )
   {
        print "Calculating missing grid point $r $c\n";
        &CalcGridPoint($grid, $fixed, $r-1,$c-1) 
   }
   else
   {
       print "Point $r $c already defined - not calculated\n";
   }
   }

for( my $c=0; $c<$nlt; $c++ ) {
  for( my $r=0; $r<$nln; $r++ ) {
     my($de,$dn) = ('*','*');
     ($de,$dn) = @{$fixed->{$r}->{$c}} if defined $fixed->{$r}->{$c};
     ($de,$dn) = @{$grid->[$r]->[$c]} if defined $grid->[$r]->[$c];
     printf DMPOUT "V%d,%d: %s %s\n",$r+1,$c+1,$de,$dn;
     }
   }

sub CalcGridPoint {
  my( $grid, $fixed, $r, $c ) = @_;
  my $ls = new LeastSquares(4);
  my $lt = $ltmin+$c*$dlt;
  my $clt = cos($lt*0.01745329);
  my $n = 0;
  for my $range (1..$maxrange)
  {
  for( my $r1 = $r-$range; $r1 <= $r+$range; $r1++ ) {
    next if $r1 < 0;
    for( my $c1 = $c-$range; $c1 <= $c+$range; $c1++ ) {
       next if $c1 < 0;
       next if abs($c1 -$c) < $range && abs($r1-$r) < $range;
       next if ! defined $grid->[$r1]->[$c1];
       $n++;
       my $x = ($r1-$r)*$dln*$clt; 
       my $y = ($c1-$c)*$dlt;
       my($dx,$dy) = @{$grid->[$r1]->[$c1]};
       $ls->add( $dx, [1,0,$x,-$y] );
       $ls->add( $dy, [0,1,$y,$x] );
       }
    }
    last if $n >= $minusept;
    }
    eval {
     die "Insufficient valid data points within range\n" if $n < $minusept;
     my $soln = $ls->solve;
     my $v1 = sprintf("%.6f",$soln->[0]);
     my $v2 = sprintf("%.6f",$soln->[1]);
     $fixed->{$r}->{$c} = [$v1,$v2];
     print "Values calculated successfully ($v1,$v2)\n";
     };
    if( $@ )
    {
        print "Unable to calculate points\n$@\n";
    }
  }
   
