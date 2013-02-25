use strict;
use FindBin;
use lib $FindBin::Bin;
use GridFile;
use Data::Dumper;

@ARGV == 2 || die;
my $grid = new GridFile($ARGV[0]) || die "Cannot open grid file\n";
$grid->WriteToFile($ARGV[1]."_velocity.txt",format=>'CSV',header=>'col row longitude latitude dE dN',skipmissingvalues=>1,separator=>"\t",pointids=>1);
open(my $output, ">$ARGV[1]"."_strain.txt") || die "Cannot open strain file ".$ARGV[1]."_strain.txt\n";
print $output "longitude\tlatitude\tdilatation\tshear\tshape\n";
my ($nx,$ny) = $grid->GridSize();


my $dtor = atan2(1,1)/45.0;
my $roe = 6371000;
my $ppb = 1.0e9;


my $r1 = $grid->RowData(0);
for my $r (1..$ny-1)
{
    my $r0 = $r1;
    $r1 = $grid->RowData($r);
    for my $c (1..$#$r0)
    {
        my $p00 = $r0->[$c-1];
        my $p10 = $r0->[$c];
        my $p01 = $r1->[$c-1];
        my $p11 = $r1->[$c];
        # Skip points with no real data
        next if $p00->[2] == 0 && $p00->[3] == 0;
        next if $p01->[2] == 0 && $p01->[3] == 0;
        next if $p10->[2] == 0 && $p10->[3] == 0;
        next if $p11->[2] == 0 && $p11->[3] == 0;
        
        # print "p00 ",join(" ",@$p00),"\n";
        # print "p01 ",join(" ",@$p01),"\n";
        # print "p10 ",join(" ",@$p10),"\n";
        # print "p11 ",join(" ",@$p11),"\n";
        # Midpoint of cell

        my ($x0, $y0, $x1, $y1) = ($p00->[0], $p00->[1], $p11->[0], $p11->[1] );

        my $x = ($x0 + $x1)/2;
        my $y = ($y0 + $y1)/2;

        my $shape = "POLYGON(($x0 $y0, $x0 $y1, $x1 $y1, $x1 $y0, $x0 $y0))";

        # Approx distance across grid cells in x an y directions

        my $coslat = cos($dtor*$y);
        my $dx = ($p10->[0] - $p00->[0])*$dtor*$roe*$coslat;
        my $dy = ($p01->[1] - $p00->[1])*$dtor*$roe;

        # Differentials of deformation with respect to x and y

        my $uxx = $ppb*(($p10->[2]+$p11->[2])-($p00->[2]+$p01->[2]))/(2*$dx);
        my $uyx = $ppb*(($p10->[3]+$p11->[3])-($p00->[3]+$p01->[3]))/(2*$dx);

        my $uxy = $ppb*(($p01->[2]+$p11->[2])-($p00->[2]+$p10->[2]))/(2*$dy);
        my $uyy = $ppb*(($p01->[3]+$p11->[3])-($p00->[3]+$p10->[3]))/(2*$dy);

        # print "$dx $dy $dtor $coslat $uxx $uxy $uyx $uyy\n";

        # Invariants of strain...

        my $dil = ($uxx+$uyy)/2;
        my $shear = sqrt((($uxx-$uyy)*($uxx-$uyy))/4 + abs($uxy*$uyx));

        print $output join("\t",$x,$y,$dil,$shear,$shape),"\n";
    }
}

