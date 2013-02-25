use strict;
use FindBin;
use lib $FindBin::Bin;
use Getopt::Std;

# Script to build a grid binary file using the IGNS velocity model software;
#

my $exefile = $FindBin::Bin.'\\gns_velocity\\gns_velocity_linz.exe';
die "Need GNS velocity program $exefile\n" if ! -x $exefile;

my $syntax = <<EOD;

Syntax: perl build_velocity_grid.pl  grid_def_file  output_file

The grid file header requires the necessary fields of a grid dump file,
that is:
  XMIN:   Min longitude
  XMAX:   Max longitude
  NGRDX:  Number of longitude grid points (number of cells + 1)
  YMIN:   Min latitude
  YMAX:   Max latitude
  NGRDY:  Number of latitude grid points
  H1:     Three lines of descriptive text header
  H2:
  H3:
  CRDSYS: Code for the reference coordinate system

Optional parameters:
  SOLNFILE:  The gns solution file
  EULER:     The Euler pole lat, long and rotation rate for the
             output reference frame

EOD

my %opts;
getopts('k',\%opts);
my $keepfiles = $opts{k};

@ARGV == 2 || die $syntax;

my( $cfgfile, $dmpfile) = @ARGV;

open(CFG,"<$cfgfile") || die "Cannot open configuration file $cfgfile\n";
open(GRD,">$dmpfile") || die "Cannot open output file $dmpfile\n";

my $param;

while( <CFG> ) {
   next if ! ~ /^(\w+)\:\s+(.*?)\s*$/;
   $param->{$1} = $2;
   }
close CFG;

foreach (qw/XMIN XMAX NGRDX YMIN YMAX NGRDY HEADER0 CRDSYS VRES/) {
   die "Parameter $_ missing from $cfgfile\n" if ! exists $param->{$_};
   }

foreach (qw/XMIN XMAX YMIN YMAX VRES/ ) {
   die "Invalid $_ parameter\n" if $param->{$_} !~ /^\-?\d+\.?\d*(e[+-]?\d+)?$/;
   }

foreach (qw/NGRDX NGRDY/ ) {
   die "Invalid $_ parameter\n" if $param->{$_} !~ /^\d+$/;
   }

$param->{CRDSYS} = uc( $param->{CRDSYS} );
die "Invalid CRDSYS parameter\n" if $param->{CRDSYS} !~ /^\w+$/;

die "Invalid HEADER0 parameter\n" if $param->{HEADER0} =~ /^\s*$/;

my ($format,$xmn, $xmx, $nx, $ymn, $ymx, $ny, $m1,$m2, $m3, $csystem, $vres,$solnfile, $euler) =
   @$param{ qw/FORMAT XMIN XMAX NGRDX YMIN YMAX NGRDY HEADER0 HEADER1 HEADER2 CRDSYS VRES SOLNFILE EULER/ };
$format ||= 'GRID1L';
$solnfile ||= 'solution.gns';
$euler ||= '0 0 0';

die "GNS solution file $solnfile does not exist\n" if ! -r $solnfile;
die "Invalid Euler parameters $euler\n" 
    if $euler !~ /^\s*\-?\d+(:?\.\d+)?\s+\-?\d+(?:\.\d+)?\s+\-?\d+(?:\.\d+)?\s*$/;
my @eulercoeff = split(' ',$euler);

# Print the grid file header and space for a pointer to the index..

print GRD <<"EOD";
FORMAT: $format
HEADER0: $m1
HEADER1: $m2
HEADER2: $m3
CRDSYS:  $csystem
NGRDX: $nx
NGRDY: $ny
XMIN: $xmn
XMAX: $xmx
YMIN: $ymn
YMAX: $ymx
VRES: $vres
NDIM: 2
LATLON: 1
VALUES: REAL
EOD

my $datfile = 'tmp_bvg_lat_long.dat';
my $outfile = 'tmp_bvg_velocity.out';

my $latinc = ($ymx-$ymn)/($ny-1);
my $loninc = ($xmx-$xmn)/($nx-1);

my $ngp = $nx*$ny;
my $ngpc = 0;
my $vmax = 0.0;

open(DAT,">$datfile") || die;
print DAT "$ngp\n";

for (my ($r,$lat) = (1, $ymn); $r <= $ny; $r++, $lat += $latinc ) 
{
    for( my ($c, $lon )  = (1, $xmn); $c <= $nx; $c++, $lon += $loninc ) 
    {
        $ngpc++;
        print DAT "$ngpc $lat $lon\n";
    }
}
close(DAT);

system( $exefile, $solnfile, $datfile, $outfile, @eulercoeff );

unlink($datfile) if ! $keepfiles;


open(DAT,"$outfile") || die "GNS velocity output file $outfile not created\n";
my $rec='';
while($rec = <DAT>)
{
    last if $rec =~ /^\s*\d+\s+/;
}

$ngpc = 0;
my $nbad = 0;
for (my ($r,$lat) = (1, $ymn); $r <= $ny; $r++, $lat += $latinc ) 
{
    for( my ($c, $lon )  = (1, $xmn); $c <= $nx; $c++, $lon += $loninc ) 
    {
        $ngpc++;
        my @recdata = split(' ',$rec);
        my $id = $recdata[0];
        die "Data sequence error $ngpc > $id\n" if $id < $ngpc && $rec;
        my ($de, $dn) = ('*','*');
        if( $id == $ngpc )
        {
            $rec = <DAT>;
            $de = $recdata[3]/1000.0;
            $dn = $recdata[4]/1000.0;
            my $dea = abs($de);
            my $dna = abs($dn);
            $vmax = $dea if $dea > $vmax;
            $vmax = $dna if $dna > $vmax;
        }
        else
        {
            $nbad++;
        }
	print GRD "V$c,$r: $de $dn\n";
   }
   }
close(DAT);
unlink($outfile) if ! $keepfiles;

print "$ngp grid points calculated - $nbad out of range\n";
print "\n";
print "Maximum grid value is $vmax\n";
$vmax /= 32000;
print "Recommended VRES value is $vmax\n";

close(GRD);

