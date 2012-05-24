#!/usr/bin/perl
#
#  Script to calculate a fault model on a grid of points.  
#
#  The extents of the grid are defined by the bounding box of points 
#  in WKT shapes in a text WKT file, assumed to be latitude/longtiude.  
#
#  The resolution is defined in metres and is used to determine the grid
#  spacing.

use strict;
use POSIX;

my $prog = $0;
$prog =~ s/[^\\\/]*$/calc_okada.exe/;
die "Cannot find Okada program\n" if ! -x $prog;

my @progprm = ();
if( $ARGV[0] eq '-p')
{
    push(@progprm, shift(@ARGV ));
    push(@progprm, shift(@ARGV ));
    push(@progprm, shift(@ARGV ));
}

@ARGV==4 || die "Syntax: model_file wkt_boundary_file resolution output_grid_file\n";

my( $modelfile, $bdyfile, $resolution, $gridfile ) = @ARGV;

unlink($gridfile);

if( $resolution =~ /^res\:(.*)$/ )
{
    my $resf = $1;
    open(RF,"<$resf") || die "Cannot open resolution file $resf\n";
    $resolution = 0;
    while(my $inrec = <RF>)
    {
        my @f = split(' ',$inrec);
        my $res = $f[2]+0;
        next if $res < 10;
        $resolution = $res if $resolution==0 || $resolution > $res;
    }
    close(RF);
    print "Read grid resolution of $resolution from $resf\n";
}
$resolution += 0;
die "Invalid grid resolution $resolution\n" if $resolution < 10;

open(F,"<$modelfile") || die "Invalid model file $modelfile\n";
close(F);
open(W,"<$bdyfile") || die "Invalid wkt boundary file $bdyfile\n";
my $wkt = join('',<W>);
close(W);

printf "Reading lat/lon range from $bdyfile\n";

my $npt = 0;
my ($lnmin,$lnmax,$ltmin,$ltmax);
while( $wkt =~ /(\([\d\-\.\,\s]+\))/g)
{
    my $pts = $1;
    while( $pts =~ /[(\,]\s*([\d\.\-]+)\s+([\d\.\-]+)/g )
    {
        my($lon,$lat) = ($1,$2);
        # print "Point: $lon $lat\n";
        if( $npt )
        {
            if( $lon < $lnmin ) { $lnmin=$lon;}  elsif( $lon > $lnmax ) { $lnmax = $lon;}
            if( $lat < $ltmin ) { $ltmin=$lat;} elsif( $lat > $ltmax ) {$ltmax = $lat;}
        }
        else
        {
            $lnmin = $lnmax = $lon;
            $ltmin = $ltmax = $lat;
        }
        $npt++;
    }
}
printf("File contain $npt coordinate values - using bounding box\n");
printf("Range: $lnmin $lnmax $ltmin $ltmax\n");

my $dtor = atan2(1,1)/45.0;

my $maxlat = abs($ltmax) > abs($ltmin) ? abs($ltmax) : abs($ltmax);
my $ltres = ($resolution/6400000.0)/$dtor;
my $lnres = $ltres/cos($maxlat*$dtor);

# Round to 2 sig fig
$ltres = sprintf("%.1e",$ltres)+0;
$lnres = sprintf("%.1e",$lnres)+0;

printf "Resolution calculated: Longitude $lnres  Latitude $ltres\n";

$lnmin = floor($lnmin/$lnres)*$lnres;
$lnmax = ceil($lnmax/$lnres)*$lnres;

$ltmin = floor($ltmin/$ltres)*$ltres;
$ltmax = ceil($ltmax/$ltres)*$ltres;

my $nlon = int(($lnmax-$lnmin)/$lnres + 0.5);
my $nlat = int(($ltmax-$ltmin)/$ltres + 0.5);

my $tempfile = $gridfile.".range";
open(F,">$tempfile");
print F "grid $lnmin $ltmin $lnmax $ltmax $nlon $nlat\n";
close(F);

print "Calculating grid: $lnmin $ltmin $lnmax $ltmax $nlon $nlat\n";

system($prog,@progprm, $modelfile,$tempfile,$gridfile);

print "Grid generated in $gridfile\n" if -r $gridfile;
print "Grid not built!\n" if ! -r $gridfile;

unlink($tempfile);

exit(1) if ! -r $gridfile;
exit(0)

