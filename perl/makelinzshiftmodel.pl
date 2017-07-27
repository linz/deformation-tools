#!/usr/bin/perl
#===============================================================================
#
# PROGRAM:             %M%
#
# VERSION:             %I%
#
# WHAT STRING:         %W%
#
# DESCRIPTION:         Script converts an ASCII definition of a shift
#                      model to a binary representation of the data
#
# PARAMETERS:          def_file   The name of the file defining the shift model
#                      bin_file   The name of the binary file generated.
#                      
# DEPENDENCIES:        
#
# MODIFICATION HISTORY
# NAME                 DATE        DESCRIPTION
# ===================  ==========  ============================================
# Chris Crook           3/02/2004  Created
#===============================================================================
#
#

use strict;

use FindBin;
use lib $FindBin::RealBin.'/perllib';

use Getopt::Std;
use Packer;

my $prgdir = $FindBin::Bin;

my $syntax = <<EOD;

Syntax: [options] source_file binary_file

Options are:
   -f format (LINZSHIFT1B, LINZSHIFT1L, LINZSHIFT2B, LINZSHIFT2L)

EOD

my %opts;
getopts('f:',\%opts);
my $forceformat = $opts{f};

my $endianness = {
   LINZSHIFT1B =>
     { sig => "LINZ shift model v1.0B\r\n\x1A",
       gridparam => '-f GRID2B',
       trigparam => '-f TRIG2B',
       bigendian => 1
       },
   LINZSHIFT1L =>
     { sig => "LINZ shift model v1.0L\r\n\x1A",
       gridparam => '-f GRID2L',
       trigparam => '-f TRIG2L',
       bigendian => 0
       },
   LINZSHIFT2B =>
     { sig => "LINZ shift model v2.0B\r\n\x1A",
       gridparam => '-f GRID2B',
       trigparam => '-f TRIG2B',
       bigendian => 1
       },
   LINZSHIFT2L =>
     { sig => "LINZ shift model v2.0L\r\n\x1A",
       gridparam => '-f GRID2L',
       trigparam => '-f TRIG2L',
       bigendian => 0
       },
   };

my @temp_files;
my $ntmpfile = 1000;


my %model_param = qw/
    FORMAT         (LINZSHIFT1B|LINZSHIFT1L|LINZSHIFT2B|LINZSHIFT2L)
    VERSION_NUMBER \d+\.\d+
    COORDSYS       \w+
    DESCRIPTION    .*
    /;

my %comp_param = qw/
    MODEL_TYPE         (trig|grid)
    NEGATIVE           (yes|no)
    DESCRIPTION        .*
    /;

my %comp_param2 = qw/
    MODEL_TYPE         (trig|grid)
    DESCRIPTION        .*
    FACTOR           \-?\d+(?:\.\d+)?
    GROUP_ID       \d+
    /;

my $def_model;

eval {
   my $model = &LoadDefinition( $ARGV[0] );
   &CheckSyntax($model);
   my $format = uc($forceformat || $model->{FORMAT});
   $model->{FORMAT} = $format;
   my $endian = $endianness->{$format};
   die "Invalid model format $format specified\n" if ! $endian;
   &BuildAllComponents($model,$endian);
   
   &WriteModelBinary($ARGV[1],$model,$endian);
};
if( $@ ) {
   print STDERR $@;
   }
unlink @temp_files;



sub LoadDefinition {
   
   my ($deffile) = @_;
   
   open(DEF,"<$deffile") 
      || die "Cannot open shift model definition file $deffile\n";
   
   my $curmod;
   my $curseq;
   my $curcomp;
   my $curobj;
   my $objname;
   
   while(<DEF>) {
      next if /^\s*($|\#)/;   # Skip blank lines, comments....
      chomp;
      my ($rectype, $value ) = split(' ',$_,2);
      $rectype = uc($rectype);
      
      if( $rectype =~ /^(SHIFT_MODEL_(MODEL|COMPONENT))$/ ) {
         if( $rectype eq 'SHIFT_MODEL_MODEL' ) {
            die "Only one SHIFT_MODEL_MODEL can be defined\n" if $curmod;
            $curmod = { type=>$rectype, name=>$value, params=>\%model_param,
                        sequences=>[], filename=>$deffile };
            $curobj = $curmod;
            }
         else {
            die "Must define SHIFT_MODEL_MODEL before $rectype\n" if ! $curmod;
            my $cparams=$curmod->{FORMAT} =~ /2/ ? \%comp_param2 : \%comp_param;
            $curcomp = { type=>$rectype, source=>$value, params=>$cparams, 
                        components=>[], sequence=>$curseq,
                        GROUP_ID=>0, FACTOR=>1.0 };
            push(@{$curmod->{components}}, $curcomp );
            $curobj = $curcomp;
            }
         }
   
      elsif( ! $curobj ) {
         die "SHIFT_MODEL_MODEL must be defined first in file\n";
         }
    
      elsif( ! $curobj->{params}->{$rectype} ) {
         die "Parameter $rectype is not valid in ".$curobj->{type}."\n";
         }
   
      else {
         if( $rectype eq 'DESCRIPTION' ) {
             $value = '';
             while(<DEF>) {
               last if /^\s*end_description\s*$/i;
               $value .= $_;
               }
             }
          $curobj->{$rectype} = $value;
          }
      }
   close(DEF);
   return $curmod;
   }

sub ModelFile
{
    my($model,$filename)=@_;
    return $filename if -e $filename;
    my $modeldir=$model->{filename};
    return $filename if $modeldir !~ /[\\\/]/;
    $modeldir =~ s/[^\\\/]*$//;
    return $modeldir.$filename;
}
   
sub TempFileName {
   my $tmpfilename;
   do {
     $tmpfilename = sprintf("deftmp%d",$ntmpfile++)
     }
   while -r $tmpfilename;

   push(@temp_files,$tmpfilename);
   return $tmpfilename;
   }


sub BuildGridComponent {
   my( $comp, $filename, $params ) = @_;
   my $origname = $filename;
   if( -T $filename ) {
     my $tmpfile = &TempFileName();
     my $command = "perl $prgdir/makegrid.pl $params $filename $tmpfile";
     system($command);
     die "Cannot create grid file from $filename\n$command\n" if ! -r $tmpfile;
     $filename = $tmpfile;
     }

   use GridFile;
   my $grid;
   eval {
      $grid = new GridFile $filename;
      };
   if( $@ ) {
      die "Cannot create or use grid file $origname\n";
      };
   if( ! $grid ) {
      die "Invalid grid in $filename\n" if ! $grid;
      };

   my $coordsys = $grid->CrdSysCode();
   my $ndim = $grid->Dimension();
   my ($xmin,$ymin,$xmax,$ymax) = $grid->Range();
   my $filesize = -s $filename;

   $comp->{sourcefile} = 
     { name=>$filename,
       origname=>$origname,
       size=>$filesize,
       coordsys=>$coordsys,
       ndim=>$ndim,
       range=>[$ymin,$ymax,$xmin,$xmax]
       };
   }

sub BuildTrigComponent {
   my( $comp, $filename, $params ) = @_;
   my $origname = $filename;
   if( -T $filename ) {
     my $tmpfile = &TempFileName();
     my $command = "perl $prgdir/maketrig.pl $params $filename $tmpfile";
     system($command);
     die "Cannot create grid file from $filename\n$command\n" if ! -r $tmpfile;
     $filename = $tmpfile;
     }

   use TrigFile;
   my $trig;
   eval {
      $trig = new TrigFile $filename;
      };
   if( $@ ) {
      die "Cannot create or use trig file $origname\n";
      };
   if( ! $trig ) {
      die "Invalid trig in $filename\n" if ! $trig;
      };

   my $coordsys = $trig->CrdSysCode();
   my $ndim = $trig->Dimension();
   my ($xmin,$ymin,$xmax,$ymax) = $trig->Range();
   my $filesize = -s $filename;

   $comp->{sourcefile} = 
     { name=>$filename,
       origname=>$origname,
       size=>$filesize,
       coordsys=>$coordsys,
       ndim=>$ndim,
       range=>[$ymin,$ymax,$xmin,$xmax]
       };
   }

sub BuildComponent {
   my ($model,$comp,$endian) = @_;
   my $source = $comp->{source};
   my ($filename, $params) = split(' ',$source,2);
   $filename=ModelFile($model,$filename);
   die "Cannot find component source file $filename\n" if ! -r $filename;

   my $type = $comp->{MODEL_TYPE};
   if( $type eq 'grid' ) {
      &BuildGridComponent($comp,$filename,$params.' '.$endian->{gridparam});
      }
   elsif( $type eq 'trig' ) {
      &BuildTrigComponent($comp,$filename,$params.' '.$endian->{trigparam});
      }
   else {
      die "Invalid component MODEL_TYPE $type\n";
      }
   }

sub CheckObject {
   my ($obj) = @_;
   my $nerror = 0;
   my $params = $obj->{params};
   foreach my $k (sort keys %$params) {
      my $pattern = $params->{$k};
      my $value = $obj->{$k};
      my $result;
      $result = $value =~ /^$pattern$/s;
      if( ! $result ) {
        print $value eq '' ? "Missing value" : "Invalid value $value",
              " for $k in ",$obj->{type},"\n";
        $nerror++;
        }
      }
   return $nerror;
   }

sub ExpandRange {
   my ($oldrange,$newrange) = @_;
   if( ! $oldrange ) {
       $oldrange = [@$newrange];
       }
   else {
       if( $oldrange->[0] > $newrange->[0] ) { $oldrange->[0] = $newrange->[0];}
       if( $oldrange->[1] < $newrange->[1] ) { $oldrange->[1] = $newrange->[1];}
       if( $oldrange->[2] > $newrange->[2] ) { $oldrange->[2] = $newrange->[2];}
       if( $oldrange->[3] < $newrange->[3] ) { $oldrange->[3] = $newrange->[3];}
       }
   return $oldrange;
   }
 

sub BuildAllComponents {
   my ($model,$endian) = @_;
   
      foreach my $c (@{$model->{components}}) {
         &BuildComponent($model,$c,$endian);
         $model->{range} = &ExpandRange( $model->{range}, $c->{sourcefile}->{range} );
         }
   }

sub CheckSyntax {
   my ($model) = @_;
   
   my $nerror = &CheckObject( $model );
 
   my $ncomponents = scalar(@{$model->{components}});
   if( $ncomponents < 1 ) {
      print "Shift model has no components\n";
      $nerror++;
      }
   $model->{ncomponents} = $ncomponents;
   foreach my $c (@{$model->{components}}) {
      $nerror += &CheckObject($c);
      }
 
   die "Failed with invalid syntax\n" if $nerror;
   }

   
sub WriteModelBinary {
   my($binfile,$model,$endian) = @_;
   my $pack = new Packer( $endian->{bigendian} );

   open(OUT,">$binfile") || die "Cannot create output file $binfile\n";
   binmode(OUT);

   # Print out file signature...

   print OUT $endian->{sig};
   my $indexloc = tell(OUT);
   print OUT $pack->long(0);

   # Print out component data ...

   my $version=1;
   $version = $1 if $model->{FORMAT} =~ /(\d+)[BL]$/;

   foreach my $c (@{$model->{components}}) {
       my $sf = $c->{sourcefile};
       open(SF,"<$sf->{name}") || 
          die "Cannot open component souce file $sf->{name}\n";
       binmode(SF);
       $sf->{loc} = tell(OUT);
       my $buf;
       sysread(SF,$buf,$sf->{size});
       print OUT $buf;
       close(SF);
     }

   # Print out model, components ... 

   my $loc = tell(OUT);

   print OUT $pack->string( @{$model}{qw/name VERSION COORDSYS/});
   print OUT $pack->string4( @{$model}{qw/DESCRIPTION/});
   print OUT $pack->double( @{$model->{range}} );
   print OUT $pack->short( $model->{ncomponents} );
                          
   foreach my $c (@{$model->{components}}) {
       print OUT $pack->string( $c->{DESCRIPTION} );
       print OUT $pack->double( @{$c->{sourcefile}->{range}} );
       print OUT $pack->short( lc($c->{MODEL_TYPE}) eq 'trig' ? 1 : 0 );
       print OUT $pack->short( $c->{sourcefile}->{ndim} );
       print OUT $pack->short( lc($c->{NEGATIVE}) eq 'no' ? 0 : 1 ) if $version <= 1;
       print OUT $pack->short( $c->{GROUP_ID} || 0 ) if $version > 1;
       print OUT $pack->double( $c->{FACTOR} || 0.0 ) if $version > 1;
       print OUT $pack->long( $c->{sourcefile}->{loc} );
     }

  seek(OUT,$indexloc,0);
  print  OUT $pack->long($loc);

  close(OUT);
  }

__END__

# Example of an input file ...
#
# The model consists of one or more components.
# Each component in the sequence can be either a triangulated model, or 
# a grid model defining horizontal, vertical, or 3d deformation.  
# The source file referenced can be a pre-built binary file,
# or more usefully an ascii source file from which the binary can be built.
# The latter is preferable, as it will ensure that all components are built
# with the same endian-ness.
#
# The format can be one LINZSHIFT1B,LINZSHIFT1L,LINZSHIFT2B,LINZSHIFT2L
# 
# If it is a version 2 format then the component can include a GROUP_ID
# and FACTOR field in place of NEGATIVE.
# Where there are consecutive components with the same id only the first
# matching component will be used

# An example follows:

SHIFT_MODEL_MODEL Darfield earthquake 2000.0 shift
FORMAT LINZSHIFT1B
VERSION_NUMBER 1.0
COORDSYS NZGD2000
DESCRIPTION
This is the description of the model
This is a first try
END_DESCRIPTION

SHIFT_MODEL_COMPONENT darfield_far_field.grd
MODEL_TYPE grid
NEGATIVE yes|no
DESCRIPTION
Far field deformation
END_DESCRIPTION

SHIFT_MODEL_COMPONENT darfield_near_field.trg
MODEL_TYPE trig
NEGATIVE yes|no
DESCRIPTION
Near field deformation
END_DESCRIPTION
