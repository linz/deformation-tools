#===============================================================================
# Module:             CSV_File.pm
#
# Description:       Defines packages: 
#                      CSV_File
#
#                    This module provides a interface to a CSV format file 
#                    (actually not strictly comma separated ... its 
#                    definable.  It is not terribly quick!  It does support
#                    quoting, quotes that span lines, and either double 
#                    quoting or an escape character to represent quotes within
#                    quoted strings.  It works by building fancy (but complex/
#                    slow) regexps.  
#
#                    The regexps can be made more efficient by using the 
#                    (?> independent subexpression regexp, but this limits
#                    the perl versions that can be used.
#
#                    For many applications this would be much more efficiently
#                    replaced with the Text:CSV_XS module from CPAN
#
#                    This module is for reading files only.
#
#  $Id: CSV_File.pm,v 1.2 2001/02/09 03:37:27 cms Exp $
#
#  $Log: CSV_File.pm,v $
#  Revision 1.2  2001/02/09 03:37:27  cms
#  Added option to explicitly close the CSV_File.
#
#  Revision 1.1  2000/03/29 03:13:12  cms
#  Building lib module
#
#  Revision 1.1  1999/09/10 13:28:30  ccrook
#  Initial revision
#
#
#===============================================================================

use strict;


#===============================================================================
# Note: This was written using FileHandle, but did not compile with use 
# strict refs failed on references to main::STDOUT.  This may be fixed
# in more recent versions of the standard perl modules.
#
# use FileHandle;

#===============================================================================
#
#   Class:       CSV_File
#
#   Description: Defines the following routines:
#                  $csv = new CSV_File($file, $def)
#                  $csv->DESTROY()
#                  $csv->NextRecord()
#                  $csv->columns()
#                  $csv->column(@i)
#                  $csv->colno(@fields)
#                  $csv->fields()
#                  $csv->field(@fields)
#                  $csv->setfield($field,$value)
#                  $csv->fieldno(@i)
#                  $csv->IsBlank()
#                  $csv->OpenOutputFile( $filename, $options )
#                  $csv->WriteRecord()
#
#   Output options:
#     def=>$def
#     fields=>\@fields
#     noquote_fields=>\@noquotefields

#
#===============================================================================

package CSV_File ;
use vars qw/$csvdef/;

my $fileno = 0;

$csvdef = { delim=>',', quote=>'"', doublequote=>1, noheaders=>0, literal=>'' };



#===============================================================================
#
#   Method:       new
#
#   Description:  $object = new CSV_File($file, $def)
#
#   Parameters:   $file       The name of the file to opend
#                 $def        The definition of the file.. a hash with elements
#                               delim   The file delimiter
#                               quote   The quote character
#                               literal The quote escape character
#                               doublequote  If true then repeated quotes are
#                                            intepreted as a literal quote
#                               noheaders If true then the first line is
#                                         treated as data rather than column
#                                         headers.
#                               columns
#                               spanlines If true than the data within quotes 
#                                         can span more than one line.
#
#   Returns:     $object     The blessed hash reference. 
#
#   Globals:
#
#===============================================================================

sub new {
  my ( $class, $file, $def ) = @_;

  my $fh = "CSV_FILE_HANDLE_$fileno"; $fileno++;
  
  no strict "refs";
  open($fh,$file) || die "Cannot open file $file\n";
  use strict "refs";

  my ($delim, $quote, $literal, $doublequote, 
      $noheaders, $columns, $spanlines)
    = @{$def}{'delim','quote','literal','doublequote',
              'noheaders','columns','spanlines'};
  $delim = "," if ! $delim;
  my $self = { fh=>$fh, name=>$file, columns=>$columns, 
               delim=>$delim, quote=>$quote, 
               literal=>$literal, doublequote=>$doublequote,
               spanlines=>$spanlines,
               recno=>0 };
  $delim = quotemeta($delim);
  $literal = quotemeta($literal);
  $quote = quotemeta($quote);
  my $fieldpattern = "[^".$delim.$quote.$literal."]";
  $fieldpattern = $fieldpattern."|$literal." if $literal;
  $fieldpattern = "(?:$fieldpattern)*";
  my $unclosed_quote;
  if( $quote ) {
     my $fp = "[^$quote$literal]";
     $fp = $fp."|$literal." if $literal;
     $fp = $fp."|$quote$quote" if $doublequote;
     $fp = "(?:$fp)*";
     $fieldpattern = "(?:$quote($fp)$quote|($fieldpattern))";
     $unclosed_quote = "(?:$quote$fp)" if $spanlines;
     }
  else {
     $fieldpattern = "($fieldpattern)";
     }
  $fieldpattern = "$delim$fieldpattern";
  $self->{fieldpattern} = $fieldpattern;
  $self->{unclosed_quote} = "$delim$unclosed_quote" if $unclosed_quote;

  bless $self, $class;
  $noheaders = 1 if $columns;
  if( ! $noheaders ) {
    $self->NextRecord || die "Header record missing in $file\n";
    $columns = $self->{columns} = $self->{data};
    $self->{recno} = 0;
    }
  my @ucolumns = map {uc($_)} @$columns;
  $self->{ucolumns} = \@ucolumns;
  # Create a subref to allow lookup of data from column names
  my $getfieldfunc;
  my $setfieldfunc;
  my $getcolnofunc;
  if( $columns ) {
     my %lookup;
     foreach (0..$#$columns) { $lookup{uc($columns->[$_])} = $_; }
     $getfieldfunc = sub { 
                          my $fn = uc($_[0]); 
                          return exists $lookup{$fn} ?
                                    $self->{data}->[$lookup{$fn}] :
                                    undef;
                          };
     $setfieldfunc = sub { 
                          my $fn = uc($_[0]); 
                          die "Invalid field $_[0] requested\n"
                              if ! exists $lookup{$fn};
                          $self->{data}->[$lookup{$fn}]  = $_[1];
                      };
     $getcolnofunc = sub {
                          my $fn = uc($_[0]); 
                          return exists $lookup{$fn} ?
                                    $lookup{$fn} :
                                    undef;
                          };
     }
  else {
     $getfieldfunc = sub { undef };
     $setfieldfunc = sub { undef };
     $getcolnofunc = sub { undef };
     }
  $self->{getfieldfunc} = $getfieldfunc;
  $self->{setfieldfunc} = $setfieldfunc;
  $self->{getcolnofunc} = $getcolnofunc;
  $self->{outputfile} = undef;
  $self->{outputfh} = undef;
  $self->{outputdef} = undef;
  $self->{outputcols} = undef;
  $self->{outputfields} = undef;
  $self->{noquotefields} = undef;
  
  return $self;
  }


sub filename {
  return $_[0]->{name};
  }


sub recno {
  return $_[0]->{recno};
  }


sub lineno {
  return $_[0]->{lineno};
  }

#===============================================================================
#
#   Method:       DESTROY
#
#   Description:  csv->DESTROY()
#
#   Parameters:   $self       
#
#   Returns:      
#
#   Globals:
#
#===============================================================================

sub DESTROY {
  my ($self) = @_;
  $self->close();
  }


#===============================================================================
#
#   Method:       close
#
#   Description:  Explicitely closes the file.  
#                 csv->close()
#
#   Parameters:   $self       
#
#   Returns:      
#
#   Globals:
#
#===============================================================================

sub close {
  my ($self) = @_;
  my $fh = $self->{fh};
  no strict "refs";
  close($fh) if $fh;
  use strict "refs";
  $self->{fh} = undef;
  $fh = $self->{outputfh};
  no strict "refs";
  close($fh) if $fh;
  use strict "refs";
  $self->{outputfh} = undef;
  }

#===============================================================================
#
#   Method:       OpenOutputFile
#
#   Description:  Creates an output file to which records can be copied.
#                   $csv->SetOutput( $filename, option=>$value, .. )
#
#   Options can include:
#      def=>$def        definition of CSV file parameters
#      fields=>$fields  array ref of output field names.  Can include
#                       "outputfieldname=inputfieldname" to rename a field
#      noquote=>$fields array ref of output field names which are not quoted by
#                       by default
#
#
#   Returns:      none
#
#   Globals:
#
#===============================================================================
 

sub OpenOutputFile
{
    my ($self,$filename,@options) = @_;
    my $fh = $self->{outputfh};
    no strict "refs";
    CORE::close($fh) if $fh;
    use strict "refs";
    $self->{outputfh} = undef;

    open($fh,">",$filename) || die "Cannot open $filename\n";

    my $options = ref($options[0]) ? $options[0] : {@options};
    $self->{outputfile} = $filename;
    $self->{outputfh} = $fh;
    $self->{outputdef} = $options->{def} || $self;
    $self->{outputfields} = undef;
    $self->{noquote} = $options->{noquote};
    $self->{noquotefields} = [];
    if( ! $self->{noheaders} )
    {
        my $fields = $options->{fields};
        $fields = $self->columns if ! $fields;
        $fields = [split(' ',$fields)] if ! ref($fields);
        my @outputnames = ();
        my @inputnames = ();
        foreach my $f (@$fields)
        {
            my($of,$if) = split(/\=/,$f,2);
            $if = $of if ! $if;
            push(@outputnames,$of);
            push(@inputnames,$if);
        }

        $self->{outputfields} = \@outputnames;
        $self->{outputcols} = [$self->colno(@inputnames)];

        $fields = $self->{noquote};
        $fields = [split(' ',$fields)] if ! ref($fields);
        $self->{noquote} = $fields;
        my $nq = $self->{noquotefields};
        foreach my $f (@$fields)
        {
            $f = uc($f);
            foreach my $i (0..$#outputnames)
            {
                $nq->[$i] = 1 if uc($outputnames[$i]) eq $f;
            }
        }
    }

    my $def = $self->{outputdef};
    my ($delim, $quote, $literal, $doublequote) 
        = @{$def}{'delim','quote','literal','doublequote'};

    my $needquote = quotemeta($delim);
    $needquote .= "|" . quotemeta($quote) if $quote;
    $self->{needquotere} = qr/$needquote/;

    my $quotefield = "";
    $quotefield .= "\$d =~ s/".quotemeta($quote)."/".quotemeta($quote).quotemeta($quote)."/g;" 
       if $doublequote;
    $quotefield .= "\$d =~ s/".quotemeta($quote)."/".quotemeta($literal).quotemeta($quote)."/g;" 
       if $literal;

    my $subcode = 'sub { my($d)=@_;'.$quotefield.'return $quote.$d.$quote; }';
    $self->{quotesub} = eval $subcode;
    
    if( ! $self->{noheaders} && ! $self->{outputdef}->{noheaders} )
    {
        my $td = $self->{data};
        my $tnq = $self->{noquotefields};

        $self->{data} = $self->{outputfields};
        my $quoteheader = exists $options->{quoteheader} ? $options->{quoteheader} : 1;
        my $nqf=[];
        if( ! $quoteheader )
        {
            $nqf = [(1) x scalar(@{$self->{data}})];
        }
        $self->{noquotefields} = $nqf;
        $self->WriteRecord;

        $self->{data} = $td;
        $self->{noquotefields} = $tnq;
    }
} 

sub WriteRecord
{
    my($self) = @_;
    my $fh = $self->{outputfh};
    die "Cannot WriteRecord - output file not defined\n" if ! $fh;
    my $data = $self->{data};
    my $cols = $self->{outputcols};
    $cols = [0..$#$data] if ! $cols;

    my $noquote = $self->{noquotefields};
    my $qre =$self->{needquotere};
    my $qsub = $self->{quotesub};

    foreach my $i (@$cols)
    {
        my $d = $data->[$i];
        my $nq = $noquote->[$i];
        $nq = 0 if $nq && $d =~ $qre;
        $d = &$qsub($d) if ! $nq;
        print $fh $self->{outputdef}->{delim} if $i;
        print $fh $d;
    }
    print $fh "\n";
}


#===============================================================================
#
#   Method:       NextRecord
#
#   Description:  Reads the next record from the file.  Returns true if
#                 it exists.
#                   $gotrecord = $csv->NextRecord
#
#   Parameters:   
#
#   Returns:      $gotrecord   True if another record is available.
#
#   Globals:
#
#===============================================================================


sub NextRecord {
  my( $self) = @_;
  $self->{fields} = undef;
  $self->{data} = undef;
  my $fh = $self->{fh};
  no strict "refs";
  my $record = <$fh>;
  use strict "refs";
  my $input_line_number = $.; 

  return 0 if ! $record;
  $record =~ s/\r?\n$//;

  my $fldpat = $self->{fieldpattern};
  my $delim = $self->{delim};
  my $qliteral = quotemeta($self->{literal});  # Quoted literal character
  my $quote = $self->{quote};
  my $qquote = quotemeta($quote);
  my $doublequote = $self->{doublequote};
  my $rc = $delim.$record;
  my $unclosed_quote = $self->{unclosed_quote};
  if( $unclosed_quote) {
     while( $rc =~ /^(?:$fldpat)*(?:$unclosed_quote)$/ ) {
        no strict "refs";
        my $nextline = <$fh>;
        use strict "refs";
        last if ! $nextline;
        $nextline =~ s/\r?\n$//;
        $rc .= "\n".$nextline;
        }
     }

  die "Bad record in ".$self->{name}." line $input_line_number\n$record\n" 
      if $rc !~ /^(?:$fldpat)*$/;
  my @data=();
  while( $rc =~ /\G(?:$fldpat)/g ) {
    my $d = $1;
    if( $quote ne '' ) {
      if( $d eq '' ) {
        $d = $2;
        }
      else {
        $d =~ s/$qquote$qquote/$quote/g if $doublequote;
        }
      }
    $d =~ s/$qliteral(.)/$1/g if $qliteral;
    push(@data,$d);
    }
  $self->{data} = \@data;
  $self->{recno}++;
  $self->{lineno} = $input_line_number;
  return 1;
  }


#===============================================================================
#
#   Method:       columns
#
#   Description:  Returns the list of column names.
#                  @columns = $csv->columns()
#
#   Parameters:   $self       
#
#   Returns:      
#
#   Globals:
#
#===============================================================================


sub columns {
  my($self) = @_;
  return $self->{columns};
  }


#===============================================================================
#
#   Method:       column
#
#   Description:  Returns the name for the specified column number
#                  $data = $csv->column($i)
#                  @data = $csv->column(@i)
#
#   Parameters:   @i          A list of field numbers (the first field is 0)
#
#   Returns:      @data       The names of the fields
#
#   Globals:
#
#===============================================================================



sub column {
  my ($self,@i) = @_;
  my @values = @{$self->{columns}}[@i];
  return wantarray ? @values : $values[0];
  }


#===============================================================================
#
#   Method:       colno
#
#   Description:  Returns the column number for a specified name
#                  $i = $csv->colno($fieldname)
#                  @i = $csv->colno(@fieldnames)
#
#   Parameters:   @fieldnames     The names of the fields
#
#   Returns:      @i              The column numbers
#
#   Globals:
#
#===============================================================================



sub colno {
  my ($self,@fields) = @_;
  my $getcolnofunc = $self->{getcolnofunc};
  my @values = map { &$getcolnofunc($_)} @fields;
  return wantarray ? @values : $values[0];
  }

#===============================================================================
#
#   Method:       hascolumns 
#
#   Description:  Tests that the file has the specified columns
#                  if( $csv->hascolumns(@fieldnames) ) { ... }
#
#   Parameters:   @fieldnames     The names of the fields
#
#   Returns:      $ok             True (1) if the file has the columns, 0 ow.
#
#   Globals:
#
#===============================================================================


sub hascolumns {
  my($self,@fields) = @_;
  foreach my $v ($self->colno(@fields))
  {
      return 0 if ! defined($v);
  }
  return 1;
  }

#===============================================================================
#
#   Method:       fields
#
#   Description:  Returns the current record as an array reference
#                  $data = $csv->fields()
#
#   Parameters:   
#
#   Returns:      $data      An array reference to an array of field values
#
#   Globals:
#
#===============================================================================


sub fields {
  my($self) = @_;
  return $self->{data};
  }


#===============================================================================
#
#   Method:       field
#
#   Description:  Return data for specific fields
#                  @data = $csv->field(@fields)
#
#   Parameters:   @fields     An array of field names
#
#   Returns:      @data       An array of corresponding field values
#
#   Globals:
#
#===============================================================================



sub field {
  my ($self,@fields) = @_;
  my $getfieldfunc = $self->{getfieldfunc};
  my @values = map { &$getfieldfunc($_)} @fields;
  return wantarray ? @values : $values[0];
  }


#===============================================================================
#
#   Method:       setfield
#
#   Description:  Sets the value of a data field in a CSV file record, replacing
#                 the value read from the file
#                  $csv->setfield($field,$value)
#
#   Parameters:   $fields     The field name
#                 $value      The value to store in the field
#
#   Globals:
#
#===============================================================================

sub setfield {
  my( $self,%values ) = @_;
  my $setfieldfunc = $self->{setfieldfunc};
  foreach my $k (keys %values)
  {
     &$setfieldfunc($k,$values{$k});
  }
  }


#===============================================================================
#
#   Method:       fieldno
#
#   Description:  Returns data for fields specified by field number
#                  @data = $csv->fieldno(@i)
#
#   Parameters:   @i          The field numbers
#
#   Returns:      @data       The corresponding data
#
#   Globals:
#
#===============================================================================



sub fieldno {
  my ($self,@i) = @_;
  my @values = @{$self->{data}}[@i];
  return wantarray ? @values : $values[0];
  }


#===============================================================================
#
#   Method:       IsBlank
#
#   Description:  Returns true if the entire record is blanks
#                  $blank = $csv->IsBLank
#
#   Parameters:   
#
#   Returns:      $blank   1 if any field contains data, 0 otherwise
#
#   Globals:
#
#===============================================================================


sub IsBlank {
  my ($self) = @_;
  foreach (@{$self->{data}} ) {
    return 0 if $_ ne '';
    }
  return 1;
  }

1;

