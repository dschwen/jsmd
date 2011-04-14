#!/usr/bin/perl

use File::Temp qw(tempfile);

my $source = 0;
my $tmp_fh;

while(<STDIN>) {
  if( $source == 0 ) {
    print $_;
  
    if( /<script/ && !/<\/script>/ ) {
      $source = 1;
      $tmp_fh = new File::Temp( UNLINK => 1 );
    }
  } else {
    if( /<\/script>/ ) {
      $source = 0;
      print `yui-compressor --type js  $tmp_fh`;
      print "\n$_";
    } else {
      print $tmp_fh $_;
    }
  }
  $source .= $_;
}
