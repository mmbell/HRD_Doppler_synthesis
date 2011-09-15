#!/usr/bin/perl -w
use POSIX;

# Cedric to Grads converter
# By Michael Bell
# Thanks to Paul Hein & Brenda Dolan at Colorado State University
# and Stacy Brodzik of University of Washington for much of the code used here
#

######## Change the unsigned short to a signed short ########

sub make_signed_short {
   my $ref = $_[0];
   my $i=0;
   foreach (@$ref) {
      if ( $$ref[$i] >= 32768 ) {
         $$ref[$i] = $$ref[$i] - 65536;
      }
      $i++;
   }
}

sub swapwords {
   my $ref = $_[0];
   my $lastindex = $_[1];
   my $i=0;
   my $tmp;
   for ($i=0; $i < $lastindex; $i=$i+2) {
      $tmp = $$ref[$i];
      $$ref[$i] = $$ref[$i+1];
      $$ref[$i+1] = $tmp;
   }
}

sub swapdatawords {
   my $ref = $_[0];
   my $xysize = $_[1];
   my $zsize = $_[2];
   my $i=0;
   my $tmp;
   for ($k=0; $k < $zsize; $k++) {
     for ($i=0; $i < $xysize; $i=$i+2) {
       $ii = $xysize*$k + $i;
       if ($i+1 >= $xysize || $ii+1 >= $#$ref) { last; }
       $tmp = $$ref[$ii];
       $$ref[$ii] = $$ref[$ii+1];
       $$ref[$ii+1] = $tmp;
     }
   }
}

%monhash = (
    1 => 'JAN',
    2 => 'FEB',
    3 => 'MAR',
    4 => 'APR',
    5 => 'MAY',
    6 => 'JUN',
    7 => 'JUL',
    8 => 'AUG',
    9 => 'SEP',
    10 => 'OCT',
    11 => 'NOV',
    12 => 'DEC',
    );

$pi = acos(-1);

######### The Start Of Main ############

###### Get filename from command line ######

$machine = "";
if ($#ARGV != 1 && $#ARGV != 0) {die "Usage:  ced2grads [sun | intel] filename\n";}
if ($#ARGV == 0) {
   $infile = $ARGV[0];
} else {
   $machine = $ARGV[0];
   $infile = $ARGV[1];
}

# Open input file and change name for output file

open(IN, "< $infile") || die "Unable to open $infile\n";
binmode IN;
$outfile = $infile;
if ($outfile =~ /ced$/) {
   $outfile =~ s/\.ced/\.dat/;
} elsif ($outfile =~ /mud$/) {
   $outfile =~ s/\.mud/\.dat/;
} else {
   die "INPUT FILE NAME ERROR: $infile\n";
}
$ctlfile = $infile;
if ($ctlfile =~ /ced$/) {
    $ctlfile =~ s/\.ced/\.ctl/;
} elsif ($ctlfile =~ /mud$/) {
    $ctlfile =~ s/\.mud/\.ctl/;
} else {
    die "INPUT FILE NAME ERROR: $infile\n";
}

print "$infile will be converted to $ctlfile, $outfile\n";

####### Read in ced file headers pulling off needed info ######

# Make certain it is a ced file

   $buf="";
   $knt = read(IN, $buf, 1540);
   if($knt != 1540) {die "Error at ced header -  knt = $knt";}
   @ahead = unpack("a*",$buf);
   if (substr($ahead[0],0,4) ne "CED1") {die " Not a Cedric file\n";}

# determine correct byte ordering

if ($machine eq "intel") {
   $s = "v";
} elsif ($machine eq "sun") {
   $s = "n";
} else {
   $s = "v";
   @shead = unpack("V*",$buf);
   if ($shead[1] != 1) {
      @shead = unpack("N*",$buf);
      $s = "n";
      if ($shead[1] != 1) { 
         print "WARNING Byte Ordering Flag NOT Found!\n";
         print "Assuming Sun/Network (Big Endian) byte order\n";
         print "If errors use: ced2nc intel file.ced\n";

#        print "If ced file was created on Sun then convert on Sun, else\n";
#        print "if created on Intel Linux then convert on Intel Linux.\n";
#        $t1 = 1;
#        $t1p=pack("s",$t1);
#        $t1v=unpack("v",$t1p);
#        if ($t1v == 1) {
#           print "Assuming Intel/Vax (Little Endian) byte order was used to create this file\n";
#           $s = "v";
#        } else {
#           print "Assuming Sun/Network (Big Endian) byte order was used to create this file\n";
#           $s = "n";
#        }
      }

   }
}

# read in volume (i.e. mud) header

   $knt = read(IN, $buf, 1020);
   if($knt != 1020) {die "Error at mud header -  knt = $knt";}
   @head = unpack("$s*",$buf);
# Swap bytes here
if ($s eq "v") {
    &make_signed_short(\@head);
    &swapwords(\@head,$#head);
}

# get the id and title
   if($head[14] == 0){
   $buf = pack("$s*",($head[12],$head[13]));
   }
   else{
   $buf = pack("$s*",($head[12],$head[13],$head[14]));
   }

   $ident = unpack("a*",$buf);
   $datet=$head[115]+1900;
   $date=substr($datet,2,2);
   $title = "Cartesian data from ".$ident." at ".$date."/".$head[116]."/".
            $head[117]." ".$head[118].":".$head[119].":".$head[120];
   $intmissing = $head[66];
   #$datestr = $datet."-".$head[116]."-".
   #  $head[117]." ".$head[118].":".$head[119].":".$head[120];
if ($head[118] < 10) {
    $hr = "0".$head[118];
} else {
    $hr = $head[118];
}
if ($head[119] < 10) {
    $min = "0".$head[119];
} else {
    $min = $head[119];
}
if ($head[117] < 10) {
    $day = "0".$head[117];
} else {
    $day = $head[117];
}
$timestr = $hr.":".$min."Z".$day.$monhash{$head[116]}.$datet;
print "$title\n";

# get grid x,y,z info

$dimsiz[0] = $head[161];
$dimnam[0] = "longitude";
$xfirst = $head[159]/$head[67];
$xdelta = $head[162]/1000.0;

$dimsiz[1] = $head[166];
$dimnam[1] = "latitude";
$yfirst = $head[164]/$head[67];
$ydelta = $head[167]/1000.0;

$dimsiz[2] = $head[171];
$dimnam[2] = "altitude";
$zfirst = $head[169]/1000.0;
$zdelta = $head[172]/1000.0;

$dimsiz[3] = 1;
$dimnam[3] = "time";


# get grid info
$lat_origin = $head[32] + $head[33]/60 + $head[34]/3600/100;
$lon_origin = ($head[35] + $head[36]/60 + $head[37]/3600/100);
#$lat_origin = 24.7;
#$lon_origin = -78.5 + 360;
print "Lat, Lon: $lat_origin, $lon_origin\n";
print "dimsizes 0=$dimsiz[0] 1=$dimsiz[1] 2=$dimsiz[2] 3=$dimsiz[3]\n";

# get variable names/attributes

   $nvars = $head[174];
print "$nvars, $zfirst\n";
   $idx = 175;
   for ($i=0; $i<$nvars; $i++) {
      $buf = pack("$s*",($head[$idx],$head[$idx+1],$head[$idx+2],$head[$idx+3]));
      $varnam[$i] = unpack("a*",$buf);
      $j = index($varnam[$i]," ");
      $varnam[$i] = substr($varnam[$i],0,$j);
      $scale[$i] = $head[$idx+4];
print "$varnam[$i]\n";
      $idx = $idx + 5;
   }
# if same name need to change
   for ($i=0; $i<$nvars; $i++) {
     for ($ii=0; $ii<$nvars; $ii++) {
       if (($varnam[$i] eq $varnam[$ii]) && ($i != $ii)) {
          $varnam[$ii] = $varnam[$ii].1;
          print "$varnam[$i] already exists so changing name to $varnam[$ii]\n";
       }
     }
   }

   @head = ();

###### Write out information to the grads file ######

print "Writing out Grads variable and attribute information to ctl\n";
open(OUT, ">$ctlfile") or die("Couldn't open $ctlfile $!\n");
print OUT "dset ^$outfile\n";
print OUT "title ELDORA data\n";
print OUT "undef $intmissing\n";
$latrad = $lat_origin * $pi/180.0;
$fac_lat = 111.13209 - 0.56605 * cos(2.0 * $latrad)
    + 0.00012 * cos(4.0 * $latrad) - 0.000002 * cos(6.0 * $latrad);
$fac_lon = 111.41513 * cos($latrad)
    - 0.09455 * cos(3.0 * $latrad) + 0.00012 * cos(5.0 * $latrad);
$lonstart = $lon_origin + $xfirst/$fac_lon;
$londelta = $xdelta/$fac_lon;
print OUT "xdef $dimsiz[0] levels\n";
for $x (0 .. $dimsiz[0]-1) {
  $xloc = $lon_origin + ($xfirst+($x-1)*$xdelta)/$fac_lon;
  print OUT " $xloc\n";
}
$latstart = $lat_origin + $yfirst/$fac_lat;
$latdelta = $ydelta/$fac_lon;
print OUT "ydef $dimsiz[1] levels\n";
for $y (0 .. $dimsiz[1]-1) {
  $yloc = $lat_origin + ($yfirst+($y-1)*$ydelta)/$fac_lat;
  print OUT " $yloc\n";
}
print OUT "zdef $dimsiz[2] levels\n";
for $i (0 .. $dimsiz[2]-1) {
    $alt[$i] = $zfirst + ($zdelta * $i);
    print OUT "$alt[$i] ";
} print OUT "\n";

print OUT "tdef 1 linear $timestr 1MN\n";
print OUT "vars ".$nvars."\n";
for ($i=0; $i<$nvars; $i++) {
    print OUT "$varnam[$i] $dimsiz[2] 99 variable\n";
}
print OUT "endvars\n";
close OUT;

# Read in data values
print "Reading in Cedric data and writing out Grads data\n";

$size = $dimsiz[0]*$dimsiz[1];
for ($lev=0; $lev<$dimsiz[2]; $lev++) {

 $knt = read(IN, $buf, 20);
 if($knt != 20) {die "Error at Level header $lev -  knt = $knt";}

 for ($i=0; $i<$nvars; $i++) {

   $knt = read(IN, $buf, 2*$size);
   if($knt != 2*$size) {die "Error at ced data $i -  knt = $knt";}
   if ($lev == 0) {
     $tmpfile[$i] = "ced2grads".$i.".".$$;
     $out[$i] = $tmpfile[$i]."handle";
     open( $out[$i], "> $tmpfile[$i]") || die " unable to open $tmpfile[$i]\n";
     binmode $out[$i];
   }
   print { $out[$i] } $buf;
 }
}
$buf = "";
for ($i=0; $i<$nvars; $i++) { close($out[$i]); }

print "Data in tmp files - ready to unpack, scale and write\n";

##Write out data values

open(DAT, ">$outfile") or die("Couldn't open $outfile $!\n");
binmode DAT;
for ($i=0; $i<$nvars; $i++) {


     open( $out[$i], "$tmpfile[$i]") || die "Unable to open $tmpfile[$i]\n";
     binmode $out[$i];
     @ls = stat $tmpfile[$i];
     $sizels = $ls[7];
     $knt = read($out[$i], $buf, $sizels);
     if($knt != $sizels) {die "Error at tmp ced data $i -  knt = $knt";}
     close($out[$i]);
     system("rm $tmpfile[$i]");

     @data = unpack("$s*",$buf);
     &make_signed_short(\@data);
     &swapdatawords(\@data,$size,$dimsiz[2]);
     # &swapdatawords(\@data,$dimsiz[0],$dimsiz[1],$dimsiz[2]);

     for ($j=0; $j<$size*$dimsiz[2]; $j++) {
        if ($data[$j] == $intmissing) {
           $data[$j] = $intmissing;
        } else {
           $data[$j] = $data[$j]/$scale[$i];
        }
	# Pack the float data
	$buf = pack("f*", $data[$j]);
	#print "$j, $data[$j]\n";
	print DAT $buf;
	#$stop = <STDIN>;

	$buf="";
     }

}

##Close output file
close(DAT);

print "Done!\n";
exit(0);

###############The End####################
