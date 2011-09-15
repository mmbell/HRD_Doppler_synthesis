#!/usr/bin/perl -w

# Cedric to Zebra netcdf converter
# By Michael Bell & Stacy Brodzik
# Thanks to Paul Hein & Brenda Dolan at Colorado State University
# for much of the code used here
#

use NetCDF;

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


######### The Start Of Main ############

###### Get filename from command line ######

$machine = "";
if ($#ARGV != 1 && $#ARGV != 0) {die "Useage:  ced2nc4swap [sun | intel] filename\n";}
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
   $outfile =~ s/\.ced/\.nc/;
} elsif ($outfile =~ /mud$/) {
   $outfile =~ s/\.mud/\.nc/;
} else {
   die "INPUT FILE NAME ERROR: $infile\n";
}
print "$infile  ";
print "$outfile\n";

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
         print "WARNING Byte Ordering Flag NOT Founded!\n";
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
   print"\nDATE: $date";
   $title = "Cartesian data from ".$ident." at ".$date."/".$head[116]."/".
            $head[117]." ".$head[118].":".$head[119].":".$head[120];
   $intmissing = $head[66];
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
if ($head[116] < 10) {
    $mon = "0".$head[116];
} else {
    $mon = $head[116];
}
if ($head[120] < 10) {
    $sec = "0".$head[120];
} else {
    $sec = $head[120];
}
$linuxstr = $datet."-".$mon."-".$day." ".$hr.":".$min.":".$sec;
$macstr = $mon.$day.$hr.$min.$datet.".".$sec;
#   $datestr = $datet."-".$head[116]."-".
#     $head[117]." ".$head[118].":".$head[119].":".$head[120];

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
$lon_origin = ($head[35] + $head[36]/60 + $head[37]/3600/100) + 360.;
#$lat_origin = 24.7;
#$lon_origin = -78.5 + 360;
print "Lat, Lon: $lat_origin, $lon_origin\n";
print "$head[32] $head[33]  $head[34] : $head[35] $head[36] $head[37]\n";
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

   $ndims = 4;
   $missing = -999.0;
   @head = ();

###### Write out information to the netcdf file ######

print "Writing out NetCDF variable and attribute information\n";

##Create generic netcdf file - output

$ncid2 = NetCDF::create($outfile, NetCDF::CLOBBER);
die "Couldn't open/create netCDF file\n" if $ncid2 < 0;

##Write out dimension info

for ($i=0; $i<$ndims; $i++) {
   $status =  NetCDF::dimdef($ncid2, $dimnam[$i], $dimsiz[$i]);
   die "Couldn't write dimension info\n" if $status < 0;
}

##Write out variable info (Note number of var dims parameter is not in the perl version)

for ($i=0; $i<$nvars; $i++) {
   @vdims1 = ();
   for ($j=0; $j<4; $j++) {
       $vdims1[$j] = 3 - $j;
   }
   $status =  NetCDF::vardef($ncid2, $varnam[$i], NetCDF::FLOAT,
                             \@vdims1);
   die "Couldn't write variable info\n" if $status < 0;

   ##Write out attribute info (Note length parameter is not in the perl version)
   if ($varnam[$i] =~ /D/) {
       # Probably reflectivity
       $unit = "dBZ";
   } elsif ($varnam[$i] =~ /V/) {
       $unit = "m/s";
   } else {
       $unit = "Unknown";
   }

   $status =  NetCDF::attput($ncid2, $i, "units",NetCDF::CHAR,
			     \$unit);
   die "Couldn't write attribute info\n" if $status < 0;

##Write out attribute info (Note length parameter is not in the perl version)

     $status =  NetCDF::attput($ncid2, $i, "missing_value",NetCDF::FLOAT,
                               \$missing);
     die "Couldn't write attribute info\n" if $status < 0;
}


$lonunits = "degrees_E";
$status =  NetCDF::vardef($ncid2, 'longitude', NetCDF::FLOAT, \( 0 ));
die "Couldn't write variable info\n" if $status < 0;
$status =  NetCDF::attput($ncid2, $nvars, 'units', NetCDF::CHAR, \$lonunits);
die "Couldn't write attribute info\n" if $status < 0;

$latunits = "degrees_N";
$status =  NetCDF::vardef($ncid2, 'latitude', NetCDF::FLOAT, \( 1 ));
die "Couldn't write variable info\n" if $status < 0;
$status =  NetCDF::attput($ncid2, $nvars+1, 'units', NetCDF::CHAR, \$latunits);
die "Couldn't write attribute info\n" if $status < 0;

$altunits = "km";
$status =  NetCDF::vardef($ncid2, 'altitude', NetCDF::FLOAT, \( 2 ));
die "Couldn't write variable info\n" if $status < 0;
$status =  NetCDF::attput($ncid2, $nvars+2, 'units', NetCDF::CHAR, \$altunits);
die "Couldn't write attribute info\n" if $status < 0;

$timeunits = "seconds since 1970-01-01 0:00:00 0:00";
$status =  NetCDF::vardef($ncid2, 'time', NetCDF::LONG, \( 3 ));
die "Couldn't write variable info\n" if $status < 0;
$status =  NetCDF::attput($ncid2, $nvars+3, 'units', NetCDF::CHAR, \$timeunits);
die "Couldn't write attribute info\n" if $status < 0;

$status =  NetCDF::endef($ncid2);
die "Couldn't leave define mode\n" if $status < 0;

##########Read in and Write out data values############### 

print "Reading in Cedric data and writing out NetCDF data\n";

##Write out coordinate values

##Create coordinate variables and attributes
for $i (0 .. $dimsiz[0]) {
    for $j (0 .. $dimsiz[1]) {
	$lat[$j] = $lat_origin + ($yfirst +$j * $ydelta)/111.2; 
        $lon[$i] = $lon_origin + ($xfirst +$i * $xdelta)/(111.2 * cos($lat[$j]/57.296));
    }
}
for $i (0 .. $dimsiz[2] ) {
  $alt[$i] = $zfirst + ($zdelta * $i);
}
$status =  NetCDF::varput($ncid2, $nvars, (0), ($dimsiz[0]), \@lon);
die "Couldn't write out values\n" if $status < 0;
$status =  NetCDF::varput($ncid2, $nvars+1, (0), ($dimsiz[1]), \@lat);
die "Couldn't write out values\n" if $status < 0;
$status =  NetCDF::varput($ncid2, $nvars+2, (0), ($dimsiz[2]), \@alt);
die "Couldn't write out values\n" if $status < 0;

# Convert time to seconds and record as time
$platform = `uname`;
if ($platform =~ /Linux/) {
	$inttime = `date -u -d '$linuxstr' +%s`;
print "datestr=$linuxstr and inttime=$inttime\n";
} elsif ($platform =~ /Darwin/) {
	$inttime = `date -j -u +%s "$macstr"`;
print "datestr=$macstr and inttime=$inttime\n";
} else {
	print "Not sure about $platform, trying Linux date command...\n";
	$inttime = `date -u -d '$linuxstr' +%s`;
}
chomp($inttime);
@unixtime = ( $inttime );
$status =  NetCDF::varput($ncid2, $nvars+3, (0), ($dimsiz[3]), \@unixtime);
die "Couldn't write out values\n" if $status < 0;


# Read in data values

$size = $dimsiz[0]*$dimsiz[1];
for ($lev=0; $lev<$dimsiz[2]; $lev++) {

 $knt = read(IN, $buf, 20);
 if($knt != 20) {die "Error at Level header $lev -  knt = $knt";}

 for ($i=0; $i<$nvars; $i++) {

   $knt = read(IN, $buf, 2*$size);
   if($knt != 2*$size) {die "Error at ced data $i -  knt = $knt";}
   if ($lev == 0) {
     $tmpfile[$i] = "ced2nc".$i.".".$$;
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

     for ($j=0; $j<4; $j++) {
         $start[$j] = 0;
         $count[$j] = $dimsiz[3 - $j];
     }

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
           $data[$j] = $missing;
        } else {
           $data[$j] = $data[$j]/$scale[$i];
        }
     }

     $status =  NetCDF::varput($ncid2, $i, \@start, \@count, \@data);
     die "Couldn't write out values\n" if $status < 0;
}

##Close output file

$status =  NetCDF::close($ncid2);
die "Couldn't close outfile\n" if $status < 0;

print "Done!\n";
exit(0);

###############The End####################
