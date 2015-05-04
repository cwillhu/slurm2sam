#!/usr/local/bin/perl -w
#
#Author:
# Chris Williams
# Harvard Informatics And Scientific Applications
# http://informatics.fas.harvard.edu

use warnings FATAL => "all";
use strict;
use POSIX;
use Getopt::Long;
use Data::Dumper;
use List::MoreUtils qw(uniq);
use FindBin qw($Bin);
use lib "$Bin/../lib/perl5";
use lib "$Bin/../lib/site_perl";
use Time::Piece;
use Time::Seconds;
use Date::Parse;
use Number::Format qw(:subs :vars);
use constant { true => 1, false => 0 };

$KILO_SUFFIX = "KB";
$MEGA_SUFFIX = "MB";
$GIGA_SUFFIX = "GB";

my $jFile = "";
my $outname = "SlurmHistory";
my $help = 0;

Getopt::Long::GetOptions(
  'jobfile=s' => \$jFile,
  'name:s'    => \$outname,
  'help'      => \$help,
  ) or die "Incorrect input! Use -h for usage.\n";

if ($help) {
  print "\nUsage: slurm2sam.pl -jobfile <jobfile> [Options]\n";
  print "Options:\n";
  print "  -jobfile    Job history file, as created by slurmhistory.sh\n";
  print "  -name       Output name\n";
  print "  -help       Display usage information.\n";
  exit 0;
}

($jFile eq "") && die "Error: A job history file must be provided. Use -h for usage.\n";

my $samFile       = $outname . ".jobs.sam";
my $hours_gffFile = $outname . ".hours.gff3";
my $days_gffFile  = $outname . ".days.gff3";
my $fastaFile     = $outname . ".genome.fasta";
my $statsFile     = $outname . ".general_stats.NOT_FOR_IGV.txt";
my $jref = {};
my ($minstart, $maxend);

#read job history file
($jref, $minstart, $maxend) = readJobsFile($jFile);

#write out a few stats
writeStats($statsFile, $jref, $minstart, $maxend);

#extend beginning and ending time points, and have them occur at midnight (looks better in IGV).
my $minstartObj = Time::Piece->strptime($minstart, "%s") - 3*ONE_DAY; 
my $mindayObj = Time::Piece->strptime($minstartObj->year." ".$minstartObj->month." ".$minstartObj->mday, "%Y %b %d");
my $maxendObj = Time::Piece->strptime($maxend, "%s") + 3*ONE_DAY;
my $maxdayObj = Time::Piece->strptime($maxendObj->year." ".$maxendObj->month." ".$maxendObj->mday, "%Y %b %d");
$minstart = $mindayObj->epoch;
$maxend = $maxdayObj->epoch;

#write the input files for IGV
writeSam($samFile, $jref, $minstart);
writeGff($days_gffFile, $hours_gffFile, $minstart, $maxend);
writeFasta($fastaFile, $minstart, $maxend);

exit 0;

sub writeStats { 
  my $statsFile = shift;
  my $jref = shift;
  my $minstart = shift;
  my $maxend = shift;
  my %stateCount;

  my @users = uniq sort map { $_->{"User"} } values %$jref;
  my @partitions = uniq sort map { $_->{"Partition"} } values %$jref;
  $stateCount{$_->{"State"}}++ for values %$jref;
  my $njobs = Number::Format::format_number(scalar keys %$jref);

  print $#users+1, " unique users\n";
  print "$njobs jobs\n";

  open(STATS, ">$statsFile") or die "Can't open file $statsFile\n";
  print STATS "\nJob History Statistics\n";
  print STATS "----------------------\n\n";
  print STATS $#partitions+1, " partition(s): ", join(", ", @partitions), "\n";
  print STATS $#users+1, " unique users\n";
  print STATS "$njobs jobs\n\n";
  print STATS "Earliest time point: ", Time::Piece->strptime($minstart, "%s")->cdate, "\n";
  print STATS "Latest time point: ", Time::Piece->strptime($maxend, "%s")->cdate, "\n\n";
  print STATS "Job State Counts:\n";
  printf STATS " %-11s %d\n", "$_:", $stateCount{$_} for sort { $stateCount{$b} <=> $stateCount{$a} } keys %stateCount;
  close(STATS);
}

sub readJobsFile {
  my $jFile = shift;
  open(IN, "<$jFile") or die "Can't open file $jFile\n";
  chomp(my $header = <IN>);

  my @types = split('\|',$header);
  my $minstart = Date::Parse::str2time("2100", "GMT");
  my $maxend = Date::Parse::str2time("1995", "GMT");
  my %jobs;
  my $i=1;
  my $skip = "";
  while (my $line = <IN>) {
    chomp($line);
    my @vals = split('\|',$line);

    my %temp;
    @temp{@types} = @vals;
    #print Dumper(\%temp), "\n\n";

    my $id = $temp{"JobID"}; #get jobid root (ie, without trailing ".0", ".batch", etc.)
    $id =~ s/\..*//;

    if ($temp{"Start"} eq "Unknown" || $temp{"End"} eq "Unknown") {
      $skip = $id;
      next;
    }
  
    if ($id eq $skip) {
      next;
    }
 
    if (!defined($jobs{$id})) { #we're at a new job
      $jobs{$id}{"User"} = $temp{"User"};
      $jobs{$id}{"JobName"} = $temp{"JobName"};
      $jobs{$id}{"Partition"} = $temp{"Partition"};
      $jobs{$id}{"AllocCPUS"} = $temp{"AllocCPUS"};
      $jobs{$id}{"NNodes"} = $temp{"NNodes"};
      $jobs{$id}{"ReqMem"} = reqmem2num($temp{"ReqMem"}, $temp{"NNodes"}, $temp{"AllocCPUS"});
      $jobs{$id}{"MaxRSS"} = bytes2num($temp{"MaxRSS"});
      $jobs{$id}{"State"} = $temp{"State"} =~ s/ by \d+//r;
      $jobs{$id}{"Start"} = Date::Parse::str2time($temp{"Start"}, "GMT");
      $jobs{$id}{"End"} = Date::Parse::str2time($temp{"End"}, "GMT");
      $jobs{$id}{"Suspended"} = $temp{"Suspended"};
      $jobs{$id}{"TotalCPU"} = $temp{"TotalCPU"};
      $jobs{$id}{"CPUTime"} = $temp{"CPUTime"};
      $jobs{$id}{"Submit"} = Date::Parse::str2time($temp{"Submit"}, "GMT");
    } else {  #we're at a step of the current job. Update stats if necessary.
      my $step_maxrss = bytes2num($temp{"MaxRSS"});
      my $step_start = Date::Parse::str2time($temp{"Start"}, "GMT");
      my $step_end = Date::Parse::str2time($temp{"End"}, "GMT");
      my $step_totalcpu = stringToSeconds($temp{"TotalCPU"});
      my $step_cputime = stringToSeconds($temp{"CPUTime"});
      if ($step_maxrss > $jobs{$id}{"MaxRSS"}) { $jobs{$id}{"MaxRSS"} = $step_maxrss; }
      if ($step_start < $jobs{$id}{"Start"}) { $jobs{$id}{"Start"} = $step_start; }
      if ($step_end > $jobs{$id}{"End"}) { $jobs{$id}{"End"} = $step_end; }
      if ($step_totalcpu > stringToSeconds($jobs{$id}{"TotalCPU"})) { $jobs{$id}{"TotalCPU"} = $temp{"TotalCPU"}; }
      if ($step_cputime > stringToSeconds($jobs{$id}{"CPUTime"})) { $jobs{$id}{"CPUTime"} = $temp{"CPUTime"}; }
      if ($jobs{$id}{"AllocCPUS"} < $temp{"AllocCPUS"}) { $jobs{$id}{"AllocCPUS"} = $temp{"AllocCPUS"}; }
    }
    ($minstart, $maxend) = update_endpoints($minstart, $maxend, $jobs{$id}{"Start"},$jobs{$id}{"End"}); 

    $i++;

    if ($i == 50000) {
      last;
    }
  }
  close(IN);
  return (\%jobs, $minstart, $maxend);
} # end readJobsFile

sub writeSam {  #write out job stats in SAM format
  my $samFile = shift;
  my $jref = shift;
  my $minstart = shift;

  open(SAM, ">$samFile") or die "Can't open file $samFile\n";
  print(SAM "\@HD VN:1.5 SO:sorted\n");
  foreach my $id  (    
      sort { $jref->{$a}{"Start"} <=> $jref->{$b}{"Start"} } keys %$jref
      ) {
    my $name = $jref->{$id}{"JobName"};
    my $user = $jref->{$id}{"User"};
    my $state = $jref->{$id}{"State"};
    my $submit = Time::Piece->strptime($jref->{$id}{"Submit"}, "%s")->cdate;
    my $suspended = $jref->{$id}{"Suspended"};
    my $partition = $jref->{$id}{"Partition"};
    my $nnodes = $jref->{$id}{"NNodes"};

    my $startObj = Time::Piece->strptime($jref->{$id}{"Start"}, "%s");
    my $endObj = Time::Piece->strptime($jref->{$id}{"End"}, "%s");
    my $startString = $startObj->cdate;
    my $endString = $endObj->cdate;
    my $startIndex = floor( ($jref->{$id}{"Start"} - $minstart)/60 ) + 1;  #gff index must be 1-based; shift start and end to the right by 1.
    my $endIndex = floor( ($jref->{$id}{"End"} - $minstart)/60 ) + 1;
    my $length = $endIndex - $startIndex + 1;
    my $duration = Time::Seconds->new($length*60);
    my $durationFormatted = $duration->pretty;

    my $reqmem = $jref->{$id}{"ReqMem"};
    my $reqmemFormatted = ($reqmem == -1) ? "NA" : Number::Format::format_bytes($reqmem);
    my $maxrss = $jref->{$id}{"MaxRSS"};
    my $memFormatted = ($maxrss == -1) ? "NA" : Number::Format::format_bytes($maxrss);
    my $memUsageTag = getMemUsageTag($maxrss,$reqmem);

    my $alloccpus = $jref->{$id}{"AllocCPUS"};
    my $coresUsed = calcCores($jref->{$id}{"TotalCPU"},$jref->{$id}{"CPUTime"});
    my $coresUsedFormatted = sprintf("%.2f",$coresUsed);
    my $coresUsageTag = getCoresUsageTag($coresUsed, $alloccpus);

    #Query template name, flag, refname, 1-based leftmost pos, mapping qual, cigar, rnext pos, pnext pos, observed template len, seqment seq, qual, optional
    print(SAM join("\t","$id","0","1","$startIndex","255","${length}M","*","0","0","*","*",
      "A1:Z:AveCoresUsed: $coresUsedFormatted ($alloccpus allocated), ReqMem: $reqmemFormatted, MaxRAM: $memFormatted, JobID: $id, Time (Start to End): $durationFormatted (Suspended time: $suspended)",
      "A2:Z:Dates: Submit: $submit, Start: $startString, End: $endString",
      "CA:Z:Allocated Cores: $alloccpus", 
      "CU:Z:CPU Tag: $coresUsageTag",     # has custom color map
      "ME:Z:Memory Tag: $memUsageTag",    # has custom color map
      "NO:Z:Nodes: $nnodes",              
      "PA:Z:Partition: $partition",    
      "ST:Z:State: $state",               
      "US:Z:User: $user\n"));
  }
  close(SAM);
  
} # end writeSam

sub writeGff {  #write out timeline annotations in GFF format
  my $dayFile = shift;
  my $hourFile = shift;
  my $minstart = shift;
  my $maxend = shift;

  my $numSeconds = $maxend - $minstart + 1;
  my $numDays = floor($numSeconds / (60 * 60 * 24)) + 1;  #add one to round up and ensure numDays is at least 1.
  my $minstartObj = Time::Piece->strptime($minstart, "%s");  #the earliest timepoint
  my $maxendObj = Time::Piece->strptime($maxend, "%s");  #the latest timepoint

  open(GFFD, ">$dayFile") or die "Can't open file $dayFile\n";
  open(GFFH, ">$hourFile") or die "Can't open file $hourFile\n";
  for (my $d=1; $d <= $numDays+1; $d++) {
    #write day annotation for this day
    my $daystartObj = $minstartObj + ($d-1)*ONE_DAY;
    my $dayendObj = $minstartObj + $d*ONE_DAY;
    my $indexstart = floor(($daystartObj->epoch-$minstart)/60) + 1; #convert to minutes; add one to round up and start counting from 1.
    my $indexend = floor(($dayendObj->epoch-$minstart)/60) + 1;
    my $color = ($d % 2 == 0) ? "#814BAD" : "#B81113"; #color days alternately purple, brown.
    my $daystr = $daystartObj->mon."/".$daystartObj->mday;
    print(GFFD join("\t","1",$daystr,"gene",$indexstart,$indexend,"0",".","0","ID=".$daystr.";color=$color\n")); 

    #write hour annotations for this day
    for (my $h=1; $h <= 24; $h++) { 
      my $hourstartObj = $daystartObj + ($h-1)*ONE_HOUR;
      my $hourendObj = $daystartObj + $h*ONE_HOUR;
      my $indexstart = floor(($hourstartObj->epoch-$minstart)/60) + 1; 
      my $indexend = floor(($hourendObj->epoch-$minstart)/60) + 1;
      my $color = ($h % 2 == 0) ? "#4B98AD" : "#4BAD58"; #color hours alternately blue, green.
      my $hourstr = $daystr." ".$hourstartObj->hour.":00";  #since we started at midnight on the first day, the minutes is always 00.
      print(GFFH join("\t","1",$hourstr,"gene",$indexstart,$indexend,"0",".","0","ID=".$hourstr.";color=$color\n")); 
    }
  }
  close(GFFD);
  close(GFFH);
} # end writeGff

sub writeFasta {  #make "genome" fasta file for this job set
  open(GEN, ">$fastaFile") or die "Can't open file $fastaFile\n";
  print(GEN ">chr1\n");
  my $minutes = ($maxend - $minstart + 1)/60;  
  my $numlines = floor($minutes/60) + 1;  #60 nucleotides per line. Add 1 to round up.
  for (my $k=1; $k <= $numlines; $k++) {
    for (my $j=1; $j <= 60; $j++) {
      print(GEN "C");  #choice of base is irrelevant 
    }
    print(GEN "\n");
  }
  close(GEN);
} # end writeFasta


sub update_endpoints {
  my $minstart = shift;
  my $maxend = shift;
  my $start = shift;
  my $end = shift;

  if ($start < $minstart) {
    $minstart = $start;
  }
  if ($end > $maxend) {
    $maxend = $end;
  }
  return ($minstart, $maxend);
}

sub reqmem2num {
  my $reqmem = shift;
  my $nnodes = shift;
  my $ntasks = shift;

  if ($reqmem eq "") {
    $reqmem = -1;
  } elsif ($reqmem =~ m/^([.0-9]+)n$/) {
    $reqmem = $1;
  } elsif ($reqmem =~ m/^([.0-9]+)Kn$/) {
    $reqmem = $1*1024;
  } elsif ($reqmem =~ m/^([.0-9]+)Mn$/) {
    $reqmem = $1*1024**2;
  } elsif ($reqmem =~ m/^([.0-9]+)Gn$/) {
    $reqmem = $1*1024**3;
  } elsif ($reqmem =~ m/^([.0-9]+)c$/) {
    $reqmem = $1*$ntasks/$nnodes;
  } elsif ($reqmem =~ m/^([.0-9]+)Kc$/) {
    $reqmem = $1*1024*$ntasks/$nnodes;
  } elsif ($reqmem =~ m/^([.0-9]+)Mc$/) {
    $reqmem = $1*1024**2*$ntasks/$nnodes;
  } elsif ($reqmem =~ m/^([.0-9]+)Gc$/) {
    $reqmem = $1*1024**3*$ntasks/$nnodes;
  }
  return($reqmem);
}

sub bytes2num {  #better than using Number::Format
  my $mem = shift;

  if    ($mem eq "") { $mem = -1; }
  elsif ($mem =~ m/^([.0-9]+)KB?$/i) { $mem = $1*1024; }
  elsif ($mem =~ m/^([.0-9]+)MB?$/i) { $mem = $1*1024**2; }
  elsif ($mem =~ m/^([.0-9]+)GB?$/i) { $mem = $1*1024**3; }
  elsif ($mem =~ m/^([.0-9]+)TB?$/i) { $mem = $1*1024**4; }

  return($mem);
}

sub getMemUsageTag {
  my $RAMused = shift;
  my $RAMreserved = shift;

  my $percUsed = 0;
  my $usageTag = "";

  if ($RAMreserved == 0 || $RAMreserved == -1) {
    $usageTag = "No memory reservation reported.";
  } elsif ($RAMused == 0 || $RAMused == -1) {
    $usageTag = "No memory usage reported.";
  } else {
    my $delta = $RAMreserved - $RAMused;
    $percUsed = 100*$RAMused/$RAMreserved;
    if (abs($delta) > 1024**3 ) { #if at least 1GB of RAM was unused, or used w/o reservation...
      if ($percUsed > 200) {
        $usageTag = "Over 200% of memory reservation was used. (Usage was at least 1 GB over reservation)";
      } elsif ($percUsed > 100 && $percUsed <= 200) {
        $usageTag = "Over 100%, up to 200%, of memory reservation was used. (Usage was at least 1 GB over reservation)";
      } elsif ($percUsed > 80 ) {
        $usageTag = "Between 80 and 100 % of reserved memory was used. (More than 1GB went unused)";
      } elsif ($percUsed > 60 ) {
        $usageTag = "Between 60 and 80 % of reserved memory was used. (More than 1GB went unused)";
      } elsif ($percUsed > 40 ) {
        $usageTag = "Between 40 and 60 % of reserved memory was used. (More than 1GB went unused)";
      } elsif ($percUsed > 20 ) {
        $usageTag = "Between 20 and 40 % of reserved memory was used. (More than 1GB went unused)";
      } else {
        $usageTag = "Between 0 and 20 % of reserved memory was used. (More than 1GB went unused)";
      }
    } else {
      $usageTag = "Less than 1GB memory went unused, or was used without reservation.";
    }
  }
  return ($usageTag);
}

sub getCoresUsageTag {
  my $coresUsed = shift;
  my $coresAllocated = shift;
  my $usageTag = "";

  if ($coresAllocated > 0) {
    my $percUsed = 100*$coresUsed/$coresAllocated;
    if ($percUsed == 0) {
      $usageTag = "No CPU usage reported."
    } elsif ($percUsed < 20) {
      $usageTag = "Less than 20% of CPU reservation used.";
    } elsif ($percUsed < 40 ) {
      $usageTag = "Between 20 and 40 % of CPU reservation used.";
    } elsif ($percUsed < 60) {
      $usageTag = "Between 40 and 60 % of CPU reservation used.";
    } elsif ($percUsed < 80) {
      $usageTag = "Between 60 and 80 % of CPU reservation used.";
    } elsif ($percUsed <= 100) {
      $usageTag = "Between 80 and 100 % of CPU reservation used.";
    } elsif ($percUsed > 100 && $percUsed <= 200) {
      $usageTag = "Over 100%, up to 200%, of CPU reservation used.";
    } else {
      $usageTag = "Over 200% of CPU reservation used.";
    }
  } else {
    $usageTag = "No CPU allocation reported.";
  }

  return $usageTag;
}

sub calcCores {
  my $totalcpu = shift;
  my $cputime = shift;

  my $totalsec = stringToSeconds($totalcpu);
  my $cpusec = stringToSeconds($cputime);

  if ($cpusec == 0) { return 0; } 
  
  return $totalsec/$cpusec;
}

sub stringToSeconds {
  my $timeStr = shift;
  my $days = 0;
  my $hours = 0;
  my $minutes = 0;
  my $seconds = 0;

  $timeStr =~ m/^(([0-9]+)-)?([0-9]{1,2})?:([0-9]{1,2})?:([0-9]{1,2})$/;
  if ($2) { $days=$2; }
  if ($3) { $hours=$3; }
  if ($4) { $minutes=$4; }
  if ($5) { $seconds=$5; }

  return($days*24*60*60 + $hours*60*60 + $minutes*60 + $seconds);
}
