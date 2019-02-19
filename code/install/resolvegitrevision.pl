#!/usr/bin/env perl

# Note:
#branch="release/5.1.0"
#branch="release/vla/1.0.0"
#println branch

#if (branch =~ ".*release/[A-Za-z].*") {
#    println "A: " + branch
#}

#if (branch =~ /.*release\/\d+\.\d+\.\d+/) {
#    println "B: " + branch
#}

$casaBranchHint=$ARGV[0];
$headGrep="-mas-";
$tagid="mas";

$debug=0;
print "headGrep: $headGrep\n" if $debug;
print "tagId: $tagid\n" if $debug;

# Get the hint from an environment variable. This is used for detached head builds
# Default grep is master


if ($casaBranchHint)  {
    print "Branch hint provided: $casaBranchHint\n" if $debug;
    if ($casaBranchHint =~ "^(feature|bugfix).*CAS.*") {
        @splat = split '/', $casaBranchHint;
        print "Splat: @splat\n" if $debug;
        $headGrep = "@splat[0]-@splat[1]";
    } elsif ($casaBranchHint =~ "^CAS-.*" ) {
        $headGrep = $casaBranchHint;
    } elsif ($casaBranchHint =~ ".*release/\d+\.\d+\.\d+") {
        $headGrep = $casaBranchHint;
        $tagid="rel";
    } elsif ($casaBranchHint =~ "^bambooprtest.*" ) {
        $headGrep = "bambooprtest";
    }
} else {
    print "No branch hint available $casaBranchHint\n" if $debug;

}
print "headGrep: $headGrep\n" if $debug;

# Check where the current "HEAD" points to.
$branch = `git rev-parse --abbrev-ref HEAD`;
chomp($branch);
print "Current branch: $branch\n" if $debug;

# Figure version for detached HEAD
if ( $branch eq "HEAD") {
  print "Detached head.\n" if $debug;
  $headTag = `git tag --points-at HEAD | grep -- $headGrep`;
  print "Head tag: $headTag\n" if $debug;
  # I don't understand what the purpose of the commented out section was
  #if (headTag) {
      #      if [[ -z "${headTag// }" ]]; then
      #         # Get the nearest tag and add Desc
      #         headCommit=`git rev-parse HEAD`
      #         headTag=`git tag --points-at HEAD | grep $headGrep | xargs`
      #     fi
  #}
  $casaVersionDesc = $headTag;
  print("Version Description: $casaVersionDesc\n") if $debug;
  # $CASAFORKPOINTHINT is the fork point commit
  # You can obtain this by executing  "git merge-base --fork-point master"
  # while in the branch, but before detaching the HEAD
  $forkpoint = `git merge-base master $branch`;
  print("Default Fork Point: $forkpoint\n") if $debug;
  if ($casaBranchHint =~ ".*release.*") {
        $forkpoint=`git merge-base $casaBranchHint $branch`;
  }
  print("Fork Point: $forkpoint\n") if $debug;
  #
  $headTag=`git describe --abbrev=0 --tags --match='[0-9]*.[0-9]*.[0-9]*-mas-[0-9]*' \$(git rev-parse $forkpoint)`;
  if ($tagid eq "rel") {
     $headTag=`git describe --abbrev=0 --tags --match='[0-9]*.[0-9]*.[0-9]*-rel-[0-9]*' \$(git rev-parse $forkpoint)`;
  }
  chomp($headTag);
  print("headTag: $headTag\n") if $debug;

  @splat = split '-',$headTag;
  if ($casaVersionDesc =~".*-$tagid-.*") {
      print "$splat[-1]\n";
  }
  else {
      print "$splat[-1];$casaVersionDesc\n";
  }
}
# Version from branch
else {
     if ($branch eq "master") { $tagMatcher = "mas"; }
     elsif ($branch =~ ".*release/\d+\.\d+\.\d+") { $tagMatcher = "rel"; }
     else {
         $tagMatcher = $branch;
     }

     if ($tagMatcher =~ "/") {
         $tagMatcher =~ s/\//-/;
     }
     print "tagMatcher: $tagMatcher\n" if $debug;
     $branchTag=`git tag --points-at HEAD | grep -- -$tagMatcher-`;
     $headCommit=`git rev-parse HEAD`;
     $casaVersionDesc="";
     $needsId=0;
     if (!$branchTag) {
         $headCommit=`git rev-parse HEAD`;
         #$branchTag=`git describe --abbrev=0 --match='[0-9]*.[0-9]*.[0-9]*-$tagMatcher-[0-9]*'`;
         $needsId=1
     }
     chomp ($branchTag);
     $casaVersionDesc = "ID $headCommit";
     chomp($casaVersionDesc);
     print "needsId: $needsId\n" if $debug;
     print "branchTag: $branchTag\n" if $debug;
     print "casaVersionDesc: $casaVersionDesc\n" if $debug;
     @splat = split '-', $branchTag;
     if ($needsId) {
         if ($tagMatcher eq "rel" || $tagMatcher eq "mas") {
             print "$splat[-1];$casaVersionDesc\n"
         }
         else {
             print ";$casaVersionDesc\n"
         }
     }
     else {
         if ($tagMatcher eq "rel" || $tagMatcher eq "mas") {
             print "$splat[-1];\n"
         } else {
             print ";$branchTag\n"
         }
     }
}
