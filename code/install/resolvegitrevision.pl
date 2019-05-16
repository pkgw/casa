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

$casaBranchHint = $ENV{'CASABRANCHHINT'};
$tagMatcher="-mas-";
$casaVariant="";
#$tagid="mas";

$debug=$ARGV[0];
print "tagMatcher: $tagMatcher\n" if $debug;
#print "tagId: $tagid\n" if $debug;

# Get the hint from an environment variable. This is used for detached head builds
# Default grep is master

# The "default" Casa release that is not tied to any instrument
$defaultReleaseMatch = '.*release/\d+\.\d+\.\d';
$instrumentReleaseMatch = 'release/.*/\d+\.\d+\.\d';

if ($casaBranchHint)  {
    print "Branch hint provided: $casaBranchHint\n" if $debug;
    if ($casaBranchHint =~ "^(feature|bugfix).*CAS.*") {
        @splat = split '/', $casaBranchHint;
        print "Splat: @splat\n" if $debug;
        $tagMatcher = "@splat[0]-@splat[1]";
    } elsif ($casaBranchHint =~ "^CAS-.*" ) {
        $tagMatcher = $casaBranchHint;
    } elsif ($casaBranchHint =~ m/$defaultReleaseMatch/) {
        @splat = split '/', $casaBranchHint;
        $tagMatcher = "@splat[1]-rel-";
    } elsif ($casaBranchHint =~ m/$instrumentReleaseMatch/) {
        $tagMatcher = $casaBranchHint;
        $tagMatcher =~s/\//-/g;
        #$tagid="rel";
    } elsif ($casaBranchHint =~ "^bambooprtest.*" ) {
        $tagMatcher = "bambooprtest";
    }
} else {
    print "No branch hint available $casaBranchHint\n" if $debug;
}
print "tagMatcher: $tagMatcher\n" if $debug;

# Check where the current "HEAD" points to.
$branch = `git rev-parse --abbrev-ref HEAD`;
chomp($branch);
print "Current branch: $branch\n" if $debug;

$branchTag = "";
# Figure version for detached HEAD
if ( $branch eq "HEAD") {
  print "Detached head.\n" if $debug;
  $branchTag = `git tag --points-at HEAD | grep -- $tagMatcher`;
  print "Head tag: $branchTag\n" if $debug;
}
else {
     if ($branch eq "master") { $tagMatcher = "-mas-"; }
     elsif ($branch =~ m/$defaultReleaseMatch/) { $tagMatcher = "-rel-"; }
     else {
         $tagMatcher = $branch;
     }

     if ($tagMatcher =~ "/") {
         $tagMatcher =~ s/\//-/g;
     }
     print "tagMatcher: $tagMatcher\n" if $debug;
     $branchTag=`git tag --points-at HEAD | grep -- $tagMatcher`;
     $headCommit=`git rev-parse HEAD`;
     $casaVersionDesc="";
}
##

$needsId=0;
if (!$branchTag) {
     $headCommit=`git rev-parse HEAD`;
     #$branchTag=`git describe --abbrev=0 --match='[0-9]*.[0-9]*.[0-9]*-$tagMatcher-[0-9]*'`;
     $needsId=1;
}
chomp ($branchTag);
$casaVersionDesc = "ID $headCommit";
chomp($casaVersionDesc);
print "needsId: $needsId\n" if $debug;
print "branchTag: $branchTag\n" if $debug;
print "casaVersionDesc: $casaVersionDesc\n" if $debug;
@splat = split '-', $branchTag;
if ($needsId) {
    # if ($tagMatcher eq "-rel-" || $tagMatcher eq "-mas-" || $tagMatcher =~ /^release-/) {
    #     if ($tagMatcher =~ /^release-/) {
    #         @mComp = split '-', $tagMatcher;
    #         $casaVariant = uc $mComp[1];
    #     }
    #     print "$splat[-1];$casaVersionDesc;$casaVariant\n"
    # }
    # else {
        print "0;$casaVersionDesc;$casaVariant\n"
    #}
}
else {
    if ($tagMatcher =~ "-rel-" || $tagMatcher eq "-mas-" || $tagMatcher =~ /^release-/) {
        if ($tagMatcher =~ /^release-/) {
            @mComp = split '-', $tagMatcher;
            $casaVariant = uc $mComp[1];
        }
        print "$splat[-1];;$casaVariant\n"
    } else {
        print "0;$branchTag;$casaVariant\n"
    }
}
