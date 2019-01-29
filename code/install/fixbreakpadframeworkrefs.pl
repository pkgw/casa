use Cwd;
use Cwd 'abs_path';
$builddir=$ARGV[0];
@files = `find $builddir/breakpad/breakpad-distro/src/client/mac/build/Release/ -type f`;
#print @files;
# Fix library references
foreach $file ( @files ) {
      $type = `file $file`;
      if ($type =~ m|Mach|) {
        print "$file\n";
        @references =`otool -L $file`;
        foreach $reference ( @references ) {
          print "  $reference\n";
          if ($reference =~ m|\@executable_path/../Frameworks|) {
            print "Changing \@executable_path to \@loader_path for $file\n";
            $estring = "\@executable_path/../Frameworks/";
            $lstring = "$builddir/breakpad/breakpad-distro/src/client/mac/build/Release/";
            my ($oldref, $reftail) = split / /, $reference, 2;
            $newref = $oldref;
            $newref =~ s/$estring/$lstring/g;
            print "Old ref: $oldref\n";
            print "New ref: $newref\n";
            system( "/usr/bin/install_name_tool -change $oldref $newref $file " ) == 0 
                or die "/usr/bin/install_name_tool failed to update rpath. Fix your configuration or disable Breakpad build";
          }
          elsif ($reference =~ m|\@executable_path/../Resources|) {
            print "Changing \@executable_path to \@loader_path for $file\n";
            $estring = "\@executable_path/../Resources/";
            $lstring = "$builddir/breakpad/breakpad-distro/src/client/mac/build/Release/";
            my ($oldref, $reftail) = split / /, $reference, 2;
            $newref = $oldref;
            $newref =~ s/$estring/$lstring/g;
            print "Old ref: $oldref\n";
            print "New ref: $newref\n";
            system( "/usr/bin/install_name_tool -change $oldref $newref $file " ) == 0
                or die "install_name_tool failed to update rpath. Fix your configuration or disable Breakpad build";
          }
        }
      }
}

# Fix id
foreach $file ( @files ) {
    $type = `file $file`;
    if ($type =~ m|Mach|) {
        print "$file\n";
        @references =`otool -L $file`;
        foreach $reference ( @references ) {
            if ($reference =~ m|\@executable_path/../|) {
            print "Changing id for $file\n";
            $estring = "\@executable_path/../Frameworks/";
            $lstring = "$builddir/breakpad/breakpad-distro/src/client/mac/build/Release/";
            my ($oldref, $reftail) = split / /, $reference, 2;
            $newref = $oldref;
            $newref =~ s/$estring/$lstring/g;
            print "Old ref: $oldref\n";
            print "New ref: $newref\n";
            print "File: $file\n";
            `/usr/bin/install_name_tool -id $newref $file `;
          }
        }
    }
}
