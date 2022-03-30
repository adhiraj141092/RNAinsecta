use strict;
use warnings;

while (<>) {
    if (/>\w+/) {
        
        print "$1"; # just the seq id
    }
    
}
