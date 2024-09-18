#!/usr/bin/perl

#
# Rotational Translational Block method for large scale normal mode analysis
#
# Copyright (c) 2000-2005 Florence Tama, Yves-Henri Sanejouand
# Developed at Paul Sabatier University, Toulouse, France
#
# References:
# 
# RTB method
# 
# Durand, P, Trinquier, G & Sanejouand, YH. New Approach for Determining
# Low-Frequency Normal-Modes in Macromolecules. Biopolymers 34: 759-71
# (1994).
# 
# Tama, F, Gadea, FX, Marques, O & Sanejouand, YH. Building-block
# approach for determining low-frequency normal modes of
# macromolecules. Proteins 41: 1-7 (2000).
# 
# Elastic network model
# 
# Tirion, MM. Large amplitude elastic motions in proteins from a single-
# parameter, atomic analysis. Phys Rev Lett 77: 1905-8 (1996). 
#

# input
# 0 => pdb file

open( A, "$ARGV[0]" ) or die ("error :$!");
$icomp   = 0;
$i       = 0;
$allatom = 0;
$compres = 0;

while (<A>) {
    if (/ATOM/) {
        $i++;
        chomp($_);
        $key[$i]      = substr( $_, 0,  6 );
        $serial[$i]   = substr( $_, 6,  5 );
        $name[$i]     = substr( $_, 12, 4 );
        $altloc[$i]   = substr( $_, 16, 1 );
        $resnam[$i]   = substr( $_, 17, 3 );
        $chainid[$i]  = substr( $_, 21, 1 );
        $resseq[$i]   = substr( $_, 22, 4 );
        $icode[$i]    = substr( $_, 26, 1 );
        $x[$i]        = substr( $_, 30, 8 );
        $y[$i]        = substr( $_, 38, 8 );
        $z[$i]        = substr( $_, 46, 8 );
        $occup[$i]    = substr( $_, 54, 6 );
        $tempfact[$i] = substr( $_, 60, 6 );
        $segid[$i]    = substr( $_, 72, 4 );

        if ( $name[$i] ne ' CA ' ) {
            $allatom = 1;
        }
        if ( $resseq[$i] != $compres ) {
            $numbres++;
            $compres = $resseq[$i];
        }
    }
}

$natom = $i;
$i     = 0;
$nresb = 0;

### set number of residue per block ###

if ( $allatom == 0 ) {
    if ( $numbres <= 4000 ) {
        $nresb = 3;
    }
    if ( $numbres > 4000 and $numbres <= 8000 ) {
        $nresb = 5;
    }
    if ( $numbres > 8000 and $numbres <= 11000 ) {
        $nresb = 7;
    }
    if ( $numbres > 11000 ) {
        $nresb = 10;
    }
}

if ( $allatom == 1 ) {
    if ( $numbres <= 1000 ) {
        $nresb = 1;
    }
    if ( $numbres > 1000 and $numbres <= 4000 ) {
        $nresb = 3;
    }
    if ( $numbres > 4000 and $numbres <= 8000 ) {
        $nresb = 5;
    }
    if ( $numbres > 8000 and $numbres <= 11000 ) {
        $nresb = 7;
    }
    if ( $numbres > 11000 ) {
        $nresb = 10;
    }
}

warn " $numbres residues, $nresb residues per blocks\n";

######

$chain      = $chainid[1];
$chainid[0] = $chain;
$segib      = 1;
$segid[0]   = $segib;
$comp       = 0;
$resseq[0]  = $resseq[1] - 1;
$test       = $resseq[1] + $nresb;

### If CA only ########
if ( $allatom == 0 ) {
    warn "CA atoms only\n";

    for ( $i = 1 ; $i <= $natom ; $i++ ) {

        $segid[$i] = $segid[ $i - 1 ];
        $diff = $resseq[$i] - $resseq[ $i - 1 ];
        if ( $chainid[$i] ne $chain ) {
            $comp = $segid[ $i - 1 ] - $segid[ $i - 3 ];
            if ( $comp != 0 ) {
                $segid[$i]       = $segid[ $i - 1 ];
                $segid[ $i - 1 ] = $segid[ $i - 3 ];
                $segid[ $i - 2 ] = $segid[ $i - 3 ];
                $test            = $resseq[$i] + $nresb;
                $chain           = $chainid[$i];
            }
            else {
                $segib++;
                $segid[$i] = $segib;
                $chain     = $chainid[$i];
                $test      = $resseq[$i] + $nresb;
            }
        }
        elsif ( $diff != 1 ) {
            if ( $chainid[$i] eq $chainid[ $i - 1 ] ) {
                $comp = $segid[ $i - 1 ] - $segid[ $i - 3 ];
                if ( $comp != 0 ) {
                    $segid[$i]       = $segid[ $i - 1 ];
                    $segid[ $i - 1 ] = $segid[ $i - 3 ];
                    $segid[ $i - 2 ] = $segid[ $i - 3 ];
                    $test            = $resseq[$i] + $nresb;
                }
                else {
                    $segib++;
                    $segid[$i] = $segib;
                    $test = $resseq[$i] + $nresb;
                }
            }
        }
        elsif ( $resseq[$i] == $test ) {
            $segib++;
            $segid[$i] = $segib;
            $test = $resseq[$i] + $nresb;
        }
    }

    # check last block is at least made of 3 atoms

    $difflast = $segid[$natom] - $segid[ $natom - 2 ];
    if ( $difflast != 0 ) {
        $segid[$natom] = $segid[ $natom - 2 ];
        $segid[ $natom - 1 ] = $segid[ $natom - 2 ];
    }
    $chainid[$natom] = ' ';
}
### END if CA only

### If all atoms ###

if ( $allatom == 1 ) {
    warn "all atoms\n";
    for ( $i = 1 ; $i <= $natom ; $i++ ) {
        $segid[$i] = $segib;
        if ( $resseq[$i] != $resseq[ $i - 1 ] ) {
            $diffres = $resseq[$i] - $resseq[ $i - 1 ];
            $diff    = $resseq[$i] - $test;
            if ( $diffres < 0 ) {
                $segib++;
                $segid[$i] = $segib;
                $test = $resseq[$i] + $nresb;
            }
            if ( $diffres >= 2 ) {
                $segib++;
                $segid[$i] = $segib;
                $test = $resseq[$i] + $nresb;
            }
            elsif ( $diff == 0 ) {
                $segib++;
                $segid[$i] = $segib;
                $test = $resseq[$i] + $nresb;
            }
        }
    }
    $chainid[$natom] = ' ';
}

for ( $i = 1 ; $i <= $natom ; $i++ ) {
    printf(
        "ATOM  %5d %4s%1s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n",
        $serial[$i],  $name[$i],   $altloc[$i], $resnam[$i],
        $chainid[$i], $resseq[$i], $icode[$i],  $x[$i],
        $y[$i],       $z[$i],      $occup[$i],  $tempfact[$i],
        $segid[$i]
    );
}
