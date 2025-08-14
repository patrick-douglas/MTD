#!/usr/bin/env perl

# Copyright 2013-2021, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Reads an assembly_summary.txt file, which indicates taxids and FTP paths for
# genome/protein data.  Performs the download of the complete genomes from
# that file, decompresses, and explicitly assigns taxonomy as needed.

use strict;
use warnings;
use File::Basename;
use File::Copy;
use Getopt::Std;
use List::Util qw/max/;

my $PROG = basename $0;

# Specify the original directory where your downloaded .gz files are stored
my $local_download_dir = "/media/me/4TB_BACKUP_LBN/Compressed/MTD/Kraken2DB_micro/library/bacteria/all/";

# Specify the new directory where you want to copy the files
my $new_download_dir = "$ENV{HOME}/MTD/kraken2DB_micro/library/bacteria/all";

# Create the new directory if it doesn't exist
mkdir $new_download_dir unless -d $new_download_dir;

# Copy all files from the original directory to the new directory
opendir(my $dh, $local_download_dir) || die "Can't opendir $local_download_dir: $!";
while (my $file = readdir($dh)) {
    next if $file =~ /^\./;  # skip hidden files/directories
    my $source = "$local_download_dir/$file";
    my $destination = "$new_download_dir/$file";
    copy($source, $destination) or die "Copy failed: $!";
}
closedir($dh);

# Update $local_download_dir to the new directory
$local_download_dir = $new_download_dir;

# Manifest hash maps filenames (keys) to taxids (values)
my %manifest;
while (<>) {
    next if /^#/;
    chomp;
    my @fields = split /\t/;
    my ($taxid, $asm_level, $ftp_path) = @fields[5, 11, 19];
    # Possible TODO - make the list here configurable by user-supplied flags
    next unless grep {$asm_level eq $_} ("Complete Genome", "Chromosome");
    next if $ftp_path eq "na";  # Skip if no provided path

    my $filename = basename($ftp_path);
    my $local_path = "$local_download_dir/$filename.gz";
    
    if (-e $local_path) {
        $manifest{$local_path} = $taxid;
    } else {
        my $alt_local_path = "$local_download_dir/${filename}_genomic.fna.gz";
        if (-e $alt_local_path) {
            $manifest{$alt_local_path} = $taxid;
        } else {
            print STDERR "$PROG: Local file $local_path or $alt_local_path not found. Skipping.\n";
        }
    }
}

open MANIFEST, ">", "manifest.txt"
    or die "$PROG: can't write manifest: $!\n";
print MANIFEST "$_\n" for keys %manifest;
close MANIFEST;

print STDERR "Manifest verification complete:\n";
print STDERR "  Found: " . scalar(keys %manifest) . "\n";
#print STDERR "  Not found: " . (58075 - scalar(keys %manifest)) . "\n";
#print STDERR "  Replaced: 0\n";

print STDERR "Step 1/2: Processing locally downloaded files\n";
my $output_file = "library.fna";
open OUT, ">", $output_file
    or die "$PROG: can't write $output_file: $!\n";
my $projects_added = 0;
my $sequences_added = 0;
my $ch_added = 0;
my $ch = "bp";
my $max_out_chars = 0;
for my $in_filename (keys %manifest) {
    my $taxid = $manifest{$in_filename};
    open IN, "gunzip -c $in_filename |" or die "$PROG: can't read $in_filename: $!\n";
    while (<IN>) {
        if (/^>/) {
            s/^>/>kraken:taxid|$taxid|/;
            $sequences_added++;
        } else {
            $ch_added += length($_) - 1;
        }
        print OUT;
    }
    close IN;
    # Remove the unlink command to preserve the original files
    # unlink $in_filename;
    $projects_added++;
    my $out_line = progress_line($projects_added, scalar keys %manifest, $sequences_added, $ch_added) . "...";
    $max_out_chars = max(length($out_line), $max_out_chars);
    my $space_line = " " x $max_out_chars;
    print STDERR "\r$space_line\r$out_line" if -t STDERR;
}
close OUT;
print STDERR " done.\n" if -t STDERR;

print STDERR "All files processed, cleaning up...\n";
print STDERR " done, library complete.\n";

sub progress_line {
    my ($projs, $total_projs, $seqs, $chs) = @_;
    my $line = "Processed ";
    $line .= ($projs == $total_projs) ? "$projs" : "$projs/$total_projs";
    $line .= " project" . ($total_projs > 1 ? 's' : '') . " ";
    $line .= "($seqs sequence" . ($seqs > 1 ? 's' : '') . ", ";
    my $prefix;
    my @prefixes = qw/k M G T P E/;
    while (@prefixes && $chs >= 1000) {
        $prefix = shift @prefixes;
        $chs /= 1000;
    }
    if (defined $prefix) {
        $line .= sprintf '%.2f %s%s)', $chs, $prefix, $ch;
    } else {
        $line .= "$chs $ch)";
    }
    return substr($line, 0, 79);
}

