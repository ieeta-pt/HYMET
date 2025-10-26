#!/usr/bin/perl
use strict;
use warnings;

# Assume the current directory is the root of the tool (HYMET)
my $base_path = '.';  # Current directory (HYMET)

sub ensure_tool {
    my ($tool) = @_;
    my $path = `command -v $tool 2>/dev/null`;
    chomp $path;
    die "Error: required tool '$tool' not found in PATH.\n" unless $path;
}

sub run_or_die {
    my ($cmd_ref, $message) = @_;
    my $rc = system(@{$cmd_ref});
    if ($rc != 0) {
        my $exit = $rc >> 8;
        die "$message (exit code $exit).\n";
    }
}

ensure_tool($_) for qw(wget unzip python3);

# Define the necessary directories
my $taxonomy_files_dir = "$base_path/taxonomy_files";
my $scripts_dir = "$base_path/scripts";

# Create directories if they don't exist
mkdir $taxonomy_files_dir unless -d $taxonomy_files_dir;
mkdir $scripts_dir unless -d $scripts_dir;

# URL for taxonomy files
my $taxdmp_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip";

# Download taxonomy files
print "Downloading taxonomy files...\n";
my $taxdmp_zip = "$taxonomy_files_dir/taxdmp.zip";
run_or_die([
    'wget', '-q', '-O', $taxdmp_zip, $taxdmp_url
], "Failed to download taxonomy archive from $taxdmp_url");

# Unzip the downloaded file
print "Unzipping taxonomy files...\n";
run_or_die([
    'unzip', '-qo', $taxdmp_zip, '-d', $taxonomy_files_dir
], "Failed to extract taxonomy archive to $taxonomy_files_dir");

# Execute the Python script
print "Executing Python script...\n";
run_or_die([
    'python3', "$scripts_dir/taxonomy_hierarchy.py"
], "Failed to generate taxonomy hierarchy");

print "Configuration completed.\n";
