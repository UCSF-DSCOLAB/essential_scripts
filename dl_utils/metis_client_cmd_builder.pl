#! /usr/bin/perl

use strict;
use warnings;
use POSIX qw(strftime);

sub print_help {
    print <<'END_HELP';
Usage: perl metis_client_cmd_builder.pl [OPTIONS]

Options:
  --client-path=PATH    Specify the path to the metis_client
  --project=NAME        Name of the project to work on. E.g., "xhlt2"
  --file-ext=EXT        File extension(s) to use for matching files. Multiple file extensions can 
                        be separated by spaces (must enclose by quotes). E.g., ".fastq.gz .fq.gz"
  --cmd=COMMAND         Define the metis_client command to execute for each matching file
  --cmd-suffix=SUFFIX   Suffix to append to the command. E.g., "." to append the "put" command. The 
                        same values is applied to each command.
  --help                Show this help message and exit

Example:
  perl path/to/metis_client_cmd_builder.pl --client-path=/c4/home/$USER/metis/bin/metis_client --cmd=restrict --file-ext=".vcf.gz .vcg.gz.tbi" --project=demo
  perl path/to/metis_client_cmd_builder.pl --client-path=/c4/home/$USER/metis/bin/metis_client --cmd=get --file-ext=".fastq.gz .fq.gz" --project=demo --cmd-suffix='.'
END_HELP
    exit;
}

sub parse_args {
    my @args = @_;
    my %options;

    # Define allowed and required options
    my %allowed  = map { $_ => 1 } qw(client-path cmd cmd-suffix file-ext project);
    my %required = map { $_ => 1 } qw(client-path cmd file-ext project);

    foreach my $arg (@args) {
        if ($arg eq '--help' or $arg eq '-h') {
            print_help();
        }
        elsif ($arg =~ /^--(\w[\w-]*)=(.*)$/) {
            my ($key, $value) = ($1, $2);
            if ($allowed{$key}) {
                $options{$key} = $value;
            } else {
                warn "Unknown option ignored: --$key\n";
            }
        } else {
            warn "Invalid argument format (expected --key=value): $arg\n";
        }
    }

    # Check for missing required options
    my @missing = grep { !exists $options{$_} } keys %required;
    if (@missing) {
        die "Missing required option(s): " . join(', ', @missing) . "\nUse --help to see usage.\n";
    }

    return %options;
}


sub random_string {
    my $length = shift || 10;
    my @chars = ('A'..'Z', 'a'..'z', 0..9);
    my $str = '';
    $str .= $chars[int rand @chars] for 1..$length;
    return $str;
}

sub make_timestamp {
    return strftime("%Y-%m-%d_%H.%M.%S", localtime);
}

my $timestamp = make_timestamp();

#my $rand_str = random_string(10);

my %options=parse_args(@ARGV);
my $p = $options{'project'};
my $client_path = $options{'client-path'};
my $cmd = $options{'cmd'};
my $cmd_suffix = "";
$cmd_suffix = $options{'cmd-suffix'} if(exists $options{'cmd-suffix'});
my $file_ext = $options{'file-ext'};
my @file_exts = split(" ", $file_ext);

print "###############################################\n";
print "Input options:\n";
for my $i (keys %options) {
  print "  $i=$options{$i}\n";
}
my $metis_client_cmds_file = "$p-$timestamp.cmds";
print "Output:\n";
print "  metis_client command file=$metis_client_cmds_file\n";
print "###############################################\n";

# Check if the janus token is set
print "Checking Janus token and input parameters..... ";
my @token_check = `$client_path metis://$p $cmd 2>&1`;
if($token_check[0] =~ /No environment variable TOKEN is set/ or $token_check[0] =~ /No project is selected/ or $token_check[0] =~ /Your token is expired/ or $token_check[0] =~ /Invalid command/ or $token_check[0] =~ /No such file or directory/) {
  print "\nError:\n";
  print @token_check;
  exit 1;
}
print "Active Janus token is found and input parameters are valid..... Done.\n";

print "Finding files that match the input file extension(s)..... ";
open(O, ">$metis_client_cmds_file");
my @files = ();
for my $f_ext (@file_exts) {
  open(I, "$client_path metis://$p/ find name=~%$f_ext |");
  my @tmp_files = <I>;
  push @files, @tmp_files;
  close(I);
}
print "Found ", $#files+1, " files matching '$file_ext'..... Done.\n";

print "Building metis_client commands..... ";
my $dir = "";
for my $f (@files) {
  chomp $f;
  next if($f =~ /No results/);
  if($f =~ /:/){ 
    $f=~s/://g;
    $dir = $f;
  }
  else { 
    my $file="/$dir/$f";
    print O "$cmd '$file' $cmd_suffix\n";
  }
}
close(O);
print "Done.\n";
print "-----------------------------------------------\nRun the following command to run metis_client commands:\n";
print "  $client_path metis://$p < $p-$timestamp.cmds > $p-$timestamp.log 2> $p-$timestamp.err\n";


