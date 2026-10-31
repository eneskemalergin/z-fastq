#!/usr/bin/env bash
# Check the Markdown git tracks or would track: the wiki's layout, and that every link, image, and heading anchor it names exists.
#
#   docs_check.sh
#
# Wiki pages are top-level wiki/*.md files with Home, _Sidebar, and _Footer present; wiki images live in
# wiki/assets/. A wiki link names a page without .md; a repository link names a tracked path relative to its file.
# Links to this repository's wiki, blob, and raw URLs are checked the same way, so the README's wiki links and the
# wiki's benchmark figures cannot point at a page, heading, or file that does not exist. Text inside code is skipped.
# External sites are not fetched.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

repository=${GITHUB_REPOSITORY:-eneskemalergin/z-fastq}

git ls-files -z --cached --others --exclude-standard | perl -0 /dev/fd/3 "$repository" 3<<'PERL'
use strict;
use warnings;
use utf8;
use File::Basename qw(dirname);
use File::Spec::Functions qw(canonpath);

binmode(STDERR, ':encoding(UTF-8)');
my ($repository) = @ARGV;
my $repo_url = quotemeta("https://github.com/$repository");
my $raw_url = quotemeta("https://raw.githubusercontent.com/$repository");
my @files = map { s/\0\z//r } <STDIN>;
my %tracked = map { $_ => 1 } @files;
my %dirs;
for my $file (@files) {
    my $dir = $file;
    $dirs{$dir} = 1 while ($dir =~ s{/[^/]*\z}{});
}

my $failed = 0;
sub fail { print STDERR "$_[0]\n"; $failed = 1; }

my %pages;
for my $file (grep { m{^wiki/} } @files) {
    my $name = substr($file, length 'wiki/');
    next if $name =~ m{^assets/[^/]+\z};
    if ($name !~ m{^[^/]+\.md\z}) {
        fail("wiki files are top-level Markdown pages or wiki/assets/ images: $file");
        next;
    }
    fail("wiki page name has a character a wiki cannot hold: $file") if $name =~ m{[\\:*?"<>|]};
    fail("wiki page is empty: $file") if !-s $file;
    (my $page = $name) =~ s/\.md\z//;
    if (exists $pages{lc $page}) {
        fail("wiki pages differ only in case: $pages{lc $page} and $file");
    } else {
        $pages{lc $page} = $file;
    }
}
for my $required (qw(Home _Sidebar _Footer)) {
    fail("required wiki page is missing: wiki/$required.md") if !exists $pages{lc $required};
}

my (%text, %anchors);
sub markdown {
    my ($file) = @_;
    return $text{$file} if exists $text{$file};
    open(my $fh, '<:encoding(UTF-8)', $file) or die "cannot read $file: $!\n";
    local $/ = "\n";
    my ($fence, @lines) = (undef);
    while (my $line = <$fh>) {
        if (defined $fence) {
            $fence = undef if $line =~ /^\s*\Q$fence\E\s*$/;
            push @lines, "\n";
            next;
        }
        if ($line =~ /^\s*(`{3,}|~{3,})/) {
            $fence = $1;
            push @lines, "\n";
            next;
        }
        push @lines, $line;
    }
    close($fh);
    my %seen;
    my %found;
    for my $line (@lines) {
        $found{$1} = 1 while $line =~ /<a\s+(?:name|id)="([^"]+)"/g;
        next if $line !~ /^#{1,6}\s+(.*?)\s*#*\s*$/;
        my $slug = $1;
        $slug =~ s/!?\[([^\]]*)\]\([^)]*\)/$1/g;
        $slug =~ s/<[^>]+>|`//g;
        $slug = lc $slug;
        $slug =~ s/[^\p{L}\p{N}\s_-]//g;
        $slug =~ s/\s/-/g;
        my $n = $seen{$slug}++;
        $found{$n ? "$slug-$n" : $slug} = 1;
    }
    $anchors{$file} = \%found;
    return $text{$file} = join('', map { s/(`+).*?\1//gr } @lines);
}

sub check_anchor {
    my ($from, $file, $anchor) = @_;
    return if !defined $anchor || $anchor eq '' || $file !~ /\.md\z/;
    markdown($file);
    fail("$from links to a missing heading: $file#$anchor") if !$anchors{$file}{lc $anchor};
}

sub check_wiki_page {
    my ($from, $target) = @_;
    my ($page, $anchor) = split(/#/, $target, 2);
    $page =~ s{^\./}{};
    return check_anchor($from, $from, $anchor) if $page eq '';
    fail("$from links to a wiki page with .md: $target") if $page =~ /\.md\z/;
    $page =~ s/\.md\z//;
    my $file = $pages{lc $page};
    return fail("$from links to a missing wiki page: $target") if !defined $file;
    check_anchor($from, $file, $anchor);
}

sub check_path {
    my ($from, $path, $anchor) = @_;
    $path =~ s{/\z}{};
    return fail("$from links to an untracked path: $path") if !$tracked{$path} && !$dirs{$path};
    check_anchor($from, $path, $anchor);
}

for my $file (grep { /\.md\z/ } @files) {
    my $text = markdown($file);
    my $in_wiki = $file =~ m{^wiki/};
    my @targets;
    push @targets, $1 while $text =~ /\]\(\s*<?([^)\s>]+)>?(?:\s+"[^"]*")?\s*\)/g;
    push @targets, $1 while $text =~ /\b(?:src|href)="([^"]+)"/g;
    while ($text =~ /\bsrcset="([^"]+)"/g) {
        push @targets, map { (split ' ', $_)[0] } split(/,/, $1);
    }
    for my $target (@targets) {
        $target =~ s/&amp;/&/g;
        if ($target =~ m{^$repo_url/wiki(?:/([^?#]*))?(?:#(.*))?\z}) {
            check_wiki_page($file, (defined $1 && $1 ne '' ? $1 : 'Home') . (defined $2 ? "#$2" : ''));
        } elsif ($target =~ m{^$repo_url/(?:blob|tree)/main/([^?#]+)(?:#(.*))?\z} ||
            $target =~ m{^$raw_url/main/([^?#]+)\z}) {
            check_path($file, $1, $2);
        } elsif ($target =~ m{^[a-z][a-z0-9+.-]*:}i) {
            next;
        } elsif ($in_wiki && $target !~ m{^assets/}) {
            check_wiki_page($file, $target);
        } else {
            my ($path, $anchor) = split(/#/, $target, 2);
            if ($path eq '') {
                check_anchor($file, $file, $anchor);
                next;
            }
            $path = canonpath(dirname($file) . "/$path");
            1 while $path =~ s{(?:^|/)(?!\.\.)[^/]+/\.\.(?=/|\z)}{};
            $path =~ s{^/}{};
            check_path($file, $path, $anchor);
        }
    }
}

exit 1 if $failed;
printf "checked links in %d Markdown files and %d wiki pages\n", scalar(grep { /\.md\z/ } @files), scalar(keys %pages);
PERL
