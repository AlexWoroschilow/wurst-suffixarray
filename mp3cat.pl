#!/usr/bin/perl

# use modules
use MP3::Tag;
use File::Find;

# set up counter to track number of files processed 
$count = 0;

# set up array of search directories
@dirs = (".", "/tmp", "/usr/local/music/", "/mnt/cdrom");

# print page and table header
print <<STARTHTML;
<html>
<head></head>
<body>
<table border=1 cellspacing=0 cellpadding=5> <tr> <td>Filename</td> <td>Artist</td> <td>Song</td> <td>Album</td> <td>Genre</td> <td>Year</td> </tr> 
STARTHTML

# look for files in each directory
find(\&displayMP3Info, @dirs);

# this function is called every time a file is found 
sub displayMP3Info {
	# if the file has an MP3 extension
	if (/\.mp3$/) 
	{
		# increment counter
		$count++;

		# create new MP3-Tag object
		$mp3 = MP3::Tag->new($_);

		# get tag information
		$mp3->get_tags();

		# check to see if an ID3v1 tag exists
		# if it does, print track information
		if (exists $mp3->{ID3v1})
		{
			print "<tr>\n";
			print "<td><a href=\"$File::Find::name\">$_</a></td>\n";
			print "<td>" . $mp3->{ID3v1}->artist . "</td>\n";
			print "<td>" . $mp3->{ID3v1}->title . "</td>\n";
			print "<td>" . $mp3->{ID3v1}->album . "</td>\n";
			print "<td>" . $mp3->{ID3v1}->genre . "</td>\n";
			print "<td>" . $mp3->{ID3v1}->year . "</td>\n";
			print "</tr>\n";
		}

		# clean up
		$mp3->close();
	}
}

