use strict;
use warnings;

my $project_dir = shift;

my @tmp=split /\//,$project_dir;

my $pro=$tmp[-1];

#my $back_dir="/RRM/RNA/TJPROJ/ncRNA/sRNA/"$pro;

print "rsync -I  -c -a --stats --progress $project_dir /RRM/RNA/ncRNA/$pro\n";
print "echo -e \"`du -sh $project_dir`\\t`du -sh /RRM/RNA/ncRNA/$pro`\" >> project_size.txt\n";
print "perl /PUBLIC/source/RNA/smallRNA/version3/dir_process/change_step2_3.pl /ifs/TJPROJ3/RNA/WORK/project_status/$pro*";
print "\n";
print "\n";
