use strict;
use warnings;

#config
my $IFS=",";
my $removeINDELS=0;

my $usage="$0 file1 file2\n";

if (scalar @ARGV != 2 || ! -f $ARGV[0] || ! -f $ARGV[1])
{
	die "ERROR\n$usage";
}

my $file1=$ARGV[0];
my $file2=$ARGV[1];

mopen(my $FH1,$file1);
my @content1=<$FH1>;
close($FH1);

mopen(my $FH2,$file2);
my @content2=<$FH2>;
close($FH2);

shift @content1;
shift @content2;

my %vars1=map{join("\t",(split($IFS,$_))[0,1])=>[(split($IFS,$_))[3,4]]} @content1;
my %vars2=map{join("\t",(split($IFS,$_))[0,1])=>[(split($IFS,$_))[3,4]]} @content2;
my %commonvars;

if($removeINDELS==1)
{
	removeINDELS(\%commonvars);
	removeINDELS(\%vars1);
}

foreach my $var (keys %vars1)
{
	if (exists $vars2{$var})
	{
		$commonvars{$var}=$vars2{$var};
		delete $vars1{$var};
		delete $vars2{$var};
	}
}


my $nprivate1=scalar keys %vars1;
my $nprivate2=scalar keys %vars2;
my $ncommon=scalar keys %commonvars;

print("Vars1 in 2:".$ncommon/($nprivate1+$ncommon).", vars2 in 1:".$ncommon/($nprivate2+$ncommon).", private1: $nprivate1, private2: $nprivate2, common: $ncommon\n");

##Quick and dirty solution

foreach my $var (keys %commonvars)
{
	print("$var\n");
}
exit;

###DEBUG
#foreach my $var (keys %vars1)
#{
#	print("$var\n");
#}
#foreach my $var (keys %vars2)
#{
#	print("$var\n");
#}


# Open with error message
###################################################################################
sub mopen
{
    open($_[0],$_[1]) or die "ERROR: impossible to open the file $_[1] in ".($_[1]=~m/^>/?"write":"read")."mode\n";
    return $_[0];
}

# Remove indels
sub removeINDELS
{
	my ($refhash)=@_;
	foreach my $key (keys %$refhash)
	{
		if(length @{$refhash->{$key}}[0] != length @{$refhash->{$key}}[1])
		{
			delete $refhash->{$key};
		}
	}
}
