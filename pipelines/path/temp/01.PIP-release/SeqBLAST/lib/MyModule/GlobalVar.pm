package MyModule::GlobalVar;

use warnings;
use strict;
use File::Spec;
use File::Basename;
use Cwd qw(abs_path);
use Config::IniFiles;
use Exporter;
use vars qw(@ISA @EXPORT);


@ISA = qw(Exporter);
@EXPORT = qw($PIPELINE_PATH $BIN_PATH $DB_PATH $CONFIG_FILE $PERL $WATCHDOG $SGE $MAX_AVAILEBLE_MEM $PERL 
             $PYTHON $YhrunBatch $MaxThread $NT $MAKEBLASTDB $BLASTN $ACC2TAX $ID2TAX $SEQID2SPEID
             $MKDIR $UNPIGZ $ZCAT $AWK $SED $LN $WHOAMI $CHMOD $SBATCH $GREP $SQUEUE $MV $PIGZ $BLASTDB);


# Path
my $modulePath = dirname(File::Spec->rel2abs(__FILE__));
our $PIPELINE_PATH  = abs_path("$modulePath/../..");
our $BIN_PATH       = "$PIPELINE_PATH/bin";
our $DB_PATH        = "$PIPELINE_PATH/db";
our $CONFIG_FILE    = "$PIPELINE_PATH/conf/configure.ini";


# Config file
&getCfg($CONFIG_FILE);


#===============================================================================
#   Subroutine
#===============================================================================
sub getCfg {
    my ($cfg_file, ) = @_;
    tie my %ini, "Config::IniFiles", (-file, => $cfg_file);

    our $WATCHDOG          = $ini{"software"}{"watchDog"};
    our $SGE               = $ini{"software"}{"SGE"};
    our $YhrunBatch        = $ini{"software"}{"YhrunBatch"};
    our $PERL              = $ini{"software"}{"perl"};
    our $PYTHON            = $ini{"software"}{"python"};
    our $MAKEBLASTDB       = $ini{"software"}{"makeblastdb"};
    our $BLASTN            = $ini{"software"}{"blastn"};
    our $MKDIR             = $ini{"software"}{"mkdir"};
    our $UNPIGZ            = $ini{"software"}{"unpigz"};
	our $PIGZ              = $ini{"software"}{"pigz"};
    our $ZCAT              = $ini{"software"}{"zcat"};
    our $AWK               = $ini{"software"}{"awk"};
    our $SED               = $ini{"software"}{"sed"};
    our $LN                = $ini{"software"}{"ln"};
    our $WHOAMI            = $ini{"software"}{"whoami"};
    our $CHMOD             = $ini{"software"}{"chmod"};
    our $SBATCH            = $ini{"software"}{"sbatch"};
    our $SQUEUE            = $ini{"software"}{"squeue"};
    our $GREP              = $ini{"software"}{"grep"};
    our $MV                = $ini{"software"}{"mv"};
    our $NT                = $ini{"database"}{"NT"};
    our $ACC2TAX           = $ini{"database"}{"acc2tax"};
    our $ID2TAX            = $ini{"database"}{"id2tax"};
    our $SEQID2SPEID       = $ini{"database"}{"seqID2speID"};
	our $BLASTDB           = $ini{"database"}{"BLASTDB"};
    our $MAX_AVAILEBLE_MEM = $ini{"resource"}{"MaxMem"};
    our $MaxThread         = $ini{"resource"}{"MaxThread"};
}


1;
