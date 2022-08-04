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
@EXPORT = qw($PIPELINE_PATH $BIN_PATH $DB_PATH $TOOL_PATH $CONFIG_FILE $DB_LIST_KRAKEN $DB_HIGHLIGHT 
             $REF_H_CONFIG $PYTHON3 $MAX_AVAILEBLE_MEM $REF_H_DBPATH_KRAKEN $REF_H_DBMEM_KRAKEN $PERL
             $MaxThread $BWA $PIGZ $SAMTOOLS $SRC $TXT2EXCEL $BLASTN $BLASTDB $KRAKENDB $ID2TAX $SEQID2SPEID $SEQTAXID2TAX);


# Path
my $modulePath = dirname(File::Spec->rel2abs(__FILE__));
our $PIPELINE_PATH  = abs_path("$modulePath/../../..");
our $BIN_PATH       = "$PIPELINE_PATH/bin";
our $DB_PATH        = "$PIPELINE_PATH/db";
our $SRC            = "$PIPELINE_PATH/src";
our $TOOL_PATH      = "$PIPELINE_PATH/tools";
our $CONFIG_FILE    = "$PIPELINE_PATH/conf/configure.ini";
#our $DB_LIST_KRAKEN = "$DB_PATH/db_Kraken/db.list";
#our $DB_HIGHLIGHT   = "$DB_PATH/db_Kraken/highlight.list";


# Config file
&getCfg($CONFIG_FILE);


# Kraken db info
#&getKrakenDbInfo($DB_LIST_KRAKEN);


#===============================================================================
#   Subroutine
#===============================================================================
sub getCfg {
    my ($cfg_file, ) = @_;
    tie my %ini, "Config::IniFiles", (-file, => $cfg_file);

	our $REF_H_CONFIG      = \%ini;
	our $PYTHON3           = $ini{"software"}{"python"};
	our $PERL              = $ini{"software"}{"perl"};
	our $BWA               = $ini{"software"}{"bwa"};
	our $PIGZ              = $ini{"software"}{"pigz"};
	our $SAMTOOLS          = $ini{"software"}{"samtools"};
	our $BLASTN            = $ini{"software"}{"blastn"};
	our $TXT2EXCEL         = $ini{"software"}{"txt2excel"};
	our $BLASTDB           = $ini{"database"}{"BLASTDB"};
	our $ID2TAX            = $ini{"database"}{"id2tax"};
	our $SEQID2SPEID       = $ini{"database"}{"seqID2speID"};
	#our $SEQTAXID2TAX      = $ini{"database"}{"seqTaxID2Tax"};
	our $KRAKENDB          = $ini{"database"}{"KrakenDB"};
	our $MAX_AVAILEBLE_MEM = $ini{"resource"}{"MaxMem"};
	our $MaxThread         = $ini{"resource"}{"MaxThread"};
}

#sub getKrakenDbInfo {
#    my ($dbList, ) = @_;
#    $dbList = abs_path($dbList);
#    my $dbDir = dirname($dbList);
    
#    my %h_dbPath = ();
#    my %h_dbMem = ();
#    open IN, "$dbList" or die $!;
#    while (<IN>) {
#        chomp;
        
#        my($dbName, $dbPath, $dbMem) = split /\s+/;
#        $dbPath = "$dbDir/$dbPath"; 
#        $h_dbPath{$dbName} = $dbPath;
#        $h_dbMem{$dbName}    = $dbMem;
#    }
#    close IN;
    
#    our $REF_H_DBPATH_KRAKEN = \%h_dbPath;
#    our $REF_H_DBMEM_KRAKEN  = \%h_dbMem;
#}


1;
