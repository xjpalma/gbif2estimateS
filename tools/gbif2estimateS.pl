#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Text::CSV_XS;
use Geo::Hash;
use Geo::Hash::XS;
use Geo::Hash::Grid;
use Bio::Community;
use Bio::Community::Member;
use Bio::Community::Meta;
use Bio::Community::Tools::Accumulator;
use Data::Dumper;

$| = 1;

my ($file_input, $file_output, $help, $verbose, $mverbose, $ini_year, $end_year, $precision, $separator, $species);

my $start_time = time();

######## options/help ########
GetOptions(
    'help|h'                    => \$help,
    'verbose|v'                 => \$verbose,
    'mverbose|vv'               => \$mverbose,
    'input|i=s'                 => \$file_input,
    'output|o=s'                => \$file_output,
    'iniyear|iy=i'              => \$ini_year,
    'endyear|ey=i'              => \$end_year,
    'level|l=i'                 => \$precision,
    'sep|s=s'                   => \$separator,
    "species|n=s"               => \$species
) or die "\n Invalid commmand line options.\n";

my $USAGE =<<USAGE;
    Usage:
         $0 [-h|--help]
            [-v|--verbose]
            [-vv|--mverbose]
            [-i|--input]
            [-o|--output]
            [-iy|--iniyear]
            [-ey|--endyear]
            [-l|--level]
            [-s]|--sep]
            [-n]|--species]

    Options:
          -h,  --help          show this help
          -v,  --verbose       verbose mode
          -vv, --mverbose      more verbose mode
          -i,  --input         input file
          -o,  --output        output file
          -iy, --iniyear       start year
          -ey, --endyear       end year
          -l,  --level         level precision [1-12]
          -s,  --sep           csv separator (default is TAB)
          -n,  --species       species list allowed Ex: -n 2356,4586,7456,...

USAGE

## check arguments
if($help) {
    print "$USAGE\n";
    exit 1;
}
if(!$file_input){
    print "  $0: Argument required. \n";
    print "$USAGE\n";
    exit 1
}
if(!$file_output){
    $file_output = 'output.dat';
}
if(!$ini_year){
    $ini_year = 0;
}
if(!$end_year){
    $end_year = 9999;
}
if(!$precision){
    $precision = 5;
}
if(!$separator){
    $separator = "\t";
}
my @species_list = ();
if($species){
    @species_list = split(/,/,$species);
}

my $csv = Text::CSV_XS->new({ sep_char => "\t", quote_char => undef, auto_diag => 1});

##
# Open the file for read data:
open my $csv_fh, "<:encoding(UTF-8)", "$file_input"  or die "$file_input: $!";

# Open the file for write data:
open my $OUTPUT, "+>:encoding(UTF-8)", "localization.dat" or die "localization.dat: $!";
print $OUTPUT "latitude;longitude;species\n";

# Set the character which will be used to indicate the end of a line.
local $/ = "\n";

###################################
print "1. build grid\n";
## Cell dimensions vary with latitude
## Length  Cell width  Cell height
##  1     ≤ 5000km     ×    5000km
##  2     ≤ 1250km     ×    625km
##  3     ≤ 156km      ×    156km
##  4     ≤ 39.1km     ×    19.5km
##  5     ≤ 4.89km     ×    4.89km
##  6     ≤ 1.22km     ×    0.61km
##  7     ≤ 153m       ×    153m
##  8     ≤ 38.2m      ×    19.1m
##  9     ≤ 4.77m      ×    4.77m
##  10    ≤ 1.19m      ×    0.596m
##  11    ≤ 149mm      ×    149mm
##  12    ≤ 37.2mm     ×    18.6mm


my $south_west_latitude  = 36;
my $south_west_longidude = -10;
my $north_east_latitude  = 44;
my $north_east_longitude = 3.5;


#open OUTGRID, '+>', "grid.dat";
#print OUTGRID "lat;lon;geohash\n";
my $grid = Geo::Hash::Grid->new(
        sw_lat    => $south_west_latitude,
        sw_lon    => $south_west_longidude,
        ne_lat    => $north_east_latitude,
        ne_lon    => $north_east_longitude,
        precision => $precision,
    );
my $hash_count        = $grid->count;
my $geohash_array_ref = $grid->hashes;

print "2. init big data\n";
my %big_data = ();
foreach my $geohash (@$geohash_array_ref) {
    my $gh = Geo::Hash->new;
    my ($lat,$lon) = $gh->decode( $geohash );

    ## filter geohash localizations
    #if(($lon >= -10 and $lon <= 0.5 and $lat >= 36 and $lat <= 40.5) or   ############################################################
    #   ($lon >= -10 and $lon <= 3.5 and $lat >= 40.5 and $lat <= 44)){
    if(($lon >= -10 and $lon <= 3.5 and $lat >= 36 and $lat <= 44)){
        $big_data{$geohash}{'sum'} = 0;
        $big_data{$geohash}{'species'} = ();
        #print OUTGRID "$lat;$lon;$geohash\n";
    }
    else{
        #print "skip  $geohash ($lat, $lon)\n" if($verbose);
        next;
    }
}
#close OUTGRID;

###################################
my @fields = ();
my ($ncol_gbifid, $ncol_lat, $ncol_lon, $ncol_species, $ncol_year, $ncol_month, $ncol_day);

my $lat_min=99;
my $lat_max=-99;
my $lon_min=999;
my $lon_max=-999;

# Loop through each line:
print "3. read each gbif observation (be patient...)\n";
my $nobsrepeat = 0;
my $nobsrejected = 0;
my $obsaccepted = 0;
my $nobs_outside_grid_or_time = 0;


## bind columns
my @cols = @{$csv->getline ($csv_fh)};
my $row = {};
$csv->bind_columns (\@{$row}{@cols});

my $nlines = 0;
while ($csv->getline($csv_fh)) {
    $nlines++;
    #gbifID  abstract    accessRights    accrualMethod   accrualPeriodicity  accrualPolicy   alternative audience    available   bibliographicCitation   conformsTo  contributor coverage
    #created creator date    dateAccepted    dateCopyrighted dateSubmitted   description educationLevel  extent  format  hasFormat   hasPart hasVersion  identifier  instructionalMethod isFormatOf  isPartOf
    #isReferencedBy  isReplacedBy    isRequiredBy    isVersionOf issued  language    license mediator    medium  modified    provenance  publisher   references  relation    replaces    requires    rights
    #rightsHolder    source  spatial subject tableOfContents temporal    title   type    valid   institutionID   collectionID    datasetID   institutionCode collectionCode  datasetName ownerInstitutionCode
    #basisOfRecord   informationWithheld dataGeneralizations dynamicProperties   occurrenceID    catalogNumber   recordNumber    recordedBy  individualCount organismQuantity    organismQuantityType
    #sex lifeStage   reproductiveCondition   behavior    establishmentMeansoccurrenceStatus  preparations    disposition associatedReferences    associatedSequences associatedTaxa  otherCatalogNumbers
    #occurrenceRemarks   organismIDorganismName  organismScope   associatedOccurrences   associatedOrganisms previousIdentifications organismRemarks materialSampleID    eventID parentEventID
    #fieldNumber eventDate   eventTime   startDayOfYear  endDayOfYear    year    month   day verbatimEventDate   habitat samplingProtocol    samplingEffort  sampleSizeValue sampleSizeUnit
    #fieldNotes  eventRemarks    locationID  higherGeographyID   higherGeography continent   waterBody   islandGroupisland   countryCode stateProvince   county  municipality    locality
    #verbatimLocality    verbatimElevation   verbatimDepth   minimumDistanceAboveSurfaceInMeters maximumDistanceAboveSurfaceInMeters locationAccordingTo locationRemarks decimalLatitude
    #decimalLongitude    coordinateUncertaintyInMeters   coordinatePrecision pointRadiusSpatialFit   verbatimCoordinateSystem    verbatimSRS footprintWKT    footprintSRS
    #footprintSpatialFit georeferencedBy georeferencedDate   georeferenceProtocol    georeferenceSources georeferenceVerificationStatus  georeferenceRemarks geologicalContextID
    #earliestEonOrLowestEonothemlatestEonOrHighestEonothem   earliestEraOrLowestErathem  latestEraOrHighestErathem   earliestPeriodOrLowestSystem    latestPeriodOrHighestSystem
    #earliestEpochOrLowestSeries latestEpochOrHighestSeries  earliestAgeOrLowestStage    latestAgeOrHighestStage lowestBiostratigraphicZone  highestBiostratigraphicZonelithostratigraphicTerms
    #group   formation   member  bed identificationID    identificationQualifier typeStatus  identifiedBy    dateIdentified  identificationReferences    identificationVerificationStatus
    #identificationRemarks   taxonID scientificNameID    acceptedNameUsageID parentNameUsageID   originalNameUsageID nameAccordingToID   namePublishedInID   taxonConceptID  scientificName
    #acceptedNameUsage   parentNameUsage originalNameUsage   nameAccordingTo namePublishedIn namePublishedInYear higherClassification    kingdom phylum  class   order   family  genus
    #subgenus    specificEpithet infraspecificEpithet    taxonRank   verbatimTaxonRank   vernacularName  nomenclaturalCode   taxonomicStatus nomenclaturalStatus taxonRemarks
    #datasetKey  publishingCountry   lastInterpreted elevation   elevationAccuracy   depth   depthAccuracy   distanceAboveSurface    distanceAboveSurfaceAccuracy    issue   mediaType
    #hasCoordinate   hasGeospatialIssues taxonKey    acceptedTaxonKey    kingdomKey  phylumKey   classKey    orderKey    familyKey   genusKey    subgenusKey speciesKey  species genericName acceptedScientificName
    #typifiedName    protocol    lastParsed  lastCrawled repatriated
    my $gbifid      = $row->{gbifID};
    my $lat         = $row->{decimalLatitude};
    my $lon         = $row->{decimalLongitude};
    my $speciesKey  = $row->{speciesKey};
    my $year        = $row->{year};
    my $month       = $row->{month};
    my $day         = $row->{day} ;
    my $date        = $year.'-'.$month.'-'.$day;

    if($date !~ /^\d{4}-\d+-\d+$/){
        $nobsrejected ++;
        next;
    }

    #print "$gbifid  $speciesKey, $lat, $lon, $date \n";

    ## filter zones and year
    if(
        #(($lon >= -10 and $lon <= 0.5 and $lat >= 36 and $lat <= 40.5) or ($lon >= -10 and $lon <= 3.5 and $lat >= 40.5 and $lat <= 44)) and
        ($lon >= -10 and $lon <= 3.5 and $lat >= 36 and $lat <= 44) and
        $year >= $ini_year and $year <= $end_year and
        (!@species_list or grep { $_ eq $speciesKey} @species_list)
        ){


        #if($nlines > 10){
        #    exit;
        #}


        #print "$gbifid  $lat  $lon  $speciesKey  $date \n";

        ## min max latitude
        $lat_min = $lat if($lat < $lat_min);
        $lat_max = $lat if($lat > $lat_max);

        ## min max longitude
        $lon_min = $lon if($lon < $lon_min);
        $lon_max = $lon if($lon > $lon_max);


        my $gh = Geo::Hash::XS->new();
        my $geohash = $gh->encode($lat, $lon, $precision);

        if(exists $big_data{$geohash}){

            ## filter: only consider one observation per day, per grid
            ## create key to save only one obs per day, per grid
            if(exists $big_data{$geohash}{dates}{$date}{$speciesKey}){
                #print "xxxx $geohash  $date   $speciesKey \n";
                $nobsrepeat ++;
                $nobsrejected ++;
                next;
            }
            else{
                #print "goto $geohash  $date   $speciesKey \n";
                $big_data{$geohash}{dates}{$date}{$speciesKey} = 1;

                $big_data{$geohash}{sum} += 1;
                $big_data{$geohash}{species}{$speciesKey}{id} = $speciesKey;
                if(exists $big_data{$geohash}{species}{$speciesKey}){
                    $big_data{$geohash}{species}{$speciesKey}{count} += 1;
                }else{
                    $big_data{$geohash}{species}{$speciesKey}{count} = 1;
                }

                $obsaccepted ++;
            }

        }else{
            #print " New geohash: $geohash  ->  $lat, $lon, $precision\n";
            #print " EXIT\n";
            #exit;
        }

        print $OUTPUT "$lat;$lon;$speciesKey\n";
    }
    else{
        $nobsrejected ++;
        $nobs_outside_grid_or_time ++;
        #print "not accepted: $gbifid  $lat  $lon  $speciesKey  $date \n";
    }
}

close $csv_fh;
close $OUTPUT;

#print Dumper(\%big_data);
print "\n";
print "\tobservations outside grid, time or selected species: $nobs_outside_grid_or_time\n";
print "\tobservations repeated: $nobsrepeat\n";
print "\tobservations rejected: $nobsrejected\n";
print "\tobservations accepted: $obsaccepted\n";
print "\tobservations total: $nlines\n\n";

## ============================================
print "4. process big data and output results\n";

my $size_geo = keys %big_data;

open OUTPUT, '+>:encoding(UTF-8)', $file_output;

open OUTGEO, '+>', "geohash.dat";
print OUTGEO "lat;lon;geohash;has_species\n";

my @community_list = ();
my $community_accepted = 0;
foreach my $geohash (sort keys %big_data){

    my $gh = Geo::Hash->new;
    my ( $lat, $lon ) = $gh->decode( $geohash );

    #print "  = Community $geohash ($lat, $lon)\n" if($verbose);
    my $community = Bio::Community->new( -name => $geohash );

    my @community_member_list = ();

    ## Add members to community
    while (my ($specie_key, $specie_val) = each %{ $big_data{$geohash}{'species'} } ){
        my $member = Bio::Community::Member->new( -id => $big_data{$geohash}{'species'}{$specie_key}{'id'} );
        my $member_count = $big_data{$geohash}{'species'}{$specie_key}{'count'};
        $community->add_member( $member, $member_count );

        for( my $i=1; $i <= $member_count; $i++ ){
            push @community_member_list,  $member->id();
        }
    }


    my $members_count = $community->get_members_count;
    my $richness = $community->get_richness;
    #print "   There are $members_count members in the community\n" if($verbose);
    #print "   The total diversity is $richness species\n" if($verbose);

    if($mverbose){
        while (my $member = $community->next_member) {
            my $member_id     = $member->id;
            my $member_count  = $community->get_count($member);
            my $member_rel_ab = $community->get_rel_ab($member);
            my $member_rank = $community->get_rank($member);
            print "     The relative abundance of member $member_id is $member_rel_ab % ($member_count counts). Rank is $member_rank\n";
        }
    }

    my $members = $community->get_all_members();
    #print Dumper($members);

    ## print header
    if($members_count > 2){

        push @community_list, $community;

        $community_accepted++;

        print OUTPUT "\"Cell $geohash ($lat, $lon)\"\t*SampleSet*\t1\t1\t1\r\n";
        print OUTPUT "$richness\t$members_count\r\n";
        for( my $i=1; $i<=$members_count; $i++ ){
            my $sample_id = sprintf("%05d", $i);
            print OUTPUT "\t$sample_id";
        }
        print OUTPUT "\r\n";

        ## set matrix
        my %matrix = (); ## matrix data for random get
        for( my $i=1; $i<=$members_count; $i++ ){
            ## get random member
            my $rand_index = 0 + int(rand($#community_member_list - 0));
            my $rand_member = $community_member_list[$rand_index];

            ## remove getted member from community
            splice @community_member_list, $rand_index, 1;
            #print "$geohash   $members_count  ".@community_member_list."  --> $rand_index (0 - ".@community_member_list.") --> $rand_member\n" if($verbose);

            while (my $member = $community->next_member) {
                my $member_id     = $member->id;
                if($member_id eq $rand_member){
                    push( @{ $matrix{$member_id} }, 1);
                }else{
                    push( @{ $matrix{$member_id} }, 0);
                }
            }
        }

        ## print matrix
        foreach my $member_id (sort {lc $a cmp lc $b} keys %matrix) {
            print OUTPUT "$member_id";
            for( my $i=0; $i<$members_count; $i++ ){
                print OUTPUT "\t".$matrix{$member_id}[$i];
            }
            print OUTPUT "\r\n";
        }

        #last if($community_accepted == 100);
    }
    print OUTGEO "$lat;$lon;$geohash;$members_count\n";

}
close OUTPUT;
close OUTGEO;

open INPUT,  '<',  $file_output      or die "Can't read old file: $!";
open OUTPUT, '>', "$file_output.new" or die "Can't write new file: $!";
print OUTPUT "*MultipleSampleSets*\t$community_accepted\t\"PT Community with more then 2 members";
if($ini_year != 0){
    print OUTPUT "; start year: $ini_year";
}
if($end_year != 9999){
    print OUTPUT ", end year: $end_year";
}
print OUTPUT "\"\r\n";

while( <INPUT> ){
        print OUTPUT $_;
}
close OUTPUT;
close INPUT;

unlink $file_output;
rename "$file_output.new", $file_output;


## Meta Community
my $meta = Bio::Community::Meta->new( -name => 'PT' );
$meta->add_communities( \@community_list );
print "\n== Metacommunities with more then 2 individuals:\n";
print "\t".$meta->get_communities_count." communities\n";
print "\t".$meta->get_richness." species\n";
print "\t".$meta->get_members_count." individuals\n";

my $end_time = time();
my $run_time = $end_time - $start_time;
print "\n\nTime elapsed: $run_time seconds\n";
