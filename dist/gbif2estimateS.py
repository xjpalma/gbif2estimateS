#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on Thu Oct 31 11:09:38 2019

@author: Jorge palma
         MARETEC/Instituto Superior Técnico
         Universidade de Lisboa
'''

import sys
import os
import gc
import argparse
import traceback
import time
import datetime
import pandas as pd
import json
import random
import string

## python-geohash: https://pypi.org/project/python-geohash/
import geohash
## pyshp: https://pythonhosted.org/Python%20Shapefile%20Library/
import shapefile
## https://shapely.readthedocs.io
from shapely.geometry import Point, Polygon, shape

sys.tracebacklimit=0

'''dev'''
sys.tracebacklimit=1


output_default_file = 'output.dat'
species_loc_file = 'speciesloc.dat'
grid_file = 'grid.dat'


class GeohashMaker(object):
    def __init__(self, precision, shapefile, geojson, bbox):
        self.precision = precision
        self.bbox = bbox
        self.shapefile = shapefile
        self.geojson = geojson

    def create_grid(self):
        if self.bbox:
            return self._buil_cell_tiles_from_bbox(self.bbox)
        elif self.shapefile:
            return self._build_cell_tiles_from_shapefile(self.shapefile)
        elif self.geojson:
            return self._build_cell_tiles_from_geojson(self.geojson)

    def is_in_grid(self, coordinates):
        if self.bbox:
            return self._is_coordinates_in_bounding_box(coordinates, self.bbox)
        elif self.shapefile:
            return self._is_coordinates_in_shapefile(coordinates, self.shapefile)
        elif self.geojson:
            return self._is_coordinates_in_geojson(coordinates, self.geojson)

    def build_cell(self, coordinates):
        return geohash.encode(coordinates[1], coordinates[0], precision=self.precision)

    def get_precision(self):
        return self.precision

    '''bounding_box'''
    def _is_geohash_in_bounding_box(self, current_geohash, bbox_coordinates):
        '''Checks if the box of a geohash is inside the bounding box
        :param current_geohash: a geohash
        :param bbox_coordinates: bounding box coordinates, [lon1, lat1, lon2, lat2]
        :return: true if the center of the geohash is in the bounding box
        '''

        # decode return [latitude, longitude]
        (latitude, longitude) = geohash.decode(current_geohash)

        geohash_in_bounding_box = (bbox_coordinates[0] < longitude < bbox_coordinates[2]) and \
                                  (bbox_coordinates[1] < latitude < bbox_coordinates[3])

        return geohash_in_bounding_box

    def _is_coordinates_in_bounding_box(self, coordinates, bbox_coordinates):
        '''Checks if coordinates is inside the bounding box
        :param coordinates: [lon, lat]
        :param bbox_coordinates: bounding box coordinates, [lon1, lat1, lon2, lat2]
        :return: true if the coordinate is in the bounding box
        '''

        coordinates_in_bounding_box = (bbox_coordinates[0] < coordinates[0] < bbox_coordinates[2]) and \
                                      (bbox_coordinates[1] < coordinates[1] < bbox_coordinates[3])

        return coordinates_in_bounding_box

    def _buil_cell_tiles_from_bbox(self, bbox_coordinates):
        '''Computes all geohash tile in the given bounding box
        :param bbox_coordinates: the bounding box coordinates of the geohashes
        :return: a list of geohashes
        '''

        checked_geohashes = set()
        geohash_stack = set()
        geohashes = []

        '''get center of bounding box, assuming the earth is flat'''
        center_longitude = (bbox_coordinates[0] + bbox_coordinates[2]) / 2
        center_latitude = (bbox_coordinates[1] + bbox_coordinates[3]) / 2

        center_geohash = self.build_cell([center_longitude, center_latitude])

        geohashes.append(center_geohash)
        geohash_stack.add(center_geohash)
        checked_geohashes.add(center_geohash)

        while len(geohash_stack) > 0:
            current_geohash = geohash_stack.pop()
            neighbors = geohash.neighbors(current_geohash)
            for neighbor in neighbors:
                if neighbor not in checked_geohashes and self._is_geohash_in_bounding_box(neighbor, bbox_coordinates):
                    geohashes.append(neighbor)
                    geohash_stack.add(neighbor)
                    checked_geohashes.add(neighbor)

        geohashes.sort()

        return geohashes

    '''shapefile'''
    def _is_coordinates_in_shapefile(self, coordinates, shpfile):

        ''' open shapefile'''
        sf = shapefile.Reader(shpfile)

        '''get features'''
        shapes = sf.shapes()
        first_shp = shapes[0]

        ''' get points coordinates for each point in the shape '''
        points = first_shp.points
        polygon = Polygon(points)

        point = Point(coordinates[0], coordinates[1])
        return polygon.contains(point)

    def _build_cell_tiles_from_shapefile(self, shpfile):
        '''Computes all geohash tiles in the given shapefile
        :param shapefile: shapefile
        :return: a list of geohashes
        '''

        ''' open shapefile'''
        sf = shapefile.Reader(shpfile)

        '''get features'''
        shapes = sf.shapes()
        if len(shapes) > 1:
            print("More than one feature was found. Only first will be selected.")
            input("Press Enter to continue...")
        '''only use first feature'''
        first_shp = shapes[0]

        ''' get shape type. only if shapetype is polygon'''
        shape_type = first_shp.shapeType
        if shape_type != 5:
            handle_error(msg='Shapefile feature be a polygon')

        ''' get points coordinates for each point in the shape '''
        points = first_shp.points
        polygon = Polygon(points)

        checked_geohashes = set()
        geohash_stack = set()
        geohashes = []

        '''get center of bounding box, assuming the earth is flat'''
        center_latitude = polygon.centroid.coords[0][1]
        center_longitude = polygon.centroid.coords[0][0]

        center_geohash = self.build_cell([center_longitude, center_latitude])

        geohashes.append(center_geohash)
        geohash_stack.add(center_geohash)
        checked_geohashes.add(center_geohash)

        while len(geohash_stack) > 0:
            current_geohash = geohash_stack.pop()
            neighbors = geohash.neighbors(current_geohash)
            for neighbor in neighbors:
                point = Point(geohash.decode(neighbor)[::-1])
                if neighbor not in checked_geohashes and polygon.contains(point):
                    geohashes.append(neighbor)
                    geohash_stack.add(neighbor)
                    checked_geohashes.add(neighbor)

        geohashes.sort()

        return geohashes

    '''geojson'''
    def _is_coordinates_in_geojson(self, coordinates, jsonfile):
        '''Checks if coordinates is inside the polygon
        :param coordinates: [lon, lat]
        :geojson file with polygon
        :return: true if the coordinate is in polygon
        '''

        with open(jsonfile) as f:
            try:
                data = json.load(f)
                polygon = shape(data["geometry"])
                point = Point(coordinates[0], coordinates[1])
                return polygon.contains(point)

            except ValueError as e:
                    handle_error(msg='Invalid GEOJSON format')

    def _build_cell_tiles_from_geojson(self, jsonfile):
        '''Computes all geohash tiles in the given geojson file
        :param jsonfile: geojson (polygon)
        :return: a list of geohashes
        '''
        with open(jsonfile) as f:
            try:
                data = json.load(f)

                polygon = shape(data["geometry"])
                geom_type = polygon.geom_type

                if geom_type != 'Polygon':
                    handle_error('SyntaxError', 'Invalid GEOJSON format: Must be a Polygon type')

                checked_geohashes = set()
                geohash_stack = set()
                geohashes = []

                '''get center of bounding box, assuming the earth is flat'''
                center_longitude = polygon.centroid.coords[0][0]
                center_latitude = polygon.centroid.coords[0][1]

                center_geohash = self.build_cell([center_longitude, center_latitude])

                geohashes.append(center_geohash)
                geohash_stack.add(center_geohash)
                checked_geohashes.add(center_geohash)

                while len(geohash_stack) > 0:
                    current_geohash = geohash_stack.pop()
                    neighbors = geohash.neighbors(current_geohash)
                    for neighbor in neighbors:
                        point = Point(geohash.decode(neighbor)[::-1])
                        if neighbor not in checked_geohashes and polygon.contains(point):
                            geohashes.append(neighbor)
                            geohash_stack.add(neighbor)
                            checked_geohashes.add(neighbor)

                geohashes.sort()

                return geohashes

            except ValueError as e:
                handle_error(msg='Invalid GEOJSON format')


class BioCommunity(object):
    units_count = 0
    members = {}

    def __init__(self, name):
        self.name = name
        self.reset()

    def add_member(self, member, count=1):
        member_id = member.getid()
        try:
            self.members[member_id] = {'member': member, 'count': count}
            self.units_count += count
            return True
        except KeyError:
            return False

    def remove_member(self, member, count=1):
        member_id = member.getid()
        try:
            if self.members[member_id][count] <= count:
                self.units_count -= self.members['member_id'][count]
                del self.members[member_id]
            else:
                self.members[member_id][count] -= count
            return True
        except KeyError:
            return False

    def add_communities(self, communities):
        for community in communities:
            members = community.get_all_members()
            for speciekey in members:
                member = BioMember(speciekey)
                self.add_member(member, members[speciekey]['count'])

    def get_units_count(self):
        return self.units_count

    def get_all_members(self):
        return self.members

    def set_name(self, name):
        self.name = name
        return True

    def get_name(self):
        return self.name

    ''' Report the community richness or number of different types of members.
    This is a form of alpha diversity. '''

    def get_richness(self):
        return len(self.members)

    def reset(self):
        '''Re-initialize the community'''
        self.units_count = 0
        self.members = {}
        return True


class BioMember(object):
    def __init__(self, id=''):
        self.id = id
        if not id:
            self.id = self.randomStringDigits()

    def getid(self):
        return self.id

    def setid(self, id):
        self.id = id

    def randomStringDigits(stringLength=6):
        '''Generate a random string of letters and digits '''
        lettersAndDigits = string.ascii_letters + string.digits
        return ''.join(random.choice(lettersAndDigits) for i in range(stringLength))


def get_parser():
    ''' Get parser object '''
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='read gbif and make input file to EstimateS')

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()

    parser.add_argument('-v', dest='verbose', help='verbose', action='store_true')
    parser.add_argument('-vv', dest='vverbose', help='more verbose', action='store_true')

    ## Create io files group
    subparser_io = parser.add_argument_group(title='IO group')
    subparser_io.add_argument('-i', dest='input', help='csv gbif results', required=True)
    subparser_io.add_argument('-s', dest='separator', help='csv separator', default='\t', required=False)
    subparser_io.add_argument('-o', dest='output', help='output file', default=output_default_file, required=False)

    ## Create time group
    subparser_time = parser.add_argument_group(title='time group')
    subparser_time.add_argument('-str', dest='strdate', help="the Start Date format YYYYMMDD",
                                type=lambda d: datetime.datetime.strptime(d, '%Y%m%d'), required=False)
    subparser_time.add_argument('-end', dest='enddate', help="the End Date format YYYYMMDD",
                                type=lambda d: datetime.datetime.strptime(d, '%Y%m%d'), required=False)

    ## Create grid group
    subparser_grid = parser.add_argument_group(title='grid group')
    subparser_grid.add_argument('-g', dest='grid_type', choices=['geohash'], default='geohash', required=False)
    subparser_grid.add_argument('-p', dest='precision', type=int, help='grid precision', default=5, required=False)

    subparser_grid_exclusive = subparser_grid.add_mutually_exclusive_group(required=True)
    subparser_grid_exclusive.add_argument('-shp', dest='shapefile', help='shapefile with polygon', required=False)
    subparser_grid_exclusive.add_argument('-geojson', dest='geojson', help='geojson file with polygon', required=False)
    subparser_grid_exclusive.add_argument('-bbox', dest='bbox', nargs='+', type=float, help='bounding box: x1 y1 x2 y2', required=False)

    ## Create species group
    subparser_specie = parser.add_argument_group(title='specie group')
    subparser_specie.add_argument('-n', dest='species', nargs='+', default=[], help='species allowed', required=False)

    args = parser.parse_args()

    if args.vverbose:
        args.verbose = True

    if not os.path.isfile(args.input):
        raise IOError('No such file {}'.format(args.input))

    args.outdir = os.path.dirname(args.output)
    outfile = os.path.basename(args.output)
    ## verify if is a path and create it
    if args.outdir:
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
        args.outdir = args.outdir + '/'
        ## verify if is a path with filename
        if not outfile:
            args.output = args.outdir + '/output.dat'

    if not args.strdate:
        args.strdate = datetime.datetime.strptime('1900-01-01', '%Y-%m-%d')

    if not args.enddate:
        args.enddate = datetime.datetime.strptime('2100-01-01', '%Y-%m-%d')

    if args.shapefile:
        if not os.path.isfile(args.shapefile):
            handle_error('FileNotFoundError', 'Shapefile not found')

    if args.geojson:
        if not os.path.isfile(args.geojson):
            handle_error('FileNotFoundError', 'JSON file not found')

    return args


def handle_error(error='', msg=''):
    formatted_lines = traceback.format_exc().splitlines()
    print()

    if error:
        exec('raise ' + error + '(\'' + msg + '\')')
    elif msg:
        print(msg)
    else:
        print()

    sys.exit(1)


def test_csv_header(csv_columns, must_exist_columns):
    not_exist_columns = set()
    for elem in must_exist_columns:
        if not elem in csv_columns:
            not_exist_columns.add(elem)

    if not_exist_columns:
        print("Error: Missing columns in input file.\n       Input file doesn\'t have " + \
              '"{}"'.format('", "'.join(not_exist_columns)) + \
              ' columns')

        return False

    del not_exist_columns, csv_columns, must_exist_columns

    return True;


def progress(count, total, status=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 0)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('{:} {:}% ...{:}\r'.format(bar, percents, status))
    sys.stdout.flush()


def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)


def elapsed_time(start_time):
    end_time = time.time()

    if isinstance(start_time, list):
        print()
        for st in start_time:
            secs = (end_time - st)
            print("Time elapsed: {}".format(str(datetime.timedelta(seconds=secs))))
        print()

    else:
        secs = (end_time - start_time)
        print()
        print("Time elapsed: {}".format(str(datetime.timedelta(seconds=secs))))
        print()


if __name__ == "__main__":
    start_time = time.time()

    args = get_parser()

    '''#### build grid'''
    print('1. build grid')
    grid = []
    if args.grid_type == 'geohash':
        grid_maker = GeohashMaker(args.precision, args.shapefile, args.geojson, args.bbox)
        grid = grid_maker.create_grid()

    else:
        handle_error(msg='Error: only accept geohash grid type')

    '''#### init big_data variable'''
    print('2. init big data')
    big_data = {}
    for cell in grid:
        big_data[cell] = {}
        '''how many species in cell'''
        big_data[cell]['sum'] = 0
        '''list of species in cell'''
        big_data[cell]['species'] = {}
        '''used to consider only one observation (specie and time) in cell'''
        big_data[cell]['dates'] = {}

    '''create localization.dat file'''
    f = open(args.outdir + species_loc_file, 'w+')
    f.write("latitude;longitude;species\n")

    '''#### read csv file'''
    print('3. read each gbif observation (be patient...)')
    nobs_accepted = 0
    nobs_rejected = 0
    nobs_repeated = 0
    nobs_outside_grid_or_time = 0
    nobs_wrong_format = 0
    nobs = 0
    usecols = ['gbifID', 'decimalLatitude', 'decimalLongitude', 'speciesKey', 'year', 'month', 'day']
    chunksize = 10 ** 5
    filesize = os.path.getsize(args.input)
    linesize = 820
    for df in pd.read_csv(args.input, sep=args.separator, chunksize=chunksize, engine='c', low_memory=False, usecols=usecols, skip_blank_lines=True):
        s_time = time.time()

        nlines = len(df.index)
        nobs += nlines

        ''' verify if all columns exist in header csv'''
        csv_columns = df.columns.tolist()
        test_csv_header(csv_columns, usecols)

        '''
        gbifID  abstract    accessRights    accrualMethod   accrualPeriodicity  accrualPolicy   alternative audience    available   bibliographicCitation   conformsTo  contributor coverage
        created creator date    dateAccepted    dateCopyrighted dateSubmitted   description educationLevel  extent  format  hasFormat   hasPart hasVersion  identifier  instructionalMethod isFormatOf  isPartOf
        isReferencedBy  isReplacedBy    isRequiredBy    isVersionOf issued  language    license mediator    medium  modified    provenance  publisher   references  relation    replaces    requires    rights
        rightsHolder    source  spatial subject tableOfContents temporal    title   type    valid   institutionID   collectionID    datasetID   institutionCode collectionCode  datasetName ownerInstitutionCode
        basisOfRecord   informationWithheld dataGeneralizations dynamicProperties   occurrenceID    catalogNumber   recordNumber    recordedBy  individualCount organismQuantity    organismQuantityType
        sex lifeStage   reproductiveCondition   behavior    establishmentMeansoccurrenceStatus  preparations    disposition associatedReferences    associatedSequences associatedTaxa  otherCatalogNumbers
        occurrenceRemarks   organismIDorganismName  organismScope   associatedOccurrences   associatedOrganisms previousIdentifications organismRemarks materialSampleID    eventID parentEventID
        fieldNumber eventDate   eventTime   startDayOfYear  endDayOfYear    year    month   day verbatimEventDate   habitat samplingProtocol    samplingEffort  sampleSizeValue sampleSizeUnit
        fieldNotes  eventRemarks    locationID  higherGeographyID   higherGeography continent   waterBody   islandGroupisland   countryCode stateProvince   county  municipality    locality
        verbatimLocality    verbatimElevation   verbatimDepth   minimumDistanceAboveSurfaceInMeters maximumDistanceAboveSurfaceInMeters locationAccordingTo locationRemarks decimalLatitude
        decimalLongitude    coordinateUncertaintyInMeters   coordinatePrecision pointRadiusSpatialFit   verbatimCoordinateSystem    verbatimSRS footprintWKT    footprintSRS
        footprintSpatialFit georeferencedBy georeferencedDate   georeferenceProtocol    georeferenceSources georeferenceVerificationStatus  georeferenceRemarks geologicalContextID
        earliestEonOrLowestEonothemlatestEonOrHighestEonothem   earliestEraOrLowestErathem  latestEraOrHighestErathem   earliestPeriodOrLowestSystem    latestPeriodOrHighestSystem
        earliestEpochOrLowestSeries latestEpochOrHighestSeries  earliestAgeOrLowestStage    latestAgeOrHighestStage lowestBiostratigraphicZone  highestBiostratigraphicZonelithostratigraphicTerms
        group   formation   member  bed identificationID    identificationQualifier typeStatus  identifiedBy    dateIdentified  identificationReferences    identificationVerificationStatus
        identificationRemarks   taxonID scientificNameID    acceptedNameUsageID parentNameUsageID   originalNameUsageID nameAccordingToID   namePublishedInID   taxonConceptID  scientificName
        acceptedNameUsage   parentNameUsage originalNameUsage   nameAccordingTo namePublishedIn namePublishedInYear higherClassification    kingdom phylum  class   order   family  genus
        subgenus    specificEpithet infraspecificEpithet    taxonRank   verbatimTaxonRank   vernacularName  nomenclaturalCode   taxonomicStatus nomenclaturalStatus taxonRemarks
        datasetKey  publishingCountry   lastInterpreted elevation   elevationAccuracy   depth   depthAccuracy   distanceAboveSurface    distanceAboveSurfaceAccuracy    issue   mediaType
        hasCoordinate   hasGeospatialIssues taxonKey    acceptedTaxonKey    kingdomKey  phylumKey   classKey    orderKey    familyKey   genusKey    subgenusKey speciesKey  species genericName acceptedScientificName
        typifiedName    protocol    lastParsed  lastCrawled repatriated
        '''

        for index, row in df.iterrows():

            if args.verbose:
                if nlines < chunksize:
                    progress(index, nlines)
                else:
                    progress(index*linesize, filesize)

            '''get values'''
            try:
                gbifid = str(row['gbifID'])
                speciekey = int(float(row['speciesKey']))
                lon = round(float(row['decimalLongitude']), 6)
                lat = round(float(row['decimalLatitude']), 6)
                year = row['year']
                month = row['month']
                day = row['day']
                date_obj = datetime.datetime(int(year), int(month), int(day), 0, 0)
                date = date_obj.strftime("%Y-%m-%d")

                #print(index, gbifid, speciekey, lon, lat, date)

            except Exception as exception:
                # traceback.print_exc()
                nobs_wrong_format += 1
                nobs_rejected += 1
                continue
            else:
                '''test if observation is in domain, in time and in species list'''
                if grid_maker.is_in_grid([lon, lat]) and date_obj >= args.strdate and date_obj <= args.enddate and (not args.species or speciekey in args.species):

                    cell = grid_maker.build_cell([lon, lat])

                    if not cell in grid:
                        nobs_outside_grid_or_time += 1
                        #print(cell + '  ' + str(lon) + '  ' + str(lat) + '   is not in grid')
                        #handle_error(msg=cell + '  ' + str(lon) + '  ' + str(lat) + '   is not in grid')
                        continue

                    '''
                    filter: only consider one observation per day, per grid
                    create key to save only one obs per day, per grid
                    '''
                    try:
                        if speciekey in big_data[cell]['dates'][date]:
                            # print('repeated: ' + cell + '  ' + date + '  ' + speciekey)
                            nobs_repeated += 1
                            nobs_rejected += 1
                            continue

                        else:
                             # print('accepted: ' + cell + '  ' + date + '  ' + speciekey)
                            big_data[cell]['dates'][date].append(speciekey)

                            big_data[cell]['sum'] += 1

                            if speciekey in big_data[cell]['species']:
                                big_data[cell]['species'][speciekey]['count'] += 1
                            else:
                                big_data[cell]['species'][speciekey] = {'speciekey': speciekey, 'count': 1}

                            nobs_accepted += 1

                            f.write("{0:.6f};{1:.6f};{2:d}\r\n".format(lat, lon, speciekey))

                    except KeyError:
                        # print('accepted: ' + cell + '  ' + date + '  ' + speciekey)
                        big_data[cell]['dates'][date] = []
                        big_data[cell]['dates'][date].append(speciekey)

                        big_data[cell]['sum'] += 1

                        if speciekey in big_data[cell]['species']:
                            big_data[cell]['species'][speciekey]['count'] += 1
                        else:
                            big_data[cell]['species'][speciekey] = {'speciekey': speciekey, 'count': 1}

                        nobs_accepted += 1

                        f.write("{0:.6f};{1:.6f};{2:d}\r\n".format(lat, lon, speciekey))

                else:
                    # print('out of domain or out of time period: ' + cell + '  ' + str(lon) + '  ' + str(lat) + '  ' + date + '  ' + str(speciekey))
                    nobs_outside_grid_or_time += 1
                    nobs_rejected += 1

        if args.vverbose: elapsed_time([s_time, start_time])

        del df
        gc.collect()

    f.close()

    print()
    print('\tobservations outside grid, time or selected species: {0}'.format(nobs_outside_grid_or_time))
    print('\tobservations wrong format (no date): {0}'.format(nobs_wrong_format))
    print('\tobservations repeated: {0}'.format(nobs_repeated))
    print('\tobservations rejected: {0}'.format(nobs_rejected))
    print('\tobservations accepted: {0}'.format(nobs_accepted))
    print('\tobservations total: {0}'.format(nobs))
    print()

    ''' delete unecessary variables '''
    for c in big_data:
        del big_data[c]['dates']


    print('4. process big data and output results')

    '''open output files'''
    fout = open(args.output, 'w+')
    fgeo = open(args.outdir + grid_file, 'w+')
    fgeo.write("lat;lon;geohash;has_species\r\n")

    community_list = []
    community_accepted = 0
    units_count = 0
    for cell in big_data:
        (lat, lon) = coordinates = geohash.decode(cell)

        '''create new community'''
        community = BioCommunity(cell)
        community_unit_list = []
        for speciekey, v in big_data[cell]['species'].items():
            member = BioMember(speciekey)
            member_count = v['count']

            '''add member to community'''
            community.add_member(member, member_count)

            '''add member to member list'''
            for i in range(member_count):
                community_unit_list.append(member.getid())

        community_units_count = community.get_units_count()
        richness = community.get_richness()

        '''only consider community with more than 2 observations'''
        if community_units_count > 2:
            community_accepted += 1
            units_count += community_units_count

            if args.vverbose: print("    There are {mc:} units in community {cell:} ({lat:}, {lon:}. The total diversity is {rich:} species)".format(mc=community_units_count, cell=cell, lat=lat, lon=lon, rich=community.get_richness()))

            '''add community to list'''
            community_list.append(community)

            '''print header
            "Cell eydsh (37.28759765625, -7.53662109375)" * SampleSet * 1   1   1
             8   8
                   00001   00002   00003   00004   00005   00006   00007   00008
            '''
            fout.write("\"Cell {cell:} ({lat:}, {lon:})\"\t*SampleSet*\t1\t1\t1\r\n".format(cell=cell, lat=lat, lon=lon))
            fout.write("{r:}\t{uc:}\r\n".format(r=richness, uc=community_units_count))

            for i in range(1, community_units_count + 1):
                fout.write("\t{0:05d}".format(i))
            fout.write("\r\n")

            '''set matrix data for random get'''
            matrix = {}
            members = community.get_all_members()

            '''init matrix'''
            for speciekey in members:
                matrix[speciekey] = []

            for i in range(community_units_count):
                '''get random member'''
                random_member = community_unit_list.pop(random.randrange(len(community_unit_list)))

                for speciekey in members:
                    if speciekey == random_member:
                        matrix[speciekey].append(1)
                    else:
                        matrix[speciekey].append(0)
            '''
            print matrix
            2474051    0       0       1       0       0       0       0       0
            2492606    0       0       0       0       0       0       0       1
            2492867    0       0       0       0       1       0       0       0
            2495000    0       0       0       0       0       1       0       0
            2498415    0       0       0       0       0       0       1       0
            5229493    1       0       0       0       0       0       0       0
            6092830    0       0       0       1       0       0       0       0
            9515886    0       1       0       0       0       0       0       0
            '''

            for speciekey in sorted(matrix):
                fout.write("{0:d}".format(speciekey))
                for i in range(community_units_count):
                    fout.write("\t{0:}".format(int(matrix[speciekey][i])))
                fout.write("\r\n")


        fgeo.write("{lat:};{lon:};{cell:};{uc:}\r\n".format(lat=lat, lon=lon, cell= cell, uc=community_units_count))

    fout.close()
    fgeo.close()

    '''add first line in output file'''
    first_line = "*MultipleSampleSets*\t{0:}\t\"PT Community with more then 2 members".format(community_accepted)
    if args.strdate.year != 1900:
        first_line += '; start year: ' + str(args.strdate.year)

    if args.enddate.year != 2100:
        first_line += '; end year: ' + str(args.enddate.year)

    first_line += "\"\r\n"

    line_prepender(args.output, first_line)


    '''print stat'''
    meta = BioCommunity('tmp')
    meta.add_communities(community_list)
    print("\n== Metacommunities with more then 2 individuals:")
    print("\t{0:} communities".format(len(community_list)))
    print("\t{0:} species".format(meta.get_richness()))
    print("\t{0:} individuals".format(units_count))


    '''the end'''
    elapsed_time(start_time)
