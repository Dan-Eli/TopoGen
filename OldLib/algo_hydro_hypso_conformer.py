#!/usr/local/bin/python
# -=- encoding: utf-8 -=-
#####################################################################################################################################

"""
    This algorithm corrects hypsographic features so that they are no more in conflict with hydrography.

    
    *Usage*:
        from algo_hydro_hypso_conformer import AlgoHydroHypsoConformer

    Limits and constraints
        This application do not, for the moment, correct any error that either touches or crosses waterbody chordal axis or
        touches the virtual shorelines on both sides of watercourses. This may change in the future.

"""
__revision__ = "--REVISION-- : $Id: algo_hydro_hypso_conformer.py 312 2011-05-09 19:35:40Z langlois $"

#####################################################################################################################################

from collections import Iterable

from shapely.geometry import Point, MultiPoint, LineString, MultiLineString, Polygon
from shapely.ops import linemerge

from lib_genmetal import MA_Point, MA_LineString, MA_Polygon, GenUtil, Holder, \
                         GenStatistics, Parameters, SpatialContainer, Algorithm, \
                         PolygonModifier, ChordalAxisTransformer, PolygonMerger, \
                         InternalError

########################################################

# Public constant

# Properties definition

# Define attribute name
_UUID_MASTER = '_uuid_master'
_UUID_RELATED = '_uuid_related'
_UUID_FEAT = '_uuid_feat'
_GEN_FTYPE = 'gen_ftype'
_BDG_ID = '_bdg_id'
_MAX_TOLERANCE_WB = '_max_tolerance_wb'
_MAX_TOLERANCE_SLW = '_max_tolerance_slw'
_TEST_SIMPLE_LINE = '_test_simple_line'
_TEST_CROSSING_LINE = '_test_crossing_line'
_TEST_SIDEDNESS = '_test_sidedness'
_SEGMENT_LENGTH = '_segment_length'
_SEGMENT_LENGTH_RANK = '_segment_length_rank'

# Define attribute value
_BOTTLENECK = 'bottleneck'
_TO_BE_DELETED = 'to_be_deleted'
_TO_BE_BUFFERED = 'to_be_buffered'

# Define statistics constant
SLW_CORRECTED = 'SLW_C'
SLW_UNCORRECTED = 'SLW_UC'
WB_CORRECTED = 'WB_C'
WB_UNCORRECTED = 'WB_UC'
TOTAL_IN_LINES = 'TOTAL_IN_LINES'
TOTAL_IN_POLYGONS = 'TOTAL_IN_POLYGONS'
TOTAL_CC = 'TOTAL_CC'
TOTAL_CA = 'TOTAL_CA'
TOTAL_B_SLW = 'TOTAL_B_SLW'
TOTAL_SLW = 'TOTAL_SLW'
TOTAL_S_SLW = 'TOTAL_S_SLW'
TOTAL_S_WB = 'TOTAL_S_WB'
TOTAL_S_WBX = 'TOTAL_S_WBX'
TOTAL_WB = 'TOTAL_WB'
TOTAL_B_WB = 'TOTAL_B_WB'

#Internal key word constants
_ALGO = "genHydroHypsoConformer"

_SEARCH_TOLERANCE = .001 # Search tolerance in meter for the RTree search

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
class Statistics(GenStatistics):
    """Class that contains the statistics for this algorithm

    Attributes
        stat_names: Name of the statistics for the SpikeStatistics class. These name are
                    used by the Statistics class

    """

    # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    def __init__(self):
        """Initialize the attributes of an object of the class"""

        GenStatistics.__init__(self)
        self.stats_names = ((_ALGO,
                             SLW_CORRECTED, SLW_UNCORRECTED, WB_CORRECTED, WB_UNCORRECTED,
                             TOTAL_IN_LINES, TOTAL_IN_POLYGONS,
                             TOTAL_CC, TOTAL_CA, TOTAL_B_SLW, TOTAL_SLW, TOTAL_S_SLW,
                             TOTAL_S_WB, TOTAL_S_WBX, TOTAL_WB, TOTAL_B_WB,
                             GenUtil.INVALID, GenUtil.CROSSING_LINE, GenUtil.SIDEDNESS, GenUtil.SIMPLE_LINE))

    # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    def get_stats (self, type=GenStatistics.SUMMARY):
        """Extract the current statistics and build  a list of string that forms the statistical message"

        Parameters:
            type: Give the form of statistics to extract. Can take 2 values.
                SUMMARY: Summary information
                DETAILED: Detailed information

        """

        str_out = []
        str_out.append( "Hydro Hypso COnformer Statistics" )

        str_out.append( "----------------------------" )
        str_out.append("Number of _SLW_CORRECTED : " + str(self.get_stats_name_count_total(_SLW_CORRECTED)))
        str_out.append("Number of _SLW_UNCORRECTED : " + str(self.get_stats_name_count_total(_SLW_UNCORRECTED)))
        str_out.append("Number of _WB_CORRECTED : " + str(self.get_stats_name_count_total(_WB_CORRECTED)))
        str_out.append("Number of _WB_UNCORRECTED : " + str(self.get_stats_name_count_total(_WB_UNCORRECTED)))
        str_out.append( "----------------------------" )

        return str_out


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
class AlgoHydroHypsoConformer(Algorithm):
    """This is the main class for this algorithm

    """

    # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    def __init__(self, test_simple_line=True, test_crossing_line=True, test_sidedness=False, debug=False):
        """Initialize the attributes of an object of the class DPAlgorithm

        Parameters:
            test_simple_line: Enable (True) or disable (False) simple line constraint test
            test_crossing_line: Enable (True) or disable (False) crossing line constraint test
            debug: Flag to enable (True) or disable (False) the output of extra information

        Return value:
            None

        """

        Algorithm.__init__(self)
        
        self.params = Parameters()
        self.params.test_simple_line = test_simple_line
        self.params.test_crossing_line = test_crossing_line
        self.params.test_sidedness = test_sidedness
        self.params.debug = debug

        self.stats = Statistics()
        self.stats.add_iteration()
        self.polygon_id = -1


    # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    def check_features(self):
        """
        Check if the features passed in parameters are of the valid class type and have the required attributes

        Parameters: None

        Return value: None

        """

        # Check each feature in turn
        properties_to_check = [_UUID_MASTER, _GEN_FTYPE]
        
        GenUtil.print_debug(self.params, '\nSTART CHECKING features:')
        GenUtil.print_debug(self.params, '  - checking for properties: %s' % (properties_to_check))
        
        for feature in self.features:
            if isinstance(feature, MA_LineString):
                GenUtil.check_feature_integrity(feature, MA_LineString, properties_to_check)
                
            if isinstance(feature, MA_Polygon):
                GenUtil.check_feature_integrity(feature, MA_Polygon, properties_to_check)

        GenUtil.print_debug(self.params, 'END CHECKING features')


    # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    def _load_feature(self, feature, container, stats_item=None):
        """Load a feature into a spatial container and optionnally update statistics accordingly.
        
        *Note*: Due to a still unknown situation, FME sometimes pass aggregate features which are badly 
                managed by the infrastructure. Care is here take to kick out these features before processing.
        """
        
        # HANDLE the FME aggregate case 
        if hasattr(feature, "_sci_id"):
            GenUtil.print_debug(self.params, '   feature is already loaded : %s' % (str(feature.ma_properties)))
        else:
            container.add_feature(feature)
##            if (stats_item):
##                stats_item = GenUtil._make_iterable(stats_item)
##                [self.stats.add_stats(s) for s in stats_item]

        return


    # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    def _load_features(self):
        """Overload the Algorithm load_features method

        *Parameters*: *None*

        *Returns*: *None*

        """

        # Create the spatial containers that will maintain all spatial features
        sc_CC = SpatialContainer()
        sc_CA = SpatialContainer()
        sc_SHR = SpatialContainer()
        sc_SLW = SpatialContainer()
        sc_B_SLW = SpatialContainer()
        sc_S_SLW = SpatialContainer()
        sc_S_SLWX = SpatialContainer()
        sc_WB = SpatialContainer()
        sc_B_WB = SpatialContainer()
        sc_S_WB = SpatialContainer()
        sc_S_WBX = SpatialContainer()

        GenUtil.print_debug (self.params, '\nSTART LOADING features into spatial containers:')
        GenUtil.print_debug (self.params, '  number of features to load: %d' % len(self.features))
        
        # Load each feature in the appropriate spatial container depending of its type
        for feature in self.features:
            GenUtil.print_debug(self.params, '  feature.ma_properties: %s' % str(feature.ma_properties))
            
            if (feature.ma_properties[_GEN_FTYPE] == 'CC'):
                self._load_feature(feature, sc_CC, ['TOTAL_IN_LINES', 'TOTAL_CC'])

            elif (feature.ma_properties[_GEN_FTYPE] == 'CA'):
                self._load_feature(feature, sc_CA, ['TOTAL_IN_LINES', 'TOTAL_CA'])

#            elif (feature.ma_properties[_GEN_FTYPE] in ['B_SLW','B_SLW_100','B_SLW_75','B_SLW_50','B_SLW_25']):
            elif (feature.ma_properties[_GEN_FTYPE] in ['B_SLW_100','B_SLW_75','B_SLW_50','B_SLW_25']):
                self._load_feature(feature, sc_B_SLW, ['TOTAL_IN_LINES', 'TOTAL_B_SLW'])

            elif (feature.ma_properties[_GEN_FTYPE] == 'SLW'):
                self._load_feature(feature, sc_SLW, ['TOTAL_IN_LINES', 'TOTAL_SLW'])

            elif (feature.ma_properties[_GEN_FTYPE] == 'S_SLW'):
                self._load_feature(feature, sc_S_SLW, ['TOTAL_IN_LINES', 'TOTAL_S_SLW'])

            elif (feature.ma_properties[_GEN_FTYPE] == 'S_SLWX'):
                self._load_feature(feature, sc_S_SLWX, ['TOTAL_IN_LINES', 'TOTAL_S_SLWX'])

            elif (feature.ma_properties[_GEN_FTYPE] == 'S_WB'):
                self._load_feature(feature, sc_S_WB, ['TOTAL_IN_LINES', 'TOTAL_S_WB'])

            elif (feature.ma_properties[_GEN_FTYPE] == 'S_WBX'):
                self._load_feature(feature, sc_S_WBX, ['TOTAL_IN_LINES', 'TOTAL_S_WBX'])

            elif (feature.ma_properties[_GEN_FTYPE] == 'WB'):
                self._load_feature(feature, sc_WB, ['TOTAL_IN_POLYGONS', 'TOTAL_WB'])

            elif (feature.ma_properties[_GEN_FTYPE] in ['B_WB_100','B_WB_75','B_WB_50','B_WB_25']):
                self._load_feature(feature, sc_B_WB, ['TOTAL_IN_POLYGONS', 'TOTAL_B_WB'])

            elif (feature.ma_properties[_GEN_FTYPE] == 'B_SHR'):
                self._load_feature(feature, sc_SHR, ['TOTAL_IN_POLYGONS', 'TOTAL_B_SLW'])

            else:
                print "  ERROR: Missing spatial container for gen_ftype = %s" % (feature.ma_properties[_GEN_FTYPE])

        GenUtil.print_debug (self.params, 'END LOADING features into spatial containers')

        #Reset the main features container in order to reload the transformed features at the end
        self.features = []

        # make sure to return all spatial containers to work with
        return (sc_CC, sc_CA, sc_SHR, sc_SLW, sc_B_SLW, sc_S_SLW, sc_S_SLWX, sc_WB, sc_B_WB, sc_S_WB, sc_S_WBX)


    # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    def _cut_line_using_segment(self, line, segment):
        """Cut a line in 3 parts based on a line segment.
        """
        distance_start = line.project(Point(segment.coords[0]))
        distance_end   = line.project(Point(segment.coords[-1]))
                        
        if (distance_end > distance_start):
            delta = distance_end - distance_start
        else:
            delta = distance_start - distance_end
            distance_start, distance_end = distance_end, distance_start     #swap segment ends to be in line with the line
        
#        first_cut = self._cut(line, distance_start)
#        last_cut  = self._cut(first_cut[1], delta)
        first_cut = GenUtil.cut_line_distance(line, distance_start)
        last_cut  = GenUtil.cut_line_distance(first_cut[1], delta)
        
        return [first_cut[0], last_cut[0], last_cut[1]]
       
    
    # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    def _cut(self, line, distance):
        # Cuts a line in two at a distance from its starting point
        if distance <= 0.0 or distance >= line.length:
            return [LineString(line)]
        coords = list(line.coords)
        for i, p in enumerate(coords):
            pd = line.project(Point(p))
            if pd == distance:
                return [
                    LineString(coords[:i+1]),
                    LineString(coords[i:])]
            if pd > distance:
                cp = line.interpolate(distance)
                return [
                    LineString(coords[:i] + [(cp.x, cp.y)]),
                    LineString([(cp.x, cp.y)] + coords[i:])]    
    
   
    # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    def _process_contour_over_waterbody (self, segments, sc_cc, sc_b_wb):
        """
        """

        GenUtil.print_debug (self.params, '\nSTART PROCESSING contour over waterbodies:')
        GenUtil.print_debug (self.params, '  number of segments to process: %d' % (len(segments)))       

        # SORT segments list based on their length rank (as set by FME)
        # This is done to make sure that segments are processed in the right order (shortest to longest)
        sorted(segments, key=lambda by_item: by_item.ma_properties[_SEGMENT_LENGTH_RANK])
                 
        # PROCESS each segment in turn
        for segment in GenUtil.make_iterable(segments):

            GenUtil.print_debug(self.params, '  SEGMENT.ma_properties: %s' % str(segment.ma_properties))    
            
            # extract segment identifiers for associated contour and waterbody
            _uuid_master = segment.ma_properties[_UUID_MASTER]      # contour id
            _uuid_related = segment.ma_properties[_UUID_RELATED]    # waterbody id
            
            GenUtil.print_debug(self.params, '    using ids :: _uuid_master: %s | _uuid_related: %s' % (_uuid_master, _uuid_related))            

            # get contour from spatial container
            feat_list = sc_cc.get_features(filter="feature.ma_properties['%s']=='%s'" % (_UUID_MASTER, _uuid_master))
            if (len(feat_list) == 1):        
                contour = feat_list[0]
            else:
                raise InternalError ("%s :: Waterbody : Contours container returns %d feature(s), should be 1" % (_ALGO, len(feat_list)))
            
            # get first waterbody's buffer from spatial container and make sure that it is a polyline 
            everything_ok = True
            feat_list = sc_b_wb.get_features(filter="feature.ma_properties['%s']=='%s' and feature.ma_properties['%s']=='%s'" % (_UUID_MASTER,_uuid_related,_GEN_FTYPE,'B_WB_100'))
            if (len(feat_list) == 1):
                buffer_P = feat_list[0]        
                buffer_L = LineString(buffer_P.exterior)
            else:
                everything_ok = False
                           
            if (everything_ok):
                # VALIDATE if the segment to correct goes in but also goes out the buffer.
                # If not, this is impossible to correct correctly so do not try it.
                segment_in_out = buffer_L.intersection(segment)
                #
                if ((segment_in_out.geom_type == "MultiPoint") and (len(segment_in_out) == 2)):
    
                    # SPLIT buffer line based on segment to correct
                    chunks = GenUtil.make_iterable(buffer_P.intersection(contour))

                    for chunk in chunks:
                        if (segment.intersection(chunk)):                                       # we found the segment to correct
                            parts_buffer_L = self._cut_line_using_segment(buffer_L, chunk)      # so, cut buffer line into chunks to work with
                            parts_contour = self._cut_line_using_segment(contour, chunk)        # so for contour line
                            
                            mid_chunk = parts_buffer_L[1]
                            ends_chunks = linemerge([parts_buffer_L[0], parts_buffer_L[2]])
                            
                            length_mid_chunk = mid_chunk.length
                            length_ends_chunks = ends_chunks.length
                            
                            point_origin = Point(buffer_L.coords[0])
                            dist_mid_chunk_to_origin = mid_chunk.distance(point_origin)
                            dist_ends_chunks_to_origin = ends_chunks.distance(point_origin)
                            
                            if (length_mid_chunk < length_ends_chunks):
                                id_shortest = 2     # refers to chunks[1]
                            else:
                                id_shortest = 13    # refers to chunks[0] and chunks[2]
                            
                            if ((id_shortest == 13) and (dist_ends_chunks_to_origin < dist_mid_chunk_to_origin)):
                                replacement_part = ends_chunks
                            else:
                                replacement_part = mid_chunk       
                            
                            # VERIFY constraints on the proposed replacement buffer part
                            # ... first set new segment properties to those of the segment to modify
                            replacement_part.ma_properties = segment.ma_properties
                            # ... then copy the original contour id (needed by the constraints test)
                            replacement_part.ma_properties['contour_sci_id'] = contour._sci_id
                            # ... finally, execute the constraints test
                            replacement_part_is_valid = self.check_constraints(replacement_part)                           
                            
                            if (replacement_part_is_valid):
                                # $%?&! problem with shapely which do not offer rounding functionality
                                # so, round first and last coordinates to 8 digits to ensure afterwards line merging
                                #
                                # ... round coordinates of unchanged first contour part                               
                                p1 = list(parts_contour[0].coords)                        
                                p1fX = round(p1[0][0],8)
                                p1fY = round(p1[0][1],8)
                                p1lX = round(p1[-1][0],8)
                                p1lY = round(p1[-1][1],8)
                                p1f = tuple([p1fX,p1fY])
                                p1l = tuple([p1lX,p1lY])
                                p1[0]=p1f
                                p1[-1]=p1l
                                #
                                # ... round coordinates of replacement_part                            
                                p2 = list(replacement_part.coords)
                                p2fX = round(p2[0][0],8)
                                p2fY = round(p2[0][1],8)
                                p2lX = round(p2[-1][0],8)
                                p2lY = round(p2[-1][1],8)
                                p2f = tuple([p2fX,p2fY])
                                p2l = tuple([p2lX,p2lY])
                                p2[0]=p2f
                                p2[-1]=p2l
                                #
                                #  ... round coordinates of unchanged last contour part                     
                                p3 = list(parts_contour[2].coords)
                                p3fX = round(p3[0][0],8)
                                p3fY = round(p3[0][1],8)
                                p3lX = round(p3[-1][0],8)
                                p3lY = round(p3[-1][1],8)
                                p3f = tuple([p3fX,p3fY])
                                p3l = tuple([p3lX,p3lY])
                                p3[0]=p3f
                                p3[-1]=p3l
                                
                                # test if the length of replacement part is greater than 0. If not, do not assemble nor update.
                                if (LineString(p2).length >= 0.0005):
                                    # assemble new contour parts for line merging
                                    contour_parts = []
                                    contour_parts.append(LineString(p1))
                                    contour_parts.append(LineString(p2))
                                    contour_parts.append(LineString(p3))
                                                            
                                    # merge parts into a single line                        
                                    new_contour = linemerge(contour_parts)
                                    
                                    # finally ... make sure this is a single line and ...
                                    # ... if so, update contour with his new coordinates
                                    if isinstance(new_contour, LineString):
                                        contour.update_coords(list(new_contour.coords), sc_cc)
            else:
                GenUtil.print_debug(self.params, '[%s] :: Something wrong while processing waterbodies:' % (_ALGO))
                GenUtil.print_debug(self.params, '  Shorelines container returns %d feature(s), should be 1' % (len(feat_list)))
                GenUtil.print_debug(self.params, '  segment properties: %s' % str(segment.ma_properties))
                GenUtil.print_debug(self.params, '  contour properties: %s' % str(contour.ma_properties))                    

        GenUtil.print_debug (self.params, 'END PROCESSING contour over waterbodies')   
                                                    
        return
    

    # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    def process_contour_over_watercourse (self, segments, sc_contours, sc_buffers):
        """
        """

        GenUtil.print_debug (self.params, '\nSTART PROCESSING contour over watercourses:')
        GenUtil.print_debug (self.params, '  number of segments to process: %d' % (len(segments)))       

        # SORT segments list based on their length rank (as set by FME)
        # This is done to make sure that segments are processed in the right order (shortest to longest)
        sorted(segments, key=lambda by_item: by_item.ma_properties[_SEGMENT_LENGTH_RANK])
                
        # PROCESS each segment in turn
        for segment in GenUtil.make_iterable(segments):
            
            GenUtil.print_debug(self.params, '  SEGMENT.ma_properties: %s' % str(segment.ma_properties))

            # extract segment identifiers for associated contour and single line watercourse
            _uuid_master = segment.ma_properties[_UUID_MASTER]      # contour id
            _uuid_related = segment.ma_properties[_UUID_RELATED]    # watercourse id
            _uuid_feat = segment.ma_properties[_UUID_FEAT]          # feature id
            
            GenUtil.print_debug(self.params, '    using ids :: _uuid_master: %s | _uuid_related: %s' % (_uuid_master, _uuid_related))            
            
            # get contour from spatial container
            feat_list = sc_contours.get_features(filter="feature.ma_properties['%s']=='%s'" % (_UUID_MASTER, _uuid_master))
            if (len(feat_list) == 1):        
                contour = feat_list[0]
            else:
                raise InternalError ("[%s] :: Watercourse : Contours container returns %d feature(s), should be 1" % (_ALGO, len(feat_list)))
            
            # get first waterbody's buffer from spatial container and make sure that it is a polyline 
            everything_ok = True
            feat_list = sc_buffers.get_features(filter="feature.ma_properties['%s']=='%s' and feature.ma_properties['%s']=='%s'" % (_UUID_MASTER,_uuid_related,_GEN_FTYPE,'B_SLW_100'))
            if (len(feat_list) == 1):
                buffer_P = feat_list[0]        
                buffer_L = LineString(buffer_P.exterior)
            else:
                everything_ok = False
                                    
            if (everything_ok):
                # VALIDATE if the segment to correct goes in but also goes out the buffer.
                # If not, this is impossible to correct correctly so do not try it.
                segment_in_out = buffer_L.intersection(segment)
                #
                if ((segment_in_out.geom_type == "MultiPoint") and (len(segment_in_out) == 2)):
    
                    # SPLIT buffer line based on segment to correct
                    chunks = GenUtil.make_iterable(buffer_P.intersection(contour))
    
                    for chunk in chunks:
                        if (segment.intersection(chunk)):                                       # we found the segment to correct
                            parts_buffer_L = self._cut_line_using_segment(buffer_L, chunk)      # so, cut buffer line into chunks to work with
                            parts_contour = self._cut_line_using_segment(contour, chunk)        # so for contour line
                            
                            mid_chunk = parts_buffer_L[1]
                            ends_chunks = linemerge([parts_buffer_L[0], parts_buffer_L[2]])
                            
                            length_mid_chunk = mid_chunk.length
                            length_ends_chunks = ends_chunks.length
                            
                            point_origin = Point(buffer_L.coords[0])
                            dist_mid_chunk_to_origin = mid_chunk.distance(point_origin)
                            dist_ends_chunks_to_origin = ends_chunks.distance(point_origin)
                            
                            if (length_mid_chunk < length_ends_chunks):
                                id_shortest = 2     # refers to chunks[1]
                            else:
                                id_shortest = 13    # refers to chunks[0] and chunks[2]
                            
                            if ((id_shortest == 13) and (dist_ends_chunks_to_origin < dist_mid_chunk_to_origin)):
                                replacement_part = ends_chunks
                            else:
                                replacement_part = mid_chunk       
    
                            # VERIFY constraints on the proposed replacement buffer part
                            # ... first set new segment properties to those of the segment to modify
                            replacement_part.ma_properties = segment.ma_properties
                            # ... then copy the original contour id (needed by the constraints test)
                            replacement_part.ma_properties['contour_sci_id'] = contour._sci_id
                            # ... finally, execute the constraints test
                            replacement_part_is_valid = self.check_constraints(replacement_part)                           
                            
                            if (replacement_part_is_valid):
                                # $%?&! problem with shapely which do not offer rounding functionality
                                # so, round first and last coordinates to 7 digits to ensure afterwards line merging
                                #
                                # ... round coordinates of unchanged first contour part                               
                                p1 = list(parts_contour[0].coords)                        
                                p1fX = round(p1[0][0],7)
                                p1fY = round(p1[0][1],7)
                                p1lX = round(p1[-1][0],7)
                                p1lY = round(p1[-1][1],7)
                                p1f = tuple([p1fX,p1fY])
                                p1l = tuple([p1lX,p1lY])
                                p1[0]=p1f
                                p1[-1]=p1l
                                #
                                # ... round coordinates of replacement_part                            
                                p2 = list(replacement_part.coords)
                                p2fX = round(p2[0][0],7)
                                p2fY = round(p2[0][1],7)
                                p2lX = round(p2[-1][0],7)
                                p2lY = round(p2[-1][1],7)
                                p2f = tuple([p2fX,p2fY])
                                p2l = tuple([p2lX,p2lY])
                                p2[0]=p2f
                                p2[-1]=p2l
                                #
                                #  ... round coordinates of unchanged last contour part                     
                                p3 = list(parts_contour[2].coords)
                                p3fX = round(p3[0][0],7)
                                p3fY = round(p3[0][1],7)
                                p3lX = round(p3[-1][0],7)
                                p3lY = round(p3[-1][1],7)
                                p3f = tuple([p3fX,p3fY])
                                p3l = tuple([p3lX,p3lY])
                                p3[0]=p3f
                                p3[-1]=p3l
                                
                                if (LineString(p2).length >= 0.0005):
                                    # assemble new contour parts for line merging
                                    contour_parts = []
                                    contour_parts.append(LineString(p1))
                                    contour_parts.append(LineString(p2))
                                    contour_parts.append(LineString(p3))
                                                            
                                    # merge parts into a single line                        
                                    new_contour = linemerge(contour_parts)
                                    
                                    # finally ... make sure this is a single line and ...
                                    # ... if so, update contour with his new coordinates
                                    if isinstance(new_contour, LineString):
                                        contour.update_coords(list(new_contour.coords), sc_contours)
            else:
                GenUtil.print_debug(self.params, '[%s] :: Something wrong while processing watercourses:' % (_ALGO))
                GenUtil.print_debug(self.params, '  Shorelines container returns %d feature(s), should be 1' % (len(feat_list)))
                GenUtil.print_debug(self.params, '  segment properties: %s' % str(segment.ma_properties))
                GenUtil.print_debug(self.params, '  contour properties: %s' % str(contour.ma_properties))                    

        GenUtil.print_debug (self.params, 'END PROCESSING contour over watercourses')       
                                                    
        return


    # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    def check_constraints(self, new_seg):
        """
        Check if segment is topologically valid.
        
        *Parameters*:
            - old_seg : (LineString) the segment to verify
            - new_seg : (LineString) the proposed new segment that must be verified
            
        *Returns*:
            - Flag indicating if segment is valid (True) or no (False)
        """

        mother_id = new_seg.ma_properties['contour_sci_id']

        line_simple_line = new_seg
        
        # use a 'buffer' LineString around the new segment to make sure to trap any other feature. 
        line_crossing_line = LineString(list(new_seg.buffer(0.1,2).exterior.coords))
                
        polygon_sidedness = None
        conflict_type = GenUtil.test_constraints (self, None, line_simple_line, line_crossing_line, 
                                                  polygon_sidedness, self.sc_CC, mother_id)

        if (conflict_type is not None):
            GenUtil.print_debug (self.params, '  conflict_type : %s' % (conflict_type))
            
        return (conflict_type is None)
    
        
    # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    def process(self):
        """Main routine of the algorithm

        """

        GenUtil.print_debug (self.params, "START of %s algorithm)" % (_ALGO))
        GenUtil.print_debug (self.params, "Parameters description:")
        GenUtil.print_debug (self.params, "  - Check constraints:")
        GenUtil.print_debug (self.params, "    - sidedness      : %s" % (self.params.test_sidedness))
        GenUtil.print_debug (self.params, "    - crossing line  : %s" % (self.params.test_crossing_line))
        GenUtil.print_debug (self.params, "    - simple line    : %s" % (self.params.test_simple_line))
        GenUtil.print_debug (self.params, "  - Debug status     : %s" % (self.params.debug))

        # Check the feature's class and attributes
        self.check_features()

        # Load the shapely features into the spatial container
        (self.sc_CC, self.sc_CA, self.sc_SHR, self.sc_SLW, self.sc_B_SLW, self.sc_S_SLW, self.sc_S_SLWX, self.sc_WB, self.sc_B_WB, self.sc_S_WB, self.sc_S_WBX) = self._load_features()
       
        GenUtil.print_debug (self.params, '\nSTART Spatial Containers counts:')
        GenUtil.print_debug (self.params, '  sc_CC    : %d' % len(self.sc_CC.get_features()))
        GenUtil.print_debug (self.params, '  sc_CA    : %d' % len(self.sc_CA.get_features()))
        GenUtil.print_debug (self.params, '  sc_SHR   : %d' % len(self.sc_SHR.get_features()))
        GenUtil.print_debug (self.params, '  sc_SLW   : %d' % len(self.sc_SLW.get_features()))
        GenUtil.print_debug (self.params, '  sc_B_SLW : %d' % len(self.sc_B_SLW.get_features()))
        GenUtil.print_debug (self.params, '  sc_S_SLW : %d' % len(self.sc_S_SLW.get_features()))
        GenUtil.print_debug (self.params, '  sc_S_SLWX: %d' % len(self.sc_S_SLWX.get_features()))
        GenUtil.print_debug (self.params, '  sc_WB    : %d' % len(self.sc_WB.get_features()))
        GenUtil.print_debug (self.params, '  sc_B_WB  : %d' % len(self.sc_B_WB.get_features()))
        GenUtil.print_debug (self.params, '  sc_S_WB  : %d' % len(self.sc_S_WB.get_features()))
        GenUtil.print_debug (self.params, '  sc_S_WBX : %d' % len(self.sc_S_WBX.get_features()))
        GenUtil.print_debug (self.params, 'END Spatial Containers counts')

        #- MAIN PROCESSING IS HERE -----------------------
        # First, PROCESS contour segments over waterbodies
        segments = self.sc_S_WB.get_features()      # get list of segments to correct
        self._process_contour_over_waterbody(segments, self.sc_CC, self.sc_B_WB)

        # Then, PROCESS contour segments over single line watercourses
        segments = self.sc_S_SLW.get_features()      # get list of segments to correct
        self.process_contour_over_watercourse(segments, self.sc_CC, self.sc_B_SLW)

        self.extract_features_out(self.sc_CC.get_features())        
        GenUtil.print_debug (self.params, "End of %s" % (_ALGO))