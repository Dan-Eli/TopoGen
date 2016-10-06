#!/usr/local/bin/python
# -=- encoding: utf-8 -=-
#####################################################################################################################################

"""
    This algorithm allows to detect bottlenecks into a polygon.

    In order to detect bottlenecks the algorithm receive in input one polygon and the triangles created on the polygon.  The algorithm
    is not modifying any input geometry but it will add a property ("bottleneck") to each input triangle.  If a triangle is part of a
    bottleneck it will set the value to 1 otherwise it will assign the value 0.
    
    The algorithm has 2 parameters:
        - width: The width of the bottleneck to detect
        - perimeter: In order to avoid the detection of false bottleneck which are usually triangles in the extremity of the polygon, 
                     we verify the perimeter by one of the areas form by clipping the original area into 2 sub polygon is lower 
                     than than the value of this parameter.  By default the value of this parameter is equal to the parameter 
                     "width" multiply by 4.
                     
    The polygon entities entering this algorithm must each have 2 properties:
        code: Identify the type of polygon. Can take 2 values
              T: Polygone of type triangle
              P: Polygon
        link_id: The identifier that links one and only one polygon to [1..n] triangles.

    Usage:
        import detect_bottleneck

    Limits and constraints
        A program calculating the triangles must have been executed before executing this algorithm 


"""
__revision__ = "--REVISION-- : $Id: algo_detect_bottleneck.py 312 2011-05-09 19:35:40Z dpilon $"

#####################################################################################################################################

from collections import Iterable
from copy import deepcopy

from shapely.geometry import Point, MultiPoint, LineString, MultiLineString, Polygon
from shapely.ops import cascaded_union
from shapely.ops import linemerge

from lib_genmetal import MA_Point, MA_LineString, MA_Polygon, GenUtil, Holder, \
                         GenStatistics, Parameters, SpatialContainer, Algorithm, \
                         PolygonModifier, ChordalAxisTransformer, PolygonMerger

########################################################

# Public constant

# Properties definition

# Define attribute name
_ACTION = 'action'
_BUFFER = 'buffer'
_CODE = 'code'
_CENTER_LINE = 'center_line'
_LENGTH = 'length'
_POLYGON = 'polygon'
_CAT_ID = 'cat_id'
_WEIGHTED_WIDTH = 'weighted_width'
_WIDTH = 'width'

# Define attribute value
_BOTTLENECK = 'bottleneck'
_TO_BE_DELETED = 'to_be_deleted'
_TO_BE_BUFFERED = 'to_be_buffered'

# Define statistics constant
_CONFLICT = "conflict"
_ELIMINATED = "eliminated"
_WIDEN = "widen"

#Internal key word constants
_ALGO = 'Detect bottleneck'


_SEARCH_TOLERANCE = .001 # Search tolerance in meter for the RTree search



                
class Statistics(GenStatistics):
    """Class that contains the statistics for this algorithm

    Attributes
        stat_names: Name of the statistics for the SpikeStatistics class. These name are
                    used by the Statistics class

    """

    def __init__(self):
        """Initialize the attributes of an object of the class"""

        GenStatistics.__init__(self)
        self.stats_names = ((_ALGO,_WIDEN, _ELIMINATED, _CONFLICT))

    def get_stats (self, type=GenStatistics.SUMMARY):
        """Extract the current statistics and build  a list of string that forms the statistical message"

        Parameters:
            type: Give the form of statistics to extract. Can take 2 values.
                SUMMARY: Summary information
                DETAILED: Detailed information

        """

        str_out = []
        str_out.append( "Correct bottleneck statistics" )

        str_out.append( "----------------------------" )
        str_out.append("Number of widen bottlenecks : " + str(self.get_stats_name_count_total( _WIDEN)))
        str_out.append("Number of eliminated bottlenecks : " + str(self.get_stats_name_count_total( _ELIMINATED)))
        str_out.append("Number of conflictual bottlenecks : " + str(self.get_stats_name_count_total( _CONFLICT)))
        str_out.append( "----------------------------" )

        return str_out
    
class AlgoCorrectBottleneck(Algorithm):
    """This is the main class for this algorithm

    """

    def __init__(self, critical_width=20, minimal_width=25, minimal_length=75, test_sidedness=True, test_crossing_line=True, 
                 width_factor=3. ):
        """Initialize the attributes of an object of the class DPAlgorithm

        Parameters:
            critical_width: First width tolerance for detecting bottleneck
            minimal_width: Second width tolerance for detecting bottleneck. Minimal width always greater than critical width
            minimal_length: Treshold length of bottlenecks between buffer and delete bottleneck
            test_sidedness: Flag to enable or disable the topology sidedness testing
            test_crossing_line: Flag to enable or disable the topology crossing testing 
            width_factor: The width factor is multiplied by the with in order to obtain the false bottleneck perimeter length. This length is used
                          to detect False bottleneck. False bottleneck are triangles in the extremity  of a polygon.  
                          By default the width factor is equal 4. 
            debug: Flag to enable (True) or diable (False) the output of extra information

        Return value:
            None

        """

        Algorithm.__init__(self)
        
        self.params = Parameters()
        
        self.params.critical_width = critical_width
        self.params.minimal_width = minimal_width
        self.params.minimal_length = minimal_length
        self.params.test_sidedness = test_sidedness
        self.params.test_crossing_line = test_crossing_line
        
        self.params.perimeter_tol = minimal_width * width_factor
        self.params.buffer_width = minimal_width/2.
        
        self.params.debug = False
        
        self.stats = Statistics()
        self.stats.add_iteration()
        self.polygon_id = -1

    def _check_features(self):
        """
        Check if the features passed in parameters are of the good class type and have the good attributes
        
        Parameters: None
        
        Return value: None
        
        """
        
        # Check the line string
        properties_to_check = [_CAT_ID]
        for feature in self.features:
            if isinstance(feature, MA_LineString):
                GenUtil.check_feature_integrity(feature, MA_LineString, properties_to_check)
            else:
                GenUtil.check_feature_integrity(feature, MA_Polygon, properties_to_check)
                
    def _classify_cat_id(self):
        """Classify all the input feature according to the CAT_ID
        
        For the same LINK_ID we must find one and only one polygon (MA_Polygon) and one or more triangle (MA_LineString) 
        otherwise there is an error
        
        Parameters: None
        
        Return value:
            Dictionary of link ID's
        
        """
        
        cat_ids = {}
        
        # Pass through each feature and create
        for feature in self.features:
            cat_id = feature.ma_properties[_CAT_ID]
            if not cat_ids.has_key(cat_id):
                # This LINK_ID was never used create a holder
                holder = Holder(polygon=[], triangles=[], cat_id=cat_id)
                cat_ids[cat_id] = holder
            # Add the feature
            if (isinstance(feature, MA_LineString)):
                cat_ids[cat_id].triangles.append(feature)
            else:
                cat_ids[cat_id].polygon.append(feature)
                
        # Pass through each item of the dictionary to check the condition
        for holder in cat_ids.itervalues():
            if (len(holder.polygon) == 1):
                holder.polygon = holder.polygon[0]
            else:
                raise IntegrityError ("There is no polygon for CAT_ID: %i" %cat_id)
            
            if (len(holder.triangles) == 0 ): 
                raise IntegrityError ("There are no triangles for CAT_ID: %i" %cat_id)
            
        return cat_ids

    def _build_bounding_box(self, tolerance, coord):
        """Adjust the bounding box
        
        Parameters: 
            line: LineString to use
            coord: (x,y) tuple
            
         Return value: bounding box tuple (xmin, ymin, xmax, ymax)  
            
        """
        
        xmin = coord[0] - tolerance
        ymin = coord[1] - tolerance
        xmax = coord[0] + tolerance
        ymax = coord[1] + tolerance
        
        return (xmin, ymin, xmax, ymax) 
   
    def _classify_bottleneck_areas(self, bottleneck_areas):
        """Classify the bottleneck area as "to be buffered" or "to be deleted"
        
           To be deleted a bottleneck must must the following criteria:
              - The weightd width of the bottleneck must be below the critical width
              - The length of the bottleneck must be over the minimal length
                
           To be buffered must meet the following criteria
               - All the bottleneck that do not meet the "to be deleted" 2 criterias
                 are considered to be buffered
                 
           - The weighted width is a better approximation of the bottlenecks than if we were taking the simple width of the 
             triangle
           - For one triangle the weigthed width is the width multiply by the length of the triangle.
           - For the bottleneck area the total weigthed width is the sum of all the weighted width of the triangles forming
             this bottleneck.
           
        """
        
        # Loop over each bottleneck area to classify them
        for bottleneck_area in bottleneck_areas:
            # Calculate the total weighted width of the bottleneck
            lst_weighted_width = bottleneck_area.ma_properties[_WEIGHTED_WIDTH]
            sum_weighted_width = sum(lst_weighted_width)

            # Calculate the length of the bottleneck
            lst_lst_center_line = bottleneck_area.ma_properties[_CENTER_LINE] # It's a list of list0
            
            # Tranform a list of list into a simple list of center lines
            center_lines = []
            for lst_center_line in lst_lst_center_line:
                for center_line in lst_center_line:
                    center_lines.append(center_line)
            
            # Merge all the certer line as one multi line string
            center_line = linemerge(center_lines) 
            center_lines = GenUtil.make_iterable(center_line)
            
            
            # Creates the buffer that will be used if the line needs to be buffered
            lst_buffer = []
            sum_length = 0
            buffer_width = self.params.minimal_width/2.
            min_buffer_width = buffer_width * 0.75
            for center_line in center_lines:
                sum_length += center_line.length
                # Buffer creation but we never creates a buffer for the full length of the center line
                # becasue the resulting in nice to see... it is to much buffereing
                if center_line.length < self.params.minimal_width:
                    # If the length of the center line is below the minimal width
                    # We create a buffer around the mid point of the center line
                    mid_point = center_line.interpolate(0.5, normalized=True)
                    buffer = mid_point.buffer(buffer_width, resolution=5)
                else:
                    # We chop a portion of the line at the start and end of the line
                    first_split = buffer_width
                    distance_split = center_line.length - buffer_width
                    lst_lines = GenUtil.cut_line_distance(center_line, distance_split) # Chop the start
                    line = lst_lines[0]
                    distance_split = buffer_width
                    lst_lines = GenUtil.cut_line_distance(line, distance_split) # Chop the end
                    line = lst_lines[1]
                    # Buffer the remaining of the line
                    buffer = line.buffer(buffer_width, resolution=5)
                lst_buffer.append(buffer)
                # The min buffer is used to reduce the soften effect of the main buffer of the resulting buffered polygon
                min_buffer = center_line.buffer(min_buffer_width, resolution=5)
                lst_buffer.append(min_buffer)

            # Creates the final buffer
            buffer = cascaded_union(lst_buffer)

            # Add the properies of the bottleneck
            bottleneck_area.ma_properties[_LENGTH] = sum_length
            mean_weighted_width = sum_weighted_width / sum_length
            bottleneck_area.ma_properties[_WEIGHTED_WIDTH] = mean_weighted_width
            bottleneck_area.ma_properties[_BUFFER] = buffer
            
            # Apply the rules to determiner if the bottle is to be deleted or be buffered
            if (sum_length < self.params.minimal_length):
                # The bottleneck is to be buffered
                action = _TO_BE_BUFFERED 
            else:
                if (mean_weighted_width < self.params.critical_width):
                    # The bottleneck is to be deleted
                    action = _TO_BE_DELETED
                else:
                    # The bottleneck is to be buffered
                    action = _TO_BE_BUFFERED
                                    
            bottleneck_area.ma_properties[_ACTION] = action
                    
    def attribute_to_keep(self, ma_properties, lst_attribute):
        """Keep only the attributes in the list of attributes all the others are deleted
        
        Parameters: 
            ma_properties: Dictionary of attributes
            lst_attributes: List of attribute names to keep
        
        Return value: None
        
        """
        
        del_property = []
        
        # creates a list of properties to delete
        for property in ma_properties:
            if property not in lst_attribute:
                del_property.append(property)
                
        # Delete the properties
        for property in del_property:
            del ma_properties[property]
        
    def _output_features(self, total_bottleneck_areas, total_multi_center_line, lst_modified_polygon ):
        """Output features by placing them into self.features 
        
        Parameters: 
            - total_bottleneck_areas: List of MA_Polygon representing bottleneck
            - total_multi_center_line: List of LineString and MultiLineString representing the center of each polygon
            - lst_modified_polygon: List of MA_Polygon representing the modified original polygons
        
        Return values: None
        
        """
        
        # Add an attribute code and output all the bottleneck area
        for bottleneck_area in total_bottleneck_areas:
            bottleneck_area.ma_properties[_CODE] = _BOTTLENECK
            self.attribute_to_keep(bottleneck_area.ma_properties, [_CODE, _ACTION, _WEIGHTED_WIDTH, _LENGTH])
            self.features.append(bottleneck_area)
                
###        print "**** Len total_center_line: ", len(total_multi_center_line)
        if (len(total_multi_center_line) >= 1):
###            multi_center_line = linemerge(total_center_line)
###            lst_center_line = GenUtil.make_iterable(multi_center_line)
            # Add an attribute code and output each center line
            for multi_center_line in total_multi_center_line:
                iter_multi_center_line = GenUtil.make_iterable(multi_center_line)
###                print "iter multi center line: ", iter_multi_center_line
                for center_line in iter_multi_center_line:
###                    print "center line: ", center_line
                    ma_center_line = MA_LineString(list(center_line.coords), dual=None)
                    ma_center_line.ma_properties[_CODE] = _CENTER_LINE
                    self.attribute_to_keep(ma_center_line.ma_properties, [_CODE])
                    self.features.append(ma_center_line)
                    
        # Add an attribute code and output all the modified polygons
        for modified_polygon in lst_modified_polygon:
            modified_polygon.ma_properties[_CODE] = _POLYGON
            self.features.append(modified_polygon)

    
    def process(self):
        """Main routine of the algorithm
        
        """

        GenUtil.print_debug (self.params, "Start of %s algorithm" %_ALGO)

        
        # Check attributes and geometry type
        self._check_features()
        
        cat_ids = self._classify_cat_id()
        self.features = []

        if (self.params.debug):
            #Test if print is needed before scanning the s_container for nothing... waste of time...
            nbr_lines = len(self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING")) 
            nbr_points = len(self.s_container.get_features(filter="feature.feature_type==GenUtil.POINT")) 
            GenUtil.print_debug (self.params, "Number of lines imported: %s"  %(nbr_lines) )
            GenUtil.print_debug (self.params, "Number of points imported: %s"  %(nbr_points) )
        
        # Creates a list of tuple containing the cat_id and the polygon
        lst_tuple_cat_id_polygon = []
        for holder in cat_ids.itervalues():
            tuple_cat_id_polygon = (holder.cat_id, holder.polygon)
            lst_tuple_cat_id_polygon.append(tuple_cat_id_polygon)
        
        # The object PolygonModifier is used to modify a list of polygon and check the topology 
        polygon_modifier = PolygonModifier(lst_tuple_cat_id_polygon,self.params.test_sidedness, 
                                           self.params.test_crossing_line)
                      

        total_bottleneck_areas = []
        total_multi_center_line = []
        
        # Process each cat_id or polygon independently
        i_count = 0
        for holder in cat_ids.itervalues():
###            print "Processing id: ", i_count
            i_count += 1
            cat_id = holder.cat_id
            chordal_axis = ChordalAxisTransformer(holder.polygon, holder.triangles, self.params.minimal_width)
            holder.polygon = None
            holder.triangles = None
            triangles = chordal_axis.get_triangles()
            
            # From the triangles extract the center line
            lst_center_line = []
            for triangle in triangles:
                lst_center_line += triangle.ma_properties[_CENTER_LINE]
###            print "lst_center_line: ", len(lst_center_line)
            if ( len(lst_center_line) >= 1 ):
                multi_center_line = linemerge(lst_center_line)
                total_multi_center_line.append(multi_center_line)
                
            # From the triangles keep only the triangles that form bottleneck regions
            triangles = [triangle for triangle in triangles if triangle.ma_properties[_CODE] == ChordalAxisTransformer.BOTTLENECK]
            
            pol_triangles = []
            for triangle in triangles:
                # Create a polygon from the triangle which is a MA_LineString
                pol_triangle = MA_Polygon(triangle.coords_dual[0:3])
                # Add some properties to the MA_Polygon
                pol_triangle.ma_properties[_WIDTH] = triangle.ma_properties[_WIDTH]
                center_line = triangle.ma_properties[_CENTER_LINE]
                pol_triangle.ma_properties[_CENTER_LINE] = center_line
                lst_length = [line.length for line in center_line]
                length = sum(lst_length)
                # Add the weighted width of the triangle
                pol_triangle.ma_properties[_WEIGHTED_WIDTH] = pol_triangle.ma_properties[_WIDTH] * length 
                pol_triangles.append(pol_triangle)
            
            # Merge the bottleneck polygons and create a list of the attribute merged
            polygon_merger = PolygonMerger(pol_triangles)
            bottleneck_areas = polygon_merger.getMerged([_WIDTH, _WEIGHTED_WIDTH, _CENTER_LINE])
                        
            # Classify the bottleneck as to be delete or to be buffered
            self._classify_bottleneck_areas(bottleneck_areas)
            
            if (len(bottleneck_areas) != 0):
                # Apply the bottleneck to the original polygon
                polygon_modifier.modify(cat_id, bottleneck_areas)            
                # Keep a copy of the bottleneck
                total_bottleneck_areas += bottleneck_areas
            
        # Extract the polygons from the polygon modifier
        lst_modified_polygon = polygon_modifier.get_polygons()
        for i in xrange(polygon_modifier.nbr_conflict): 
            self.stats.add_stats(_CONFLICT)
        for i in xrange(polygon_modifier.nbr_eliminated): 
            self.stats.add_stats(_ELIMINATED)
        for i in xrange(polygon_modifier.nbr_widen): 
            self.stats.add_stats(_WIDEN)
            
        polygon_modifier = None
        
        # Output features
        self._output_features(total_bottleneck_areas, total_multi_center_line, lst_modified_polygon )
                 
        GenUtil.print_debug (self.params, "End of %s" %(_ALGO))