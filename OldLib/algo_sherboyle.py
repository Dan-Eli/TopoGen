#!/usr/local/bin/python
# -=- encoding: utf-8 -=-

#####################################################################################################################################

"""
This algorithm implements an adapted version of the Boyle algorithm to smooth line features.
    
Usage:
    import algo_sherboyle OR
    from algo_sherboyle import AlgoSherBoyle

Limits and constraints
    None

__revision__ = "--REVISION-- : $Id: algo_sherboyle.py 128 2010-12-21 16:53:31Z langlois $"

"""

#####################################################################################################################################

from shapely.geometry import Point, LineString, Polygon

from lib_genmetal import MA_Point, MA_LineString, MA_Polygon, InvalidParameterError, GenUtil,\
                         GenStatistics, PointErrorPosition, LineStringErrorPosition, IterationResults, Parameters,\
                         SpatialContainer, Algorithm

########################################################

#Internal key word constants
_ALGO = 'SherBoyle'
_FORWARD_LOOK = 'forward_look'


class Statistics(GenStatistics):
    """Maintains statistics for the current algorithm

    Attributes
        stat_names: Dictionary of statistics counters names maintained by Statistics class. These counters
                    are then incremented by calls made to the Statistics class

    """

    def __init__(self):
        """Initialize the attributes of an object of the class"""

        GenStatistics.__init__(self)
        self.stats_names = ((_ALGO,'TOTAL_PRE','TOTAL_POST','SMOOTHED','UNTOUCHED',GenUtil.INVALID,GenUtil.CROSSING_LINE,GenUtil.SIDEDNESS,GenUtil.SIMPLE_LINE,))
        
    def get_stats(self, type=GenStatistics.SUMMARY):
        """Extract the current statistics and create a list of strings that builds the statistical message"

        Parameters:
            type: Always set to a summary of statistics

        """

        title = "%s Algorithm Statistics" % (_ALGO)
        str_out = []
        str_out.append("=" * len(title))
        str_out.append("%s Algorithm Statistics" % (_ALGO))
        str_out.append("-" * len(title))

        str_out.append("Features IN  : " + str(self.get_stats_name_count_total('TOTAL_PRE')))
        str_out.append("Features OUT : " + str(self.get_stats_name_count_total('TOTAL_POST')))
        str_out.append("-" * len(title))
        str_out.append("   BOYLED Features    :  " + str(self.get_stats_name_count_total('SMOOTHED')))
        str_out.append("   UNTOUCHED Features :  " + str(self.get_stats_name_count_total('UNTOUCHED')))
        str_out.append("   SIDEDNESS          :  " + str(self.get_stats_name_count_total(GenUtil.SIDEDNESS)))
        str_out.append("   SIMPLE_LINE        :  " + str(self.get_stats_name_count_total(GenUtil.SIMPLE_LINE)))
        str_out.append("   CROSSING_LINE      :  " + str(self.get_stats_name_count_total(GenUtil.CROSSING_LINE)))
        str_out.append("=" * len(title))

        return str_out


class AlgoSherBoyle(Algorithm):
    """
    This is the main class for the algorithm

    Attributes:
        - params: Parameters for the algorithm
    """

    def __init__(self, move_ends = False,
                       test_crossing_line=False,
                       test_simple_line=False, 
                       test_sidedness=False, 
                       keep_iterative_results=False, 
                       debug=False):
        """Initialize the attributes of an object of this class

        Parameters:
            test_crossing_line: Flag to enable(True)/disable(False) checking of the line crossing constraint
            test_simple_line: Flag to enable(True)/disable(False) checking of the simple line constraint
            test_sidedness: Flag to enable(True)/disable(False) checking of the sidedness constraint
            keep_iterative_results: Flag to enable(True)/disable(False) keeping of the iterative results
            debug: Flag to enable(True)/disable(False) debug output

        Return value:
            None

        """

        Algorithm.__init__(self)
        
        self.params = Parameters()
        self.params.move_ends = move_ends
        self.params.test_crossing_line = test_crossing_line
        self.params.test_simple_line = test_simple_line 
        self.params.test_sidedness = test_sidedness 
        self.params.keep_iterative_results = keep_iterative_results 
        self.params.debug = debug
        
        self.stats = Statistics()


    def check_features(self):
        """
        Check if the features passed in parameters are of the good class type and have the good attributes
        
        Parameters: None
        
        Return value: None
        
        """
        
        # Check the line string
        class_type = MA_LineString
        properties_to_check = ["forward_look"]
        for feature in self.features:
            if isinstance(feature, MA_LineString):
                GenUtil.check_feature_integrity(feature, MA_LineString, properties_to_check)        

    
    def _set_line_attributes (self):
        """This routine sets the attributes to the line

        The routine checks:
            - if the first/last vertices are the same the line is closed; otherwise it is open
            - if the line has less than 3 vertices it is at its simplest form; otherwise it is not

        Parameters:
            s_container: Spatial container containing the spatial features and the spatial index

        Return value: None

        """

        # Select all and only lines
        for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING"):
            self.stats.add_stats('TOTAL_PRE')
            line.is_simplest = False
            line.is_closed = False
            line.is_boyled = False
            
            # Verify line closure
            if (len(line.coords_dual) >= 3):
                if (GenUtil.distance(line.coords_dual[0], line.coords_dual[-1]) <= GenUtil.ZERO):
                    line.is_closed = True
                    
            # Verify if line is too short to work with. Line should be 'boyled' only if it has at least forward_look + 2 vertices
            if (len(line.coords_dual) < line.ma_properties[_FORWARD_LOOK] + 2):
                line.is_simplest = True
                
            if (line.is_simplest):
                self.stats.add_stats('UNTOUCHED')
                
        return
           

    def sherboyle_this_feature(self, line):
        """To be completed"""
                       
        move_ends = self.params.move_ends
        test_sidedness = self.params.test_sidedness
        test_simple_line = self.params.test_simple_line
        test_crossing_line = self.params.test_crossing_line
        
        # Extract feature coordinates
        featureCoords = line.coords_dual

        # Calculate ratio value for linear interpolation
        ratio = 1.0/float(line.ma_properties[_FORWARD_LOOK])

        # Duplicate current feature coords into new structure to accelerate further calculation speed
        current_coords = list(featureCoords)
        nb_coords = len(current_coords)
        
        # Prepare structure to maintain new coords and append first vertex to it
        new_coords = []
        new_coords.append(current_coords[0])

        for i in range(0,nb_coords-2):

            # capture coordinates of current vertex
            vertex_coords = new_coords[i]
            vertex_x = vertex_coords[0]
            vertex_y = vertex_coords[1]

            # identify forward_look vertex and capture its coordinates
            la_vertex = i + line.ma_properties[_FORWARD_LOOK]
            if (la_vertex > nb_coords - 1):
                la_vertex = nb_coords - 1
            la_vertex_coords = current_coords[la_vertex]
            la_vertex_x = la_vertex_coords[0]
            la_vertex_y = la_vertex_coords[1]

            # interpolate new coordinates
            new_x = ((la_vertex_x - vertex_x) * ratio) + vertex_x
            new_y = ((la_vertex_y - vertex_y) * ratio) + vertex_y

            # insert new coordinates into new structure
            new_coords.append((new_x, new_y))

        # Append last line vertex
        new_coords.append(current_coords[-1])
        
        # Test if new line respects user specified constraints and reject modification if not     
        line_simple_line = LineString(new_coords)
        line_crossing_line = LineString(new_coords)
        polygon_sideness = GenUtil.calculate_sidedness_polygon(LineString(new_coords), LineString(line.coords_dual))

        conflict_type = GenUtil.test_constraints (self, None, line_simple_line, line_crossing_line, 
                                                  polygon_sideness, self.s_container, line._sci_id)
        
        if (conflict_type is None):
            self.stats.add_stats(_ALGO)
            line.is_boyled = True
            line.update_coords(new_coords, self.s_container)
        else:
            line.is_simplest = False
        
        return
    

    def process(self):

        GenUtil.print_debug (self.params, "Start of %s simplification algorithm" % (_ALGO))
        GenUtil.print_debug (self.params, "Parameter description:")
        GenUtil.print_debug (self.params, "  - Simple line constraint: %s" % (self.params.test_simple_line))
        GenUtil.print_debug (self.params, "  - Crossing line constraint: %s" % (self.params.test_crossing_line))
        GenUtil.print_debug (self.params, "  - Sidedness constraint: %s" % (self.params.test_sidedness))
        GenUtil.print_debug (self.params, "  - Keep iterative results: %s" % (self.params.keep_iterative_results))

        # Start an iteration 'container' with its own statistics
        self.stats.add_iteration()
        
        # Check the feature's class and attributes
        self.check_features()
        
        # Load the shapely features into the spatial container
        self.s_container = self.load_features()

        if (self.params.debug):
            nbr_lines = len(self.s_container.get_features(filter="feature.feature_type=GenUtil.LINE_STRING")) 
            nbr_points = len(self.s_container.get_features(filter="feature.feature_type==GenUtil.POINT")) 
            GenUtil.print_debug (self.params, "Number of lines imported: %s" % (nbr_lines))
            GenUtil.print_debug (self.params, "Number of points imported: %s" % (nbr_points))
        
        self._set_line_attributes()
                 
        GenUtil.print_debug (self.params, 'Start of iteration # %s' % (self.stats.get_nbr_iteration()))
            
        # Iterate through all lines which are not at their simplest state
        line_simplified = True
        
        while (line_simplified):
            
            self.stats.add_iteration()
            GenUtil.print_debug (self.params, 'Start of iteration # %s' % (self.stats.get_nbr_iteration()))
            # At each new iteration we reset all the error position (we don't have to keep the old error position) 
            self.out_error_positions = []  
            line_simplified = False
    
            for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING and not feature.is_simplest"):
                self.sherboyle_this_feature(line)
            
                if (line.is_boyled):
                    self.stats.add_stats('SMOOTHED')
                else:
                    self.stats.add_stats('UNTOUCHED')

            nb_features = len(self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING"))
            self.stats.add_stats('TOTAL_POST', nb_features)

            if (self.params.keep_iterative_results):
                self.iter_results.add_iteration()
                self.iter_results.add_features(self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING"))
                                                         
        GenUtil.print_debug (self.params, 'End of iteration # %s' % (self.stats.get_nbr_iteration()))

        self.extract_features_out(self.s_container.get_features())        
        GenUtil.print_debug (self.params, "End of %s" % (_ALGO))        
        #self.stats.print_stats()
