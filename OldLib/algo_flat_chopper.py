#!/usr/local/bin/python
# -=- encoding: utf-8 -=-

#####################################################################################################################################

"""
This algorithm implements the FlatChopper algorithm to chop line features into segments.

This algorithm differs from FME as it tries to intelligently chop the line in segments where 
both end points are situated on a straight portion of the line. User can specify a fixed number 
of segments to chop in or a theoretical number of vertices on each segment. In this latter case, 
the algorithm uses a window of vertices to find the straighter part where to chop.
    
Usage:
    import algo_flat_chopper  OR
    from algo_flat_chopper import AlgoFlatChopper

Limits and constraints
    None

__revision__ = "--REVISION-- : $Id: flatChopper.py 128 2010-12-21 16:53:31Z langlois $"

"""

#####################################################################################################################################

from lib_genmetal import MA_Point, MA_LineString, MA_Polygon, InvalidParameterError, GenUtil,\
                         GenStatistics, PointErrorPosition, LineStringErrorPosition, IterationResults, Parameters,\
                         SpatialContainer, Algorithm

########################################################

#Internal key word constants
_ALGO = 'FlatChopper'
_NB_ELEMENTS = 'nb_elements'
_NB_VERTICE = 'nb_vertice'
_IS_CHOPPED = 'isChopped'
_NB_CHOPS = 'nbChops'
_SEGMENT_INDEX_ATTRIBUTE = 'segment_index_attribute'

# If psyco package is available, use it to speed up processing (2x)
try:
    import psyco
    psyco.full()
except ImportError:
    pass

class Statistics(GenStatistics):
    """Maintains statistics for the current algorithm

    Attributes
        stat_names: Dictionary of statistics counters names maintained by Statistics class. These counters
                    are then incremented by calls made to the Statistics class

    """

    def __init__(self):
        """Initialize the attributes of an object of the class"""

        GenStatistics.__init__(self)
        self.stats_names = ((_ALGO, 'TOTAL_PRE', 'TOTAL_POST', 'UNTOUCHED', 'CHOPPED',))
        
    def get_stats (self, type=GenStatistics.SUMMARY):
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
        str_out.append("    CHOPPED Features   :  " + str(self.get_stats_name_count_total('CHOPPED')))
        str_out.append("    UNTOUCHED Features :  " + str(self.get_stats_name_count_total('UNTOUCHED')))
        str_out.append("=" * len(title))

        return str_out


class AlgoFlatChopper(Algorithm):
    """
    This is the main class for the algorithm

    Attributes:
        - params: Parameters for the algorithm
    """

    def __init__(self, 
                 keep_iter_results=False, 
                 debug=False):
        """Initialize the attributes of an object of this class

        Parameters:
            keep_iter_results: Flag to enable(True)/disable(False) keeping of the iterative results
            debug: Flag to enable(True)/disable(False) debug output

        Return value:
            None

        """

        Algorithm.__init__(self)        
        self.params = Parameters()
        self.params.keep_iter_results = keep_iter_results 
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
        properties_to_check = [_NB_ELEMENTS, _NB_VERTICE]
        for feature in self.features:
            if isinstance(feature, MA_LineString):
                GenUtil.check_feature_integrity(feature, MA_LineString, properties_to_check)
    
    def flatChop_this_feature(self, line, s_container):
        """ ffff """
        
        nb_elements = line.ma_properties[_NB_ELEMENTS]
        nb_vertice  = line.ma_properties[_NB_VERTICE]

        # Extract feature coordinates and copy in new list structure to work with
        featureCoords = line.coords_dual
        coords_list = list(featureCoords)

        nb_coords = len(coords_list)

        chops = []

        # User asked to chop original line into lines while looking for the most appropriate vertex in a
        # window of vertices each side of the theoretical chopping vertex.
            
        # Check if line has enough vertices to work with else simply pass out the feature
        if (nb_coords > (nb_vertice * 2) + 1) or (nb_vertice >= nb_elements):
                
            # Calculate bounds for window of vertices
            bound = nb_elements
            bound_low = bound - nb_vertice
            bound_high = bound + nb_vertice + 1

            # Iterate through the line and calculate new one.
            # When bound is found, related coordinates are removed from original line.
            inProgress = True
            while (inProgress):
                nb_coords=len(coords_list)

                # Exit if end of line
                if (nb_coords <= nb_elements):
                    inProgress = False
                else:
                    # Extract coordinates of theoretical line and ask to find the appropriate vertex where to chop
                    coords_sublist = coords_list[bound_low:bound_high]
                    flatter_angle, flatter_index = GenUtil.find_flatter_vertice(coords_sublist)

                    # Chopping vertex is found, extract line from original one and chop it
                    chop = bound_low + flatter_index
                    feature_coords = coords_list[0:chop+1]  # chop + 1 is to make sure to keep line continuous
                    del coords_list[0:chop]

                    # Append the coordinates to the list of new lines
                    chops.append(feature_coords)

            # If remaining number of vertices < nb_vertice / 2, then merge them to last segment else create new segment
            if (len(coords_list) < (nb_vertice/2)):
                chops[-1].extend(coords_list[1:])
            else:
                # Append the remaining coordinates to the list of new lines
                chops.append(coords_list)

        else:
            # Line is too short to be chopped using current parameters. Line is kept untouched
            chops.append(coords_list)

        # Now, create a new feature in the spatial container for each newly created segment
        GenUtil.print_debug (self.params, "Number of segmented lines : %s" % (len(chops)))
        
        if (len(chops) > 1):
            line.update_coords([(0,0),(0,0)], s_container)
            line.ma_properties[_NB_CHOPS] = len(chops)
            line.ma_properties[_SEGMENT_INDEX_ATTRIBUTE] = 1
            line.ma_properties[_IS_CHOPPED] = True
            for i in range(1,len(chops)):
                new_line = line.cloner()
                new_line.update_coords(chops[i])
                new_line.ma_properties[_SEGMENT_INDEX_ATTRIBUTE] = i + 1
                s_container.add_feature(new_line)
            line.update_coords(chops[0], s_container)            
        else:
            line.ma_properties[_SEGMENT_INDEX_ATTRIBUTE] = 0
            line.ma_properties[_IS_CHOPPED] = False
            line.ma_properties[_NB_CHOPS] = 0

        line.chops = chops

        return len(chops)

    def process(self):

        GenUtil.print_debug (self.params, "Start of %s algorithm)" % (_ALGO))
        GenUtil.print_debug (self.params, "Parameter description:")
        GenUtil.print_debug (self.params, "  - Keep iterative results: %s" % (self.params.keep_iter_results))

        # Check the feature's class and attributes
        self.check_features()
        
        # Load the shapely features into the spatial container
        self.s_container = self.load_features()

        if (self.params.debug):
            nbr_lines = len(self.s_container.get_features(filter="feature.feature_type=GenUtil.LINE_STRING")) 
            nbr_points = len(self.s_container.get_features(filter="feature.feature_type==GenUtil.POINT")) 
            GenUtil.print_debug (self.params, "Number of lines imported: %s" % (nbr_lines))
            GenUtil.print_debug (self.params, "Number of points imported: %s" % (nbr_points))
          
        self.stats.add_iteration()
        GenUtil.print_debug (self.params, 'Start of iteration # %s' % (self.stats.get_nbr_iteration()))
            
        for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING"):
            self.stats.add_stats('TOTAL_PRE')

            nbChops = self.flatChop_this_feature(line, self.s_container)
            
            if (nbChops > 1):
                self.stats.add_stats('CHOPPED')
                self.stats.add_stats('TOTAL_POST', nbChops)
            else:
                self.stats.add_stats('UNTOUCHED')
                self.stats.add_stats('TOTAL_POST')
                                                          
        GenUtil.print_debug (self.params, 'End of iteration # %s' % (self.stats.get_nbr_iteration()))

        self.extract_features_out(self.s_container.get_features())        
        GenUtil.print_debug (self.params, "End of %s" % (_ALGO))        
        #self.stats.print_stats()
