#!/usr/local/bin/python
# -=- encoding: utf-8 -=-

#####################################################################################################################################

"""
    This algorithm allows to close open lines along a neatline or an offset neatline.

    - In order to be closed a feature must meet the following criterias 
        - The feature must be a MA_LineString.
        - The line feature must be open (first/last vertice not touching)
        - The first/last vertice must be within a tolerance from the neatline
    - All the features not meeting those criteria's are outputted untouched
        
    - When a feature is closed; if the neatline needs to be offset the offset is applied to the
      neatline the open line is snapped to the neatline and closed using the shortest exterior on the neatline 
    - 
    
    The algorithm has 2 parameters:
        - offset: Buffer value in ground unit to applied to the neatline. If 0 no buffer is applied
        - tolerance: maximum distance from the neatline for the first/last vertice of an open line in order to be closed
    
    The entities entering the algorithm have each one property
    - code: Identify the type of feature. Can take 2 values
            N: MA_LineString feature corresponding to the neatline
            F: MA features to be processed; only MA_LineString are closed; all other are output untouched
            
    Usage:
        import close_line

    Limits and constraints
        Once must supply a neatline for the features 


"""
__revision__ = "--REVISION-- : $Id: algo_detect_bottleneck.py 312 2011-05-09 19:35:40Z dpilon $"

#####################################################################################################################################

from shapely.geometry import Point, LineString, Polygon
from copy import deepcopy
from shapely.ops import linemerge
from lib_genmetal import MA_Point, MA_LineString, MA_Polygon, GenUtil,\
                         GenStatistics, Parameters,\
                         SpatialContainer, Algorithm

########################################################

# Public constant
NEATLINE = "N"
FEATURE = "F"

# Properties definition
_CODE = 'code'
_TO_BE_CLOSED = "To_be_closed"
_NBR_CLOSED_LINES = "nbr_closed_lines"
_NBR_UNCLOSED_LINES = "nbr_unclosed_lines"


#Internal key word constants
_ALGO = 'Close line'

class Statistics(GenStatistics):
    """Class that contains the statistics for this algorithm

    Attributes
        stat_names: Name of the statistics for the SpikeStatistics class. These name are
                    used by the Statistics class

    """

    def __init__(self):
        """Initialize the attributes of an object of the class"""

        GenStatistics.__init__(self)
        self.stats_names = ((_ALGO,_NBR_CLOSED_LINES, _NBR_UNCLOSED_LINES))

    def get_stats (self, type=GenStatistics.SUMMARY):
        """Extract the current statistics and build  a list of string that forms the statistical message"

        Parameters:
            type: Give the form of statistics to extract. Can take 2 values.
                SUMMARY: Summary information
                DETAILED: Detailed information

        """

        str_out = []
        str_out.append( "Close line statistics" )

        str_out.append( "---------------------" )
        str_out.append("Number of line closed: "   + str(self.get_stats_name_count_total( _NBR_CLOSED_LINES)))
        str_out.append("Number of line untouched: " + str(self.get_stats_name_count_total( _NBR_UNCLOSED_LINES)))
        str_out.append( "---------------------" )

        return str_out
    
class AlgoLineCloser(Algorithm):
    """
    This is the main class for this algorithm

    Attributes:
        - params: Parameters for the algorithm
    """

    def __init__(self, buffer=0., tolerance=0. ):
        """Initialize the attributes of an object of the class DPAlgorithm

        Parameters:
            offset: buffer to apply to the neat line if 0. no buffer is applied.
            tolerance: maximum distance from the neatline for an open line in order to be closed 

        Return value:
            None

        """

        Algorithm.__init__(self)
        
        self.params = Parameters()
        
        self.params.buffer = buffer
        self.params.tolerance = tolerance
        self.params.debug = False
        
        self.stats = Statistics()
        self.stats.add_iteration()
       
    def _check_features(self):
        """Check if the features passed in parameters are of the good class type and have the good attributes
        
        Parameters: None
        
        Return value: None
        
        """
        
        # Check the line string
        properties_to_check = [_CODE]
        for feature in self.features:
            if isinstance(feature, MA_LineString):
                GenUtil.check_feature_integrity(feature, None, properties_to_check)
                
        return
    
    def _classify_features(self):
        """Classify the feature in 3 groups: neatline, line to be closed and all the others features
        
        To be classify as to be closed, the line must meet the following criterias
            - To be a open line string within a distance from the neat line
        
        Parameters: None
        
        Return value: None
        
        """

        # first extract only the neatline       
        self.neatline = [feature for feature in self.features if feature.ma_properties[_CODE] == NEATLINE]        
        if (len(self.neatline) == 1):
            # Check if only one neatline is present
            self.neatline = self.neatline[0]
            # Check if the neatline is closed
            if not GenUtil.is_line_closed(self.neatline):
                raise Exception ("The neatline should be closed...")
        else:
            raise Exception ("There must be one and only one neatline...") 
        
        # Classify the lines as to be closed or other...
        self.lines_to_close = []
        self.others = []
        for feature in self.features:
            if feature.ma_properties[_CODE] != NEATLINE:
                # Process all except the neatline
                if (isinstance(feature, MA_LineString)):
                    # Select the MA_LineString
                    if (not GenUtil.is_line_closed(feature)):
                        # Select only open line string
                        first_coord = feature.coords_dual[0]
                        last_coord = feature.coords_dual[-1]
                        if (feature.distance(Point(first_coord)) < self.params.tolerance and
                            feature.distance(Point(last_coord))  < self.params.tolerance ):
                             # Select only the line string located within a tolerance from the neat line
                             self.lines_to_close.append(feature)
                             self.stats.add_stats(_NBR_CLOSED_LINES)
                        else:
                             self.others.append(feature)
                             self.stats.add_stats(_NBR_UNCLOSED_LINES)
                    else:
                        self.others.append(feature)
                        self.stats.add_stats(_NBR_UNCLOSED_LINES)
                else:
                    self.others.append(feature)

    def _close_open_line_on_neatline(self, line_to_close, buffered_neatline):
        """Close an open line on the neatline
        
        The method do th efollowing steps
            - Snaps the line to close on the buffered neatline.
            - Insert 2 vertice on the buffered neat line
            - Cut the neatline in 2 lines at the position of the vertice inserted
            - Create to closed line using the 2 lines of the neatline
            - Choose only the shortest line as the closed one
        
        Parameters:
            - line_to_close: MA_LineString of the line to close
            - buffered_neatline: MA_LineString of the buffered neatline
            
        Retiurn value: 
            - MA_LineString of the closed line
            
        """
        
        # Create a clone of the neatline as we will edit the original line
        buffered_neatline = buffered_neatline.cloner() 
        
        # Extract the first and last vertice of the line to close
        first_coord = line_to_close.coords_dual[0]
        last_coord = line_to_close.coords_dual[-1]
        
        #Project the first/last vertice on the neatline
        dist_first_coord = buffered_neatline.project(Point(first_coord), normalized =True)
        dist_last_coord = buffered_neatline.project(Point(last_coord), normalized=True)
        
        # Extend the line to closed to make it touch the neatline
        first_point = buffered_neatline.interpolate(dist_first_coord, normalized =True)
        last_point  = buffered_neatline.interpolate(dist_last_coord, normalized =True)
        lst_coords = list(line_to_close.coords)
        if (self.params.buffer != 0.):
            extend_lst_coords = [first_point.coords[0]] + lst_coords + [last_point.coords[0]]
        else:
            #If there is no buffer on the neatline there is no need to extend the line to be closed
            extend_lst_coords = lst_coords
        extend_line_to_close = LineString(list(extend_lst_coords))
        
        # On the neatline make sure the first point of the line to close is located before the last point  
        if (dist_last_coord < dist_first_coord):
            # Flip the point
            dist_last_coord, dist_first_coord = dist_first_coord, dist_last_coord
            first_point, last_point = last_point, first_point
        first_coord = first_point.coords[0]
        last_coord = last_point.coords[0]
            
        # Add the vertex on the neat line
        GenUtil.add_vertex(buffered_neatline, dist_first_coord, 0.)
        GenUtil.add_vertex(buffered_neatline, dist_last_coord, 0.)
            
        # Extract the position (index) of the first and last point on the neat line
        nbr_coords = len(buffered_neatline.coords_dual)
        for i in range(nbr_coords):
            if GenUtil.distance(buffered_neatline.coords_dual[i], first_coord) <= GenUtil.ZERO:
                first_neatline = i # Position of the first coordinate of the line to close on the neat line
                for j in range(i,nbr_coords):
                    GenUtil.distance(buffered_neatline.coords_dual[j], first_coord)
                    if GenUtil.distance(buffered_neatline.coords_dual[j], last_coord) <= GenUtil.ZERO:
                        last_neatline = j # Position of the last coordinate of the line to close on the neat line
                        break
                break
            
        if ( (first_neatline == last_neatline) or
             (first_neatline == 0 and last_neatline == nbr_coords-1) ):
            # The first/last vertice of the line to be closed is located on the same vertice on the neatline
            closed_line = extend_line_to_close
        else:
            # Create the 2 sub neatlines bu cutting the neatline in 2  
            sub_neatline_coords_a = buffered_neatline.coords_dual[first_neatline:last_neatline+1]
            sub_neatline_a = LineString(list(sub_neatline_coords_a)) 
            sub_neatline_coords_b = buffered_neatline.coords_dual[last_neatline:] + buffered_neatline.coords_dual[1:first_neatline+1]
            sub_neatline_b = LineString(list(sub_neatline_coords_b))
            # Merge the each sub neatline with the line to be closed
            line_a = linemerge([extend_line_to_close, sub_neatline_a])
            line_b = linemerge([extend_line_to_close, sub_neatline_b])
            # Now choose the shortest line
            if (line_a.length < line_b.length):
                closed_line = line_a
            else:
                closed_line = line_b 
            
        # Force the first/last vertice to be the same
        lst_coords = list(closed_line.coords)
        lst_coords[-1] = lst_coords[0]
        # Transform the feature into a MA_LineString
        closed_line = MA_LineString(lst_coords)
        # Copy the original attributes
        closed_line.ma_properties = deepcopy(line_to_close.ma_properties)
            
        return closed_line
    
    def _close_open_lines(self):
        """This method manages the closing of all the open lines
        
        Parameters: None
        
        Return value: None
        
        """
        
        if (self.params.buffer != 0.):
            # buffer the neatline if needed
            buffered_neatline = self.neatline.buffer(self.params.buffer, resolution=1)
            buffered_neatline = MA_LineString(list(buffered_neatline.exterior.coords))
        else:
            buffered_neatline = MA_LineString(list(self.neatline.coords))
        buffered_neatline.ma_properties = deepcopy(self.neatline.ma_properties)
            
        self.closed_lines = []
        for line_to_close in self.lines_to_close:
            closed_line = self._close_open_line_on_neatline(line_to_close, buffered_neatline)
            if (not GenUtil.is_line_closed(closed_line)):
                # This should never arrive but a small check is done if a special case is not trapped...
                raise InternalError ("The line should have been closed and is not closed...")
            self.closed_lines.append(closed_line)
        
        return buffered_neatline
    
    def process(self):

        GenUtil.print_debug (self.params, "Start of %s algorithm" %_ALGO)

        self.stats.add_iteration()
        
        self._check_features()
        
        self._classify_features()
        
        buffered_neatline = self._close_open_lines()
        
        #Output the features
        self.features = []
        self.features = self.closed_lines + self.others + [buffered_neatline]
        
        GenUtil.print_debug (self.params, "End of %s" %(_ALGO))
