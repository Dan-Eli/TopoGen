#!/usr/local/bin/python
# -=- encoding: utf-8 -=-

#####################################################################################################################################

"""
    This algorithm implement the insert vertice algorithm.

    This algorithm will insert a vertice into each line segment when another vertice
    of another line stands within the predefined tolerance

    Usage:
        import insert_vertice

    Limits and constraints


"""
__revision__ = "--REVISION-- : $Id: algo_insert_vertice.py 454 2011-09-27 12:52:53Z dpilon $"

#####################################################################################################################################

import math

from shapely.geometry import Point, LineString, Polygon
from shapely.prepared import prep

from lib_genmetal import MA_Point, MA_LineString, MA_Polygon, InvalidParameterError, GenUtil,\
                         GenStatistics, PointErrorPosition, LineStringErrorPosition, IterationResults, Parameters,\
                         SpatialContainer, Algorithm, InvalidGeometryError

########################################################

#Internal key word constants
_ALGO = 'Insert vertice'
_TOLERANCE = "tolerance"

class Statistics(GenStatistics):
    """Class that contains the statistics for the DP algorithm

    Attributes
        stat_names: Name of the statistics for the SpikeStatistics class. These name are
                    used by the Statistics class

    """

    def __init__(self):
        """Initialize the attributes of an object of the class"""

        GenStatistics.__init__(self)
        self.stats_names = ((_ALGO,))

    def get_stats (self, type=GenStatistics.SUMMARY):
        """Extract the current statistics and build  a list of string that forms the statistical message"

        Parameters:
            type: Give the form of statistics to extract. Can take 2 values.
                SUMMARY: Summary information
                DETAILED: Detailed information

        """

        str_out = []
        str_out.append( "Insert Vertice algorithm Statistics" )

        str_out.append( "-----------------------------------" )
        str_out.append("Vertices added: " + str(self.get_stats_name_count_total( _ALGO)))

        return str_out
    
class AlgoInsertVertice(Algorithm):
    """
    This is the main class for the spike algorithm

    Attributes:
        - params: Parameters for the algorithm
    """

    def __init__(self, debug=False):
        """Initialize the attributes of an object of the class DPAlgorithm

        Parameters:
            debug: Flag for debug output
                True: Output debug information
                False: Do no output debug information 

        Return value:
            None

        """

        Algorithm.__init__(self)
        
        self.params = Parameters()

        self.params.debug=debug
        
        self.stats = Statistics()

    def check_features(self):
        """
        Check if the features passed in parameters are of the good class type and have the good attributes
        
        Parameters: None
        
        Return value: None
        
        """
                    
        # Check the line string
        class_type = MA_LineString
        properties_to_check = [_TOLERANCE]
        for feature in self.features:
            if isinstance(feature, MA_LineString):
                GenUtil.check_feature_integrity(feature, MA_LineString, properties_to_check)
        
    def _set_line_attributes (self):
        """This routine sets the attributes to the line

        The routine checks:
            - if the tolerance is 0 a line does not need any further procvessing

        Parameters: None

        Return value: None

        """

        # Select only the lines
        for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING"):
            line.is_simplest = False
            if ( line.ma_properties[_TOLERANCE] == 0.):
                # A line with a tolerance of 0 is at its simplest form
                line.is_simplest = True
                
        return

    def _build_bounding_box(self, tolerance, coord):
        """Adjust the bounding box
        
        Parameters: 
            line: LineString to use
            coord: (x,y) tuple
            
         Return value: bounding box tuple (xmin, ymin, xmax, ymax)  
            
        """
        
        delta = 1.1 * tolerance
        
        xmin = coord[0] - delta
        ymin = coord[1] - delta
        xmax = coord[0] + delta
        ymax = coord[1] + delta
        
        return (xmin, ymin, xmax, ymax)

    def add_vertice_on_others(self, s_container, line):
         """ Add a vertice to the near by lines if they are close from the line passes in parameter
         
         This routine is scanning each vertice of the line to see if there are any line nearby.  If a nearby line
         is found,  a vertice is inserted if the on that line there is no vertices within the tolerance
         
         Parameter:
             s_container: SpatialContainer containing the line passed in parameters 
             line: MA_LineString line to process
             
         Return value: None
         
         """
         
         # Loop over each vertice
         for coord in line.coords_dual:
             coord_bbox = self._build_bounding_box(line.ma_properties[_TOLERANCE], coord)
             target_lines = s_container.get_features(coord_bbox, remove_keys=[line._sci_id])
             for target_line in target_lines:
                 # Check if this point is within the tolerance of the target line
                 if (GenUtil.distance(coord, target_line.coords_dual[0])  <= line.ma_properties[_TOLERANCE] or
                     GenUtil.distance(coord, target_line.coords_dual[-1]) <= line.ma_properties[_TOLERANCE]):
                     # The first or last vertice in within tolerance of the coordinate to check pass to the next line...
                     pass
                 else:
                     coord_point = Point(coord)
                     # Check that the point to insert is really within the tolerance of the line
                     if (coord_point.distance(target_line) <= line.ma_properties[_TOLERANCE]):
                         point_to_insert = Point(coord)
                         # Find the relative position of the point on the target line
                         target_l_r = target_line.project(point_to_insert, normalized = True)
                         # Loop over each coord of the target line
                         if (target_l_r == 0. or target_l_r == 1.0):
                             # The target line is out of reach of the vertice to add there is nothing to do 
                             pass
                         else:
                             # Using linear referencing find the  coordinate in the target_ine which is located just before
                             # the point to insert
                             if (GenUtil.add_vertex(target_line, target_l_r, tolerance= target_line.ma_properties[_TOLERANCE],
                                                    s_container= s_container)):
                                 self.stats.add_stats(_ALGO)                                 
                     else:
                         # Line is outside the telerance process the next line
                         pass 
                     
    def add_vertice_on_itself(self, line):
        """ Add vertices on the line itself is part 
         
        This routine is breaking the line as a bunch of small 2 verices line.  It calls the routine to add verice between the line.
        It reconstructs the original line with the extra vertices added
         
        Parameter:
            Line: MA_LineString line to process
             
        Return value: None
         
        """
         
        # Only process the line that have 3 coords and more
        if (len(line.coords_dual) >= 3):
            s_container_small_lines = SpatialContainer()
            lst_small_lines = []
            # Break the line string into small two vertices line string
            for i in xrange(len(line.coords_dual)-1):
                small_line = MA_LineString(line.coords_dual[i:i+2])
                small_line.ma_properties[_TOLERANCE] = line.ma_properties[_TOLERANCE]
                lst_small_lines.append(small_line)
            
            #Add each one in the container
            s_container_small_lines.add_features(lst_small_lines)
                
            # Add vertices on the small line
            for small_line in s_container_small_lines.get_features():
                self.add_vertice_on_others(s_container_small_lines, small_line)
                
            #Rebuild the coordinates of the original line with the extra vertices id some are added
            new_coords = lst_small_lines[0].coords_dual[0]
            new_coords = [new_coords] # Create a list
            for small_line in lst_small_lines:
                new_coords.extend(small_line.coords_dual[1:])
                
#           Update the coordinates of the original line
            line.update_coords(new_coords, self.s_container) 
                                     
    def process(self):

        GenUtil.print_debug (self.params, "Start of InserVertice algorithm)")

        self.check_features()
        
        self.s_container = self.load_features()

        if (self.params.debug):
            #Test if print is needed before scanning the s_container for nothing... waste of time...
            nbr_lines = len(self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING")) 
            nbr_points = len(self.s_container.get_features(filter="feature.feature_type==GenUtil.POINT")) 
            GenUtil.print_debug (self.params, "Number of lines imported: %s"  %(nbr_lines) )
            GenUtil.print_debug (self.params, "Number of points imported: %s"  %(nbr_points) )
        
        self._set_line_attributes()
        
        self.stats.add_iteration()

        for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING and not feature.is_simplest"):
            self.add_vertice_on_itself(line)
        
        for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING and not feature.is_simplest"):
            self.add_vertice_on_others(self.s_container, line)

        self.features = [feature for feature in self.s_container.get_features()]
        
        GenUtil.print_debug (self.params, "End of %s" %(_ALGO))
