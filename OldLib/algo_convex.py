#! /usr/bin/env python
# -*- coding: UTF-8 -*-

#####################################################################################################################################

"""
    This algorithm replace the input closed lines with there convex hulls when it meets the requirements

    This algorithm calculates the convex hull around a closed line. 
    The point and open lines that do not need to be "convex" are still used to enforce topology integrity between 
    those feature that need to be "convex"
    
    Usage:
        import convex

    Limits and constraints
"""
__revision__ = "--REVISION-- : $Id: algo_convex.py 201 2011-04-04 12:45:30Z dpilon $"

#####################################################################################################################################

import math

from shapely.geometry import Point, LineString, Polygon
from shapely.prepared import prep

from lib_genmetal import MA_Point, MA_LineString, MA_Polygon, InvalidParameterError, GenUtil,\
                         GenStatistics, PointErrorPosition, LineStringErrorPosition, IterationResults, Parameters,\
                         SpatialContainer, Algorithm
                         
        
#Internal key word constants
_ALGO = 'Convex'

# Feature properties
_MAX_AREA = 'max_area'
_INCREASE_RATIO = 'increase_ratio'


# If psyco is availiable, use it to speed up processing (2x)
try:
    import psyco
    psyco.full()
except ImportError:
    pass

class Statistics(GenStatistics):
    """Class that contains the statistics for the convex algorithm
    
    Attributes
        stat_names: Name of the statistics for the ConvexStatistics class. These name are
                    used by the Statistics class

    """
    
    def __init__(self):
        """Initialize the attributes of an object of the class"""
        
        GenStatistics.__init__(self)
        self.stats_names = ((_ALGO, GenUtil.INVALID, GenUtil.CROSSING_LINE, GenUtil.SIDEDNESS, ))
        
    def get_stats (self, type=GenStatistics.SUMMARY):
        """Extract the current statistics and build  a list of string that forms the statistical message"
        
        Parameters:
            type: Give the form of statistics to extract. Can take 2 values.
                SUMMARY: Summary information
                DETAILED: Detailed information
        
        """
        
        str_out = []
        str_out.append( "%s algorithm Statistics" %(_ALGO))
        str_out.append( "-" * len(str_out[-1]))
        if (type == GenStatistics.DETAILED):
            for i in xrange((self.get_nbr_iteration())):
                str_out.append("Detailed statistics")
                str_out.append("Iteration # " + str(i))
                str_out.append("Feature convexed: " + str(self.get_stats_name_count_iter( _ALGO, i)))
                str_out.append( "--------" )
                str_out.append( "Conflicts:" )
                str_out.append( "    Crossing Line :  " + str(self.get_stats_name_count_iter( GenUtil.CROSSING_LINE, i)))
                str_out.append( "    Sidedness     :  " + str(self.get_stats_name_count_iter( GenUtil.SIDEDNESS, i)))
                str_out.append( "--------" )
        str_out.append( "Summary statistics" )
        str_out.append("Invalid geometry: " + str(self.get_stats_name_count_total( GenUtil.INVALID)))
        str_out.append("Feature convexed: " + str(self.get_stats_name_count_total( _ALGO)))
        str_out.append( "--------" )
        str_out.append( "Conflicts:" )
        str_out.append( "    Crossing Line :  " + str(self.get_stats_name_count_total( GenUtil.CROSSING_LINE)))
        str_out.append( "    Sidedness     :  " + str(self.get_stats_name_count_total( GenUtil.SIDEDNESS)))
        str_out.append( "--------" )
        str_out.append( "Number of iteration: " + str(self.get_nbr_iteration()) )
        
        return str_out

class AlgoConvex(Algorithm):
    """This is the main class for the convex algorithm 
    
    Attributes:
        - params: Parameters of the algorithm
        
    """

    def __init__(self, test_crossing_line=True, 
                       test_sidedness=True, 
                       keep_iter_results=False, 
                       debug=False):
        """Initialize the attributes of an object of the class DPAlgorithm

        Parameters:
            test_crossing_line: Flag to enable(TRUE)/disable(FLASE) checking the line crossing constraint
            test_sidedness:  to enable(TRUE)/disable(FLASE) for checking the sidedness constraint
            keep_iter_results: Flag to enable(TRUE)/disable(FLASE) keeping the iterative results
            debug: Flag to enable(TRUE)/disable(FLASE) to keep debug output

        Return value:
            None

        """

        Algorithm.__init__(self)
        
        self.params = Parameters()
        self.params.test_crossing_line=test_crossing_line 
        self.params.test_sidedness=test_sidedness 
        self.params.keep_iter_results=keep_iter_results 
        self.params.debug=debug
        
        self.stats = Statistics()
		
    def _add_line_info_convex (self):
        """This routine identifies only the closed line as candidate line to be convex
		
		Parameters: None
		    s_container: SpatialContainer object containing the features and the spatial index
		
		Return value: None 
		
		"""
        
        # Process all the lines
        for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING"):
            line.is_simplest = True
            if (len(line.coords_dual) >= 3 and 
                GenUtil.distance(line.coords_dual[0], line.coords_dual[-1]) <= GenUtil.ZERO):
                line.polygon_line = Polygon(line.coords)
                # The next lines of codes are within "try .. except" as in some cases the closed line
                # may form a bad polygon.  The command line.polygon_line.exterior will create an exception
                # in that case we just flag the line as not candidate as to be convex
                try:
                    if ( line.polygon_line.exterior is not None):
                        line.is_simplest = False
                        line.polygon_hull = line.polygon_line.convex_hull
                except ValueError:
                    pass
                
    def check_features(self):
        """
        Check if the features passed in parameters are of the good class type and have the good attributes
        
        Parameters: None
        
        Return value: None
        
        """
        
        # Check the MA_LineString
        class_type = MA_LineString
        properties_to_check = [_MAX_AREA, _INCREASE_RATIO]
        for feature in self.features:
            if (isinstance(feature, MA_LineString)):
                GenUtil.check_feature_integrity(feature, MA_LineString, properties_to_check)        
    
    def _check_convex_ratio(self, line):
        """
        Determine if the polygon need to be transformed into its convex form or not. 
        
        To be candidate to be convex, 2 criteria's must be met
            - Area size of the convex feature must be < maximum area
            - The ratio of increase must be below the increase ratio
        
        Parameter:
            line: LineString feature to check
            
        Return value
            boolean: True: The polygon must be convex
                     False: The polygon must not be convex
                     
        """
        
        polygon_hull_area = line.polygon_hull.area
        polygon_line_area = line.polygon_line.area
        if (abs(polygon_hull_area-polygon_line_area) >= GenUtil.ZERO):       
            if (polygon_hull_area <= line.ma_properties[_MAX_AREA]):
                if ( polygon_hull_area/polygon_line_area < line.ma_properties[_INCREASE_RATIO]):
                    to_convex = True                        
                else:
                    # The area to convex is over the increase ratio
                    to_convex = False
            else:
                # The area to convex is over the maximum area
                to_convex = False
        else:
            # The feature is probably already a convex one
            to_convex = False

        return (to_convex)
                        
    def convex(self, line):
        """Method to replace a line by its convex hull
        
        Parameters:
            - line: LineString to convex
            
        Return value
            Boolean: True: The line was convex
                     False: The line was not convex
            
        """
        
        convexed = False
        line.is_simplest = True
    
        if (line.polygon_line.is_valid):
            # By checking if the closed input feature is valid we remove a lot of problem 
            # and special case management
            if (self._check_convex_ratio(line)):

                line_simple_line = None # simple line constraint is not check as the convex hull is always simple
                line_crossing_line = LineString(line.polygon_hull.exterior.coords)
                #polygon_sideness = self.polygon_hull.difference(self.polygon_line)
                sidedness_polygon = GenUtil.calculate_sidedness_polygon(line.polygon_hull, line.polygon_line)
                    
                conflict_type = GenUtil.test_constraints (self, None, line_simple_line, line_crossing_line, 
                                                          sidedness_polygon, self.s_container, line._sci_id)
    
                if (conflict_type is not None):
                    line.is_simplest = False
                else:
                    self.stats.add_stats(_ALGO)
                    convexed = True
                    line.update_coords(line.polygon_hull.exterior.coords, self.s_container)
                    
            else:
                line.is_simplest=True
        else:
                line.is_simplest=True
                self.stats.add_stats(GenUtil.INVALID)
            
        return convexed

    def process(self):
        """Main routine for the convex algorithm
        
        This algorithm will convex the closed line that meet the parameters condition
		It will keep the line as a OGC simple line. It will prevent line crossing and sidedness errors.
        
        Parameters:
            None
			
		Return value
			nONE
            
        """
        
        GenUtil.print_debug (self.params, "Start of convex  algorithm)")
        GenUtil.print_debug (self.params, "Parameter description:")
        GenUtil.print_debug (self.params, "  - Crossing line constraint: %s" %(self.params.test_crossing_line))
        GenUtil.print_debug (self.params, "  - Sidedness constraint: %s" %(self.params.test_sidedness))
        GenUtil.print_debug (self.params, "  - Keep iterative results: %s" %(self.params.keep_iter_results))
        
        # Check the feature's class and attributes
        self.check_features()
        
        # Load the shapely features into the spatial container
        self.s_container = self.load_features()
     
        if (self.params.debug):
            #Test if print is needed before scanning the s_container for nothing... waste of time...
            nbr_lines = len(self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING")) 
            nbr_points = len(self.s_container.get_features(filter="feature.feature_type==GenUtil.POINT")) 
            GenUtil.print_debug (self.params, "Number of lines imported: %s"  %(nbr_lines) )
            GenUtil.print_debug (self.params, "Number of points imported: %s"  %(nbr_points) )
            
        self._add_line_info_convex()
        
        line_simplified = True
        
        # Iterate until all the line are convexed or there are no more line have to be convex
        while (line_simplified):
            # At each new iteration we create a new itearation
            self.stats.add_iteration()
            GenUtil.print_debug (self.params, 'Start of iteration # %s' %(self.stats.get_nbr_iteration()) )
            line_simplified = False
    
            for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING and not feature.is_simplest"):
                line_simplified = self.convex(line ) or line_simplified
                                    
            GenUtil.print_debug (self.params, 'Number of feature convex %s: ' %(self.stats.get_stats_name_count_iter(_ALGO)))
            GenUtil.print_debug (self.params, 'End of iteration # %s' %(self.stats.get_nbr_iteration()) )
            
            if (line_simplified):
                    # If line_simplified is True there will be an other iteration and we reset the error counter
                    self.stats.reset_stats_names([GenUtil.CROSSING_LINE, GenUtil.SIDEDNESS]) # Reset error counter
                    self.error_positions = []  # Reset error position we don't need the old one
            
            if (self.params.keep_iter_results):
                self.iter_results.add_iteration()
                self.iter_results.add_features(self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING"))
    
        self.extract_features_out(self.s_container.get_features())
        
        GenUtil.print_debug (self.params, "End of %s" % (_ALGO))
        