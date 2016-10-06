#! /usr/bin/env python
# -*- coding: UTF-8 -*-

#####################################################################################################################################

"""
    This algorithm change the angle between 3 points when the angle is below a specified threshold.

    This algorithm calculates the angle of each vertices of a line when the angle is below a certain
    threshold the vertice is moved in the direction of the bisector in order to increase the angle up 
    the threshold value. The point and lines that do not need to be unspike are still used to enforce 
    topology integrity between  those feature that need to be unspike
    
    Usage:
        import spike

    Limits and constraints
        Always works better when the line to process meet the OGC simple line.
          

"""
__revision__ = "--REVISION-- : $Id: algo_spike.py 457 2011-09-27 12:54:28Z dpilon $"

#####################################################################################################################################


import math

from shapely.geometry import Point, LineString, Polygon
from shapely.prepared import prep

from lib_genmetal import MA_Point, MA_LineString, MA_Polygon, InvalidParameterError, GenUtil,\
                         GenStatistics, PointErrorPosition, LineStringErrorPosition, IterationResults, Parameters,\
                         SpatialContainer, Algorithm
                         
########################################################

#Internal key word constants
_ALGO = 'Spike'
_MAX_LOOP_COUNTER = 200
_MINIMUM_ANGLE = "minimum_angle"
_MOVE_FIRST_LAST = "move_first_last"

class Statistics(GenStatistics):
    """Class that contains the statistics for the pike algorithm
    
    Attributes
        stat_names: Name of the statistics for the SpikeStatistics class. These name are
                    used by the Statistics class

    """
    
    def __init__(self):
        """Initialize the attributes of an object of the class"""
        
        GenStatistics.__init__(self)
        self.stats_names = ((_ALGO, GenUtil.SIMPLE_LINE, GenUtil.CROSSING_LINE, GenUtil.SIDEDNESS, ))
        
    def get_stats (self, type=GenStatistics.SUMMARY):
        """Extract the current statistics and build  a list of string that forms the statistical message"
        
        Parameters:
            type: Give the form of statistics to extract. Can take 2 values.
                SUMMARY: Summary information
                DETAILED: Detailed information
        
        """
        
        str_out = []
        str_out.append( "Spike algorithm Statistics" )
        
        str_out.append( "-" * len(str_out[-1]))
        if (type == GenStatistics.DETAILED):
            for i in xrange((self.get_nbr_iteration())):
                str_out.append("Detailed statistics")
                str_out.append("Iteration # " + str(i))
                str_out.append("Spikes removed: " + str(self.get_stats_name_count_iter( _ALGO, i)))
                str_out.append( "--------" )
                str_out.append( "Conflicts:" )
                str_out.append( "    Simple Line   :  " + str(self.get_stats_name_count_iter( GenUtil.SIMPLE_LINE, i)))
                str_out.append( "    Crossing Line :  " + str(self.get_stats_name_count_iter( GenUtil.CROSSING_LINE, i)))
                str_out.append( "    Sidedness     :  " + str(self.get_stats_name_count_iter( GenUtil.SIDEDNESS, i)))
                str_out.append( "--------" )
        str_out.append( "Summary statistics" )
        str_out.append("Spikes removed: " + str(self.get_stats_name_count_total( _ALGO)))
        str_out.append( "--------" )
        str_out.append( "Conflicts:" )
        str_out.append( "    Simple Line   :  " + str(self.get_stats_name_count_total( GenUtil.SIMPLE_LINE)))
        str_out.append( "    Crossing Line :  " + str(self.get_stats_name_count_total( GenUtil.CROSSING_LINE)))
        str_out.append( "    Sidedness     :  " + str(self.get_stats_name_count_total( GenUtil.SIDEDNESS)))
        str_out.append( "--------" )
        str_out.append( "Number of iteration: " + str(self.get_nbr_iteration()) )
        
        return str_out

class AlgoSpike(Algorithm):
    """This is the main class for the spike algorithm 
    
    Attributes:
        - requires_params: Dictionary containing the requires parameters for the spike algorithm
        - default_params: Dictionary containing the default values for the requires parameters
        
    """

    def __init__(self, test_crossing_line=True,
                       test_simple_line=True, 
                       test_sidedness=True, 
                       keep_iter_results=False, 
                       debug=False):
        """Initialize the attributes of an object of the class DPAlgorithm

        Parameters:
            test_crossing_line: Flag to enable(TRUE)/disable(FLASE) checking the line crossing constraint
            test_simple_line: Flag to enable(TRUE)/disable(FLASE) checking the simple line constraint
            test_sidedness:  to enable(TRUE)/disable(FLASE) for checking the sidedness constraint
            keep_iter_results: Flag to enable(TRUE)/disable(FLASE) keeping the iterative results
            keep_iter_results: Flag for the iterative results
            debug: Flag to enable(TRUE)/disable(FLASE) for debug output

        Return value:
            None

        """

        Algorithm.__init__(self)
        
        self.params = Parameters()
        self.params.test_crossing_line=test_crossing_line
        self.params.test_simple_line=test_simple_line 
        self.params.test_sidedness=test_sidedness 
        self.params.keep_iter_results=keep_iter_results 
        self.params.debug=debug
        
        self.stats = Statistics()
    
    def check_features(self):
        """
        Check if the features passed in parameters are of the good class type and have the good attributes
        
        Parameters: None
        
        Return value: None
        
        """
        # Check the line string
        properties_to_check = [_MINIMUM_ANGLE, _MOVE_FIRST_LAST]
        for feature in self.features:
            if isinstance(feature, MA_LineString):
                GenUtil.check_feature_integrity(feature, MA_LineString, properties_to_check)        

    
    def _set_line_attributes (self):
        """This routine sets the attributes to the line
        
        The routine checks:
            - if the first/last vertices are the same the line is closed; otherwise it is open
            - if the line has less than 3 vertices it is at its simplest form; otherwise it is not
        
        Parameters: None
        
        Return value: None
        
        """
        
        # Select only the lines
        for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING"):
            line.is_closed = False
            line.is_simplest = False
            if (len(line.coords_dual) >= 3):
                if (GenUtil.distance(line.coords_dual[0], line.coords_dual[-1]) <= GenUtil.ZERO):
                    line.is_closed = True
            else:
                line.is_simplest = True
            
            # If the minimum angle is zero it also means that the line is at simplest form
            if (line.ma_properties[_MINIMUM_ANGLE] <= GenUtil.ZERO):
                line.is_simplest = True
    
        return

    def unspike(self, line ):
        """Replace the spike that are below a certain angle
        
        Parameters:
            - line: Shapely LineString to  unspike
            
        Return values
            Status of the spike process
                True: The line was unspiked
                False: The line was not unspiked
            
        """

        line.is_simplest = True
        unspiked = False
        nb_angle = len(line.coords_dual) - 2
        #The first/last vertice formed an supplementary angle if the line is closed
        if (line.is_closed):
            # Do we really want to move the first and last
            if (line.ma_properties[_MOVE_FIRST_LAST]):
                nb_angle += 1
        for i in xrange(nb_angle):
            p1 = line.coords_dual[i]
            p2 = line.coords_dual[i+1]
            if (i+2 == len(line.coords_dual)):
                # Special case for first/last vertice
                p3 = line.coords_dual[1]
                first_last = True
            else:
                p3 = line.coords_dual[i+2]
                first_last = False
                
            # The anlge calculation is within "try..except" as somebad lines may formed invalid angle calculation
            try:
                angle = GenUtil.angle_vector(p1,p2,p3)
            except:
                #An angle of 0 is not processed...
                angle = 0.

            if (angle < line.ma_properties[_MINIMUM_ANGLE] and self._are_points_valid(p1, p2, p3, angle)):
                p2_new = self._calculate_p2_new (line.ma_properties[_MINIMUM_ANGLE], p1, p2, p3)  
 
                # Calculate new line position
                if first_last:
                    # Move the first and last vertice at the same point
                    new_line = LineString([p2_new]+list(line.coords_dual[1:-1])+[p2_new])
                else: 
                    # The vertice to move is not a first/last vertice
                    new_line = LineString(list(line.coords_dual[0:i+1]) + [p2_new] + list(line.coords_dual[i+2:]))
                
                line_simple_line = new_line
                line_crossing_line = LineString([p1, p2_new, p3])
                sidedness_polygon = GenUtil.calculate_sidedness_polygon( LineString([p1, p2, p3]),
                                                                         LineString([p1, p2_new, p3]) )
                
                conflict_type = GenUtil.test_constraints (self, None, line_simple_line, line_crossing_line, 
                                                          sidedness_polygon, self.s_container, line._sci_id)

                if (conflict_type is not None):
                    line.is_simplest = False
                else:
                    # Update with the new coordinates
                    self.stats.add_stats(_ALGO)
                    line.update_coords(new_line.coords, self.s_container)
                    unspiked = True 
            
        return unspiked
    
    def _are_points_valid(self, p1, p2, p3, angle):
        """Test if points are valid in order to unspike the angle formed by p1, p2, p3
        
        To be valid the distance between point1-2, p2-3 and point1-3 must be over
        a small delta; otherwise this means the 3 vertices are at the same place. 
        The angle formed by p1,p2,p3 must also be over a small delta otherwise
        the line fall on itself
        
        Parameters:
            p1: (x,y) tuple coordinate
            p2: (x,y) tuple coordinate
            p3: (x,y) tuple coordinate
            angle: angle in degree formed by p1,p2,p3
            
        Return value
            Flag indicating if the points area valid
                True: The points are valid
                False: The points are invalid
                
        """
        
        if ( GenUtil.distance(p1,p2) > GenUtil.ZERO and
             GenUtil.distance(p2,p3) > GenUtil.ZERO and
             GenUtil.distance(p3,p1) > GenUtil.ZERO and
             angle > GenUtil.ZERO):
            point_valid = True
        else:
            point_valid = False
            
        return point_valid
          
    def _calculate_p2_new (self, minimum_angle, p1, p2, p3):
        """ Calculate the new p2 vertices position  that conforms to the minimum_angle constraint
        
        Using a dichotomic search the method will find the angle near by the minimum_angle constraint.
        In order to prevent infinite loops... it's always possible if the algorithm a unable to converge...
        A loop counter as been installed after 100 loops it exits the loop 
        
        Parameters:
            minimum_angle: Minimum angle to reach
            p1: first x,y tuple
            p2: second x,y tuple
            p3: third x,y tuple
        
        Return value
            (x,y) tuple containing the new position of the point p2
        
        """
    
        # Calculate the distance between the each point
        dist_p1_p2 = GenUtil.distance (p1, p2)
        dist_p1_p3 = GenUtil.distance (p1, p3)
        dist_p2_p3 = GenUtil.distance (p2, p3)
    
        dist_p3_p4 = (dist_p2_p3*dist_p1_p3) / (dist_p1_p2+dist_p2_p3)
    
        p4 = GenUtil.rescale_vector ( p1, p3, 1-dist_p3_p4/dist_p1_p3) 
    
        scale_start = 0.0
        scale_end = 1.0
        scale_middle = (scale_start+scale_end)/2.0
        p2_new = GenUtil.rescale_vector ( p4, p2, scale_middle)
        angle_theta = GenUtil.angle_vector (p1, p2_new, p3)
        delta_target = angle_theta - minimum_angle
        loop_counter = 0
    
        while (delta_target  <= 0.0 or delta_target >= 3.0 ):
            if delta_target <= 0.0:
                scale_end = scale_middle
            else:
                scale_start = scale_middle
            scale_middle = (scale_start+scale_end)/2.0
            p2_new = GenUtil.rescale_vector ( p4, p2, scale_middle)
            angle_theta = GenUtil.angle_vector (p1, p2_new, p3)
            delta_target = angle_theta - minimum_angle
            if (loop_counter < _MAX_LOOP_COUNTER):
                loop_counter += 1
            else:
                # The algorithm is unable to converge... exit the loop
                break
        
        return p2_new

    def process(self):
        """Main routine for the spike_reduction algorithm
        
        This algorithm will remove spikes in line but it will prevent any topological
        errors during the line cleaning.  It will keep the line as a OGC simple line.
        It will prevent line crossing and sidedness errors.
        
        Parameters:
            interface:
                Interface object containing all input to process
            
        """
        
        GenUtil.print_debug (self.params, "Start of %s  algorithm)" %(_ALGO))
        GenUtil.print_debug (self.params, "Parameter description:")
        GenUtil.print_debug (self.params, "  - Simple line constraint: %s" %(self.params.test_simple_line))
        GenUtil.print_debug (self.params, "  - Crossing line constraint: %s" %(self.params.test_crossing_line))
        GenUtil.print_debug (self.params, "  - Sidedness constraint: %s" %(self.params.test_sidedness))
        GenUtil.print_debug (self.params, "  - Keep iterative results: %s" %(self.params.keep_iter_results))
            
        # Check the feature's class and attributes
        self.check_features()
        
        # Load the shapely features into the spatial container
        self.s_container = self.load_features ()
     
        if (self.params.debug):
            #Test if print is needed before scanning the s_container for nothing... waste of time...
            nbr_lines = len(self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING")) 
            nbr_points = len(self.s_container.get_features(filter="feature.feature_type==GenUtil.POINT")) 
            GenUtil.print_debug (self.params, "Number of lines imported: %s"  %(nbr_lines) )
            GenUtil.print_debug (self.params, "Number of points imported: %s"  %(nbr_points) )
            
        self._set_line_attributes()
        
        line_simplified = True
        
        while (line_simplified):
            # At each new iteration we create a new itearation
            self.stats.add_iteration()
            GenUtil.print_debug (self.params, 'Start of iteration # %s' %(self.stats.get_nbr_iteration()) )  
            line_simplified = False
    
            for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING and not feature.is_simplest"):
                line_simplified = self.unspike(line) or line_simplified
                                    
            GenUtil.print_debug (self.params, 'Number of spike removed %s: ' %(self.stats.get_stats_name_count_iter(_ALGO)))
            GenUtil.print_debug (self.params, 'End of iteration # %s' %(self.stats.get_nbr_iteration()) )
            
            if (line_simplified):
                    # If line_simplified is True there will be an other iteration and we reset the error counter
                    self.stats.reset_stats_names([GenUtil.SIMPLE_LINE, GenUtil.CROSSING_LINE, GenUtil.SIDEDNESS]) # Reset error counter
                    self.error_positions = []  # Reset error position we don't need the old one
            
            if (self.params.keep_iter_results):
                self.iter_results.add_iteration()
                self.iter_results.add_features(self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING"))
    
        self.features = [feature for feature in self.s_container.get_features()]
        
        GenUtil.print_debug (self.params, "End of %s" % (_ALGO))
