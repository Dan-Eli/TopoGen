#!/usr/local/bin/python
# -=- encoding: utf-8 -=-

#####################################################################################################################################

"""
    This algorithm flags the line that forms a cycle.

    This algorithm will flag lines that creates cycle. This implementation only looks for 2 lines that forms a cycle. If 3 lines
    forms a cycle it will not flag the 3 lines as a cycle. The 2 lines that forms a cycle must not contain other none closed lines.
    The input lines must be topologically clean. Once a line is used to create a cycle it cannot be used again to form another cycle.


    Usage:
        import algo_flag_cycle

    Limits and constraints


"""
__revision__ = "--REVISION-- : $Id: algo_flag_cycle.py 452 2011-09-27 12:51:58Z dpilon $"

#####################################################################################################################################

import math

from shapely.geometry import Point, LineString, Polygon
from shapely.prepared import prep

from lib_genmetal import MA_Point, MA_LineString, MA_Polygon, InvalidParameterError, GenUtil,\
                         GenStatistics, PointErrorPosition, LineStringErrorPosition, IterationResults, Parameters,\
                         SpatialContainer, Algorithm

########################################################

#Internal key word constants
_ALGO = 'Flag cycle'
_CYCLES_FLAGGED = 'Cycles flag'
_ID_CYCLE = "id_cycle"


# If psyco is availiable, use it to speed up processing (2x)
try:
    import psyco
    psyco.full()
except ImportError:
    pass

class Statistics(GenStatistics):
    """Class that contains the statistics for the DP algorithm

    Attributes
        stat_names: Name of the statistics for the SpikeStatistics class. These name are
                    used by the Statistics class

    """

    def __init__(self):
        """Initialize the attributes of an object of the class"""

        GenStatistics.__init__(self)
        self.stats_names = ((_ALGO, GenUtil.SIMPLE_LINE, GenUtil.CROSSING_LINE, GenUtil.SIDEDNESS, _CYCLES_FLAGGED))

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
                str_out.append("Cycles flagged: " + str(self.get_stats_name_count_iter( _CYCLES_FLAGGED, i)))
        str_out.append( "Summary statistics" )
        str_out.append("Cycles flagges: " + str(self.get_stats_name_count_total( _CYCLES_FLAGGED)))
        str_out.append( "Number of iteration: " + str(self.get_nbr_iteration()) )

        return str_out

class AlgoFlagCycle(Algorithm):
    """
    This is the main class for the Flag cycle algorithm

    Attributes:
        - params: Parameters for the algorithm
    """

    def __init__(self, debug=False):
        """Initialize the attributes of an object of the class FlagCycleAlgorithm

        Parameters:
            debug: Flag to enable(TRUE)/disable(FLASE) for debug output 

        Return value:
            None

        """

        Algorithm.__init__(self)
        
        self.params = Parameters()
        self.params.debug=debug
        
        self.stats = Statistics()
        
        # Reset the counter of id cycles to 0
        self.id_cycle = 0

    def _set_line_attributes (self):
        """This routine sets the attributes to the line

        The routine checks:
            - if the first/last vertices are the same the line is closed; otherwise it is open

        Parameters: None

        Return value: None

        """

        # Select only the lines
        for line in self.s_container.get_features():
            line.ma_properties[_ID_CYCLE] = 0
            if (GenUtil.distance(line.coords_dual[0], line.coords_dual[-1]) <= GenUtil.ZERO):
                line.is_closed = True
            else:
                line.is_closed = False

        return

    def flag_cycle(self, line_to_cycle):
        """
        Method called to detect if the lines forms a cycle

        Parameters:
            line_to_cycle: Line feature to find a cycle
            
        Return value:
            True: The line contains a cycle
            False: The line do not contains a cycle

        """

        cycle_flagged = False
        
        if (not line_to_cycle.is_closed):            
            first_a = line_to_cycle.coords_dual[0]
            last_a = line_to_cycle.coords_dual[-1]
            
            x = first_a[0]
            y = first_a[1]
            dist = GenUtil.ZERO*100 
            bounds = (x-dist, y-dist, x+dist, y+dist)
            lines = self.s_container.get_features(bounds, filter="not feature.is_closed", remove_keys=[line_to_cycle._sci_id])
            for line in lines:
                first_b = line.coords_dual[0]
                last_b = line.coords_dual[-1]
                # Check if the combination of first/last vertice form a cycle (a closed area
                if ( (GenUtil.distance(first_a, first_b) <= GenUtil.ZERO and
                      GenUtil.distance(last_a, last_b) <= GenUtil.ZERO) or
                     (GenUtil.distance(first_a, last_b) <= GenUtil.ZERO and
                      GenUtil.distance(last_a, first_b) <= GenUtil.ZERO) ):
                    # Check if both lines are not identical
                    reverse_coords = list(line.coords)
                    reverse_coords.reverse()
                    if (not ( line_to_cycle.almost_equals(line, 6) or
                              line_to_cycle.almost_equals(LineString(reverse_coords)) )):
                        # Check if there are unclosed lines within the polygon formed by the 2 lines
                        if  (GenUtil.distance(first_a, last_b) <= GenUtil.ZERO):
                            closed_line = list(line_to_cycle.coords_dual) + list(line.coords_dual)
                        else:
                            closed_line = list(line_to_cycle.coords_dual) + reverse_coords
                        polygon = Polygon(closed_line)
                        # remove the 2 line from the list of lines
                        potential_lines = [line_tmp for line_tmp in lines if line_tmp not in (line_to_cycle,line) ]
                        contained = False
                        for potential_line in potential_lines:
                            de_9im = polygon.relate(potential_line)
                            if ( de_9im[0] in ("T","1") ):
                                contained = True
                                break
                        if not contained:
                            self.id_cycle += 1
                            line_to_cycle.ma_properties[_ID_CYCLE] = self.id_cycle
                            line_to_cycle.is_closed = True 
                            line.ma_properties[_ID_CYCLE] = self.id_cycle
                            line.is_closed = True
                            cycle_flagged = True
                            self.stats.add_stats(_CYCLES_FLAGGED)
                            break
                        
            return cycle_flagged

    def process(self):

        GenUtil.print_debug (self.params, "Start of %s simplification algorithm)" %(_ALGO) )

        
        self.s_container = self.load_features()

        if (self.params.debug):
            #Test if print is needed before scanning the s_container for nothing... waste of time...
            nbr_lines = len(self.s_container.get_features()) 
            GenUtil.print_debug (self.params, "Number of lines imported: %s"  %(nbr_lines) )
        
        self._set_line_attributes()

        cycle_flagged = True

        while (cycle_flagged):
            # At each new iteration we create a new itearation
            self.stats.add_iteration()
            GenUtil.print_debug (self.params, 'Start of iteration # %s' %(self.stats.get_nbr_iteration()) )
            
            cycle_flagged = False

            for line in self.s_container.get_features(filter="not feature.is_closed"):

                cycle_flagged = self.flag_cycle(line) or cycle_flagged

            GenUtil.print_debug (self.params, 'Number of points removed %s: ' %(self.stats.get_stats_name_count_iter(_ALGO)))
            GenUtil.print_debug (self.params, 'End of iteration # %s' %(self.stats.get_nbr_iteration()) )

        self.features = [feature for feature in self.s_container.get_features()]
        
        GenUtil.print_debug (self.params, "End of %s" % (_ALGO))
