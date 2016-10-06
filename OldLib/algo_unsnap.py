#!/usr/local/bin/python
# -=- encoding: utf-8 -=-
#####################################################################################################################################

"""
    This algorithm implement the unsnap algorithm.

    This algorithm will move apart 2 vertices of 2 different line when they are sharing the same spatial position.  The algorithm
    This algorithm will not move the first/last vertice unless the line is closed.
    When the vertices are moved it checks that it do not create the line to self cross, cross with another line, 
    creates a sidedness problem
    
    If 2 vertices of the same line share the same position the 2 vertices are not moves apart.


    Usage:
        import algo_unsnap

    Limits and constraints


"""
__revision__ = "--REVISION-- : $Id: algo_unsnap.py 458 2011-09-27 12:55:10Z dpilon $"

#####################################################################################################################################

import math

from shapely.geometry import Point, LineString, Polygon
from shapely.prepared import prep

from lib_genmetal import MA_Point, MA_LineString, MA_Polygon, InvalidParameterError, GenUtil,\
                         GenStatistics, PointErrorPosition, LineStringErrorPosition, IterationResults, Parameters,\
                         SpatialContainer, Algorithm

########################################################

#Internal key word constants
_ALGO = 'Unsnap'
_UNSNAPPED = 'Unsnapped'

_DISTANCE = 'distance'

class Statistics(GenStatistics):
    """Class that contains the statistics for the DP algorithm

    Attributes
        stat_names: Name of the statistics for the SpikeStatistics class. These name are
                    used by the Statistics class

    """

    def __init__(self):
        """Initialize the attributes of an object of the class"""

        GenStatistics.__init__(self)
        self.stats_names = ((_ALGO, GenUtil.SIMPLE_LINE, GenUtil.CROSSING_LINE, GenUtil.SIDEDNESS, _UNSNAPPED))

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
                str_out.append("Vertices removed: " + str(self.get_stats_name_count_iter( _ALGO, i)))
                str_out.append( "--------" )
                str_out.append( "Conflicts:" )
                str_out.append( "    Simple Line   :  " + str(self.get_stats_name_count_iter( GenUtil.SIMPLE_LINE, i)))
                str_out.append( "    Crossing Line :  " + str(self.get_stats_name_count_iter( GenUtil.CROSSING_LINE, i)))
                str_out.append( "    Sidedness     :  " + str(self.get_stats_name_count_iter( GenUtil.SIDEDNESS, i)))
                str_out.append( "--------" )
        str_out.append( "Summary statistics" )
        str_out.append("Vertices unsnapped: " + str(self.get_stats_name_count_total( _UNSNAPPED)))
        str_out.append( "--------" )
        str_out.append( "Conflicts:" )
        str_out.append( "    Simple Line   :  " + str(self.get_stats_name_count_total( GenUtil.SIMPLE_LINE)))
        str_out.append( "    Crossing Line :  " + str(self.get_stats_name_count_total( GenUtil.CROSSING_LINE)))
        str_out.append( "    Sidedness     :  " + str(self.get_stats_name_count_total( GenUtil.SIDEDNESS)))
        str_out.append( "--------" )
        str_out.append( "Number of iteration: " + str(self.get_nbr_iteration()) )

        return str_out



class AlgoUnsnap(Algorithm):
    """
    This algorithm implement the unsnap algorithm.

    This algorithm will move apart 2 vertices of 2 different line when they are sharing the same spatial position.  The algorithm
    This algorithm will not move the first/last vertice unless the line is closed.
    
    If 2 vertices of the same line share the same position the 2 vertices are not moves apart.

    Attributes:
        - params: Parameters for the algorithm
    """

    def __init__(self, test_crossing_line=True,
                       test_simple_line=True, 
                       test_sidedness=True,
                       split_size=50, 
                       debug=False):
        """Initialize the attributes of an object of the class

        Parameters:
            test_crossing_line: Flag to enable(TRUE)/disable(FLASE) checking the line crossing constraint
            test_simple_line: Flag to enable(TRUE)/disable(FLASE) checking the simple line constraint
            test_sidedness:  to enable(TRUE)/disable(FLASE) for checking the sidedness constraint
            debug: Flag to enable(TRUE)/disable(FLASE) for debug output 

        Return value:
            None

        """

        Algorithm.__init__(self)
        
        self.params = Parameters()
        self.params.test_crossing_line=test_crossing_line
        self.params.test_simple_line=test_simple_line 
        self.params.test_sidedness=test_sidedness
        self.params.split_size=split_size 
        self.params.debug=debug
        
        self.stats = Statistics()
        
    def _split_lines(self):
        """This routine reloads the input line into smaller line
        
        To speed up the search process we split the input line into smaller lines.
        A the limit we could split the line into individual point.  This solution
        is working but is to memory intensive.
        
        Parameters: None
        
        Return value: SpatialContainer containing the splitted lines
        
        """
        
        # Create the spatial container that will receive all the spatial features
        s_container_split_lines = SpatialContainer()

        # Loop over each LineString
        for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING"):
                
            coords = list(line.coords_dual)
            
            # Check if we must split the line
            if (len(coords) > self.params.split_size):
                # Creates a list of list of coords of maximum size of split_size
                split_lines = [coords[x:x+self.params.split_size] for x in xrange(0, len(coords), self.params.split_size)]
            else:
                split_lines = [coords]
            
            # Managing if the last list is of size 1
            if (len(split_lines[-1]) == 1):
                # The last list element contains only one coord.  This coord is added to the previous list element 
                split_lines[-2].extend(split_lines[-1])
                # Delete the last list element
                del split_lines[-1]
                
            lst_line_split_line = []
            last_split = len(split_lines)-1
            for i, split_line in enumerate(split_lines):
                line_split_line = MA_LineString(split_line)
                
                # Add some properties the the split_line
                line_split_line.ref_line_i = (i*self.params.split_size)  # Index of the first coord in link with the original line
                line_split_line.ref_line = line        # Reference to the original line
                lst_line_split_line.append(line_split_line)
            
            # This reference is used when updating the firs/last vertice.  We need a reference to the first split
            lst_line_split_line[-1].ref_first_split_line = lst_line_split_line[0]
                
            # Add the vertice in the container
            s_container_split_lines.add_features(lst_line_split_line)
            
        return (s_container_split_lines)
        
    def _set_line_attributes (self):
        """This routine sets the attributes to the line

        The routine checks:
            - if the first/last vertices are the same the line is closed; otherwise it is open

        Parameters: None

        Return value: None

        """

        # Select only the lines
        for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING"):
            if (GenUtil.distance(line.coords_dual[0], line.coords_dual[-1]) <= GenUtil.ZERO):
                line.is_closed = True
            else:
                line.is_closed = False

        return
    
    def check_modification (self, p1, p2, p2_new, p3, split_line, i_split_line, line, i_line):
        """Performs the constrtaint validation and update the line
        
        Checks if the position of the line formed by p1, p2, p3 will break any constraint (simple line, crossing line and sidednesn. If no constraints are
        violated update the coordinates.
        
        Parameters:
            p1, p2, p3: Tuple of coordinates forming the
            split_line: The MA_LineString containing the coordinate being modified
            i_split_line: index of the coordinate being modified in the split_line object
            line: The MA_LineString containing the coordinate being modified
            i_line: index of the coordinate being modified in the split_line object
        
        Return value:
        
        """
        
        line_crossing_line = LineString([p1, p2_new, p3])
        sidedness_polygon = GenUtil.calculate_sidedness_polygon ( LineString([p1, p2, p3]), 
                                                                  LineString([p1, p2_new, p3]) )
        new_line_coords =  list(line.coords_dual)
        
        if (i_line == len(line.coords_dual)-1):
            # Processing the last vertice of a closed line; so we modify the first and the last vertice
            new_line_coords[0] = p2_new
            new_line_coords[-1] = p2_new
        else:
            new_line_coords[i_line] = p2_new
            
        line_simple_line = LineString(new_line_coords)
        conflict_type = GenUtil.test_constraints (self, None, line_simple_line, line_crossing_line, 
                                                  sidedness_polygon, self.s_container, line._sci_id)
            
        # Update the coordinates
        if conflict_type is None:
            # New coordinate list of the split line
            new_split_line_coords = list(split_line.coords_dual)
             
            if (i_line == len(line.coords_dual)-1):
                # Special case when i_line is -1 is the last vertice of the line
                # We must update the last coordinate of the current split_line
                new_split_line_coords[-1] = p2_new
                
                # Update the coordinates
                split_line.update_coords(new_split_line_coords, self.s_container_split_lines)
                
                # We must also update the first coordinate of the first split line of this line
                first_split_line = split_line.ref_first_split_line
                new_split_line_coords = list(first_split_line.coords_dual)
                new_split_line_coords[0] = p2_new
                
                # Update the coordinates of the first split_line
                first_split_line.update_coords(new_split_line_coords, self.s_container_split_lines)
                
            else:
                
                new_split_line_coords[i_split_line] = p2_new
                
                # Update the coordinate of the split line
                split_line.update_coords(new_split_line_coords, self.s_container_split_lines)
            
            # Update the coordinate of the line
            line.update_coords(new_line_coords, self.s_container)
            edit = True
            self.stats.add_stats(_UNSNAPPED)
        else:
            edit = False
            
        return edit
    
    def unsnap(self, split_line, i_split_line, first_last):
        """
        Method to unsnap vertices

        Parameters:
            split_line: Contains the coordinate to process
            i_split_line: index of the coordinate to process in split_line
            first_last: Flag to enable (True) or disable (False) the processing of a first/last vertice of a closed line
            
        Return value:
            True: The line is unsnap
            False: The line is not unsnap

        """
            
        x = split_line.coords_dual[i_split_line][0]
        y = split_line.coords_dual[i_split_line][1]
        
        line = split_line.ref_line # Reference to the original line
        dist = line.ma_properties[_DISTANCE]*1.2
        bounds = (x-dist, y-dist, x+dist, y+dist)
        neighbours = self.s_container_split_lines.get_features(bounds)
        for neighbour in neighbours:
            # The vertice must not be part of the line being checked
            if line != neighbour.ref_line:
                #Loop over the coordinates of the neighbour line to check if the vertice touch one of the vertice of the neightbour
                touch = False
                for coord in neighbour.coords_dual:
                    if (GenUtil.distance(split_line.coords_dual[i_split_line], (coord[0],coord[1]) ) <= GenUtil.ZERO):
                        touch = True
                        break
                
                # If they touch process for the unsnapping
                if (touch):
                    i_line = split_line.ref_line_i + i_split_line
                    p2 = line.coords_dual[i_line]
                    if (first_last):
                        p1 = line.coords_dual[-2]
                        p3 = line.coords_dual[1]
                    else:
                        p1 = line.coords_dual[i_line-1]
                        p3 = line.coords_dual[i_line+1]

                    p1_p3_mid = GenUtil.mid_point(p1, p3)
                    height = GenUtil.distance(p1_p3_mid, p2) 
                    if ( height >= GenUtil.ZERO):
                        scale_factor = line.ma_properties[_DISTANCE]/height
                        p2_new = GenUtil.rescale_vector(p2, p1_p3_mid, scale_factor)
                        
                        edit_line =  self.check_modification (p1, p2, p2_new, p3, split_line, i_split_line, line, i_line)
                        if not edit_line:
                            p2_new = GenUtil.rescale_vector(p2, p1_p3_mid, -scale_factor)
                            edit_line =  self.check_modification (p1, p2, p2_new, p3, split_line, i_split_line, line, i_line)
                            
                        if edit_line:
                            break

    def check_features(self):
        """
        Check if the features passed in parameters are of the good class type and have the good attributes
        
        Parameters: None
        
        Return value: None
        
        """
        
        # Check the line string
        properties_to_check = [_DISTANCE]
        for feature in self.features:
            if isinstance(feature, MA_LineString):
                GenUtil.check_feature_integrity(feature, MA_LineString, properties_to_check)        

    def _manage_split_line(self, split_line):
        """This method manage the processing of one split line
        
        The method is looping over each coordinate of the split line. 
        The first and the last vertice of each original line (not the split) are nver process unless it is the first/last vertice
        of a closed line.
        
        Parameter:
            split_line: MA_LineString object containing part of an original that must be processed
            
        Return value: None
        
        """
        line = split_line.ref_line  # Reference to the original line
        last_coord = len(line.coords_dual)-1   # Index of the last vertice in the originasl line
        # Loop over each coordinate
        for i, coord in enumerate(split_line.coords_dual):
            i_line = split_line.ref_line_i + i  # Index of the coor in the original line
            first_last = False
            
            # We nerver process or try to unsnap the first and last vertice of a line except if its a closed line
            if (i_line == 0):
                process = False
            else:
                if (i_line == last_coord):
                    if (line.is_closed):
                        process = True
                        # Processing the first/last vertice of a closed line
                        first_last = True
                    else:
                        process = False
                else:
                    process = True
            
            if (process): 
                self.unsnap(split_line, i, first_last)

    def _del_cycling_reference(self):
        """Delete the cycling reference
        
        If the cycling reference are not deleted the garbage collector is not able to reclaim the memory.
        Temporary attributes are also removed
        
        Parameter: None
            
        Return value: None
        
        """
        
        for split_line in self.s_container_split_lines.get_features():
            if hasattr(split_line, "ref_line"):   del split_line.ref_line
            if hasattr(split_line, "ref_line_i"): del split_line.ref_line_i
            if hasattr(split_line, "ref_first_split_line"): del split_line.ref_first_split_line
            
            
        for line in self.s_container.get_features():
            if hasattr(line, "is_closed"):   del line.is_closed

    def process(self):

        # Check the feature's class and attributes
        self.check_features()
        
        # Load the shapely features into the spatial container
        self.s_container = self.load_features()
        
        self._set_line_attributes()
        
        self.s_container_split_lines = self._split_lines()

        line_unsnapped = True

        while (line_unsnapped):
            # At each new iteration we create a new itearation
            self.stats.add_iteration()
            GenUtil.print_debug (self.params, 'Start of iteration # %s' %(self.stats.get_nbr_iteration()) )
            
            split_lines = self.s_container_split_lines.get_features()
            for split_line in split_lines:
                self._manage_split_line(split_line)
                
            line_unsnapped = False

            GenUtil.print_debug (self.params, 'Number of vertices moves %s: ' %(self.stats.get_stats_name_count_iter(_UNSNAPPED)))
            GenUtil.print_debug (self.params, 'End of iteration # %s' %(self.stats.get_nbr_iteration()) )
            
        # Destruction of the cyclic references and garbage collection
        self._del_cycling_reference()
        
        self.extract_features_out(self.s_container.get_features())
        
        del self.s_container_split_lines
        del self.s_container
        
        GenUtil.print_debug (self.params, "End of %s" % (_ALGO))
