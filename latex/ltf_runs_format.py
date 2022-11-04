##################################################################################
# format file associated to latex_table_formatter
# keep a generic version of this file for later use
# Documentation at the end. 
##################################################################################

# do not modify if you don't do it in the main script as well
def tablecolumn (datacol, format, alignement="c", source="file", data=[]):
    """
    Inputs : 
        * datacol : integer in 0,1,... what is the number of the column in relative order.
        * format  : format string, for instance "%.3f". 
        * alignement : string in ["l", "c", "r"]
        * source : "file" or "custom"
        * data : if source == custom, what to put in the table (ugly)
    """
    col = {"datacol" : datacol, "format" : format, "alignement" : alignement, "source" : source, "data" : data}
    return col

#########################
# source files parameters
#########################

sep="\t"
comment="#"

# syntax : [filepath, rows, cols, sep, comment]
sources = [["runs_parameters",[":"],[":"],sep,comment]]

#########################
# Table parameters
#########################

"""
    params      erreurs
    Nx, Nv, Nt, norme inf de E - Eex, norme 1 de E - Eex

"""

columns = [
    "Run", tablecolumn(0, "%s"), 
    "$N_x$", tablecolumn(1, "%d"), 
    "$N_{v_i}$", tablecolumn(2, "%d"), 
    "$N_{v_e}$", tablecolumn(3, "%d"), 
    "$\Delta t$", tablecolumn(4, "%f"), 
]

caption     = "" # Optional (ignore the warning if commented)
label       = "tab:runs"      # Optional (ignore the warning if commented)
vertical_lines_inside_blocks = False # !! weird in recursive headings
double_hlines = [] # numbering starts at 0

#########################
# output parameters 
#########################

save_as_compilable_tex = False
save_as_includable_tex = True
generate_pdf = False
print_in_terminal = True
output_filename = "ltf_runs"
verbose = False

##################################################################################
#                                                                            End #
##################################################################################

"""                        DOCUMENTATION

########################## sources syntax ##########################

Indicates where to look for data. List of raw data files and a few parameters.

* filepath : anything that can fit in open(filepath,"r").

* rows : one of the following
  - ":", "a:", ":b", "a:b" : all / starting from a / skipping last b / starting 
        from a and skipping last b rows (a,b in 0,1,...)
  - [a,b,...,o] : list of elements e, with e being
        either a number n : pick the corresponding row
        or a tuple (a,b)  : pick rows between a and b, both included (can be 
        reversed, like (5,3), or equal, like (42,42))

* cols : same specs as rows. The order generated is a "relative" ordering of    
    columns, THAT WILL BE USED TO NUMBER COLUMNS IN THE TABLE. 
    (If cols=":" everywhere, the order of the files is kept.)

* sep : column separator of the datafile.

* comment : lines that begins by this will be discarded 
    AND NOT COUNTED IN THE NUMBERING.

########################## columns syntax ##########################

columns = [
    "column name 1", tablecolumn (relative n°, "format"),
    "column name 2", [
        "subcolumn name a" : tablecolumn (relative n°, "format"),
        "subcolumn name b" : tablecolumn (relative n°, "format"),
    ],
    "column name 3", tablecolumn (relative n°, "format"),
    ...
    "column name n", [
        "subcolumn name a" : tablecolumn (relative n°, "format"),
        ...
        "subcolumn name z" : tablecolumn (relative n°, "format"),
    ],
]

A column is, by default, the concatenation of all data read in files.
The user can eventually add a column with "custom" and plug forgotten data, 
by calling tablecolumn(randomNumber, format, "custom", [row1, row2, ..., rowN])

The relative order of columns defined by cols is used here.
For instance, if file A contains variables a,b in columns 2 and 28,
and file B contains variables a,b in columns 1 and 0 respectively,
sources = [["path_to_A", rowsA, [2,28], sepA],
           ["path_to_B", rowsB, [1, 0], sepB]]
and the relative ordering is {a,b}, the relative number of a is 0
and of b is 1. 

"""
