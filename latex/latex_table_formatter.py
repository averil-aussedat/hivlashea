#######################################################################################
# Latex table formatter
#######################################################################################

import numpy as np
from os import system
from sys import argv
from importlib import import_module

funcname = "latex_table_formatter"
warn = "[%s] WARNING : " % (funcname)

####################################################################
# Récupération des infos du format file

args = argv[1:]
if (len(args)==0):
    print(warn + "Please provide a format_file. Abort")
    exit (1)
elif (len(args)==1):
    format_file_name = args[0].replace(".py","")
else:
    print ("Too many arguments : ", args, ". Was expected a format file only.")
    exit (1)

try:
    fmt = import_module(format_file_name)
except Exception as e:
    print(warn + "error opening format file '%s' : " % (format_file_name), e, ". Abort")
    exit (1)
####################################################################

####################################################################
# Vérification que le format file est bien conforme

# list of variables (just to check existence and type)
variables = [] # syntax : name, is_mandatory, default, type
variables.append (["sources", True, [], list])
variables.append (["columns", True, [], list])
variables.append (["caption", False, "", str])
variables.append (["label", False, "", str])
variables.append (["vertical_lines_inside_blocks", False, True, bool])
variables.append (["double_hlines", False, [], bool])
variables.append (["save_as_compilable_tex", False, False, bool])
variables.append (["save_as_includable_tex", False, True, bool])
variables.append (["generate_pdf", False, False, bool])
variables.append (["print_in_terminal", False, False, bool])
variables.append (["output_filename", False, "ltf_output", str])
variables.append (["verbose", False, False, bool])

for variable, is_mandatory, default, thetype in variables:
    try:
        thevar = getattr (fmt,variable)
    except Exception as e:
        if (is_mandatory):
            print(warn + variable + " is mandatory. Abort")
            exit(1)
        else:
            print("[%s] %s not found in format file : setting default " % (funcname, variable), default)
            setattr(fmt,variable,default)
        if thetype==bool: 
            try:  
               setattr(fmt,bool(thevar))
            except Exception as e:
                print(warn + "could not cast %s '" % variable, thevar, "' as bool (", e, "). Abort")
                exit(1)
        elif (type(thevar) is not thetype):
            print(warn + variable + " is not of type ", thetype, " : ", thevar, ". Abort")
            exit(1)

# columns check
def columns_format_check (col, is_ok):
    if (np.mod(len(col),2) != 0):
        is_ok = False
        print(warn + "uneven number of elements in a level of columns")
    else:
        for ic, c in enumerate(col[1::2]):
            if ((type(col[2*ic]) is not str) or (not ((type(c) is list) or (type(c) is dict)))):
                is_ok = False
                print(warn + "bad type in columns. Use strings for name and dicts/lists of columns for datas")
            if (is_ok and (type(c) is list)):
                is_ok = columns_format_check (c, is_ok)
    return is_ok

if (not columns_format_check (fmt.columns, True)):
    print(warn + "'columns' is not well formatted. Abort")
    exit (1)
if (len(fmt.columns)==0):
    print(warn + "empty 'columns'.")

# sources check
sourcesvars = [] # syntax : name, type
sourcesvars.append(["filepath", str])
sourcesvars.append(["rows", list])
sourcesvars.append(["cols", list])
sourcesvars.append(["sep", str])

try:
    for source in fmt.sources:
        for ivar, sourcevar in enumerate(sourcesvars):
            variable, thetype = sourcevar
            if (type(source[ivar]) is not thetype):
                print(warn + "sources contains invalid %s '" % variable, source[ivar], "'. Abort")
                exit(1)
except Exception as e:
    print(warn + "error while checking sources : ", e, ". Abort")
    exit(1)
####################################################################

if (fmt.verbose): 
    print("Welcome in %s." % funcname)
    print("You may need to \\usepackage{multirow} if funny 2* appear in the head.")

####################################################################
# Extraction des infos pratiques du format file
def recur_count (columns, ncols, nrows, current_row):
    for c in columns[1::2]:
        if type(c) is dict:
            ncols += 1
            nrows = np.max([nrows, current_row])              
        elif type(c) is list:
            [ncols, nrows, current_row] = recur_count (c, ncols, nrows, current_row+1)
        else: 
            print(warn + "error in extracting infos from column. Please check format.")
    return [ncols, nrows, current_row-1]

[ncols, nrows, current_row] = recur_count (fmt.columns, 0, 1, 1)

lines_head = [[] for _ in range(nrows)]
rows_dict = []

def recur_lines (columns, lines_head, rows_dict, current_row, level_width, nrows):
    for ic, c in enumerate(columns[1::2]):
        if type(c) is dict:
            # remplissage jusqu'au plancher de l'en-tête
            lines_head[current_row].append([columns[2*ic], nrows-current_row, 1])
            for irow in range(nrows-current_row-1):
                lines_head[current_row+1+irow].append([]) # laisser du blanc
            rows_dict.append (c)
            if (current_row > 0):
                level_width[-2] += 1 # ajouter la largeur de cette colonne à celle du dessus

        else: # si (comme ?) le boulot de vérification a bien été fait, type(c)=list ici
            [lines_head, rows_dict, current_row, level_width] = recur_lines (c, lines_head, rows_dict, current_row+1, level_width+[0], nrows)
            lines_head[current_row].append ([columns[2*ic], 1, level_width[-1]])
            if (current_row > 0):
                level_width[-2] += level_width[-1] # ajouter la largeur de cette colonne à celle du dessus
                level_width[-1] = 0 # passer à la colonne suivante
            else:
                level_width = [0] # passer à la colonne principale suivante : remise à 0
    return [lines_head, rows_dict, current_row-1, level_width[:-1]]

[lines_head, rows_dict, current_row, level_width] = recur_lines (fmt.columns, lines_head, rows_dict, 0, [0], nrows)
####################################################################

####################################################################
# Récupération des données
def indices_to_list (indices, varname, minn, maxx):
    """ minn included, maxx excluded """
    res = []
    for i in indices:
        if type(i) is tuple:
            [res.append(e) for e in np.arange(i[0],i[1]+1)]
        elif type(i) is int:
            res.append(i)
        elif type(i) is str:
            ii = (":"+i+":").split(":")
            try:
                ii[0] = minn   if (ii[0]=='') else int(ii[0])
                ii[1] = maxx-1 if (ii[1]=='') else int(ii[1])
            except Exception as e:
                print(warn + varname + " contains invalid element ", i, " (", e, "). Abort")
                exit(1)
            [res.append(e) for e in np.arange(ii[0],ii[1]+1)]
    res = [minn+np.mod(n-minn,maxx-minn) for n in res]
    return res

list_data_rows = [] # list of all the lines 
for filepath,rows,cols,sep,comment in fmt.sources:
    if (fmt.verbose):
        print("[%s] Reading from %s" % (funcname, filepath))
    real_row_counter = 0
    with open(filepath, "r") as datafile:
        rawlines = [rawline for rawline in datafile.readlines () if (len(rawline)>0 and rawline!=[""] and rawline[0]!=comment)]
        rawlines = [[numstr.replace("\n","") for numstr in rawline.split(sep)] for rawline in rawlines]
        rows = indices_to_list (rows, "rows", 0, len(rawlines))
        rawlines = [rawlines[irow] for irow in rows]
        for rawline in rawlines: 
            linecols = indices_to_list (cols, "cols", 0, len(rawline)) 
            rawline = [rawline[icol] for icol in linecols]
            line = []
            for f in rawline:
                try: # brutal, but very efficient
                    line.append(float(f))
                except:
                    line.append(f)
            list_data_rows.append(line)

lendatas = len(list_data_rows)

if fmt.verbose:
    print("list_data_rows : ", list_data_rows)
####################################################################

####################################################################
########## Coeur du script : construction de la table latex

output = "\\begin{table}[H]\n"
output += "\t\\centering\n"
colspec = [d["alignement"] for d in rows_dict]
if fmt.vertical_lines_inside_blocks:
    colspec = "|"+"|".join(colspec)+"|"
else:
    index = 0
    for cell in lines_head[0]:
        index += cell[2]
        colspec.insert(index, "|")
        index += 1 # shift due to |
    colspec = "|" + "".join(colspec)
output += "\t\\begin{tabular}{"+colspec+"} \\hline\n"

########## placement de l'en-tête
for iline, line in enumerate(lines_head):
    output += "\t"
    # si pas la première ligne : remplissage des lignes brisées (\cline)

    if (iline > 0):
        coldep = 1
        cline_open = False
        for iele, ele in enumerate(line):
            if (len(ele)>0): # pas vide
                if (cline_open):
                    colarr += ele[2] # continuation d'une cline ouverte : somme des largeurs
                else:
                    cline_open = True # ouverture d'une cline
                    colarr = coldep + 1
            else: # ele vide
                if (cline_open):
                    cline_open = False # fermeture d'une cline
                    output += "\\cline{%d-%d} " % (coldep, colarr-1)
                    coldep = colarr+1
                else:
                    coldep += 1 # recherche d'une future cline
        # fermeture d'une éventuelle cline laissée ouverte
        if (cline_open):
            output += "\\cline{%d-%d} " % (coldep, colarr-1)

    # remplissage des noms
    for ic, c in enumerate(line): # soit [], donc remplissage, soit [nom, lignes, colonnes]
        if (len(c)>1): 
            if (c[2] > 1): # multicolonnes
                if (ic==0):
                    output += "\\multicolumn{%d}{|c|}{%s} " % (c[2], c[0])
                else:
                    output += "\\multicolumn{%d}{c|}{%s} " % (c[2], c[0])
            elif (c[1] > 1): # multilignes
                output += "\\multirow{%s}*{%s} " % (c[1], c[0])
            else: # single output
                output += c[0] + " "
        if (ic < len(line)-1):
            output += "& "
        else:
            output += "\\\\ \n"
output += "\t\\hline \\hline \n" 

########## placement des données
for irow, therow in enumerate(list_data_rows):
    output += "\t"
    for id, d in enumerate (rows_dict):
        if (d["source"]=="file"):
            if (d["datacol"] >= len(therow)):
                print(warn + "asked for column %d but only %d in table. Abort" % (d["datacol"], len(therow)))
                exit(1)
            thevalue = therow[d["datacol"]]
        elif (d["source"]=="custom"):
            try:
                thevalue = d["data"][irow]
            except Exception as e:
                print("[%s] Error printing col %d of row %d (" % (funcname,id,irow), e, ". Abort")
                exit(1)
        else:
            print("[%s] Unknown source '%s'. Abort" % (funcname, d["source"]))
            exit(1)
        try:
            output += d["format"] % thevalue
        except Exception as e:
            print(warn + "format error (", e, "). Expected was %s, received " \
                % (d["format"]), thevalue, " of type", type(thevalue))

        if (id < ncols-1):
            output += " & "
        else:
            output += " \\\\ \\hline " + ("\\hline" if irow in fmt.double_hlines else "") + "\n"

########## end of table
output += "\t\\end{tabular}\n"
if (len(fmt.caption)>0):
    output += "\t\\caption{" + fmt.caption + "}\n"
    # + "}\n"
if (len(fmt.label)>0):
    output += "\t\\label{" + fmt.label + "}\n"
# if (len(fmt.caption)>0):
#     output += "\t}\n" # pour la bonne prise en compte du label, il doit etre inclu dans le caption (dont know why...)
output += "\\end{table}\n"
####################################################################

####################################################################
# Outputs

if (fmt.save_as_compilable_tex or fmt.generate_pdf):
    cmpfname = fmt.output_filename + ("" if ".tex" in fmt.output_filename else ".tex")
    cmpfname = cmpfname[:-len(".tex")] + "_cmp.tex" # horrible, je sais
    with open(cmpfname, "w") as thefile:
        thefile.write ("\\documentclass{article}\n")
        thefile.write ("\\usepackage[landscape]{geometry}\n")
        thefile.write ("\\usepackage{multirow}\n")
        thefile.write ("\\usepackage{float}\n")
        thefile.write ("\\begin{document}\n")
        thefile.write ("\\thispagestyle{empty}\n")
        thefile.write (output)
        thefile.write ("\\end{document}\n")

if fmt.save_as_includable_tex:
    fname = fmt.output_filename + ("" if ".tex" in fmt.output_filename else ".tex")
    with open(fname, "w") as thefile:
        thefile.write (output)
    print("Creating file ", fname)

if fmt.print_in_terminal:
    print("\n")
    print(output)
    print("\n")

if fmt.generate_pdf:
    system ("pdflatex " + cmpfname)
    system ("rm %s.aux %s.log" % (cmpfname[:-4], cmpfname[:-4]))

if fmt.verbose:
    print("Bye")

#######################################################################################
#                                                                                 Fin #
#######################################################################################

