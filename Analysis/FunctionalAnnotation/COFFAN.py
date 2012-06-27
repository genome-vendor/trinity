'''
==============================================================================
COFFAN - Comprehensive Flexible Functional Annotator 
@Created on Feb 3, 2012
@author: Tyler Weirick, Brian Couger 
@contact: tyler.weirick@gmail.com, mbcouger@gmail.com
This program is designed to integrate a wide range functional annotation data 
into a single aggregated data source in order to facilitate ascertaining 
higher biological meaning. The program was created with the knowledge that
research projects will have variable data requirements depending on the 
organism of interest or condition of analysis and thus places a high 
emphasizes on having a robust number of input options available to the user 
while remaining flexible on data integration output. With this flexibility in 
mind the program can be used for a variety or genomic research and is not 
limited to any specific scientific discipline or domain of life.  Below is a 
full description of the annotation input data that can be integrated:
'''
"""
Notes to developers: 
If you are reading this you are likely planning to customize this program to 
fit your research needs. This section will give an over view of the program 
structure and explain the functionality to simplify code modification.

-To make changes to the code
Changes with the intent to add additional parsing abilities or remove parsing 
abilities will require changes in 3 places in the code. 
1. The function which takes command line arguments.
2. The classes the hold parsed data. 
3. The section of the program where execution is done. 

1. All that is needed in this section is to add 
"""

__version__ = "$Revision: 00f8e3bb1197 $"
# $Source$

from optparse import OptionParser
from xml.etree import ElementTree as ET
from inspect import isclass
from sys import exit
from os import path
from collections import deque 

##############################################################################
#                                 Global Variables
##############################################################################

delimiter = '`'
class_name_split_str = "__"

"""
These are specifically for flat text files. 
"""
class LONGEST_KEGG():
    def __init__(self):
        self.LONGEST = 0
    def checkandsetlongest(self,check_val):
        if check_val > self.LONGEST:
            self.LONGEST = check_val
    
class LONGEST_PFAM():
    def __init__(self):
        self.LONGEST = 0
    #def longest(self):
    #    return self.LONGEST_PFAM
    
LK = LONGEST_KEGG()
LP = LONGEST_PFAM()
##############################################################################
#                                 Classes
##############################################################################

class pfaminfo():
    def __init__(self,query_name="query_name None",target_name="target_name None",
                 evalue="evalue None",target_desc="target_desc None"):      
        #@todo need to confirm which version of output this is written for
        #self.a__pfam_id              = id
        self.b__pfam_query_name      = query_name
        self.c__pfam_target_name     = target_name
        self.d__pfam_full_seq_evalue = evalue
        self.e__pfam_target_desc     = target_desc
        
class kegginfo():
    def __init__(self):
        #self.kegg_name = "None"
        self.a__kegg_identifier = "No Kegg identifier"
        self.b__kegg_pathway    = "No Kegg Pathway"
        
    def setkeggidentifier(self,kegg_identifier):
        self.a__kegg_identifier = kegg_identifier
            
    def setkeggpathway(self,kegg_pathway):
        self.b__kegg_pathway = kegg_pathway
        
class outputclass:
    """
    This class is intended to be used as a master data collection location. 
    It is also the basis for the simplification of modification of this 
    program. This class is written to be read by the output functions in a 
    way that allows the structure of the class to determine what is given as
    output. 
    
    """
    def __init__(self):
        global LK 
        global LP 
        #Fasta file from Trinity output. 
        self.a__fasta = "None"
        #BlastX XMl output
        self.b__blastx_query_len = "None"
        self.c__blastx_hit_id = "None"
        self.d__blastx_hit_def = "None"
        self.e__blastx_accesion = "None"
        self.f__blastx_length =  "None"
        self.g__blastx_hsp_evalue =  "None"
        #RSEM flat text file
        self.h__rsem_comp_id =  "rsem_comp_id None"
        self.i__rsem_raw_alignemnt_cnt = "None"
        self.j__rsem_tau_value =  " rsem_tau_value None"
        #Holds list of PFAM classes
        self.k__longest_pfam = LP
        self.k1__PFAM_MODELS_FOUND = ""
        #This list must contain a class to header printing purposes.
        self.l__pfam_entry = [ pfaminfo() ]
        
        #Holds list of KEGG classes.
        self.m__longest_kegg = LK
        self.m1__KEGG_MODELS_FOUND = ""
        #This list must contain a class to header printing purposes.
        self.n__kegg_class_list = [ kegginfo() ]  
          
        self.o__thmmm_len      = "thmmm_len None"
        self.p__thmmm_ExpAA    = "hmmm_ExpAA None"
        self.q__thmm_First60   = "None"
        self.r__thmm_Pred_Hel  = "None"
        self.s__thmmm_topology = "thmmm_topology None"
        #The value is secreted or not.
        self.t__signalp_secreted = "sp None"
        #After this is set no other attributes may be added to the class. 
        self.initialized = True
    
    def __setattr__(self, item, value):
        if (not 'initialized' in self.__dict__):
            return dict.__setattr__(self,item, value)
        elif  item in self.__dict__:
            return dict.__setattr__(self,item, value)
        else:
            print("ERROR")
            raise AttributeError(item)
    
    def updateentrynumbs(self):
        """
        This fuction is a patch to add the functionality to print the length 
        of variable length entries.
        """
        if len(self.l__pfam_entry) == 1:
            self.k1__PFAM_MODELS_FOUND = "0"
        else:
            self.k1__PFAM_MODELS_FOUND = str(len(self.l__pfam_entry))
        
        if len(self.n__kegg_class_list) == 1:
            self.m1__KEGG_MODELS_FOUND = "0"   
        else:
            self.m1__KEGG_MODELS_FOUND = str(len(self.n__kegg_class_list))
    
    
    
    
    
    def setfastaname(self,fasta_name):
            self.a__fasta = fasta_name
        
    def setblastxdata(self,blastx_query_len,blastx_hit_id,blastx_hit_def,
                      blastx_accesion,blastx_length,blastx_hsp_evalue ):       
        self.b__blastx_query_len = blastx_query_len
        self.c__blastx_hit_id = blastx_hit_id
        self.d__blastx_hit_def = blastx_hit_def
        self.e__blastx_accesion = blastx_accesion
        self.f__blastx_length =  blastx_length
        self.g__blastx_hsp_evalue =  blastx_hsp_evalue
            
    def setkeggclasslist(self,kegg_classes):
        self.n__kegg_class_list = kegg_classes
    
##############################################################################
#                             Functions
##############################################################################

def printflatfile(class_for_title):
    """
    This function will generate the title for the delimited columns of a flat 
    text file. Returns the title as a string.
    """
    global delimiter    
    #todo: Make recursive sub-class handling to allow more flexibility.
    title_out_list = []

    title_dic = class_for_title.__dict__
    title = sorted(title_dic.keys())
    print_iterations = 0
    for el in title:  
        if len(el.split(class_name_split_str)) > 1:
            el_type = type(title_dic[el])
            
            if(isclass(el_type) and el_type != type(str()) and 
            el_type != type(list()) and el_type != type(bool()) ):
                print_iterations = title_dic[el].LONGEST
            elif el_type == type([]):
                for i in range(0,print_iterations):
                    temp_title = []
                    for sub_el in title_dic[el]:
                        for term_el in sorted(sub_el.__dict__.keys()):
                            term_minus_order_tag = term_el.split("__")
                            if len(term_minus_order_tag) > 1:
                                temp_title.append(term_minus_order_tag[1]+"_"+str(i)+delimiter)
                            else:
                                temp_title.append(term_el+"_"+str(i)+delimiter)
                                
                    temp_str = ''.join(temp_title)
                    title_out_list.append(temp_str)
            else:
                #If not list or class append then title term append to list.
                term_minus_order_tag = el.split("__")
                if len(term_minus_order_tag) > 1:
                    title_out_list.append(term_minus_order_tag[1]+delimiter)
                else:
                    title_out_list.append(el+delimiter)
    title_out_list.append("\n")         
    return ''.join(title_out_list)


def printdataforflatfile(main_output_dic):
    """
    This function takes the main dictionary and will convert the information 
    of the class held in the dictionaries value section into delimited text.
    """
    global delimiter    
    title_out_list = []
    #Use the dictionary keys to iterate through the dictionary.
    for main_key in main_output_dic.keys():
        #This prevents the printing of non-desired class attributes. 
        #if True:#len(main_key.split(class_name_split_str)) > 1:
        #Get class from the main dictionary.
        output_class = main_output_dic[main_key]
        output_values_dic = output_class.__dict__
        output_values_key_list = sorted(output_values_dic)
        print_iterations=0#This handles data sets that can differ in length.
        
        #Update the values for the lengths of variable length entries. 
        output_class.updateentrynumbs()
        
        for class_key in output_values_key_list: 
            if len(class_key.split(class_name_split_str)) > 1:
                el_value = output_values_dic[class_key]
                el_type = type(el_value)
                
                if(isclass(el_type) and el_type != type(str()) and 
                el_type != type(list()) and el_type != type(bool()) ):
                    #These will specifically recognize the PL and KL classes.
                    #Since the dict was sorted the next entry will be an array 
                    #of variable length data classes.
                    print_iterations = el_value.LONGEST
                    #This is an addition to print the number of entries for variable length entries.
                   
                    
                elif el_type == type(list()):
                    # All lists in the main class should contain classes.
                    sub_class_list = el_value
                    temp_data = []
                    total_data_length = 0
                    for sub_class in sub_class_list:
                        sub_class_dict = sub_class.__dict__
                        sub_class_keys = sorted( sub_class.__dict__.keys() )
                        varibles_per_class = len(sub_class_keys)
                        for sub_class_key in sub_class_keys:
                            temp_data.append( sub_class_dict[sub_class_key]+delimiter )
                            total_data_length+=1
                    title_out_list.append(''.join(temp_data))
                    #If the class structure is followed correctly a class giving
                    #the maximum number of entries should be encountered before a  
                    #list of classes is encountered. 
                    add_x_nodata_entires = varibles_per_class*print_iterations - total_data_length
                    for i in range(0,add_x_nodata_entires):    
                        title_out_list.append("No Data"+delimiter)
                else:
                    #If not list or class append then title term append to list.
                    title_out_list.append(str(el_value)+delimiter)
        title_out_list.append("\n")     
    return ''.join(title_out_list)



def generalXMLtext(main_output_dict):
    global delimiter
    """
    This function returns an XML string of the data contained in the main 
    dictionary class.
    """    
    XML_DECLARATION = '<?xml version="1.0" encoding="UTF-8" ?>'
    MAIN_ENTRY_TAG  = "entry_info"
    output_XML_list = [XML_HEADER]
    for main_key in main_output_dict.keys():
        
        output_XML_list.append(starttag(MAIN_ENTRY_TAG)+"\n")

        for key in sorted(main_output_dict[main_key].keys()):
            """
            non-standard classes can be ignored.
            """
            key_tag = key.split(delimiter)
            output_XML_list.append(starttag()+"\n")
            
            output_XML_list.append(closetag()+"\n")
        
        output_XML_list.append(closetag(MAIN_ENTRY_TAG)+"\n")



def starttag(tag):
    return '<'+tag+'>'

def closetag(tag):
    return '</'+tag+'>'


def getargs(ver='%prog 0.2'):
    """Gets file names for input and output 
    will use the first name and append out on the end of the file
    if no out-file name is given"""
    
    desc_list = []
    start_and_break = "'''"
    read_line_bool = False
    for line in open(__file__,'r'):
        if read_line_bool:
            if not start_and_break in line:
                desc_list.append(line.replace("\n","")+"\n")
            else:
                break    
        if (start_and_break in line) and read_line_bool == False:
            read_line_bool = True
    desc = ''.join(desc_list)
    
    troubleShoot = False
    parser = OptionParser(version=ver,description=desc)
    
    #@todo: OptionParser is depreciated in Python 3.2. 
    #Need to move to the new style of parser. 

    parser.add_option("-a", "--main_fasta_name", 
                      dest="main_fasta_name", 
                      default="False",
                      help = "REQUIRED: The name of the main fasta file.")

    parser.add_option("-b", "--blastx_name", 
                      dest="blastx_name", 
                      default="False",
                      help = "Name of a blastx file.")
    
    #@todo: output from funtion
    parser.add_option("-c", "--max_evalue", 
                      dest="max_evalue",
                      default="1",
                      help = "Set this value to filter out blastx entries"+
                      ' below the given value. Can be given as a string ex:'+ 
                      ' "1e-40" or a float 1e-40'
                      )
    
    parser.add_option("-d", "--rsem_name", 
                      dest="rsem_name", 
                      default="False",
                      help = "Name of a RSEM file.")
    
    parser.add_option("-e", "--best_canidate_file_name", 
                      dest="best_canidate_file_name", 
                      default="False",
                      help = "Required for: PFAM, SignalP, THMMM")
    
    parser.add_option("-f", "--pfam_name", 
                      dest="pfam_name", 
                      default="False",
                      help = "Name of PFAM file.")
    
    parser.add_option("-g", "--kegg_file_name", 
                      dest="kegg_file_name", 
                      default="False",
                      help = "Name of KEGG file.")
    
    parser.add_option("-i", "--kegg_dict_file_name", 
                      dest="kegg_dict_file_name", 
                      default="False",
                      help = "Name of KEGG term dictionary file.")
    
    parser.add_option("-j", "--signalp_file_name", 
                      dest="signalp_file_name", 
                      default="False",
                      help = "Name of a SignalP file.")
    
    parser.add_option("-k", "--thmmm_file_name", 
                      dest="thmmm_file_name", 
                      default="False",
                      help = "Name of THMMM file.")
    
    parser.add_option("-l", "--out_file_name", 
                      dest="out_file_name", 
                      default="False",
                      help = "REQUIRED: Name for a file to output to. "+
                      "CAUTION: If a file with the same name already "+
                      "exists it will be over-written")
        
    (options, args) = parser.parse_args()
    if(troubleShoot):print(options);print(args)
            
    return options.__dict__


def getfastanames(fasta_name): 
    """
    This function is designed to read the fasta files in the output of Trinity
    and extract the names of each fasta in the file. For each name found the 
    program will exclude the first space or tab and anything that occurs after
    the space and remove the > character from the beginning of the fasta name. 
    
    Ex. of Trinity fasta name format.
    >comp295_c0_seq2 len=884 path=[0:0-628 1736:629-842 1950:843-883]
    >comp295_c0_seq1 len=1661 path=[0:0-628 629:629-1660]
    
    It will then add each 
    name as the key in the a dictionary and set the value of the dictionary 
    to None.
    """
    output_dic = {}
    
    for line in open(fasta_name,'r'):
        if line[0] == '>':
            temp_class = outputclass()
            temp_class.a__fasta  = line.split()[0][1:]
            output_dic.update( {line.split()[0][1:]:temp_class} )
    return output_dic
        
         
def getblastxinfo(blastx_name,out_dic,max_evalue):    
    #@todo: Will blast always have one hsp value?
    
    ID_tag = "Iteration_query-def"
    test_key = ""
    for event,elem in ET.iterparse(blastx_name,events=("start",)):
        
        if elem.tag == ID_tag:
            test_key = elem.text.split()[0]
            if test_key in out_dic:
                temp_class = out_dic[test_key]
            else:
                print("WARNING: BlastX entry found which is not in fasta file.")
                test_key = ""
            
        if test_key != "":
            if elem.tag == "Iteration_query-len":
                temp_class.b__blastx_query_len = elem.text
            elif elem.tag == "Hit_id":
                temp_class.c__blastx_hit_id = elem.text  
            elif elem.tag == "Hit_def":
                temp_class.d__blastx_hit_def = elem.text
            elif elem.tag == "Hit_accession":
                temp_class.e__blastx_accesion = elem.text
            elif elem.tag == "Hit_len":
                temp_class.f__blastx_length = elem.text
            elif elem.tag == "Hsp_evalue":
                
                if temp_class.g__blastx_hsp_evalue == "None":
                    if max_evalue >= float(elem.text):
                        temp_class.g__blastx_hsp_evalue = elem.text
                        out_dic[test_key] = temp_class
                    else:
                        temp_class.g__blastx_hsp_evalue = "Greater than threshold."
                        out_dic[test_key] = temp_class
                     
    return out_dic


def getrsem(rsem_name,out_dic):
    tau_value_multiplier = 1.0e6
    
    for line in open(rsem_name,'r'):
        line_list = line.split()
        if len(line_list) != 4:
            print("ERROR: There must be four columns per rsem entry") 
            print("ERROR Encountered at:",line)
            exit()
        
        if not line_list[0] in out_dic:
            print("WARNING: ",line_list[0],"Not found in RSEM but not in fasta data.")
        else:
            temp_class = out_dic[ line_list[0] ]
            temp_class.i__rsem_raw_alignemnt_cnt = line_list[1]
            temp_class.j__rsem_tau_value         = str( tau_value_multiplier*float(line_list[2]) )
            temp_class.h__rsem_comp_id           = line_list[3]
            out_dic[ line_list[0] ] = temp_class
            
    return out_dic
        
def getnamesfrombestcanidates(best_canidate_name):
    """
    This will make a dictionary data structure consisting of the *.* style 
    names as the key and the fasta style names as the value.
    EX:
    >m.2 g.2  ORF g.2 m.2 comp1136_c0_seq1:3-656(+) 
    -> {'m.2':'comp1136_c0_seq1'}
    @change: 3.4.2012 
    Added and additional output, a dictionary that stores the same key value 
    as the first dictionary, but stores a boolean that reflects the condition 
    of the first character of the fasta entry being an M or not. This is 
    meant to be used to filter out possible erroneous values for signalP data
    as sequences that do not start at the n-terminal 
    """    
    output_dic = {} 
    n_termin_dict = {}
    first_fasta_line = False
    for line in open(best_canidate_name,'r'):
        if line[0] == '>':
            #If an fasta name line is encountered split the line to get a key
            #value pair.
            key = line.split()[0][1:]
            value = line.split()[5].split(":")[0]
            output_dic.update( {key:value} )
            first_fasta_line = True
        elif first_fasta_line:
            #If the first line of a fasta entry is encountered 
            if len(line) > 0:    
                first_fasta_line = False
                if line[0] == "M"or line[0] == "m":
                    n_termin_dict.update(  {key:True} )
                else:
                    n_termin_dict.update(  {key:False} )
            else:
                print("WARNING: Empty line encountered after fasta name in best "+ 
                      "candidate file.")
        
    return output_dic,n_termin_dict
                       
def getpfam(pfam_name,best_canidate_dic,out_dic):
    global LP 
    longest_pfam_counter = 0
    last_entry_name = ""
    temp_list = []
    
    for line in open(pfam_name,'r'):
        
        if line[0] !='#':#Do not add header lines.
            line_list = line.split()
            if line_list[2] != last_entry_name:#If these are not equal then we have a new 
                if last_entry_name != "":

                    if len(temp_list) > LP.LONGEST:
                        LP.LONGEST = len(temp_list)
                    
                    out_dic[ best_canidate_dic[last_entry_name] ].l__pfam_entry  = temp_list
                    
                    temp_list = []
                last_entry_name = line_list[2]
            
            target_name    = line_list[0]
            accession      = line_list[1]
            query_name     = line_list[2]
            accession1     = line_list[3]
            e_value        = line_list[4]
            score          = line_list[5]
            bias           = line_list[6]
            E_value2       = line_list[7]
            score          = line_list[8]
            bias           = line_list[9]
            exp            = line_list[10]
            reg            = line_list[11]
            clu            = line_list[12]
            ov             = line_list[13]
            env            = line_list[14]
            dom            = line_list[15]
            rep            = line_list[16]
            inc            = line_list[17]
            desc_of_target = ' '.join(line_list[18:])
            #Could just pass values at instantiation, but this is easier to 
            #understand.
            temp_pfam_class = pfaminfo()        
            temp_pfam_class.b__pfam_query_name      = query_name
            temp_pfam_class.c__pfam_target_name     = target_name
            temp_pfam_class.d__pfam_full_seq_evalue = e_value
            temp_pfam_class.e__pfam_target_desc     = desc_of_target
            
            temp_list.append(temp_pfam_class)
            
    return out_dic 
       
def keggparser(kegg_file_name,out_dic,bc_id_dic,kegg_term_dic):        
    """
    Kegg files are flat text files containing one entry per line. Each line 
    contains a name followed by a space and zero or more identifiers delimited 
    with commas.
    """
    global LK    
    #todo: add the kegg dic
    for line in open(kegg_file_name,'r'):
        #Remove newline char and split around space
        
        line_list = line.strip().split()
        #This will leave us with the identifier and the kegg codes.  
        temp_kegg_class_list = []   
        if len(line_list) < 2:
            #An entry with no identifiers.
            temp_kegg_class_list.append(kegginfo())    
        elif len(line_list) == 2:
            #Then one or more kegg identifiers are present.
            kegg_id_list = line_list[1].split(',')
            #For computing delimited text file spacing.
            LK.checkandsetlongest(len(kegg_id_list))
            
            for kegg_id in kegg_id_list:
                temp_kegg_class = kegginfo() 
                temp_kegg_class.setkeggidentifier(kegg_id)
                temp_kegg_class.setkeggpathway(kegg_term_dic[kegg_id])
                temp_kegg_class_list.append(temp_kegg_class)
        else:
            print("Kegg error")
        #add list of kegg classes to main class 
        out_dic[bc_id_dic[line_list[0]]].setkeggclasslist(
                                          temp_kegg_class_list)   
        
    
def make_kegg_dictionary(kegg_dic_file_name):
    
    kegg_dic = {}
    for line in open(kegg_dic_file_name,'r'):
        split_line_list = line.split()
        key = split_line_list[0]
        data = ' '.join(split_line_list[1:])
        kegg_dic.update({key:data})
    return kegg_dic

def extractsignalpdata(signalp_file_name,main_dic,translator_dic,n_termin_dict):
    """
    This function reads an abridged signalp file, parses the signalp data and 
    adds the information in the (?) column to the main class. 
    """
    
    for line in open(signalp_file_name,'r'):
        #Filter out title and start of file the contains no info.
        if len(line) != 0 and line[0] != '\n' and  line[0] != '#':
            #Convert line into list by chopping in whitespace areas.
            sp_data = line.split()

            #Just in case a text editor converts leading lines to spaces.
            if len(sp_data) > 9:
                #Use translator dictionary to look up main name inorder to look 
                #up the class to add the data to.  t__signalp_qmark
                if n_termin_dict[ sp_data[0] ]:
                    main_dic[translator_dic[sp_data[0]]].t__signalp_secreted = sp_data[9]
                else:
                    main_dic[translator_dic[sp_data[0]]].t__signalp_secreted = "N-terminal not found."
    return main_dic
            
def thmmmparser(thmmm_file_name,main_dic,translator_dic):
    """
    This function reads an abridged signalp file, parses the signalp data and 
    adds the information in the (?) column to the main class. 
    """
    
    for line in open(thmmm_file_name,'r'):
        #Filter out title and start of file the contains no info.
        if len(line) != 0 and line[0] != '\n' and  line[0] != '#':
            #Convert line into list by chopping in whitespace areas.
            thmmm_data = line.split()
            #Use translator dictionary to look up main name inorder to look 
            #up the class to add the data to. 
            if len(thmmm_data) > 5:
                main_dic[translator_dic[thmmm_data[0]]].o__thmmm_len      = thmmm_data[1]
                main_dic[translator_dic[thmmm_data[0]]].p__thmmm_ExpAA    = thmmm_data[2]
                main_dic[translator_dic[thmmm_data[0]]].q__thmm_First60   = thmmm_data[3]
                main_dic[translator_dic[thmmm_data[0]]].r__thmm_Pred_Hel  = thmmm_data[4]
                main_dic[translator_dic[thmmm_data[0]]].s__thmmm_topology = thmmm_data[5]
                            
    return main_dic        
    

def fileexists(file_name):
    """
    This function check to see if a file exists. If the file exists it will 
    output the file name and path. If the file does not exist, it will output
    an error message and exit. 
    This function is being used as the logical location location to check for 
    non-existent file, after arguments are given, contains data other that 
    files. Thus to keep from writing multiple if statments I wrote this. 
    """
    if path.exists(file_name):
        return file_name
    else:
        print('ERROR: the file "'+file_name+'" does not exist at the location', 
        'given.\n Exiting.')
        exit()
        
        
##############################################################################
#                          Start Main Program
##############################################################################
#Get command line arguments and return them as a dictionary.
#Returning things as a dictionary will help reduce the probability of errors 
#as it only requires changes in the arguments section and in the section of 
#actual use to update for new values. 
opt_dic = getargs(ver='%prog 0.1')

#The main fasta file is the basic for all other data at the moment so return 
#an error if it is not given or not included. 
if opt_dic["main_fasta_name"] != "False":
    out_dic = getfastanames(fileexists(opt_dic["main_fasta_name"]))
    if opt_dic["blastx_name"] != "False":
        try:
            max_evalue = float(opt_dic["max_evalue"])
        except:
            print("evalue must be a number")
            exit()
        out_dic = getblastxinfo(fileexists(opt_dic["blastx_name"]),
                                out_dic,
                                max_evalue)

        out_dic = getrsem(fileexists(opt_dic["rsem_name"]),out_dic)

    if opt_dic["best_canidate_file_name"] != "False":
        bc_id_dic,n_terminal_dict = getnamesfrombestcanidates(
                                    fileexists(opt_dic[
                                    "best_canidate_file_name"]))
    else:
        print("Skipping best candidate data collection.")
    
    if opt_dic["pfam_name"] != "False":
        if opt_dic["best_canidate_file_name"] != "False":
            out_dic = getpfam(fileexists(opt_dic["pfam_name"]),
                              bc_id_dic,
                              out_dic)
        else:
            print("No best candidates file found, therefore, PFAM data", 
                  "collection will be skipped.")
    else:
        print("No PFAM file specified. Skipping PFAM data collection.") 
    
    if opt_dic["signalp_file_name"] != "False":
        if opt_dic["best_canidate_file_name"] != "False":
            out_dic = extractsignalpdata(opt_dic["signalp_file_name"],out_dic,bc_id_dic,n_terminal_dict)
        else:
            print("No best candidates file found, therefore, SignalP data", 
                  "collection will be skipped")
    else:
        print("No SignalP file specified. Skipping SignalP data collection.") 
    
    if opt_dic["thmmm_file_name"] != "False":
        if opt_dic["best_canidate_file_name"] != "False":
            thmmmparser(opt_dic["thmmm_file_name"] ,out_dic,bc_id_dic)
        else:
            print("No best candidates file found, therefore, THMMM data", 
                  "collection will be skipped")
    else:
        print("No THMMM file specified. Skipping THMMM data collection.") 


    if opt_dic["kegg_dict_file_name"] != "False":
        kegg_term_dic = make_kegg_dictionary(
                        fileexists(opt_dic["kegg_dict_file_name"]))
    else:
        print("No Kegg Definition file specified.", 
              "Skipping Kegg Definition data collection.") 
            
    
    if (opt_dic["kegg_file_name"] != "False" and 
        opt_dic["kegg_dict_file_name"] != "False"):
        keggparser(fileexists(opt_dic["kegg_file_name"]),
                   out_dic,
                   bc_id_dic,
                   kegg_term_dic)
    else:
        print("Skipping Kegg data collection.")

else:
    print("ERROR: a main fasta file must be specified.Exiting.")   
    exit()

title = printflatfile(outputclass()) 
data  = printdataforflatfile(out_dic) 

out_file = open(opt_dic["out_file_name"] ,"w")
out_file.write(title+data)
out_file.close()


