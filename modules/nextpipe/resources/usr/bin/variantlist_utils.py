# -*- coding: utf-8 -*-
"""
Functions for processing CLC workbench and VEP output of PanCancer Panel

@author: PatrickBasitta
"""

import pandas as pd

def adjust_region_position(CLC_variant_track_data_PAN):
    """
    
    Parameters
    ----------
    CLC_variant_track_data_PAN : DataFrame 
        
        DESCRIPTION: contains CLC_PAN_data (csv file)
    
    Function
    --------
    Rename column "Regions" to "Position"
    Add new column "End Position"
    Get chromsomal position for column "End Position"

    Returns
    -------
    CLC_variant_track_data_PAN : DataFrame 
    
        DESCRIPTION: contains the column "Position" 
        and "End Position" with respective chromosomal positions

    """
    #------------------------------
    # Rename "Region" to "Position"
    #------------------------------
    CLC_variant_track_data_PAN = CLC_variant_track_data_PAN.rename\
                            (columns={"Region": "Position"})

    #-----------------------------                    
    # Add new column "End Position"
    #-----------------------------
    CLC_variant_track_data_PAN.insert(loc=2, column="End Position", value="")
    
    #------------------------------------------------------------------------
    # Get chromosomal postition and add position value to column "Position" and 
    # "End Position" 
    #------------------------------------------------------------------------
    for i in range(len(CLC_variant_track_data_PAN["Position"])):
        
        if ".." not in CLC_variant_track_data_PAN["Position"][i] and \
            "^" not in CLC_variant_track_data_PAN["Position"][i]:
                tmp_value_snv = CLC_variant_track_data_PAN["Position"][i]
                CLC_variant_track_data_PAN.loc\
                [CLC_variant_track_data_PAN.index[i], "End Position"] = tmp_value_snv
                
        if ".." in CLC_variant_track_data_PAN["Position"][i]:
                tmp_value_del = CLC_variant_track_data_PAN["Position"][i]
                tmp_value_del_start = tmp_value_del.split("..")[0]
                tmp_value_del_end = tmp_value_del.split("..")[1]
                CLC_variant_track_data_PAN.loc\
                [CLC_variant_track_data_PAN.index[i], "Position"] =  tmp_value_del_start
                CLC_variant_track_data_PAN.loc\
                [CLC_variant_track_data_PAN.index[i], "End Position"] = tmp_value_del_end
                
        if "^" in CLC_variant_track_data_PAN["Position"][i]:
                tmp_value_ins = CLC_variant_track_data_PAN["Position"][i]
                tmp_value_ins_start = tmp_value_ins.split("^")[0]
                tmp_value_ins_end = tmp_value_ins.split("^")[1]
                CLC_variant_track_data_PAN.loc\
                [CLC_variant_track_data_PAN.index[i], "Position"] = tmp_value_ins_start
                CLC_variant_track_data_PAN.loc\
                [CLC_variant_track_data_PAN.index[i], "End Position"] = tmp_value_ins_end
                
    return CLC_variant_track_data_PAN


def convert_coltype_str_to_float(columns_to_filter, clc_data_ReferenceAllele_NO):
    """
    
    Parameters
    ----------
    columns_to_filter : List
    
        DESCRIPTION: contains column names that has to be "type" converted from
        str to float (needed for subsequent filter step)
    
    clc_data_ReferenceAllele_NO : DataFrame
        
        DESCRIPTION: contains variant data without variants equal to 
        "Reference allele" 

    Returns
    -------
    clc_data_ReferenceAllele_NO : DataFrame
        
        DESCRIPTION:  contains "type" converted values (str to float)
            

    """
    #--------------------------
    # Convert str to float
    #--------------------------
    for i in range(len(columns_to_filter)):
     
        variants_lst_sheets_temp1 = clc_data_ReferenceAllele_NO\
            [columns_to_filter[i]].apply(lambda x: x.strip("'"))
        variants_lst_sheets_temp2 = variants_lst_sheets_temp1.apply\
                                                 (lambda x: x.replace(",","."))
        variants_lst_sheets_temp3 = variants_lst_sheets_temp2.to_frame()
        variants_lst_sheets_temp4 = variants_lst_sheets_temp3.astype("float")
        clc_data_ReferenceAllele_NO.loc[:,columns_to_filter[i]] = variants_lst_sheets_temp4
        
    return clc_data_ReferenceAllele_NO


def concatenate_https_clinvar(rs):
    """
    Concatenate https_clinvar link with rs number
    Function of function "clinar_link_list"
    
    Parameters
    ----------
    rs : pandas.core.series.Series
        
    Returns
    -------
    https_clinvar_link : List
    
    """
    return "https://www.ncbi.nlm.nih.gov/clinvar/?term="+rs


def unlist_https_clinvar(link):
    """
    Unlist result of function "concatenate_https_clinvar"
    Function of function "clinar_link_list"

    Parameters
    ----------
    link : List
       
    Returns
    -------
    i : str
       
    """
    for i in link:
        return i
    
    
def clinvar_link_list(link_lst, rs_idx, rs):
    """
    Generates a clinvar link using rs numbers    
    
    Parameters
    ----------
    link_lst : empty List
        
    rs_idx : List of rs indices (obtained from pandas.core.series.Series with
                                 non-ordered indices)
        
    rs : pandas.core.series.Series
        

    Returns
    -------
    link_lst : List of https_clinvar_links

    """
    for i, idx in enumerate(rs_idx):
   
       if pd.isna(rs[rs_idx[i]]) == False:
         
            tmp_link = list(map(concatenate_https_clinvar,\
                                rs[rs_idx[i]].split(",")))
            tmp_link = unlist_https_clinvar(tmp_link)
            link_lst.append(tmp_link)

       else:
            link_lst.append(rs[rs_idx[i]])
            
    return link_lst


def adjust_chr_pos_NM_VEP(VEP_data):
    """
    
    Parameters
    ----------
    VEP_data : DataFrame
        
        DESCRIPTION: contains data obtained from ensembl-vep (txt file)

    Returns
    -------
    VEP_data : DataFrame 
        
        DESCRIPTION: contains processed/adjusted "Chromosome", "Position" and
        "Feature" column needed for merging with CLC data

    """
    #---------------------------------------                      
    # Add column "Chromosome" and "Position"
    #---------------------------------------
    VEP_data.insert(loc=1, column="End Position", value="")
    VEP_data.insert(loc=1, column="Position", value="")
    VEP_data.insert(loc=1, column="Chromosome", value="")
    
    #---------------------------------------------------
    # Get chromosome and position from column "Location"
    #---------------------------------------------------
    for i in range(len(VEP_data["Location"])):
        if len(VEP_data["Location"][i].split("-")) == 1:
            tmp_split = VEP_data["Location"][i].split("-")
            tmp_split = tmp_split[0].split(":")
            VEP_data.loc[i, "Chromosome"] =  tmp_split[0]
            VEP_data.loc[i, "Position"] =  tmp_split[1]
            VEP_data.loc[i, "End Position"] =  tmp_split[1]
            
        if len(VEP_data["Location"][i].split("-")) == 2:
            tmp_split_two = VEP_data["Location"][i].split("-")
            tmp_split = tmp_split_two[0].split(":")
            VEP_data.loc[i, "Chromosome"] =  tmp_split[0]
            VEP_data.loc[i, "Position"] =  tmp_split[1]
            VEP_data.loc[i, "End Position"] =  tmp_split_two[1]

    for i in range(len(VEP_data["Feature"])):
        tmp_NM = VEP_data["Feature"][i].split(".")
        VEP_data.loc[i, "Feature"] =  tmp_NM[0]
        
    return VEP_data
