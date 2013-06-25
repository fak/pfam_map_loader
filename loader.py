"""Script:  loader.py

Load the mapping of Pfam-A domains. The main function is loader, defined at the bottom of the document. Specify release and version of the mapping on command line. eg.: $> python loader.py chembl_15 0_1

Note on variable names: lkp is used to represent dictionaries I was too lazy to find a proper name for.

--------------------
Author:
Felix Kruger
fkrueger@ebi.ac.uk

"""
import os
import sys
import queryDevice
import yaml


def readfile(path, key_name, val_name):
    '''Read two columns from a tab-separated file into a dictionary.

    Inputs:
    path -- filepath
    key_name -- name of the column holding the key
    val_name -- name of the column holding the value

    '''
    infile = open(path, 'r')
    lines = infile.readlines()
    infile.close()
    lkp = {}
    els = lines[0].rstrip().split('\t')
    for i, el in enumerate(els):
        if el == key_name:
            key_idx = i
        if el == val_name:
            val_idx = i
    for line in lines[1:]:
        elements = line.rstrip().split('\t')
        lkp[elements[key_idx]] = elements[val_idx]
    return  lkp



def retrieve_acts(dom_string, params):
    """Run a query for act_id, tid, component_id, compd_id and domain_name.

    Inputs:
    dom_string -- A string specifying the domain names eg. "7tm_1','Pkinase','Pkinase_tyr"
    params -- dictionary holding details of the connection string

    """
    acts = queryDevice.queryDevice("""
    SELECT DISTINCT act.activity_id, ass.tid, tc.component_id, cd.compd_id, dm.domain_name
                      FROM activities act
                      JOIN assays ass
                          ON ass.assay_id = act.assay_id
                      JOIN target_dictionary td
                          ON ass.tid = td.tid
                      JOIN target_components tc
                          ON ass.tid = tc.tid
                      JOIN component_domains cd
                          ON tc.component_id = cd.component_id
                      JOIN domains dm
                          ON dm.domain_id = cd.domain_id
                     WHERE ass.assay_type IN('B','F')
                     AND td.target_type IN('PROTEIN COMPLEX', 'SINGLE PROTEIN')
                     AND act.standard_relation ='='
                     AND ass.relationship_type = 'D'
                     AND act.standard_type IN(
                       'Ki', 'Kd', 'IC50', 'EC50', 'AC50',
                       'log Ki', 'log Kd', 'log IC50', 'Log EC50', 'Log AC50'
                       'pKi', 'pKd', 'pIC50', 'pEC50', 'pAC50'
                        )
                     AND dm.domain_name IN('%(dom_string)s')
                     """ %locals() ,params['release'], params['user'], params['pword'], params['host'], params['port'])
    return acts



def map_ints(acts):
    """ Map interactions to activity ids.

    Inputs:
    acts -- output of the sql query in retrieve_acts().

    """
    lkp = {}
    for act in acts:
        (act_id, tid, component_id, compd_id, domain_name) = act
        try:
            lkp[act_id][compd_id]=domain_name
        except KeyError:
            lkp[act_id] ={}
            lkp[act_id][compd_id]=domain_name
    return lkp



def flag_conflicts(lkp):
    """Flag conflicts.

    Input:
    lkp -- a lookup dictionary of the form lkp[act_id][compd_id] = domain_name

    """
    flag_lkp  = {}
    for act_id in lkp.keys():
        if len(lkp[act_id].keys()) == 1:
            flag_lkp[act_id]=(0,0)
        elif len(lkp[act_id].keys()) > 1:
            if len(set(lkp[act_id].values())) > 1:
                flag_lkp[act_id] = (2,0)
            elif len(set(lkp[act_id].values())) == 1:
                flag_lkp[act_id] = (1,0)
    return(lkp, flag_lkp)


def apply_manual(lkp, flag_lkp, infile):
    """Apply manually resolved conflicts.

    Input:
    lkp -- a dictionary of the form lkp[act_id][compd_id] = domain_name
    flag_lkp -- a dictionary of the form flag_lkp[act_id] = (conflict_flag, manual_flag)
    path -- a filepath to the resolved_conflicts table

    """
    infile = open(path, 'r')
    lines = infile.readlines()
    infile.close()
    for line in lines[1:]:
        (act_id, compd_id, domain_name, comment) = line.rstrip().split('\t')
        try:
            lkp[act_id][compd_id][domain_name]
            flag_lkp[act_id] = (0,1)
        except KeyError:
            print 'lost entry:', line
    return(lkp, flag_lkp)


def write_table(lkp, flag_lkp, path):
    """ Write a table containing activity_id, domain_id, tid, conflict_flag, type_flag.

    Input:
    lkp -- a dictionary of the form lkp[act_id][compd_id] = domain_name
    flag_lkp -- a dictionary of the form flag_lkp[act_id] = (conflict_flag, manual_flag)
    path -- a filepath to the output file

    """
    out = open(path, 'w')
    out.write("map_id\tactivity_id\tcompd_id\tdomain_name\tconflict_flag\tmanual_flag\n")
    counter = 0
    for act_id in lkp.keys():
        compd_ids = lkp[act_id]
        (conflict_flag, manual_flag) = flag_lkp[act_id]
        for compd_id in compd_ids.keys():
            counter +=1
            domain_name = lkp[act_id][compd_id]
            out.write("%(counter)i\t%(act_id)i\t%(compd_id)i\t%(domain_name)s\t%(conflict_flag)i\t%(manual_flag)i\n"%locals())
    out.close()



def upload_sql(params):
    """ Load SQL table using connection string defined in global parameters.

    Input:
    params -- dictionary holding details of the connection string.

    """
    status = os.system("cp data/pfam_maps_v_%(version)s.tab data/pfam_maps.txt" % params)
    if status != 0:
        sys.exit("Error copying data/pfam_maps_v_%(version)s.tab to data/pfam_maps.txt" % params)
    status = os.system("mysql -u%(user)s -p%(pword)s -h%(host)s -P%(port)s -e 'DROP TABLE %(release)s.pfam_maps'" % params)
    status = os.system("mysql -u%(user)s -p%(pword)s -h%(host)s -P%(port)s -e 'CREATE TABLE pfam_maps(map_id INT NOT NULL AUTO_INCREMENT, activity_id INT, compd_id INT, domain_name VARCHAR(100), conflict_flag INT, manual_flag INT, PRIMARY KEY (map_id))' %(release)s"% params)
    if status != 0:
        sys.exit("Error creating table pfam_maps." % params)
    os.system("mysqlimport -u%(user)s -p%(pword)s -h%(host)s -P%(port)s --ignore-lines=1 --lines-terminated-by='\n' --local %(release)s data/pfam_maps.txt" % params)
    if status != 0:
        sys.exit("Error loading table pfam_maps.""" % params)


def loader(release, version):
    """Main function to load the mapping of Pfam-A domains.

    Inputs from command line:
    release -- chembl release eg "chembl_12"
    version -- version of the mapping eg 0_1

"""
    # Read config file.
    paramFile = open('mpf.yaml')
    #paramFile = open('mpf_local.yaml')
    #paramFile = open('mpf_example.yaml')
    params = yaml.safe_load(paramFile)

    # Load the list of validated domains.
    domains = readfile('data/valid_pfam_v_%(version)s.tab' % locals(), 'pfam_a', 'pfam_a')
    dom_string = "','".join(domains.keys())

    # Get activities for domains.
    acts  = retrieve_acts(dom_string, params)

    # Map interactions to activity ids.
    lkp = map_ints(acts)

   # Flag conflicts.
    (lkp, flag_lkp) = flag_conflicts(lkp)

    # Apply manually resolved conflicts.
    infile = 'data/resolved_conflicts_v_%(version)s.tab' % locals()
    (lkp, flag_lkp) = apply_manual(lkp, flag_lkp, infile)

    # Write a table containing activity_id, domain_id, tid, conflict_flag, type_flag
    outfile = 'data/pfam_maps_v_%(version)s.tab' %locals()
    write_table(lkp, flag_lkp, outfile)

    # Load SQL table.
    upload_sql(params)

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 3:  # the program name and the two arguments
        sys.exit("Must specify release (eg. chembl_15) and version (eg. 0_1).")
    release = sys.argv[1]
    version = sys.argv[2]

    loader(release, version)
