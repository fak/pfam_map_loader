"""Script:  loader.py

Load the mapping of Pfam-A domains. The main function is loader, defined at the bottom of the document. Specify release and version of the mapping on command line. eg.: $> python loader.py chembl_15 0_1

Note on variable names: lkp is used to represent dictionaries I was too lazy to find a proper name for.

--------------------
Author:
Felix Kruger
fkrueger@ebi.ac.uk
"""
import os
import subprocess
import sys
import yaml
import pgQuery

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




def retrieve_acts_psql(domains, params):
    """Run a query for act_id, tid, component_id, compd_id and domain_name. 
       This is to identify all activities associated with any given valid
       domain. These activities are then processed with the map_ints and 
       flag_conflicts function.

    Inputs:
    dom_string -- A string specifying the domain names eg. "7tm_1','Pkinase','Pkinase_tyr"
    params -- dictionary holding details of the connection string

    """
    
    acts = pgQuery.paramquery("""
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
                     AND dm.domain_name IN %(domains)s
                     """ ,locals() ,params )
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
    """Assign a set of flags to each activity: category, status, manual.

    Input:
    lkp -- a lookup dictionary of the form lkp[act_id][compd_id] = domain_name

    """
    flag_lkp  = {}
    for act_id in lkp.keys():
        if len(lkp[act_id].keys()) == 1:
            flag_lkp[act_id]=(0,0,0) # one validated domain.
        elif len(lkp[act_id].keys()) > 1:
            if len(set(lkp[act_id].values())) > 1:
                flag_lkp[act_id] = (2,1,0) # multiple validated domains.
            elif len(set(lkp[act_id].values())) == 1:
                flag_lkp[act_id] = (1,1,0) # multiple instances of one val. dom.
    return flag_lkp



def write_table(lkp, flag_lkp, manuals, params, path):
    """ Write a table containing activity_id, domain_id, tid, conflict_flag, type_flag.

    Input:
    lkp -- a dictionary of the form lkp[act_id][compd_id] = domain_name
    flag_lkp -- a dictionary of the form flag_lkp[act_id] = (conflict_flag, manual_flag)
    path -- a filepath to the output file

    """
    out = open(path, 'w')
    out.write("""activity_id\tcompd_id\tdomain_name\tcategory_flag\tstatus_flag\tmanual_flag\tcomment\ttimestamp\tsubmitter\n""")
    for act_id in set(lkp.keys()) - set(manuals.keys()): # Not processing maunal maps.
        compd_ids = lkp[act_id]
        (category_flag, status_flag, manual_flag) = flag_lkp[act_id]
        comment = params['comment']
        timestamp = params['timestamp']
        submitter = params['submitter']
        for compd_id in compd_ids.keys():
            domain_name = lkp[act_id][compd_id]
            out.write("""%(act_id)i\t%(compd_id)i\t%(domain_name)s\t%(category_flag)i\t%(status_flag)i\t%(manual_flag)i\t%(comment)s\t%(timestamp)s\t%(submitter)s\n"""%locals())
    out.close()




def upload_psql(params):
    """ Load SQL table using connection string defined in global parameters.

    Input:
    params -- dictionary holding details of the connection string.

    """
    status = subprocess.call("psql -U%(user)s  -h%(host)s -p%(port)s -d%(release)s -c 'DROP TABLE IF EXISTS pfam_maps'" % params, shell=True)
    status = subprocess.call("psql -U%(user)s  -h%(host)s -p%(port)s -d%(release)s -c 'CREATE TABLE pfam_maps(map_id INT, activity_id INT, compd_id INT, domain_name VARCHAR(100), category_flag INT, status_flag INT, manual_flag INT, comment VARCHAR(150), timestamp VARCHAR(25)), submitter VARCHAR(150)' "% params, shell=True)
    if status != 0:
        sys.exit("Error creating table pfam_maps." % params)
    params['path'] = ('/').join([subprocess.check_output('pwd', shell=True).rstrip(), 'data', 'pfam_maps.txt'])
    subprocess.call("psql -U%(user)s  -h%(host)s -p%(port)s -d%(release)s -c \"COPY pfam_maps FROM \'%(path)s\' DELIMITER \'\t\'  \"" % params, shell=True)
    if status != 0:
        sys.exit("Error loading table pfam_maps.""" % params)


def append_table(tables, outfile):
    with open(outfile, 'w') as outfile:
        prev = open(tables[0])
        prev_header = prev.readline()
        map_id = 0
        for table in tables:
            with open(table) as infile:
                header = infile.readline()
                if header != prev_header:
                    sys.exit('input tables are not same format')
                for line in infile:
                    map_id += 1
	            outfile.write('\t'.join([str(map_id),line]))


def loader():
    """Main function to load the mapping of Pfam-A domains.

    Inputs from command line:
    release -- chembl release eg "chembl_12"
    version -- version of the mapping eg 0_1

"""
    # Read config file.
    param_file = open('local.yaml')
    #param_file = open('example.yaml')
    params = yaml.safe_load(param_file)
    param_file.close()

    # Load the list of validated domains.
    domains = readfile('data/valid_pfam_v_%(version)s.tab' % params, 'pfam_a', 'pfam_a')
    dom_string = "','".join(domains.keys())
    domains = tuple(domains.keys())
    # Load a list of manually edited activities.
    manuals = readfile('data/manual_pfam_maps_v_%(version)s.tab' % params, 'activity_id', 'manual_flag')

    # Get activities for domains.
    acts  = retrieve_acts_psql(domains, params)

    # Map interactions to activity ids.
    lkp = map_ints(acts)

    # Flag conflicts.
    flag_lkp = flag_conflicts(lkp)

    # Write a table containing activity_id, domain_id, tid, conflict_flag, type_flag
    write_table(lkp, flag_lkp, manuals, params, 'data/automatic_pfam_maps_v_%(version)s.tab' %params)
    append_table(['data/manual_pfam_maps_v_%(version)s.tab' %params, 'data/automatic_pfam_maps_v_%(version)s.tab' % params], 'data/pfam_maps.txt')

    # Load SQL table.
    upload_psql(params)

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 1:  # the program name and the two arguments
        sys.exit("All parameters are specified in local.yaml or example.yaml, depending on line 173+174 ")

    loader()
