"""Script:  exporter.py

Export mapping of Pfam-A domains and create a versioned backup of the pfam_maps table eg.: $> python loader.py chembl_15 0_1


--------------------
Author:
Felix Kruger
fkrueger@ebi.ac.uk

"""
import os
import sys
import time
import pgQuery
import yaml




def retrieve_acts(params):
    """Run a query to obtain manual mappings.

    Inputs:
    params -- dictionary holding details of the connection string

    """
    acts = pgQuery.query("SELECT * from pfam_maps WHERE manual_flag = 1 " ,params)
    return acts


def write_table(acts, path):
    """ Export the manual mappings into the manual_pfam_maps_v_x_x.tab file to be fed into the next round of curation.

    Input:
    acts -- results of the manual query
    path -- a filepath to the output file

    """
    out = open(path, 'w')
    out.write("""map_id\tactivity_id\tcompd_id\tdomain_name\tcategory_flag\tstatus_flag\tmanual_flag\tcomment\ttimestamp\n""")
    for act in acts:
        # map_id = act[0] this value is generated from scratch in load.py
        act_id = act[1]
        compd_id = act[2]
        domain_name = act[3]
        category_flag = act[4]
        status_flag = act[5]
        manual_flag = act[6]
        comment = act[7]
        timestamp = act[8]
        out.write("""%(act_id)i\t%(compd_id)i\t%(domain_name)s\t%(category_flag)i\t%(status_flag)i\t%(manual_flag)i\t%(comment)s\t%(timestamp)s\n"""%locals())
    out.close()





def bkp_sql(params):
    """ Load SQL table using connection string defined in global parameters.

    Input:
    params -- dictionary holding details of the connection string.

    """
    params['timestamp'] = time.strftime('%d_%B_%Y', time.gmtime())
    status = os.system("mysqldump -u%(user)s -p%(pword)s chembl_15 pfam_maps > data/pfam_maps_%(timestamp)s.bkp" % params)
    if status != 0:
        sys.exit("There was a problem dumping the table")



def exporter():
    """Main function to load the mapping of Pfam-A domains."""
    # Read config file.
    param_file = open('local.yaml')
    #param_file = open('example.yaml')
    params = yaml.safe_load(param_file)
    param_file.close()

    # Get activities for domains.
    acts  = retrieve_acts(params)

    # Write activity on new manual_pfam_maps file
    path = 'data/manual_pfam_maps_v_%(nextversion)s.tab' %params
    write_table(acts, path)

    # Dump SQL table.
    bkp_sql(params)

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 1:  # the program name and the two arguments
        sys.exit("All parameters are specified in local.yaml or example.yaml, depending on line 173+174 ")

    exporter()
