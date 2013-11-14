"""
    Script:  coverage.py
    Identifies domains that only occur in multi-domain proteins. The main
    script is master.
    --------------------
    Felix A Kruger
    momo.sander@ebi.ac.uk
"""
####
#### import modules.
####
import queryDevice
import operator
import yaml
import time
####
#### Load parameters.
####
paramFile = open('local.yaml')
params = yaml.safe_load(paramFile)
paramFile.close()




#### Define functions.
#-----------------------------------------------------------------------------------------------------------------------


def get_el_targets(params):
    """Query the ChEMBL database for (almost) all activities that are subject to the mapping. Does not conver activities expressed in log-conversion eg pIC50 etc. This function works with chembl_15 upwards. Outputs a list of tuples [(tid, target_type, domain_count, assay_count, act_count),...]
    """
    data = queryDevice.queryDevice("""
            SELECT DISTINCT dc.tid, dc.target_type, dc.dc, COUNT(DISTINCT act.assay_id), COUNT(DISTINCT activity_id)
            FROM assays ass
            JOIN(
                      SELECT td.tid, td.target_type, COUNT(cd.domain_id) as dc
                      FROM target_dictionary td
                      JOIN target_components tc
                        ON tc.tid = td.tid
		      JOIN component_sequences cs
			ON cs.component_id = tc.component_id
                      JOIN component_domains cd
 			ON cd.component_id = cs.component_id
                      WHERE td.target_type IN('SINGLE PROTEIN', 'PROTEIN COMPLEX')
                      GROUP BY td.tid
                     ) as dc
              ON dc.tid = ass.tid
            JOIN activities act
              ON act.assay_id = ass.assay_id
            WHERE act.standard_type IN('Ki','Kd','IC50','EC50', 'AC50')
            AND ass.relationship_type = 'D'
            AND assay_type IN('B')
            AND act.standard_relation IN('=')
            AND standard_units = 'nM'
            AND standard_value <= %s
            GROUP BY dc.tid ORDER BY COUNT(activity_id)""" % (int(params['threshold']) * 1000) , params)
    print "retrieved data for ", len(data), "tids."
    return data

#-----------------------------------------------------------------------------------------------------------------------


def readfile(path, key_name, val_name):
    """Read two columns from a tab-separated file into a dictionary.
    Inputs:
    path -- filepath
    key_name -- name of the column holding the key
    val_name -- name of the column holding the value
    """
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

#-----------------------------------------------------------------------------------------------------------------------

def get_archs(el_targets, pfam_lkp):
    """Find multi-domain architectures.
    Inputs:
    el_targets -- list of eligible targets
    """
    act_lkp = {}
    arch_lkp = {}
    dom_lkp = {}
    for ent in el_targets:
        try:
            doms = pfam_lkp[ent[0]]
        except KeyError:
            print "no doms in ", ent[0]
        arch = ', '.join(sorted(doms))
        try:
            arch_lkp[arch] += 1
            act_lkp[arch] += ent[4]
        except KeyError:
            arch_lkp[arch] = 1
            act_lkp[arch] = ent[4]
        if len(doms) <= 1:
            continue
        for dom in set(doms):
            try:
                dom_lkp[dom] += 1
            except KeyError:
                dom_lkp[dom] = 1
    return(arch_lkp, dom_lkp, act_lkp)

#-----------------------------------------------------------------------------------------------------------------------

def get_doms(tids, params):
    """Get domains for a list of tids.
    Inputs:
    el_targets -- list of eligible targets
    """
    pfam_lkp = {}
    tidstr = "', '".join(str(t) for t in tids)
    data = queryDevice.queryDevice("""
            SELECT tid, domain_name
            FROM target_components tc
	    JOIN component_domains cd
	      ON cd.component_id = tc.component_id
            JOIN domains d
	      ON d.domain_id = cd.domain_id
            WHERE tc.tid IN('%s') and domain_type = 'Pfam-A'""" %tidstr, params)
    for ent in data:
        tid = ent[0]
        dom = ent[1]
        try:
            pfam_lkp[tid].append(dom)
        except KeyError:
            pfam_lkp[tid] = [dom]
    return pfam_lkp

#-----------------------------------------------------------------------------------------------------------------------


def count_valid(lkp, valid_doms):
    """Get count of architectures and activities covered by the mapping.
    """
    valz = []
    for arch in lkp.keys():
        valid = False
        doms = arch.split(', ')
        for dom in doms:
            if dom in valid_doms:
                valid = True
                break
        valz.append((lkp[arch], valid))
    valid = sum([x[0] for x in valz if x[1]])
    allz = sum([x[0] for x in valz])
    valid_archs = len([x[0] for x in valz if x[1]])
    all_archs = len(sum([x[0] for x in valz]))
    out = open('data/log.tab', 'a')
    timestamp = time.strftime('%d %B %Y %T', time.gmtime())
    comment = "only binding assays"
    release = params['release']
    threshold = params['threshold']
    out.write("%(valid)s\t%(allz)s\t%(release)s\t%(threshold)s\t%(comment)s\t%(timestamp)s\t%(valid_archs)s\t%(all_archs)s\n"
            % locals())
    out.close()

#-----------------------------------------------------------------------------------------------------------------------


def export_archs(arch_lkp, valid_doms, path):
    '''Write out multi-domain architectures in markdown tables.
    Inputs:
    arch_lkp -- dictionary of multi-domain architectures.
    '''
    sorted_archs = sorted(arch_lkp.iteritems(), key=operator.itemgetter(1), reverse = True)
    out = open('%s.md' % path ,'w')
    out.write('|architecture|count|mapped|\n')
    out.write('|:-----------|:---------|-----:|\n')
    for arch in sorted_archs:
        doms = str(arch[0]).split(', ')
        if len(doms) <= 1:
            continue
        mapped = ', '.join([x for x in doms if x in valid_doms])
        if len(mapped) == 0:
            mapped = False
        out.write("|%s|%s|%s|\n"%(arch[0], arch[1], mapped))

#-----------------------------------------------------------------------------------------------------------------------


def export_network(arch_lkp, valid_doms, path):
    '''Write out network file.
    Inputs:
    arch_lkp -- dictionary of multi-domain architectures.
    '''
    lkp = {}
    for arch in arch_lkp.keys():
        doms = arch.split(', ')
        if len(doms) <= 1:
            continue
        count = arch_lkp[arch]
        if type(doms) is str:
            continue
        for i in range(len(doms)-1):
            for j in range(i+1, len(doms)):
                dom_key = ', '.join(sorted([doms[i],doms[j]]))
                try:
                   lkp[dom_key] += count
                except KeyError:
                   lkp[dom_key] = count
    out = open('%s.tab' % path ,'w')
    out.write('dom_1\tdom_2\tcount\n')
    for link in lkp.keys():
        doms = str(link).split(', ')
        out.write("%s\t%s\t%s\n"%(doms[0], doms[1], lkp[link]))
    out.close()

#-----------------------------------------------------------------------------------------------------------------------


def export_attribs(arch_lkp, valid_doms, path):
    '''Write out network file.
    Inputs:
    arch_lkp -- dictionary of multi-domain architectures.
    '''
    out = open('%s.tab' % path ,'w')
    out.write('dom\tvalid\n')
    lkp = {}
    for arch in arch_lkp.keys():
        doms = arch.split(', ')
        if len(doms) <= 1:
            continue
        for dom in doms:
            valid = False
            if dom in valid_doms:
                valid = True
            lkp[dom] = valid
    for it in lkp.items():
        out.write("%s\t%s\n"%(it[0], it[1]))
    out.close()

#-----------------------------------------------------------------------------------------------------------------------


def export_doms(dom_lkp, valid_doms, path):
    '''Write out identified architectures in markdown tables.
    Inputs:
    dom_lkp -- dictionary of domains occuring in multi-domain architectures.
    '''
    sorted_doms = sorted(dom_lkp.iteritems(), key=operator.itemgetter(1), reverse= True)
    out = open('%s.md' % path ,'w')
    out.write('|domain |count| validated|\n')
    out.write('|:-----------|:-----|-------:|\n')
    for dom in sorted_doms:
        mapped = False
        count = dom[1]
        dom = str(dom[0])
        if dom in valid_doms:
            mapped = True
        out.write("|%s|%s|%s|\n"%(dom, count, mapped))

#-----------------------------------------------------------------------------------------------------------------------

def master(version):
    """
    Function:  master
    Run through all steps to identify mandatory muli-domain architectures.
    """
    # Load the list of validated domains.
    valid_dom_d = readfile('data/valid_pfam_v_%(version)s.tab' % locals(), 'pfam_a', 'pfam_a')
    valid_doms = valid_dom_d.keys()
    ## Load eligible targets.
    el_targets = get_el_targets(params)
    ## Get domains for tids.
    pfam_lkp = get_doms([x[0] for x in el_targets], params)
    ## Add targets with given architecture.
    (arch_lkp, dom_lkp, act_lkp) = get_archs(el_targets, pfam_lkp)
    ## Count covered acrchitectures.
    count_valid(arch_lkp, valid_doms)
    ## Count covered activities.
    count_valid(act_lkp, valid_doms)
    ##  Write multi-domain architechtures to markdown tables.
    export_archs(arch_lkp, valid_doms, 'data/multi_dom_archs_%s'% params['release'])
    ## Write domains from multi-domain architechtures to markdown tables.
    export_doms(dom_lkp, valid_doms, 'data/multi_dom_doms_%s'% params['release'])
    ## export network file.
    export_network(arch_lkp, valid_doms, 'data/multi_dom_network_%s'% params['release'])
    ## export network attribute file.
    export_attribs(arch_lkp, valid_doms, 'data/multi_dom_attributes_%s'% params['release'])


#-----------------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
        import sys

        if len(sys.argv) != 2: # the program name and one argument
                sys.exit("""Parameters are read from mpf.yaml but must specify
		 	    version for data/valid_pfam_v_%(version)s.tab""")
        version = sys.argv[1]

        master(version)
