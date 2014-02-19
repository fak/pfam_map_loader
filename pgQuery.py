"""
Function: pgQuery.py

Generic query function.
--------------
momo.sander@googlemail.com
"""
import psycopg2

def query(query, usr, pwd, db, hostname):
    """
    Processes a query without paramters.
    """
    conn = psycopg2.connect(host = params['host'], user = params['user'], password = params['pword'], database = params['release'], port = params['port'])
    curs = conn.cursor()
    curs.execute(query)
    return curs.fetchall()


def paramquery(query, param, params):
    """
    Processes a query with parameters.
    """
    conn = psycopg2.connect(host = params['host'], user = params['user'], password = params['pword'], database = params['release'], port = params['port'])

    curs = conn.cursor()
    curs.execute(query, param)
    return curs.fetchall()
