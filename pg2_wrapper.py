"""
Function: pgQuery.py

Generic query function.
--------------
momo.sander@googlemail.com
"""
import psycopg2

def sql_query(query, param, params):
    """
    Processes a query with parameters.
    """
    conn = psycopg2.connect(host = params['host'], user = params['user'], password = params['pword'], database = params['release'], port = params['port'])
    curs = conn.cursor()
    curs.execute(query, param)
    return curs.fetchall()

def sql_execute(query, param, params):
    """
    Processes a query with parameters.
    """
    conn = psycopg2.connect(host = params['host'], user = params['user'], password = params['pword'], database = params['release'], port = params['port'])
    curs = conn.cursor()
    curs.execute(query, param)
    conn.commit()
    conn.close()
    return
