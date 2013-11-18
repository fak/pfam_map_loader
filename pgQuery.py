"""
Function: pgQuery.py

Generic query function.
--------------
momo.sander@googlemail.com
"""
def query(query, usr, pwd, db, hostname):
 
  import psycopg2
  conn = psycopg2.connect(host = params['host'], user = params['user'], password = params['pword'], database = params['release'], port = params['port'])
  curs = conn.cursor()
  
  curs.execute(query)

  return curs.fetchall()


def paramquery(query, param, params):

  import psycopg2
  conn = psycopg2.connect(host = params['host'], user = params['user'], password = params['pword'], database = params['release'], port = params['port'])
  curs = conn.cursor()
  
  curs.execute(query, param)

  return curs.fetchall()
