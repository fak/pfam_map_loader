"""
  Function:  queryDevice
  --------------------
  query the ChEMBL

  momo.sander@googlemail.com
"""
import MySQLdb
def queryDevice(sqlQuery, params):
    conn = MySQLdb.connect(host= params['host'], user= params['user'], passwd= params['pword'], db = params['release'], port = params['port'])
    query = conn.cursor()
    query.execute(sqlQuery)
    queryResults = query.fetchall()
    conn.close()
    return queryResults

