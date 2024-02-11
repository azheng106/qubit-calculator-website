import os
import psycopg2
from psycopg2 import sql


DATABASE_URL = os.environ['DATABASE_URL']
conn = psycopg2.connect(DATABASE_URL, sslmode='require')

cur = conn.cursor()  # Create cursor object to execute our stuff

# Create table (if it doesn't exist)
cur.execute('''
CREATE TABLE IF NOT EXISTS entanglement_cache (
    input_parameters TEXT PRIMARY KEY,
    result TEXT
)
''')


def store_result(state, result):
    try:
        cur.execute(
            sql.SQL('INSERT INTO entanglement_cache (input_parameters, result) VALUES (%s, %s) ON CONFLICT DO NOTHING'),
            [state, result]
        )
        conn.commit()
    except psycopg2.Error as e:
        conn.rollback()

