import psycopg2
from openpyxl import load_workbook
from tqdm import tqdm
import bcrypt
from datetime import datetime
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
import os
import requests
import json
from dotenv import load_dotenv

load_dotenv()

# Retrieve environment variables
DB_HOST = os.getenv('DB_HOST')
DB_PORT = os.getenv('DB_PORT')
DB_NAME = os.getenv('DB_NAME')
DB_USER = os.getenv('DB_USER')
DB_PASSWORD = os.getenv('DB_PASSWORD')
NPMINE_WEB_APP_PASSWORD = os.getenv('NPMINE_WEB_APP_PASSWORD')
NPMINE_WEB_APP_EMAIL = os.getenv('NPMINE_WEB_APP_EMAIL')

# Establish a connection to the database using environment variables
conn = psycopg2.connect(
    host=DB_HOST,
    port=DB_PORT,
    database=DB_NAME,
    user=DB_USER,
    password=DB_PASSWORD
)

# Create a cursor object to execute SQL queries
cursor = conn.cursor()

# Define the get_or_create function for DOI
def get_or_create_doi(article_url):
    cursor.execute("SELECT id FROM doi WHERE doi = %s", (article_url,))
    doi_id = cursor.fetchone()
    if not doi_id:
        cursor.execute("INSERT INTO doi (doi) VALUES (%s) RETURNING id", (article_url,))
        doi_id = cursor.fetchone()[0]
    else:
        doi_id = doi_id[0]
    return doi_id

# Define the admin role ID, editor role ID, and user role ID
admin_role_id = 1
editor_role_id = 2
user_role_id = 3  # Add this line

# Insert roles if not already present
cursor.execute("SELECT id FROM role WHERE id = %s", (admin_role_id,))
if not cursor.fetchone():
    cursor.execute("INSERT INTO role (id, name) VALUES (%s, %s)", (admin_role_id, 'admin'))

cursor.execute("SELECT id FROM role WHERE id = %s", (editor_role_id,))
if not cursor.fetchone():
    cursor.execute("INSERT INTO role (id, name) VALUES (%s, %s)", (editor_role_id, 'editor'))

cursor.execute("SELECT id FROM role WHERE id = %s", (user_role_id,))  # Add this line
if not cursor.fetchone():  # Add this line
    cursor.execute("INSERT INTO role (id, name) VALUES (%s, %s)", (user_role_id, 'user'))  # Add this line

# Insert the "admin" user if not already present
cursor.execute("SELECT id FROM accounts WHERE username = %s", ('admin',))
admin_id = cursor.fetchone()
if not admin_id:
    # Hash the admin's password
    admin_password = NPMINE_WEB_APP_PASSWORD.encode('utf-8')
    salt = bcrypt.gensalt()
    admin_password_hash = bcrypt.hashpw(admin_password, salt).decode('utf-8')

    # Insert the admin user with the created_at field and admin role_id
    cursor.execute(
        "INSERT INTO accounts (username, email, password, created_at, updated_at, role_id) VALUES (%s, %s, %s, %s, %s, %s) RETURNING id",
        ('admin', NPMINE_WEB_APP_EMAIL, admin_password_hash, datetime.utcnow(), datetime.utcnow(), admin_role_id)
    )
    admin_id = cursor.fetchone()[0]
else:
    admin_id = admin_id[0]

conn.commit()
cursor.close()
conn.close()