"""initial commit

Revision ID: 4ce0890619df
Revises: 
Create Date: 2024-11-12 00:57:55.161765

"""
from alembic import op
import sqlalchemy as sa

class MOL(sa.types.UserDefinedType):
    def get_col_spec(self, **kw):
        return "MOL"
        
class BFP(sa.types.UserDefinedType):
    def get_col_spec(self, **kw):
        return "BFP"

# revision identifiers, used by Alembic.
revision = '4ce0890619df'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    op.execute("CREATE EXTENSION IF NOT EXISTS rdkit;")

    op.create_table('doi',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('doi', sa.String(length=150), nullable=True),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('doi')
    )
    op.create_table('role',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('name', sa.String(length=50), nullable=False),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('name')
    )
    op.create_table('accounts',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('username', sa.String(length=20), nullable=False),
    sa.Column('email', sa.String(length=120), nullable=False),
    sa.Column('name', sa.String(length=256), nullable=False),
    sa.Column('surname', sa.String(length=256), nullable=True),
    sa.Column('academic_level', sa.String(length=32), nullable=True, default=None),
    sa.Column('password', sa.String(length=64), nullable=False),
    sa.Column('created_at', sa.DateTime(), nullable=False),
    sa.Column('updated_at', sa.DateTime(), nullable=True),
    sa.Column('deleted_at', sa.DateTime(), nullable=True),
    sa.Column('role_id', sa.Integer(), nullable=False),
    sa.ForeignKeyConstraint(['role_id'], ['role.id'], ),
    sa.PrimaryKeyConstraint('id'),
    sa.UniqueConstraint('email'),
    sa.UniqueConstraint('username')
    )
    op.create_table('compounds',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('journal', sa.String(length=5000), nullable=True),
    sa.Column('compound_name', sa.String(length=1024), nullable=True),
    sa.Column('compound_image', sa.String(length=5000), nullable=True),
    sa.Column('molecule', MOL(), nullable=True),
    sa.Column('fingerprint', BFP(), nullable=True),
    sa.Column('smiles', sa.String(length=2048), nullable=True),
    sa.Column('article_url', sa.String(length=500), nullable=True),
    sa.Column('inchi_key', sa.String(length=1024), nullable=True),
    sa.Column('exact_molecular_weight', sa.Float(), nullable=True),
    sa.Column('class_results', sa.String(length=1024), nullable=True),
    sa.Column('superclass_results', sa.String(length=1024), nullable=True),
    sa.Column('pathway_results', sa.String(length=1024), nullable=True),
    sa.Column('isglycoside', sa.Boolean(), nullable=True, server_default=sa.text('true')),
    sa.Column('ispublic', sa.Boolean(), nullable=True, server_default=sa.text('true')),
    sa.Column('pubchem_id', sa.String(length=2048), nullable=True),
    sa.Column('inchi', sa.String(length=5000), nullable=True),
    sa.Column('source', sa.String(length=10), nullable=True),
    sa.Column('user_id', sa.Integer(), nullable=True),
    sa.Column('status', sa.String(length=10), nullable=False),
    sa.Column('created_at', sa.DateTime(), nullable=False),
    sa.Column('updated_at', sa.DateTime(), nullable=True),
    sa.Column('deleted_at', sa.DateTime(), nullable=True),
    sa.ForeignKeyConstraint(['user_id'], ['accounts.id'], ),
    sa.PrimaryKeyConstraint('id')
    )

    op.create_index(
        'idx_compounds_molecule_gist',
        'compounds',                    
        ['molecule'],                   
        postgresql_using='gist'         
    )

    op.create_index(
        'idx_compounds_fingerprint_gist',  
        'compounds',                       
        ['fingerprint'],                   
        postgresql_using='gist'            
    )

    op.execute(
        "UPDATE compounds SET fingerprint = morganbv_fp(molecule) WHERE molecule IS NOT NULL;"
    )

    op.execute("""
    CREATE OR REPLACE FUNCTION update_molecule_from_smiles()
    RETURNS TRIGGER AS $$
    BEGIN
        -- This check safely handles INSERTs and UPDATEs
        IF NEW.smiles IS DISTINCT FROM OLD.smiles THEN
            NEW.molecule = mol_from_smiles(NEW.smiles::cstring);
            
            -- === THIS IS THE NEW LINE ===
            -- Also update the fingerprint at the same time
            NEW.fingerprint = morganbv_fp(NEW.molecule);
        END IF;
        
        RETURN NEW;
    END;
    $$ LANGUAGE plpgsql;
    """)

    op.execute("""
    CREATE TRIGGER trg_update_molecule
    BEFORE INSERT OR UPDATE ON compounds
    FOR EACH ROW
    EXECUTE FUNCTION update_molecule_from_smiles();
    """)

    op.create_table('compound_history',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('compound_id', sa.Integer(), nullable=True),
    sa.Column('active', sa.Boolean(), nullable=True, server_default=sa.text('true')),
    sa.Column('journal', sa.String(length=5000), nullable=True),
    sa.Column('smiles', sa.String(length=5000), nullable=True),
    sa.Column('article_url', sa.String(length=500), nullable=True),
    sa.Column('inchi_key', sa.String(length=5000), nullable=True),
    sa.Column('exact_molecular_weight', sa.Float(), nullable=True),
    sa.Column('class_results', sa.String(length=1024), nullable=True),
    sa.Column('superclass_results', sa.String(length=1024), nullable=True),
    sa.Column('pathway_results', sa.String(length=1024), nullable=True),
    sa.Column('pubchem_id', sa.String(length=2048), nullable=True),
    sa.Column('inchi', sa.String(length=5000), nullable=True),
    sa.Column('source', sa.String(length=10), nullable=True),
    sa.Column('created_at', sa.DateTime(), nullable=False, server_default=sa.func.now()),
    sa.Column('updated_at', sa.DateTime(), nullable=True, server_default=sa.func.now()),
    sa.Column('deleted_at', sa.DateTime(), nullable=True),
    sa.ForeignKeyConstraint(['compound_id'], ['compounds.id'], ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('taxa',
    sa.Column('id', sa.Integer(), nullable=False),
    sa.Column('article_url', sa.String(length=5000), nullable=True),
    sa.Column('verbatim', sa.String(length=5000), nullable=True),
    sa.Column('odds', sa.Float(), nullable=True),
    sa.Column('datasourceid', sa.String(length=5000), nullable=True),
    sa.Column('taxonid', sa.Integer(), nullable=True),
    sa.Column('classificationpath', sa.String(length=5000), nullable=True),
    sa.Column('classificationrank', sa.String(length=5000), nullable=True),
    sa.Column('matchtype', sa.String(length=50), nullable=True),
    sa.Column('user_id', sa.Integer(), nullable=True),
    sa.Column('created_at', sa.DateTime(), nullable=False),
    sa.ForeignKeyConstraint(['user_id'], ['accounts.id'], ),
    sa.PrimaryKeyConstraint('id')
    )
    op.create_table('doicomp',
    sa.Column('doi_id', sa.Integer(), nullable=False),
    sa.Column('compound_id', sa.Integer(), nullable=False),
    sa.ForeignKeyConstraint(['compound_id'], ['compounds.id'], ondelete='CASCADE'),
    sa.ForeignKeyConstraint(['doi_id'], ['doi.id'], ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('doi_id', 'compound_id')
    )
    op.create_table('doitaxa',
    sa.Column('doi_id', sa.Integer(), nullable=False),
    sa.Column('taxon_id', sa.Integer(), nullable=False),
    sa.ForeignKeyConstraint(['doi_id'], ['doi.id'], ondelete='CASCADE'),
    sa.ForeignKeyConstraint(['taxon_id'], ['taxa.id'], ondelete='CASCADE'),
    sa.PrimaryKeyConstraint('doi_id', 'taxon_id')
    )
    # ### end Alembic commands ###

def downgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    op.drop_table('doitaxa')
    op.drop_table('doicomp')
    op.drop_table('taxa')
    op.drop_table('compound_history')
    op.drop_table('compounds')
    op.drop_table('accounts')
    op.drop_table('role')
    op.drop_table('doi')
    
    op.execute("DROP TRIGGER IF EXISTS trg_update_molecule ON compounds;")
    op.execute("DROP FUNCTION IF EXISTS update_molecule_from_smiles();")
    op.drop_index('idx_compounds_molecule_gist', table_name='compounds')

    op.drop_index(
        'idx_compounds_fingerprint_gist',
        table_name='compounds'
    )

    op.execute("DROP EXTENSION IF EXISTS rdkit;")
    # ### end Alembic commands ###
