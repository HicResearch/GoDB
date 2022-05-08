import sys
import re
import pymongo
sys.path.insert(1, 'webapp/app')
sys.path.insert(2, 'lib/py')

from models import DataStore
from models import _variants
from models import _filepaths
from models import _samples
from godb import GoDb
from dbconfig import DBSERVER
from dbconfig import DBNAME_TEST
from dbconfig import VARIANTS

'''
Subclass GoDb to use a test database
'''
class MyDb(GoDb):
    def __init__(self):
        try:
            connection = pymongo.MongoClient(DBSERVER)
            db = eval("connection." + DBNAME_TEST)
            variants = eval("db." + VARIANTS)
        except:
            raise Exception("Unexpected error connecting to %s @ %s" % (DBNAME_TEST, DBSERVER))
        print("Connected to db: " + DBSERVER + ", " + DBNAME_TEST)
        self.db = db
        self.variants = variants
        self.samples = db.samples
        self.filedata = db.filedata
        self.filepaths = db.filepaths
        self.genemap = db.genemap
        self.dbname = DBNAME_TEST
        self.variantbuff = []
        self.samplebuff = []
        self.filebuff = []
        self.genemapbuff = []
        self.int_fields = ["position"]
        self.flt_fields = ["all_maf", "info", "cohort_1_hwe"]
        self.p = re.compile('\d+')

'''
Subclass DataStore to use the above test database
'''
class MyDataStore(DataStore):
    def __init__(self):
        self.godb = MyDb()
        self.filepaths_coll = _filepaths(self.godb)
        self.sam_coll = _samples(self.godb)
        self.var_coll = _variants(self.filepaths_coll, self.sam_coll, self.godb)
        self.variant_totals = []
        self.sample_count = -1
        self.call_rates = {}
        self.maxrslist = 200

# global vars of test variables
PATH = 'lib/py/test/'
ASSAY = 'gd'
CHROMOSOME = '21'
RSID = 'rs181691356'

def load_vcf():
    '''
    Load test data from VCF file
    '''
    godb = MyDb()
    fullpath = PATH + '/' + ASSAY + '/' + 'chr' + CHROMOSOME + '.vcf'
    with open(fullpath,'r') as vcf:
        count = 0
        for line in vcf:
            line = line.strip()
            if (line.startswith('#')):
                pass
            else:
                godb.process_variant_detail_vcf(line, ASSAY)
                count = count + 1
    godb.flush_variant_buff()
    return(count)

def load_filepath():
    '''
    Load filepath for test file
    '''
    godb = MyDb()
    fullpath = PATH + '/' + ASSAY + '/' + 'chr' + CHROMOSOME + '.vcf'
    godb.add_filepath_detail(ASSAY, PATH, ASSAY, fullpath)

def test_load_vcf():
    '''
    Test that we can load the test data
    '''
    n = load_vcf()
    assert n == 3

def test_filepath():
    
    godb = MyDb()
    load_filepath()
    assert PATH + '/' + ASSAY == godb.get_filepath(ASSAY,CHROMOSOME)


def test_connect():
    '''
    Test connection to test database
    '''
    ds = MyDataStore()
    assert DBNAME_TEST == ds.get_db_name()

def test_get_variant_data_multi_count():
    '''
    Test that get_variant_data_multi() returns the correct
    number of records for the required rsID
    '''
    ds = MyDataStore()
    docs = ds.var_coll.get_variant_data_multi(RSID)
    assert len(docs) == 1

def test_get_variant_summary_probs():
    pass

## RUN LAST
def test_drop_db():
    '''
    drop the database once all the tests are done
    '''
    client = pymongo.MongoClient(DBSERVER)
    client.drop_database(DBNAME_TEST)
    pass