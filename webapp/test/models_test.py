import sys
import re
import pymongo
sys.path.insert(1, 'webapp/app')
sys.path.insert(2, 'lib/py')

from models import DataStore
from models import _variants
from models import _filepaths
from models import _samples
from vcfrecord import VCFrecord
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
    filename = 'chr' + CHROMOSOME + '.vcf.gz'
    godb.add_filepath_detail(ASSAY, PATH, ASSAY, [filename])

def load_samples():
    '''
    Load sample names
    '''
    godb = MyDb()
    fullpath = PATH + '/' + ASSAY + '/' + 'chr' + CHROMOSOME + '.vcf'
    with open(fullpath,'r') as vcf:
        count = 0
        for line in vcf:
            line = line.strip()
            if (line.startswith('#CHROM')):
                vcfr = VCFrecord(line)
                prf, sfx = vcfr.get_prfx_sfx()
                for idx, field in enumerate(sfx):
                    count += 1
                    godb.process_sample_detail(field, idx, ASSAY)
                break
    godb.flush_sample_buff()
    return(count)


def test_load_vcf():
    '''
    Test that we can load the test data
    '''
    n = load_vcf()
    assert n == 3

def test_filepath():
    '''
    Test filepath was loaded correctly
    '''
    godb = MyDb()
    load_filepath()
    fullpath = PATH + '/' + ASSAY + '/' + 'chr' + CHROMOSOME + '.vcf.gz'
    assert fullpath == godb.get_filepath(ASSAY,CHROMOSOME)

def test_sammples():
    '''
    Test samples were loaded correctly
    '''
    godb = MyDb()
    n = load_samples()
    assert n == 3

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
    '''
    Test variant stats
    '''
    ds = MyDataStore()
    # Set threshold = 0.9
    vars = ds.get_variant_summary_probs(RSID, 0.9)
    assert RSID == vars[0][0]["rsid"]
    assert CHROMOSOME == vars[0][0]["chromosome"]
    assert 1 == vars[0][0]["Missing"]
    assert '1.00000000' == vars[0][0]["hwe_p"]
    assert '0.750000' == vars[0][0]['a_af']
    assert 0.00092393 == vars[0][0]['ref_maf']
    # Set threshold = 0.4
    vars = ds.get_variant_summary_probs(RSID, 0.4)
    assert RSID == vars[0][0]["rsid"]
    assert CHROMOSOME == vars[0][0]["chromosome"]
    assert 0 == vars[0][0]["Missing"]
    assert '1.00000000' == vars[0][0]["hwe_p"]
    assert '0.833333' == vars[0][0]['a_af']
    assert 0.00092393 == vars[0][0]['ref_maf']

def test_get_variant_data_by_range_too_big():
    '''
    Test query range is too large
    '''
    ds = MyDataStore()
    (vars, msg) = ds.get_variant_data_by_range(21, 1, 1000000)
    assert len(vars) == 0
    assert msg.startswith('Range is too great should be 250Kb or less') is True 

def test_get_variant_data_by_range_too_small():
    '''
    Test query range is too small
    '''
    ds = MyDataStore()
    (vars, msg) = ds.get_variant_data_by_range(21, 2, 1)
    assert len(vars) == 0
    assert msg == 'Start pos is greater than End pos' 

def test_get_variant_data_by_range_nothing():
    '''
    Test nothing returned
    '''
    ds = MyDataStore()
    (vars, msg) = ds.get_variant_data_by_range(21, 1, 10)
    assert len(vars) == 0
    assert msg == 'Nothing found in range' 

def test_get_variant_data_by_range():
    '''
    Test query range 
    '''
    ds = MyDataStore()
    (vars, msg) = ds.get_variant_data_by_range('21', 9411200, 9411400)
    assert len(vars) == 2
    assert msg == '' 
    assert vars[0]['position'] == 9411245
    assert vars[0]['alleleA'] == 'C'
    assert vars[0]['alleleB'] == 'A'
    assert vars[0]['samplecount'] == 3
    assert vars[0]['info'] == 0.844373
    assert vars[0]['ref_maf'] == 0.00092393
    assert vars[0]['rsid'] == 'rs181691356'

## RUN LAST
def test_drop_db():
    '''
    drop the database once all the tests are done
    '''
    client = pymongo.MongoClient(DBSERVER)
    client.drop_database(DBNAME_TEST)
    pass