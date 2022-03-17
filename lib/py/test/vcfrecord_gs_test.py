from re import S
import sys
sys.path.insert(1, 'lib/py')
import vcfrecord

def get_vcf_line(num = 1):
    with open('lib/py/test/gs.vcf','r') as vcf:
        n = 0     
        for line in vcf:
            line = line.strip()
            if (line.startswith('#')):
                pass
            else:
                n = n + 1
                if (n == num):
                    return(line)

def test_get_chr():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(1))
    chr = vcfr.get_chr()
    assert chr == '22'

def test_get_info_af():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(1))
    af = vcfr.get_info_value('AF')
    assert af == '0.00011'

def test_get_info_info():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(1))
    info = vcfr.get_info_value('INFO')
    assert info == '1'

def test_get_info_no_info():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(2))
    info = vcfr.get_info_value('INFO')
    assert info is None

def test_get_info_ER2():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(2))
    info = vcfr.get_info_value('ER2')
    assert info == '0.15847'

def test_get_info_no_ER2():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(1))
    info = vcfr.get_info_value('ER2')
    assert info is None

def test_get_varid():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(1))
    rsid = vcfr.get_varid()
    assert rsid == '22:16482885:G:A'

def test_get_alleles():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(1))
    ref_allele, alt_allele = vcfr.get_alleles()
    assert ref_allele == 'G'
    assert alt_allele == 'A'

def test_get_posn():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(1))
    posn = vcfr.get_posn_as_int()
    assert posn == 16482885
