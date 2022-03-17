from re import S
import sys
sys.path.insert(1, 'lib/py')
import vcfrecord

def get_vcf_line(num = 1):
    with open('lib/py/test/gd.vcf','r') as vcf:
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
    assert chr == '21'

def test_get_info_af():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(1))
    af = vcfr.get_info_value('RefPanelAF')
    assert af == '0.00092393'

def test_get_info_info1():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(1))
    info = vcfr.get_info_value('INFO')
    assert info == '0.844373'

def test_get_info_info1():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(2))
    info = vcfr.get_info_value('INFO')
    assert info == '1'

def test_get_info_no_ER2():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(1))
    info = vcfr.get_info_value('ER2')
    assert info is None

def test_get_varid():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(1))
    rsid = vcfr.get_varid()
    assert rsid == 'rs181691356'

def test_get_alleles():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(1))
    ref_allele, alt_allele = vcfr.get_alleles()
    assert ref_allele == 'C'
    assert alt_allele == 'A'

def test_get_posn():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(1))
    posn = vcfr.get_posn_as_int()
    assert posn == 9411245

def test_get_calls1():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(1))
    probidx = vcfr.get_probidx()
    prfx, sfx = vcfr.get_prfx_sfx()
    geno_data = sfx[0].split(":")
    probVals = geno_data[probidx]
    call, call_idx, prob = vcfr.get_call(probVals.split(','), 0.3)
    assert call == '0/0'
    assert call_idx == 0
    assert prob == 1.0

def test_get_calls3():
    vcfr = vcfrecord.VCFrecord(get_vcf_line(3))
    probidx = vcfr.get_probidx()
    prfx, sfx = vcfr.get_prfx_sfx()
    geno_data = sfx[0].split(":")
    probVals = geno_data[probidx]
    call, call_idx, prob = vcfr.get_call(probVals.split(','), 0.3)
    assert call == '0/1'
    assert call_idx == 1
    assert prob == 0.545

