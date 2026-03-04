#!/usr/bin/env python3
"""
06. MTG2 ID 파일 형식 오류 진단 및 수정

[수정사항 - 검토보고서 기반]
1. ID 파일: FID IID (공백 구분, 2컬럼) - GRM .grm.id에서 생성
2. 표현형 파일: FID IID TRAIT (3컬럼), 결측=NA
3. 테스트 스크립트: -mod 1, -zg 사용, -nr 제거
4. GRM 파일 확인: -zg 사용 시 4컬럼 gz OK
"""
import os
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# GRM 설정
GRM_PREFIX = "Hanwoo_sample_96_geno_QC_final_grm_gcta"
GRM_GZ = f"{GRM_PREFIX}.grm.gz"
GRM_ID = f"{GRM_PREFIX}.grm.id"

def diagnose_id_file():
    """ID 파일 진단"""
    logger.info("=" * 70)
    logger.info("Diagnosing ID Files")
    logger.info("=" * 70)
    
    id_files = ['total_animals.id']
    
    for id_file in id_files:
        if os.path.exists(id_file):
            logger.info(f"\nChecking: {id_file}")
            
            with open(id_file, 'r') as f:
                lines = f.readlines()
                
            logger.info(f"  Lines: {len(lines)}")
            logger.info(f"  First 3 lines:")
            
            problems = []
            for i, line in enumerate(lines[:5]):
                raw = line.rstrip('\n')
                parts = raw.split()
                
                # Windows 줄바꿈 체크
                if '\r' in raw:
                    problems.append(f"Line {i+1}: Contains \\r (Windows line ending)")
                
                # 탭 체크
                if '\t' in raw:
                    problems.append(f"Line {i+1}: Contains TAB character")
                
                # 컬럼 수 체크
                if len(parts) != 2:
                    problems.append(f"Line {i+1}: {len(parts)} columns (expected 2)")
                
                if i < 3:
                    logger.info(f"    Line {i+1}: {len(parts)} columns - '{raw[:80]}'")
            
            if problems:
                logger.warning(f"  ❌ Problems found:")
                for p in problems:
                    logger.warning(f"    - {p}")
            else:
                logger.info(f"  ✓ Format OK (2 columns, space separated)")
        else:
            logger.warning(f"\n{id_file} not found")

def fix_id_file():
    """
    GRM .grm.id 파일에서 올바른 ID 파일 생성
    
    MTG2 -p 파일 = GRM과 동일 개체/순서
    """
    logger.info("\n" + "=" * 70)
    logger.info("Fixing ID File")
    logger.info("=" * 70)
    
    id_file = "total_animals.id"
    
    if os.path.exists(GRM_ID):
        count = 0
        with open(GRM_ID, 'r') as grm_in, open(id_file, 'w') as id_out:
            for line in grm_in:
                # Windows 줄바꿈 제거, 공백으로 구분
                parts = line.strip().replace('\r', '').split()
                if len(parts) >= 2:
                    id_out.write(f"{parts[0]} {parts[1]}\n")
                    count += 1
        
        logger.info(f"✓ Created: {id_file} ({count} individuals)")
        logger.info(f"  Source: {GRM_ID}")
    else:
        logger.error(f"❌ {GRM_ID} not found!")
        logger.info(f"  Available .grm.id files:")
        for f in os.listdir('.'):
            if f.endswith('.grm.id'):
                logger.info(f"    - {f}")

def check_pheno_files():
    """
    표현형 파일 형식 확인
    
    MTG2 -d 형식: FID IID TRAIT (3컬럼, 결측=NA)
    """
    logger.info("\n" + "=" * 70)
    logger.info("Checking Phenotype Files")
    logger.info("=" * 70)
    
    pheno_files = [f'trait{i}.pheno' for i in range(1, 5)]
    # 이전 형식 파일도 체크
    pheno_files += [f'trait{i}_matched.pheno' for i in range(1, 5)]
    
    for pheno_file in pheno_files:
        if os.path.exists(pheno_file):
            logger.info(f"\nChecking: {pheno_file}")
            
            with open(pheno_file, 'r') as f:
                lines = f.readlines()
            
            logger.info(f"  Lines: {len(lines)}")
            
            problems = []
            na_count = 0
            minus9_count = 0
            
            for i, line in enumerate(lines[:10]):
                parts = line.strip().split()
                
                if len(parts) != 3:
                    problems.append(f"Line {i+1}: {len(parts)} columns (expected 3: FID IID TRAIT)")
                
                if len(parts) >= 3:
                    if parts[2] == '-9' or parts[2] == '-9.0' or parts[2] == '-9.000000':
                        minus9_count += 1
                    if parts[2] == 'NA':
                        na_count += 1
                
                if i < 3:
                    logger.info(f"    Line {i+1}: {len(parts)} columns - {line.strip()}")
            
            # 전체 파일 결측치 확인
            for line in lines[10:]:
                parts = line.strip().split()
                if len(parts) >= 3:
                    if parts[2] == '-9' or parts[2] == '-9.0' or parts[2] == '-9.000000':
                        minus9_count += 1
                    if parts[2] == 'NA':
                        na_count += 1
            
            if minus9_count > 0:
                problems.append(f"❌ {minus9_count} lines use -9 for missing (MTG2 requires NA)")
            
            if na_count > 0:
                logger.info(f"  ✓ {na_count} lines use NA for missing (correct)")
            
            if '\r' in open(pheno_file, 'r').read():
                problems.append("Contains Windows line endings (\\r)")
            
            if problems:
                logger.warning(f"  Problems:")
                for p in problems:
                    logger.warning(f"    - {p}")
            elif not problems:
                logger.info(f"  ✓ Format OK")

def check_grm_files():
    """GRM 파일 확인"""
    logger.info("\n" + "=" * 70)
    logger.info("Checking GRM Files")
    logger.info("=" * 70)
    
    # GRM gz 파일
    if os.path.exists(GRM_GZ):
        size_mb = os.path.getsize(GRM_GZ) / (1024**2)
        logger.info(f"  ✓ {GRM_GZ} ({size_mb:.2f} MB)")
        logger.info(f"    → MTG2 사용: -zg {GRM_GZ}")
    else:
        logger.error(f"  ✗ {GRM_GZ} NOT FOUND")
        logger.info(f"  Available GRM files:")
        for f in os.listdir('.'):
            if '.grm' in f:
                logger.info(f"    - {f}")
    
    # GRM ID 파일
    if os.path.exists(GRM_ID):
        with open(GRM_ID, 'r') as f:
            n_ids = sum(1 for _ in f)
        logger.info(f"  ✓ {GRM_ID} ({n_ids} IDs)")
    else:
        logger.error(f"  ✗ {GRM_ID} NOT FOUND")
    
    # 개체 수 일치 확인
    if os.path.exists('total_animals.id') and os.path.exists(GRM_ID):
        with open('total_animals.id', 'r') as f:
            n_id = sum(1 for _ in f)
        with open(GRM_ID, 'r') as f:
            n_grm = sum(1 for _ in f)
        
        if n_id == n_grm:
            logger.info(f"  ✓ ID count matches: {n_id}")
        else:
            logger.error(f"  ✗ ID count mismatch! total_animals.id={n_id}, GRM={n_grm}")

def create_test_script():
    """단순 테스트 스크립트 생성 (매뉴얼 준수)"""
    logger.info("\n" + "=" * 70)
    logger.info("Creating Test Script")
    logger.info("=" * 70)
    
    script = f"""#!/bin/bash
# Simple test for Trait 1 (MTG2 매뉴얼 v2.22 준수)

ulimit -s unlimited
export OMP_STACKSIZE=1G

echo "Testing MTG2 with Trait 1..."

# 파일 확인
echo ""
echo "=== File Check ==="
echo "ID file:"
wc -l total_animals.id
head -3 total_animals.id

echo ""
echo "Pheno file:"
wc -l trait1.pheno
head -3 trait1.pheno

echo ""
echo "GRM file:"
ls -lh {GRM_GZ}

echo ""
echo "=== Running MTG2 ==="

mtg2 -p total_animals.id \\
     -d trait1.pheno \\
     -zg {GRM_GZ} \\
     -out trait1_test \\
     -mod 1 \\
     -sv 0.3 0.7

echo ""
if [ -f "trait1_test.vc" ]; then
    echo "✓ Success!"
    echo ""
    echo "Variance components:"
    cat trait1_test.vc
else
    echo "✗ Failed"
    echo ""
    echo "Troubleshooting:"
    echo "  1. 표현형 정규화: python 04_fix_mtg2_likelihood_nan.py"
    echo "  2. 시작값 변경: -sv 0.2 0.8 또는 -sv 0.5 0.5"
    echo "  3. 파일 형식 확인: 결측치가 NA인지 (-9 아님)"
fi
"""
    
    with open('test_mtg2.sh', 'w') as f:
        f.write(script)
    
    os.chmod('test_mtg2.sh', 0o755)
    
    logger.info("Created: test_mtg2.sh")
    logger.info("Run: bash test_mtg2.sh")

if __name__ == "__main__":
    logger.info("=" * 70)
    logger.info("MTG2 FILE FORMAT - DIAGNOSIS & FIX")
    logger.info("=" * 70)
    
    # 1. 진단
    diagnose_id_file()
    check_pheno_files()
    check_grm_files()
    
    # 2. 수정
    fix_id_file()
    
    # 3. 테스트 스크립트
    create_test_script()
    
    logger.info("\n" + "=" * 70)
    logger.info("NEXT STEPS:")
    logger.info("=" * 70)
    logger.info("""
1. Test with single trait:
   bash test_mtg2.sh

2. If successful, run all traits:
   bash run_separate_traits.sh

3. MTG2 -d 표현형 파일 규칙:
   - 3컬럼: FID IID TRAIT_VALUE
   - 결측치: NA (not -9)
   - 공백 구분 (탭도 가능)
   - Windows 줄바꿈 없어야 함
""")
