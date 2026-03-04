#!/usr/bin/env python3
"""
04. MTG2 "likelihood nan" 오류 진단 및 해결

[수정사항 - 검토보고서 기반]
1. 정규화 출력 형식: FID IID 유지 (MTG2 -d 형식)
2. 결측치: NA 문자열 사용
3. MTG2 명령어 수정: -mod 1, -zg, -nr 제거
4. 스크립트 내 명령어 형식 통일

오류 원인:
1. 표현형 데이터 스케일 문제
2. 분산 시작값이 부적절
3. GRM 행렬에 문제
4. 모델 설정 오류
"""

import numpy as np
import pandas as pd
import logging
import os

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# =============================================================================
# Step 1: 표현형 데이터 진단
# =============================================================================
def diagnose_phenotype_data(fam_file="PYNO.fam"):
    """
    표현형 데이터의 문제점 찾기
    
    PYNO.fam 구조: FID IID T1 T2 T3 T4
    """
    logger.info("=" * 70)
    logger.info("Step 1: Diagnosing Phenotype Data")
    logger.info("=" * 70)
    
    logger.info(f"\nReading: {fam_file}")
    
    try:
        data = pd.read_csv(fam_file, sep=r'\s+', header=None, dtype=str)
        logger.info(f"Shape: {data.shape} (rows x columns)")
        
        # 컬럼 구조: [0]=FID [1]=IID [2]=T1 [3]=T2 [4]=T3 [5]=T4
        if data.shape[1] < 6:
            logger.error(f"Expected at least 6 columns (FID IID T1 T2 T3 T4), got {data.shape[1]}")
            return None
        
        logger.info(f"\nColumn structure: FID(0) IID(1) T1(2) T2(3) T3(4) T4(5)")
        
        phenotypes = data.iloc[:, 2:6]  # T1~T4
        logger.info(f"Number of traits: {phenotypes.shape[1]}")
        
        logger.info("\n" + "-" * 70)
        logger.info("PHENOTYPE STATISTICS:")
        logger.info("-" * 70)
        
        problems = []
        
        for i, col in enumerate(phenotypes.columns):
            trait_data = phenotypes[col]
            
            # NA, -9, . 을 결측으로 처리
            trait_numeric = pd.to_numeric(trait_data, errors='coerce')
            valid_mask = trait_numeric.notna() & (trait_data != '-9') & (trait_data != 'NA')
            
            trait_valid = trait_numeric[valid_mask]
            n_valid = len(trait_valid)
            n_missing = len(trait_data) - n_valid
            
            if n_valid > 0:
                mean = trait_valid.mean()
                std = trait_valid.std()
                min_val = trait_valid.min()
                max_val = trait_valid.max()
                
                logger.info(f"\nTrait {i+1} (column {i+2}):")
                logger.info(f"  Valid records: {n_valid}")
                logger.info(f"  Missing records: {n_missing}")
                logger.info(f"  Mean: {mean:.4f}")
                logger.info(f"  Std Dev: {std:.4f}")
                logger.info(f"  Min: {min_val:.4f}")
                logger.info(f"  Max: {max_val:.4f}")
                logger.info(f"  Range: {max_val - min_val:.4f}")
                
                if std == 0:
                    problems.append(f"Trait {i+1}: No variation (all same value)")
                
                if std < 0.01:
                    problems.append(f"Trait {i+1}: Very low variation (std={std:.6f})")
                
                if abs(mean) > 1000:
                    problems.append(f"Trait {i+1}: Large scale (mean={mean:.2f}) → 정규화 필요")
                
                # 적절한 시작값 제안
                if std > 0:
                    total_var = std ** 2
                    suggested_vg = total_var * 0.3  # 유전력 0.3 가정
                    suggested_ve = total_var * 0.7
                    logger.info(f"  Suggested -sv: {suggested_vg:.4f} {suggested_ve:.4f}")
                    logger.info(f"  (정규화 후 -sv: 0.3 0.7)")
                
                outliers = ((trait_valid < mean - 5*std) | 
                           (trait_valid > mean + 5*std)).sum()
                if outliers > 0:
                    problems.append(f"Trait {i+1}: {outliers} extreme outliers detected")
                    logger.info(f"  ⚠️  Outliers: {outliers} records")
            else:
                problems.append(f"Trait {i+1}: No valid data")
        
        if problems:
            logger.info("\n" + "=" * 70)
            logger.info("⚠️  POTENTIAL PROBLEMS DETECTED:")
            logger.info("=" * 70)
            for prob in problems:
                logger.info(f"  - {prob}")
        else:
            logger.info("\n✓ No obvious problems detected")
        
        return data
        
    except Exception as e:
        logger.error(f"Error reading file: {e}")
        return None

# =============================================================================
# Step 2: 표현형 데이터 정규화
# =============================================================================
def normalize_phenotypes(input_file="PYNO.fam", output_file="PYNO_normalized.fam"):
    """
    표현형을 표준화 (mean=0, sd=1)
    
    입력: FID IID T1 T2 T3 T4
    출력: FID IID norm_T1 norm_T2 norm_T3 norm_T4 (동일 구조 유지)
    """
    logger.info("\n" + "=" * 70)
    logger.info("Step 2: Normalizing Phenotype Data")
    logger.info("=" * 70)
    
    data = pd.read_csv(input_file, sep=r'\s+', header=None, dtype=str)
    output_data = data.copy()
    
    # col 2~5 (T1~T4) 정규화
    for i in range(2, min(6, data.shape[1])):
        trait = pd.to_numeric(data.iloc[:, i], errors='coerce')
        # NA, -9 → 결측 처리
        valid_mask = trait.notna()
        orig_str = data.iloc[:, i]
        invalid_str = orig_str.isin(['-9', 'NA', 'na', '.', ''])
        valid_mask = valid_mask & (~invalid_str)
        
        if valid_mask.sum() > 0:
            valid_data = trait[valid_mask]
            mean = valid_data.mean()
            std = valid_data.std()
            
            if std > 0:
                # 정규화 적용
                normalized = (trait - mean) / std
                # 유효한 값만 정규화, 나머지는 NA
                output_data.iloc[:, i] = 'NA'  # 기본값 NA
                for idx in valid_mask[valid_mask].index:
                    output_data.iloc[idx, i] = f"{normalized[idx]:.6f}"
                
                logger.info(f"Trait {i-1}: Normalized (mean={mean:.4f}, std={std:.4f})")
            else:
                logger.warning(f"Trait {i-1}: Skipped (no variation)")
    
    # 공백으로 구분하여 저장 (FID IID norm_T1 norm_T2 norm_T3 norm_T4)
    output_data.to_csv(output_file, sep=' ', header=False, index=False)
    logger.info(f"\n✓ Normalized phenotypes saved to: {output_file}")
    logger.info(f"  Format: FID IID T1 T2 T3 T4 (missing=NA)")
    
    return output_file

# =============================================================================
# Step 3: MTG2 수정 스크립트 생성
# =============================================================================
def create_fixed_mtg2_scripts(grm_gz="Hanwoo_sample_96_geno_QC_final_grm_gcta.grm.gz"):
    """
    문제 해결된 MTG2 실행 스크립트 생성 (매뉴얼 v2.22 준수)
    """
    logger.info("\n" + "=" * 70)
    logger.info("Step 3: Creating Fixed MTG2 Scripts")
    logger.info("=" * 70)
    
    # 방법 1: 정규화된 데이터 + 개별 형질
    script1 = "run_mtg2_normalized.sh"
    
    content1 = f"""#!/bin/bash
# MTG2 with normalized phenotypes - 개별 형질 분석
# MTG2 매뉴얼 v2.22 준수

ulimit -s unlimited
export OMP_STACKSIZE=1G

echo "Running MTG2 with normalized phenotypes..."

# 정규화된 파일에서 개별 형질 추출
for TRAIT_NUM in 1 2 3 4; do
    echo ""
    echo "=========================================="
    echo "Analyzing Trait ${{TRAIT_NUM}}..."
    echo "=========================================="
    
    # 형질 컬럼 추출 (col 0=FID, 1=IID, 2=T1, 3=T2, ...)
    TRAIT_COL=$((TRAIT_NUM + 2))  # awk는 1-based
    awk -v col=$TRAIT_COL '{{print $1, $2, $col}}' PYNO_normalized.fam > trait${{TRAIT_NUM}}_norm.pheno
    
    echo "Phenotype file preview:"
    head -3 trait${{TRAIT_NUM}}_norm.pheno
    
    # MTG2 실행
    mtg2 -p total_animals.id \\
         -d trait${{TRAIT_NUM}}_norm.pheno \\
         -zg {grm_gz} \\
         -out trait${{TRAIT_NUM}}_result \\
         -mod 1 \\
         -sv 0.3 0.7 \\
         -bv trait${{TRAIT_NUM}}_ebv
    
    if [ -f "trait${{TRAIT_NUM}}_result.vc" ]; then
        echo "✓ Trait ${{TRAIT_NUM}} completed"
        cat trait${{TRAIT_NUM}}_result.vc
    else
        echo "✗ Trait ${{TRAIT_NUM}} failed"
    fi
done

echo ""
echo "======================================"
echo "All traits processed!"
echo "======================================"
"""
    
    with open(script1, 'w') as f:
        f.write(content1)
    os.chmod(script1, 0o755)
    logger.info(f"Created: {script1}")
    
    # 방법 2: 다양한 시작값 시도
    script2 = "run_mtg2_try_sv.sh"
    
    content2 = f"""#!/bin/bash
# MTG2 - 다양한 시작값으로 시도
# likelihood nan 발생 시 시작값 변경

ulimit -s unlimited
export OMP_STACKSIZE=1G

TRAIT_NUM=${{1:-1}}
echo "Testing Trait ${{TRAIT_NUM}} with different starting values..."

# 시작값 세트
declare -a SV_SETS=("0.3 0.7" "0.2 0.8" "0.5 0.5" "0.1 0.9" "0.4 0.6")

for SV in "${{SV_SETS[@]}}"; do
    echo ""
    echo "--- Trying -sv $SV ---"
    
    mtg2 -p total_animals.id \\
         -d trait${{TRAIT_NUM}}_norm.pheno \\
         -zg {grm_gz} \\
         -out trait${{TRAIT_NUM}}_test \\
         -mod 1 \\
         -sv $SV
    
    if [ -f "trait${{TRAIT_NUM}}_test.vc" ]; then
        echo "✓ Success with -sv $SV"
        cat trait${{TRAIT_NUM}}_test.vc
        break
    else
        echo "✗ Failed with -sv $SV"
    fi
done
"""
    
    with open(script2, 'w') as f:
        f.write(content2)
    os.chmod(script2, 0o755)
    logger.info(f"Created: {script2}")

# =============================================================================
# Step 4: 파라미터 가이드
# =============================================================================
def explain_parameters():
    """
    MTG2 파라미터 설명 (매뉴얼 v2.22 기준)
    """
    logger.info("\n" + "=" * 70)
    logger.info("MTG2 Parameter Guide (Manual v2.22)")
    logger.info("=" * 70)
    
    guide = """
올바른 명령어:
  mtg2 -p total_animals.id \\
       -d trait1.pheno \\
       -zg grm_file.grm.gz \\
       -out result \\
       -mod 1 \\
       -sv 0.3 0.7 \\
       -bv result_ebv

파라미터 설명 (매뉴얼 기준):
------------------------------
-p   : ID 파일 = PLINK .fam 형식 (FID IID ..., 첫 2컬럼만 사용)
-d   : 표현형 파일 (FID IID TRAIT). 결측치 = NA
-zg  : 단일 GRM gz 파일 (GCTA/PLINK 형식, 4컬럼 OK)
-mzg : 다중 GRM gz 파일 리스트 (리스트 안에 .grm.gz 확장자 포함)
-bg  : 단일 GRM binary 파일
-mod : 형질 수 (단일=1, 다형질=N)
-sv  : 분산 시작값 (유전분산 잔차분산)
-bv  : BLUP 출력 파일명

주의: 매뉴얼에 없는 옵션
--------------------------
-nr  : 매뉴얼에 명시되지 않음 → 사용 자제
-mod 1 0 1 3 3 1 : 잘못된 형식 → -mod 1 사용

분산 시작값 가이드:
------------------
정규화 후 (mean=0, var≈1):
  유전력 0.3 → -sv 0.3 0.7
  유전력 0.2 → -sv 0.2 0.8
  유전력 0.5 → -sv 0.5 0.5

원본 데이터 사용 시:
  총분산 = phenotypic_variance
  -sv (h2 * Vp) ((1-h2) * Vp)
"""
    
    logger.info(guide)

# =============================================================================
# 메인 실행부
# =============================================================================
if __name__ == "__main__":
    
    logger.info("=" * 70)
    logger.info("MTG2 LIKELIHOOD NAN - DIAGNOSIS & FIX")
    logger.info("=" * 70)
    
    try:
        # Step 1: 진단
        phenotype_data = diagnose_phenotype_data("PYNO.fam")
        
        if phenotype_data is not None:
            # Step 2: 정규화
            normalized_file = normalize_phenotypes("PYNO.fam", "PYNO_normalized.fam")
            
            # Step 3: 수정 스크립트 생성
            create_fixed_mtg2_scripts()
            
            # Step 4: 가이드
            explain_parameters()
            
            logger.info("\n" + "=" * 70)
            logger.info("RECOMMENDED ACTIONS:")
            logger.info("=" * 70)
            logger.info("""
1. 정규화 후 분석:
   python 04_fix_mtg2_likelihood_nan.py   # 정규화 파일 생성
   bash run_mtg2_normalized.sh            # 분석 실행

2. likelihood nan 지속 시:
   bash run_mtg2_try_sv.sh 1              # Trait 1에 다양한 시작값 시도

3. 정규화 없이 원본 사용 시:
   진단 결과의 Suggested -sv 값을 사용하세요.
""")
        
    except Exception as e:
        logger.error(f"\nError: {e}")
        import traceback
        traceback.print_exc()
