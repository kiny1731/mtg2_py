#!/usr/bin/env python3
"""
03. Multi-trait Genomic Prediction (4개 형질 개별 예측)

[수정사항 - 검토보고서 기반]
1. -mod 수정: "1 0 1 3 3 1" → "1" (형질 수)
2. GRM 옵션: -zg (단일 GRM gz) 사용으로 통일
3. 표현형 파일 형식: FID IID TRAIT (MTG2 -d 형식)
4. 결측치: NA 문자열 사용 (MTG2 매뉴얼 준수)
5. 컬럼 인덱스 통일: PYNO.fam = [FID, IID, T1, T2, T3, T4]
6. -nr 옵션 제거 (매뉴얼에 없는 옵션)

PYNO.fam 구조:
  col0=FID  col1=IID  col2=Trait1  col3=Trait2  col4=Trait3  col5=Trait4
"""

import subprocess
import logging
import os

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# 파일 설정
reference_fam = "PYNO.fam"
grm_gz_file = "Hanwoo_sample_96_geno_QC_final_grm_gcta.grm.gz"  # gz 파일 직접 지정
grm_id_file = "Hanwoo_sample_96_geno_QC_final_grm_gcta.grm.id"
output_prefix = "multi_trait_prediction"

# =============================================================================
# Step 1: 표현형 파일 확인 및 변환
# =============================================================================
def prepare_phenotype_files():
    """
    PYNO.fam 파일에서 4개 형질을 개별 파일로 추출
    
    MTG2 -d 형식: FID IID TRAIT_VALUE (결측=NA)
    """
    logger.info("=" * 70)
    logger.info("Step 1: Preparing Individual Phenotype Files")
    logger.info("=" * 70)
    
    logger.info(f"\nReading phenotype data from: {reference_fam}")
    
    # 파일 구조 확인
    with open(reference_fam, 'r') as f:
        first_line = f.readline().strip().split()
        logger.info(f"First line preview: {' '.join(first_line[:8])}...")
        logger.info(f"Total columns: {len(first_line)}")
        num_cols = len(first_line)
    
    # PYNO.fam 컬럼 구조 확인
    # 예상: FID IID T1 T2 T3 T4 (6컬럼)
    if num_cols < 6:
        logger.error(f"Expected at least 6 columns (FID IID T1 T2 T3 T4), got {num_cols}")
        return []
    
    pheno_files = []
    
    for trait_num in range(1, 5):
        # col_idx: 2=T1, 3=T2, 4=T3, 5=T4
        col_idx = trait_num + 1
        pheno_file = f"trait{trait_num}.pheno"
        
        trait_count = 0
        missing_count = 0
        
        with open(reference_fam, 'r') as fam_in, open(pheno_file, 'w') as pheno_out:
            for line in fam_in:
                parts = line.strip().split()
                if len(parts) < col_idx + 1:
                    continue
                
                fid = parts[0]
                iid = parts[1]
                trait_val = parts[col_idx]
                
                # 결측치 처리: MTG2는 NA만 인식
                if trait_val in ['-9', 'NA', 'na', '.', '']:
                    pheno_out.write(f"{fid} {iid} NA\n")
                    missing_count += 1
                else:
                    pheno_out.write(f"{fid} {iid} {trait_val}\n")
                    trait_count += 1
        
        logger.info(f"\nTrait {trait_num} (column {col_idx}):")
        logger.info(f"  Valid: {trait_count}, Missing: {missing_count}")
        logger.info(f"  Saved: {pheno_file}")
        pheno_files.append(pheno_file)
    
    return pheno_files

# =============================================================================
# Step 2: ID 파일 생성
# =============================================================================
def create_id_file():
    """
    GRM의 .grm.id 파일에서 ID 추출 (MTG2 -p 옵션용)
    -p 파일과 GRM은 같은 개체/순서여야 함
    """
    logger.info("\n" + "=" * 70)
    logger.info("Step 2: Creating ID File from GRM")
    logger.info("=" * 70)
    
    id_file = "total_animals.id"
    
    if os.path.exists(grm_id_file):
        # GRM의 .grm.id에서 직접 생성 (순서 일치 보장)
        count = 0
        with open(grm_id_file, 'r') as grm_in, open(id_file, 'w') as id_out:
            for line in grm_in:
                parts = line.strip().replace('\r', '').split()
                if len(parts) >= 2:
                    id_out.write(f"{parts[0]} {parts[1]}\n")
                    count += 1
        
        logger.info(f"ID file created from GRM: {id_file} ({count} individuals)")
    else:
        # GRM ID 파일이 없으면 .fam에서 생성
        logger.warning(f"{grm_id_file} not found, using {reference_fam}")
        count = 0
        with open(reference_fam, 'r') as fam_in, open(id_file, 'w') as id_out:
            for line in fam_in:
                parts = line.strip().split()
                if len(parts) >= 2:
                    id_out.write(f"{parts[0]} {parts[1]}\n")
                    count += 1
        
        logger.info(f"ID file created from FAM: {id_file} ({count} individuals)")
    
    return id_file

# =============================================================================
# Step 3: 예측 스크립트 생성
# =============================================================================
def create_prediction_script(pheno_files, id_file):
    """
    개별 형질 예측 스크립트 생성 (MTG2 매뉴얼 준수)
    """
    logger.info("\n" + "=" * 70)
    logger.info("Step 3: Creating MTG2 Prediction Script")
    logger.info("=" * 70)
    
    script_file = "run_separate_traits.sh"
    
    script_content = f"""#!/bin/bash
# ============================================
# MTG2 개별 형질 예측 (매뉴얼 v2.22 기반)
# ============================================

ulimit -s unlimited
export OMP_STACKSIZE=1G

echo "======================================"
echo "MTG2 4형질 개별 예측"
echo "======================================"

success_count=0

for TRAIT_NUM in 1 2 3 4; do
    echo ""
    echo "==================================="
    echo "Predicting Trait ${{TRAIT_NUM}}..."
    echo "==================================="
    
    # MTG2 실행 (매뉴얼 준수 명령어)
    #   -p  : ID 파일 (FID IID) - GRM과 동일 순서
    #   -d  : 표현형 파일 (FID IID TRAIT) - 결측=NA
    #   -zg : 단일 GRM gz 파일 (GCTA/PLINK 형식)
    #   -mod: 형질 수 (단일=1)
    #   -sv : 시작값 (유전분산 잔차분산)
    #   -bv : BLUP 출력
    
    mtg2 -p {id_file} \\
         -d trait${{TRAIT_NUM}}.pheno \\
         -zg {grm_gz_file} \\
         -out trait${{TRAIT_NUM}}_result \\
         -mod 1 \\
         -sv 0.4 0.6 \\
         -bv trait${{TRAIT_NUM}}_ebv
    
    if [ -f "trait${{TRAIT_NUM}}_result.vc" ]; then
        echo "✓ Trait ${{TRAIT_NUM}} completed!"
        echo "  Variance components:"
        cat trait${{TRAIT_NUM}}_result.vc
        ((success_count++))
    else
        echo "✗ Trait ${{TRAIT_NUM}} failed"
        echo "  Tip: likelihood nan → 표현형 정규화 또는 -sv 값 조정 필요"
    fi
done

echo ""
echo "======================================"
echo "완료: $success_count / 4 형질 성공"
echo "======================================"

if [ $success_count -eq 4 ]; then
    echo ""
    echo "결과 파일:"
    ls -lh trait*_result.vc trait*_ebv 2>/dev/null
fi
"""
    
    with open(script_file, 'w') as f:
        f.write(script_content)
    
    os.chmod(script_file, 0o755)
    
    logger.info(f"Created script: {script_file}")
    logger.info("Run: bash run_separate_traits.sh")

# =============================================================================
# 결과 통합 스크립트
# =============================================================================
def create_results_merger():
    """
    개별 형질 예측 결과를 하나의 파일로 통합
    """
    logger.info("\n" + "=" * 70)
    logger.info("Creating Results Merger Script")
    logger.info("=" * 70)
    
    script_file = "merge_predictions.py"
    
    script_content = '''#!/usr/bin/env python3
"""개별 형질 예측 결과를 통합"""

import pandas as pd
import os

def merge_trait_predictions():
    print("Merging prediction results...")
    
    # 기준 ID
    ids = pd.read_csv('total_animals.id', sep=r'\\s+', header=None,
                      names=['FID', 'IID'], dtype=str)
    result = ids.copy()
    
    for i in range(1, 5):
        ebv_file = f'trait{i}_ebv'
        if os.path.exists(ebv_file):
            # MTG2 bv 출력: trait\\n1\\nEBV\\n...
            # 첫 줄은 "trait", 이후 trait번호와 EBV가 번갈아 나옴
            ebv_data = pd.read_csv(ebv_file, sep=r'\\s+', header=0)
            result[f'Trait{i}_EBV'] = ebv_data.iloc[:, -1].values
            print(f"  ✓ Trait {i} merged ({len(ebv_data)} records)")
        else:
            print(f"  ✗ {ebv_file} not found")
    
    # 저장
    output_file = 'all_traits_predictions.txt'
    result.to_csv(output_file, sep='\\t', index=False)
    
    print(f"\\nMerged results: {output_file}")
    print(f"Total individuals: {len(result)}")
    print(f"\\nPreview:")
    print(result.head())
    print(f"\\nSummary:")
    print(result.describe())

if __name__ == "__main__":
    merge_trait_predictions()
'''
    
    with open(script_file, 'w') as f:
        f.write(script_content)
    
    os.chmod(script_file, 0o755)
    
    logger.info(f"Created merger script: {script_file}")
    logger.info("Run after predictions: python merge_predictions.py")

# =============================================================================
# 메인 실행부
# =============================================================================
if __name__ == "__main__":
    
    logger.info("=" * 70)
    logger.info("MULTI-TRAIT GENOMIC PREDICTION (4 Traits)")
    logger.info("=" * 70)
    
    logger.info(f"\nInput files:")
    logger.info(f"  Reference phenotypes: {reference_fam}")
    logger.info(f"  GRM (gz):            {grm_gz_file}")
    logger.info(f"  GRM ID:              {grm_id_file}")
    logger.info(f"\nNumber of traits: 4")
    
    try:
        # Step 1: 표현형 파일 준비
        pheno_files = prepare_phenotype_files()
        
        # Step 2: ID 파일 생성
        id_file = create_id_file()
        
        # Step 3: 예측 스크립트 생성
        create_prediction_script(pheno_files, id_file)
        
        # Step 4: 결과 통합 스크립트
        create_results_merger()
        
        logger.info("\n" + "=" * 70)
        logger.info("All scripts created successfully!")
        logger.info("=" * 70)
        
        logger.info("\nGenerated files:")
        for i in range(1, 5):
            logger.info(f"  - trait{i}.pheno (FID IID TRAIT, 결측=NA)")
        logger.info(f"  - {id_file}")
        logger.info("  - run_separate_traits.sh")
        logger.info("  - merge_predictions.py")
        
        logger.info("\n" + "=" * 70)
        logger.info("WORKFLOW:")
        logger.info("="*70)
        logger.info("\n1. bash run_separate_traits.sh")
        logger.info("2. python merge_predictions.py")
        
    except Exception as e:
        logger.error(f"\nError: {e}")
        import traceback
        traceback.print_exc()
