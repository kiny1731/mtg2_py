#!/usr/bin/env python3
"""
02. PLINK를 사용하여 Genomic Relationship Matrix (GRM) 계산

[수정사항]
- MTG2 호환성 안내 강화
- Method 2 (GCTA gz format)가 MTG2 -zg/-mzg와 직접 호환됨을 명시
"""
import subprocess
import logging
import os

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def calculate_grm_with_plink(plink_prefix):
    """
    PLINK를 사용하여 GRM 계산
    
    MTG2 호환:
    - Method 2 (.grm.gz) → MTG2 -zg 또는 -mzg로 직접 사용 가능
    - Method 3 (.grm.bin) → MTG2 -bg 또는 -mbg로 직접 사용 가능
    - Method 1 (.rel) → MTG2 -g 사용 시 3컬럼 변환 필요
    """
    
    grm_prefix = f"{plink_prefix}_grm"
    
    logger.info("=" * 60)
    logger.info("Calculating Genomic Relationship Matrix (GRM)")
    logger.info("=" * 60)
    
    try:
        # Method 1: PLINK relationship matrix (참고용)
        logger.info("\n[Method 1] Creating relationship matrix (.rel)...")
        subprocess.run([
            'plink',
            '--bfile', plink_prefix,
            '--make-rel', 'square',
            '--out', f"{grm_prefix}_plink",
            '--cow',
            '--allow-no-sex'
        ], check=True)
        
        logger.info(f"  Created: {grm_prefix}_plink.rel")
        logger.info(f"  Created: {grm_prefix}_plink.rel.id")
        logger.info(f"  ※ MTG2 -g 사용 시 3컬럼 변환 필요 (13_rel.py)")
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Method 1 failed: {e}")
    
    try:
        # Method 2: GCTA-format GRM (compressed) ← MTG2 권장
        logger.info("\n[Method 2] Creating GCTA-format GRM (compressed) ← MTG2 권장...")
        subprocess.run([
            'plink',
            '--bfile', plink_prefix,
            '--make-grm-gz',
            '--out', f"{grm_prefix}_gcta",
            '--cow',
            '--allow-no-sex'
        ], check=True)
        
        logger.info(f"  Created: {grm_prefix}_gcta.grm.gz")
        logger.info(f"  Created: {grm_prefix}_gcta.grm.id")
        logger.info(f"  ★ MTG2에서 직접 사용: -zg {grm_prefix}_gcta.grm.gz")
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Method 2 failed: {e}")
    
    try:
        # Method 3: GCTA-format GRM (binary)
        logger.info("\n[Method 3] Creating GCTA-format GRM (binary)...")
        subprocess.run([
            'plink',
            '--bfile', plink_prefix,
            '--make-grm-bin',
            '--out', f"{grm_prefix}_gcta_bin",
            '--cow',
            '--allow-no-sex'
        ], check=True)
        
        logger.info(f"  Created: {grm_prefix}_gcta_bin.grm.bin")
        logger.info(f"  Created: {grm_prefix}_gcta_bin.grm.N.bin")
        logger.info(f"  Created: {grm_prefix}_gcta_bin.grm.id")
        logger.info(f"  ★ MTG2에서 직접 사용: -bg {grm_prefix}_gcta_bin.grm.bin")
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Method 3 failed: {e}")
    
    logger.info("\n" + "=" * 60)
    logger.info("GRM calculation completed!")
    logger.info("=" * 60)
    
    # 생성된 파일 목록
    logger.info("\nGenerated GRM files:")
    grm_files = [
        f"{grm_prefix}_plink.rel",
        f"{grm_prefix}_plink.rel.id",
        f"{grm_prefix}_gcta.grm.gz",
        f"{grm_prefix}_gcta.grm.id",
        f"{grm_prefix}_gcta_bin.grm.bin",
        f"{grm_prefix}_gcta_bin.grm.N.bin",
        f"{grm_prefix}_gcta_bin.grm.id"
    ]
    
    for gf in grm_files:
        if os.path.exists(gf):
            size_mb = os.path.getsize(gf) / (1024 * 1024)
            logger.info(f"  ✓ {gf} ({size_mb:.2f} MB)")

if __name__ == "__main__":
    plink_prefix = "Hanwoo_sample_96_geno_QC_final"
    
    logger.info("Input file: {}.bed/bim/fam".format(plink_prefix))
    
    calculate_grm_with_plink(plink_prefix)
    
    logger.info("\n" + "=" * 60)
    logger.info("★ MTG2 사용 시 권장:")
    logger.info(f"  gz 형식:  -zg {plink_prefix}_grm_gcta.grm.gz")
    logger.info(f"  bin 형식: -bg {plink_prefix}_grm_gcta_bin.grm.bin")
    logger.info("=" * 60)
