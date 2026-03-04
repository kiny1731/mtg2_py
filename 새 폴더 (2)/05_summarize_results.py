#!/usr/bin/env python3
"""
05. Hanwoo Genotype Processing 결과 요약

[수정사항]
- 원본 대비 큰 변경 없음 (원본이 올바름)
"""
import os
import logging
import sys

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def summarize_results(prefix="Hanwoo_sample_96_geno_QC_final"):
    """Summary of generated files"""
    
    logger.info("=" * 70)
    logger.info("HANWOO GENOTYPE PROCESSING - RESULTS SUMMARY")
    logger.info("=" * 70)
    
    files_to_check = {
        "Binary files (QC applied)": [
            f"{prefix}.bed",
            f"{prefix}.bim",
            f"{prefix}.fam"
        ],
        "Log file": [
            f"{prefix}.log"
        ],
        "Original PED/MAP": [
            "Hanwoo_sample_96_geno_complete.ped",
            "Hanwoo_sample_96_geno.map"
        ]
    }
    
    logger.info("\n[1] Generated Files:")
    logger.info("-" * 70)
    
    for category, files in files_to_check.items():
        logger.info(f"\n{category}:")
        for f in files:
            if os.path.exists(f):
                size_mb = os.path.getsize(f) / (1024 * 1024)
                logger.info(f"  OK: {f:<50} {size_mb:>10.2f} MB")
            else:
                logger.info(f"  MISSING: {f:<50}")
    
    log_file = f"{prefix}.log"
    if os.path.exists(log_file):
        logger.info("\n[2] Quality Control Statistics:")
        logger.info("-" * 70)
        
        with open(log_file, 'r', errors='ignore') as f:
            log_content = f.read()
        
        key_stats = [
            "variants loaded from .bim file",
            "people loaded from .fam",
            "variants removed due to missing genotype data",
            "variants removed due to minor allele threshold",
            "variants removed due to Hardy-Weinberg",
            "people removed due to missing genotype data",
            "variants and",
            "people pass filters and QC"
        ]
        
        for stat in key_stats:
            for line in log_content.split('\n'):
                if stat in line:
                    logger.info(f"  {line.strip()}")
                    break

def quick_stats(prefix="Hanwoo_sample_96_geno_QC_final"):
    """Calculate quick stats from BIM/FAM"""
    
    logger.info("\n[3] Quick Population Statistics:")
    logger.info("-" * 70)
    
    bim_file = f"{prefix}.bim"
    if os.path.exists(bim_file):
        with open(bim_file, 'r') as f:
            snp_count = sum(1 for _ in f)
        logger.info(f"  Total SNPs: {snp_count:,}")
    
    fam_file = f"{prefix}.fam"
    if os.path.exists(fam_file):
        with open(fam_file, 'r') as f:
            ind_count = sum(1 for _ in f)
        logger.info(f"  Total Individuals: {ind_count:,}")

if __name__ == "__main__":
    target_prefix = "Hanwoo_sample_96_geno_QC_final"
    if len(sys.argv) > 1:
        target_prefix = sys.argv[1]
        
    summarize_results(target_prefix)
    quick_stats(target_prefix)
