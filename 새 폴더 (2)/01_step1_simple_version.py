#!/usr/bin/env python3
"""
01. Raw Genotype → PLINK PED/MAP 변환 및 QC

[수정사항]
- 기존 코드 대비 큰 변경 없음 (원본이 올바름)
- 로깅 메시지 일부 개선
"""

import os
import time
import subprocess
import logging

# -------------------------------------------------------------------------
# 초기 설정
# -------------------------------------------------------------------------
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

input_file = "combined_geno.txt"
map_path = "Hanwoo_sample_96_geno.map"
ped_output = "Hanwoo_sample_96_geno_complete.ped"
missing_genotype = "NN"

start_time = time.time()

# -------------------------------------------------------------------------
# STEP 1: Simple and Clear Conversion
# -------------------------------------------------------------------------
def convert_to_ped_simple():
    """
    가장 단순하고 명확한 변환 로직
    한 번에 모든 데이터를 메모리에 로드하여 처리
    """
    logger.info("--- Converting Raw Genotype to PLINK PED/MAP ---")
    
    # Step 1: 파일 읽기 및 개체 목록 추출
    logger.info("Reading input file...")
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    header_parts = lines[0].strip().split('\t')
    individuals = header_parts[3:]  # SNP_ID, CHR, POS 제외
    num_individuals = len(individuals)
    
    logger.info(f"Number of individuals: {num_individuals}")
    logger.info(f"Number of SNPs: {len(lines) - 1}")
    
    # Step 2: 개체별 genotype 데이터 저장소 초기화
    individual_genotypes = {ind: [] for ind in individuals}
    
    # Step 3: MAP 파일 작성 및 genotype 데이터 수집
    logger.info("Processing SNPs and creating MAP file...")
    snp_count = 0
    
    with open(map_path, 'w') as map_out:
        for line_idx in range(1, len(lines)):
            parts = lines[line_idx].strip().split('\t')
            
            if len(parts) < 4:
                logger.warning(f"Skipping line {line_idx + 1}: insufficient columns")
                continue
            
            snp_id, chr_num, pos = parts[0], parts[1], parts[2]
            genotypes = parts[3:]
            
            # MAP 파일에 기록
            map_out.write(f"{chr_num}\t{snp_id}\t0\t{pos}\n")
            
            # 각 개체의 genotype 추가
            for i, ind in enumerate(individuals):
                if i < len(genotypes):
                    geno = genotypes[i]
                    
                    if geno == missing_genotype or geno == 'NN' or len(geno) != 2:
                        individual_genotypes[ind].extend(['0', '0'])
                    else:
                        individual_genotypes[ind].extend([geno[0], geno[1]])
                else:
                    individual_genotypes[ind].extend(['0', '0'])
            
            snp_count += 1
            if snp_count % 10000 == 0:
                logger.info(f"  Processed {snp_count} SNPs...")
    
    logger.info(f"Total SNPs processed: {snp_count}")
    logger.info(f"MAP file created: {map_path}")
    
    # Step 4: PED 파일 작성
    logger.info("Writing PED file...")
    
    with open(ped_output, 'w') as ped_out:
        for idx, ind in enumerate(individuals):
            ped_line_parts = [ind, ind, '0', '0', '3', '-9']
            ped_line_parts.extend(individual_genotypes[ind])
            ped_out.write('\t'.join(ped_line_parts) + '\n')
            
            if (idx + 1) % 100 == 0:
                logger.info(f"  Written {idx + 1}/{num_individuals} individuals...")
    
    logger.info(f"PED file created: {ped_output}")
    logger.info("Conversion completed!")
    
    logger.info("\nValidating files...")
    logger.info(f"  Expected SNPs per individual: {snp_count}")
    logger.info(f"  Expected genotype fields: {snp_count * 2}")
    logger.info(f"  Total PED columns: {6 + snp_count * 2}")
    
    return ped_output, map_path

# -------------------------------------------------------------------------
# STEP 2: PLINK QC
# -------------------------------------------------------------------------
def run_plink_qc(ped_file, map_file):
    """PLINK QC 수행"""
    logger.info("\n--- Running PLINK QC ---")
    
    binary_prefix = "Hanwoo_sample_96_geno_binary"
    final_prefix = "Hanwoo_sample_96_geno_QC_final"
    
    try:
        logger.info("Step 1: Converting to binary format...")
        result = subprocess.run([
            'plink',
            '--ped', ped_file,
            '--map', map_file,
            '--make-bed',
            '--out', binary_prefix,
            '--cow',
            '--allow-no-sex'
        ], capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"PLINK stdout: {result.stdout}")
            logger.error(f"PLINK stderr: {result.stderr}")
            raise subprocess.CalledProcessError(result.returncode, result.args)
        
        logger.info("  Binary files created successfully!")
        
        logger.info("Step 2: Performing Quality Control...")
        result = subprocess.run([
            'plink',
            '--bfile', binary_prefix,
            '--maf', '0.05',
            '--mind', '0.5',
            '--geno', '0.2',
            '--hwe', '0.000001',
            '--make-bed',
            '--out', final_prefix,
            '--cow',
            '--allow-no-sex'
        ], capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"PLINK stdout: {result.stdout}")
            logger.error(f"PLINK stderr: {result.stderr}")
            raise subprocess.CalledProcessError(result.returncode, result.args)
        
        logger.info(f"  QC completed: {final_prefix}.bed/bim/fam")
        
        if os.path.exists(f"{final_prefix}.log"):
            with open(f"{final_prefix}.log", 'r') as log:
                for line in log:
                    if 'variants and' in line or 'people pass' in line:
                        logger.info(f"  {line.strip()}")
        
        return final_prefix
        
    except subprocess.CalledProcessError as e:
        logger.error(f"PLINK failed with error code {e.returncode}")
        raise

# -------------------------------------------------------------------------
# 메인 실행부
# -------------------------------------------------------------------------
if __name__ == "__main__":
    try:
        logger.info("=" * 70)
        logger.info("Hanwoo Genotype Conversion Pipeline")
        logger.info("=" * 70)
        
        # Step 1: 변환
        ped_file, map_file = convert_to_ped_simple()
        
        # Step 2: PLINK QC
        final_prefix = run_plink_qc(ped_file, map_file)
        
        logger.info("\n" + "=" * 70)
        logger.info("Pipeline completed successfully!")
        logger.info("=" * 70)
        logger.info("\nOutput files:")
        logger.info(f"  1. PED/MAP files:")
        logger.info(f"     - {ped_file}")
        logger.info(f"     - {map_file}")
        logger.info(f"  2. Binary files (QC applied):")
        logger.info(f"     - {final_prefix}.bed")
        logger.info(f"     - {final_prefix}.bim")
        logger.info(f"     - {final_prefix}.fam")
        logger.info(f"  3. Log files:")
        logger.info(f"     - {final_prefix}.log")
        logger.info("=" * 70)
        
    except Exception as e:
        logger.error("\n" + "=" * 70)
        logger.error("Pipeline FAILED!")
        logger.error("=" * 70)
        logger.error(f"Error: {e}")
        import traceback
        traceback.print_exc()
        raise
    
    finally:
        elapsed = time.time() - start_time
        logger.info(f"\nTotal runtime: {elapsed:.2f} seconds ({elapsed/60:.2f} minutes)")
