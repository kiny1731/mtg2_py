import os
import re
import logging
import subprocess
from pathlib import Path

# -------------------------------------------------------------------------
# 1. 초기 설정 및 로깅
# -------------------------------------------------------------------------
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def run_pipeline(input_geno_path: str, plink_input_prefix: str):
    """
    input_geno_path: 원본 유전자형 파일 (final report 파일 wide 포맷으로 가공한 것)
    plink_input_prefix: main 파일 실행 후 생성되는 map/ped 파일
    """
    

    input_path = Path(input_geno_path)
    work_dir = input_path.parent
    base_name = input_path.stem
    
    today = base_name.split('_')[0] if '_' in base_name else "result"
    

    final_output_text = f"{today}_Hanwoo_Ref_final_text"
    
    # -------------------------------------------------------------------------
    # STEP 1: PLINK 실행 
    # -------------------------------------------------------------------------
    logger.info("--- [STEP 1] PLINK QC Pipeline Starting ---")
    try:
        # QC 기준 설정 
        steps = [
            ['plink', '--file', plink_input_prefix, '--make-bed', '--out', f"{plink_input_prefix}_binary", '--cow'],
            ['plink', '--bfile', f"{plink_input_prefix}_binary", '--maf', '0.05', '--mind', '0.5', '--geno', '0.2', 
             '--hwe', '0.000001', '--recode', '--out', f"{plink_input_prefix}_second", '--cow'],
            ['plink', '--file', f"{plink_input_prefix}_second", '--make-bed', '--out', f"{plink_input_prefix}_final", '--cow'],
            ['plink', '--bfile', f"{plink_input_prefix}_final", '--recode', '--out', final_output_text, '--cow']
        ]
        
        for cmd in steps:
            logger.info(f"Running: {' '.join(cmd)}")
            subprocess.run(cmd, check=True, capture_output=True)
            
    except subprocess.CalledProcessError as e:
        logger.error(f"PLINK error at step: {e.cmd}\n{e.stderr.decode()}")
        return

    # -------------------------------------------------------------------------
    # STEP 2: Pruned ID/SNP 추출 (Pruned_SNP_ID.py)
    # -------------------------------------------------------------------------
    logger.info("--- [STEP 2] Comparing Pruned IDs and SNPs ---")
    

    ped_path = f"{final_output_text}.ped"
    map_path = f"{final_output_text}.map"
    all_id_path = f"{plink_input_prefix}_id.txt"
    all_map_path = f"{plink_input_prefix}.map"

    prune_in_ids = set()
    prune_in_snps = set()

    # 살아남은 데이터 읽기
    with open(ped_path, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2: prune_in_ids.add(parts[1])
            
    with open(map_path, 'r') as f:
        for line in f:
            parts = line.split('\t')
            if len(parts) >= 2: prune_in_snps.add(parts[1])

    # 제거된 리스트 생성 및 파일 저장
    prune_out_ids = []
    with open(all_id_path, 'r') as f:
        for line in f:
            iid = line.strip()
            if iid not in prune_in_ids: prune_out_ids.append(iid)

    prune_out_snps = []
    with open(all_map_path, 'r') as f:
        for line in f:
            parts = line.split('\t')
            if len(parts) >= 2 and parts[1] not in prune_in_snps:
                prune_out_snps.append(parts[1])

    # -------------------------------------------------------------------------
    # STEP 3: 최종 포맷팅 (formatting_EB_foward.py)
    # -------------------------------------------------------------------------
    logger.info("--- [STEP 3] Formatting Final Genotype File ---")
    
    final_out_path = work_dir / f"{base_name}_originalChip_done.txt"
    
    # 제거 대상 인덱싱 (SNP)
    all_snps_list = []
    with open(all_map_path, 'r') as f:
        all_snps_list = [line.split('\t')[1] for line in f if len(line.split('\t')) >= 2]
    
    remove_snp_idx = {i for i, snp in enumerate(all_snps_list) if snp in prune_out_snps}
    
    # 제거 대상 인덱싱 (ID)
    all_ids_list = []
    with open(all_id_path, 'r') as f:
        all_ids_list = [line.strip() for line in f]
    
    remove_id_idx = {i for i, iid in enumerate(all_ids_list) if iid in prune_out_ids}

    # 스트리밍 필터링 작업
    with open(input_geno_path, 'r') as f_in, open(final_out_path, 'w') as f_out:
        # 헤더 처리
        header = f_in.readline().strip().split('\t')
        new_header = header[:3] + [header[i] for i in range(3, len(header)) if (i-3) not in remove_id_idx]
        f_out.write("\t".join(new_header) + "\n")

        # 본문 처리
        for idx, line in enumerate(f_in):
            if idx in remove_snp_idx: continue
            
            parts = line.strip().split('\t')
            if not parts: continue
            
            row_out = parts[:3] + [parts[i] for i in range(3, len(parts)) if (i-3) not in remove_id_idx]
            f_out.write("\t".join(row_out) + "\n")

    logger.info(f"✨ 실행 완료! 최종 파일: {final_out_path}")

# 실행
if __name__ == "__main__":
    run_pipeline(
        input_geno_path="/media/work/yujisuk/SungYeon/plink_test/260211_sample_input_geno.txt",
        plink_input_prefix="outputtest"
    )