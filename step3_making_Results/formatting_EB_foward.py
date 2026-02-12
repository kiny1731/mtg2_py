import os
import re

# -------------------------------------------------------------------------
# 1. 입력 설정 (파일 경로)
# -------------------------------------------------------------------------
map_file = "qc 후 생성된.map"
id_list_file = "step1에서 생성된 id 리스트.txt"
origin_geno = "geno_file_making.py로 만들어진 파일.txt"
prune_snp_file = "prune_out_SNP.txt"
prune_id_file = "prune_out_ID.txt"



# 출력 파일 경로 생성 (원본 map 경로 기반)
base_path = re.match(r"(.*)\.\w+", map_file).group(1)
final_out_path = f"{base_path}_originalChip_input_done_character.txt"
removed_snp_idx_path = f"{base_path}_제거될SNPindex.txt"

# -------------------------------------------------------------------------
# 2. 리스트 로드 및 제거 대상 인덱싱
# -------------------------------------------------------------------------

# 전체 SNP 리스트 및 제거 여부 확인
all_snps = []
with open(map_file, 'r') as f:
    for line in f:
        # chr, snp, dist, pos
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            all_snps.append(parts[1])

# 전체 개체(ID) 리스트
all_ids = []
with open(id_list_file, 'r') as f:
    for line in f:
        all_ids.append(line.strip())

# 제거할 SNP/ID 리스트
with open(prune_snp_file, 'r') as f:
    prune_snps = set(line.strip() for line in f)

with open(prune_id_file, 'r') as f:
    prune_ids = set(line.strip() for line in f)

# 제거할 인덱스 계산 (Perl의 col_bind::count_index 대체)
# SNP 인덱스 (0부터 시작)
remove_snp_indices = set(i for i, snp in enumerate(all_snps) if snp in prune_snps)
# ID 인덱스 (0부터 시작, 헤더의 4번째 열부터가 ID이므로 실제 데이터 인덱스 기준)
remove_id_indices = set(i for i, id_val in enumerate(all_ids) if id_val in prune_ids)

# 제거될 SNP 인덱스 저장
with open(removed_snp_idx_path, 'w') as f:
    for idx in sorted(remove_snp_indices):
        f.write(f"{idx}\n")

# -------------------------------------------------------------------------
# 3. 데이터 필터링 및 최종 파일 생성 (Streaming 처리)
# -------------------------------------------------------------------------
print("필터링 작업 시작 (SNP 및 ID 제거)...")

with open(origin_geno, 'r') as f_in, open(final_out_path, 'w') as f_out:
    # 1. 헤더 처리
    header_line = f_in.readline()
    if header_line:
        header_parts = header_line.strip().split('\t')
        # 앞의 3컬럼(snp, chr, pos)은 유지
        new_header = header_parts[:3]
        
        # ID 컬럼 필터링 (인덱스 3번부터가 ID)
        for i in range(3, len(header_parts)):
            id_idx = i - 3
            if id_idx not in remove_id_indices:
                new_header.append(header_parts[i])
        
        f_out.write("\t".join(new_header) + "\n")

    # 2. 본문 처리 (SNP 필터링 및 ID 필터링)
    snp_ptr = 0 # 현재 읽고 있는 행이 몇 번째 SNP인지 추적
    for line in f_in:
        # 제거 대상 SNP인 경우 건너뜀
        if snp_ptr in remove_snp_indices:
            snp_ptr += 1
            continue
        
        parts = line.strip().split('\t')
        if not parts: continue
        
        # 앞의 3컬럼 유지
        row_out = parts[:3]
        
        # 유전자형 데이터 필터링 (인덱스 3번부터)
        for i in range(3, len(parts)):
            geno_idx = i - 3
            if geno_idx not in remove_id_indices:
                row_out.append(parts[i])
        
        f_out.write("\t".join(row_out) + "\n")
        snp_ptr += 1

print(f"작업이 완료되었습니다.\n최종 파일: {final_out_path}")