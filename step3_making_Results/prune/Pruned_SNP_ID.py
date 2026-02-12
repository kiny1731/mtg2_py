import os
import re

# -------------------------------------------------------------------------
# 1. 입력 설정 (경로를 본인의 환경에 맞게 수정하세요)
# -------------------------------------------------------------------------
# QC 후 결과 파일 (PLINK output)
ped_path ="아웃풋.ped"
map_path = "아웃풋.map"

# 원본 데이터 파일 (QC 전)
all_id_path ="qc전 id 리스트.txt"
all_snp_path = "qc 전.map"



# 출력 경로 설정 (ped 파일이 있는 디렉토리에 생성)
output_dir = os.path.dirname(ped_path) + "/"

# -------------------------------------------------------------------------
# 2. 개체(ID) 비교 및 리스트 생성
# -------------------------------------------------------------------------
print("개체(ID) 리스트 비교 중...")

prune_in_ids = set()
in_id_file = output_dir + "prune_in_ID.txt"
out_id_file = output_dir + "prune_out_ID.txt"

# QC 후 살아남은 ID 추출
with open(ped_path, 'r') as f_in, open(in_id_file, 'w') as f_out:
    for line in f_in:
        parts = line.strip().split() # 공백(스페이스/탭) 기준 분할
        if len(parts) >= 2:
            fid, iid = parts[0], parts[1]
            f_out.write(f"{iid}\n")
            prune_in_ids.add(iid)

# 원본과 비교하여 제거된 ID 추출
with open(all_id_path, 'r') as f_all, open(out_id_file, 'w') as f_out:
    for line in f_all:
        iid = line.strip()
        if iid not in prune_in_ids:
            f_out.write(f"{iid}\n")

# -------------------------------------------------------------------------
# 3. SNP 리스트 비교 및 생성
# -------------------------------------------------------------------------
print("SNP 리스트 비교 중...")

prune_in_snps = set()
in_snp_file = output_dir + "prune_in_SNP.txt"
out_snp_file = output_dir + "prune_out_SNP.txt"

# QC 후 살아남은 SNP 추출
with open(map_path, 'r') as f_in, open(in_snp_file, 'w') as f_out:
    for line in f_in:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            snp_name = parts[1]
            f_out.write(f"{snp_name}\n")
            prune_in_snps.add(snp_name)

# 원본과 비교하여 제거된 SNP 추출
with open(all_snp_path, 'r') as f_all, open(out_snp_file, 'w') as f_out:
    for line in f_all:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            snp_name = parts[1]
            if snp_name not in prune_in_snps:
                f_out.write(f"{snp_name}\n")

print("-" * 30)
print(f"작업 완료! 결과 파일 위치: {output_dir}")
print(f"1. 살아남은 ID: prune_in_ID.txt")
print(f"2. 제거된 ID  : prune_out_ID.txt")
print(f"3. 살아남은 SNP: prune_in_SNP.txt")
print(f"4. 제거된 SNP  : prune_out_SNP.txt")