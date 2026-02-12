import glob
import os
import datetime
from collections import defaultdict

# -------------------------------------------------------------------------
# 1. 파일 확인 및 날짜 설정
# -------------------------------------------------------------------------
"""file_list = glob.glob("*FinalReport.txt")

if len(file_list) != 1:
    print("입력 파일이 1개가 아닙니다. 확인해 주세요.")
    exit()

file_name = file_list[0]"""
file_name = "2022년도_한우참조집단_6,168ea_Analysis_FinalReport.txt"

# 오늘 날짜 기반 파일명 생성 (Perl의 $today 로직 반영)
today_full = datetime.datetime.now().strftime("%Y%m%d")
today = today_full[2:]  # Perl 코드의 split(//,$file_today,3) 결과인 뒷부분 6자리

map_file = f"{today}_input_map.txt"
sample_file = f"{today}_input_sample_list.txt"
geno_file = f"{today}_input_geno.txt"

# -------------------------------------------------------------------------
# 2. FinalReport 읽기 (Map, Sample 정보 추출)
# -------------------------------------------------------------------------
sample_check = []
snp_check = {}
snp_infor = {}

try:
    with open(file_name, 'r') as f:
        # 상단 10줄 헤더 스킵
        for _ in range(10):
            f.readline()
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4: continue
            
            snp, chr_num, pos, sample = parts[0], parts[1], parts[2], parts[3]
            
            # 샘플 리스트 수집 (순서 유지)
            if sample not in sample_check:
                sample_check.append(sample)
            
            # SNP 정보 저장
            if snp not in snp_check:
                snp_check[snp] = True
                snp_infor[snp] = {'chr': chr_num, 'pos': int(pos)}

    # 중간 파일들 저장 (Perl 로직 유지)
    with open(sample_file, 'w') as sf:
        for s in sample_check:
            sf.write(f"{s}\n")
            
    with open(map_file, 'w') as mf:
        for snp, info in snp_infor.items():
            mf.write(f"{snp}\t{info['chr']}\t{info['pos']}\n")

# -------------------------------------------------------------------------
# 3. Geno 파일 생성 (데이터 재구조화)
# -------------------------------------------------------------------------
    # 데이터를 구조적으로 저장하기 위한 딕셔너리
    # 구조: snp_geno_data[chr][pos][snp] = [geno1, geno2, ...]
    snp_geno_data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    
    with open(file_name, 'r') as f:
        for _ in range(10): f.readline()
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 7: continue
            
            snp, chr_num, pos = parts[0], parts[1], parts[2]
            top_a, top_b = parts[4], parts[5]
            
            # 유전자형 포맷팅
            geno = f"{top_a}{top_b}"
            if geno == "--":
                geno = "NN"
            
            try:
                c_int = int(chr_num)
                p_int = int(pos)
                snp_geno_data[c_int][p_int][snp].append(geno)
            except ValueError:
                continue # 숫자가 아닌 염색체 번호 등은 제외

    # 최종 출력
    with open(geno_file, 'w') as out:
        # 헤더: snp chr pos Sample1 Sample2 ...
        samples_header = "\t".join(sample_check)
        out.write(f"snp\tchr\tpos\t{samples_header}\n")
        
        # 염색체 순서(1~29) 및 위치 순서대로 정렬하여 출력
        sorted_chrs = sorted([c for c in snp_geno_data.keys() if 1 <= c <= 29])
        
        for c in sorted_chrs:
            sorted_positions = sorted(snp_geno_data[c].keys())
            for p in sorted_positions:
                for snp in snp_geno_data[c][p]:
                    genos = "\t".join(snp_geno_data[c][p][snp])
                    out.write(f"{snp}\t{c}\t{p}\t{genos}\n")

    # 중간 파일 삭제
    if os.path.exists(map_file): os.remove(map_file)
    if os.path.exists(sample_file): os.remove(sample_file)
    
    print(f"완료! 생성된 파일: {geno_file}")

except Exception as e:
    print(f"오류 발생: {e}")
