import os
import time

# -------------------------------------------------------------------------
# 1. 설정 (사용자 입력 부분)
# -------------------------------------------------------------------------
missing_genotype = "NN"
input_file = "인풋파일 경로"
map_path = "아웃풋.map"
ped_base = "아웃풋"
ani_id_path = "아웃풋_id.txt"
num_missing_threshold = 0

start_time = time.time()

# -------------------------------------------------------------------------
# 2. 전처리 및 초기화
# -------------------------------------------------------------------------
def get_first_chr(path):
    with open(path, 'r') as f:
        f.readline() # header skip
        line = f.readline()
        if not line: return None
        return line.strip().split('\t')[1]

pre_chr = get_first_chr(input_file)
mk_id_dict = {} # 개체별 유전자형 저장용
stored_ped_files = []

# -------------------------------------------------------------------------
# 3. 메인 변환 로직
# -------------------------------------------------------------------------
try:
    with open(input_file, 'r') as in_f, \
         open(map_path, 'w') as map_out, \
         open(ani_id_path, 'w') as ani_out:
        
        # 헤더 처리 (개체명 추출)
        header = in_f.readline().strip().split('\t')
        individuals = header[3:]
        print(f"개체의 마리 수 :: {len(individuals)}")
        
        # 개체별 유전자형 리스트 초기화
        for ind in individuals:
            mk_id_dict[ind] = []
            ani_out.write(f"{ind}\n")

        # 첫 번째 염색체용 ped 파일 열기
        current_ped_path = f"{ped_base}{pre_chr}_input.ped"
        ped_out = open(current_ped_path, 'w')
        stored_ped_files.append(current_ped_path)

        snp_count = 0
        for line in in_f:
            snp_count += 1
            parts = line.strip().split('\t')
            snp_id, chr_num, pos = parts[0], parts[1], parts[2]
            genotypes = parts[3:]

            # 염색체 번호가 바뀌면 이전 데이터 쓰고 파일 교체
            if chr_num != pre_chr:
                for ind in individuals:
                    ped_out.write("\t".join(mk_id_dict[ind]) + "\n")
                    mk_id_dict[ind] = [] # 메모리 해제
                
                ped_out.close()
                pre_chr = chr_num
                current_ped_path = f"{ped_base}{chr_num}_input.ped"
                ped_out = open(current_ped_path, 'w')
                stored_ped_files.append(current_ped_path)

            # 유전자형 처리 및 2-allele 검증
            allele_check = set()
            for i, geno in enumerate(genotypes):
                if geno == missing_genotype:
                    formatted_geno = "0\t0"
                else:
                    # 'AA' -> 'A\tA'
                    a1, a2 = geno[0], geno[1]
                    formatted_geno = f"{a1}\t{a2}"
                    allele_check.add(a1)
                    allele_check.add(a2)
                
                mk_id_dict[individuals[i]].append(formatted_geno)

            if len(allele_check) > 2:
                print(f"에러 :: SNP id ({snp_id})에서 allele 종류가 2개 이상입니다.")
                # exit(1) 혹은 필요에 따라 처리

            # MAP 파일 작성
            chr_label = "XY" if chr_num == "PseudoX" else chr_num
            map_out.write(f"{chr_label}\t{snp_id}\t0\t{pos}\n")

        # 마지막 남은 데이터 쓰기
        for ind in individuals:
            ped_out.write("\t".join(mk_id_dict[ind]) + "\n")
        ped_out.close()

    print(f"총 snp의 수 :: {snp_count}")

    # -------------------------------------------------------------------------
    # 4. PLINK 포맷 조립 (cbind 기능)
    # -------------------------------------------------------------------------
    print("파일 합치기 및 최종 ped 생성 시작...")
    
    # ID/Sex/Pheno 정보를 포함한 머리부분 생성
    id_info_path = f"{ped_base}_id_sex_phe.ped"
    with open(id_info_path, 'w') as f:
        for ind in individuals:
            # FamilyID IndividualID PaternalID MaternalID Sex Phenotype
            f.write(f"{ind}\t{ind}\t0\t0\t3\t-9\n")

    # 모든 염색체 ped 파일들을 옆으로 병합 (Column Bind)
    final_ped_path = f"{ped_base}_complete.ped"
    
    # 파일을 한 번에 열어서 병합
    file_handles = [open(p, 'r') for p in [id_info_path] + stored_ped_files]
    with open(final_ped_path, 'w') as out_f:
        for _ in range(len(individuals)):
            row_parts = [h.readline().strip() for h in file_handles]
            out_f.write("\t".join(row_parts) + "\n")
            
    for h in file_handles:
        h.close()

    print(f"작업 완료: {final_ped_path}")

except Exception as e:
    print(f"오류 발생: {e}")

# 시간 계산
end_time = time.time()
elapsed = end_time - start_time
print(f"소요 시간: {elapsed:.2f} seconds")