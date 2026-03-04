#!/usr/bin/env python3
"""
MTG2 데이터 진단 도구

[수정사항 - 검토보고서 기반]
1. -mzg/-zg 옵션 구분 안내 개선
2. grm_list.txt 형식 검증 추가
3. 표현형 파일 결측치 검증 (NA vs -9)
4. 전체 파이프라인 일관성 검증
"""

import gzip
import os

# 설정
id_file = "total_animals.id"
grm_prefix = "Hanwoo_sample_96_geno_QC_final_grm_gcta"
grm_gz_file = f"{grm_prefix}.grm.gz"
grm_id_file = f"{grm_prefix}.grm.id"

def diagnose():
    print("=" * 60)
    print("MTG2 Data Diagnosis Report (매뉴얼 v2.22 기반)")
    print("=" * 60)

    # [진단 1] ID 파일
    print("\n[1] ID 파일 확인")
    print("-" * 40)
    
    if not os.path.exists(id_file):
        print(f"❌ {id_file} 없음")
        print(f"   생성: awk '{{print $1, $2}}' {grm_id_file} | tr -d '\\r' > {id_file}")
        return
    
    with open(id_file, 'r') as f:
        id_lines = [line.strip() for line in f if line.strip()]
    
    n = len(id_lines)
    print(f"  개체 수: {n}")
    
    # ID 파일 형식 검증
    id_problems = []
    for i, line in enumerate(id_lines[:5]):
        parts = line.split()
        if len(parts) != 2:
            id_problems.append(f"Line {i+1}: {len(parts)} columns (expected 2)")
        if '\t' in line:
            id_problems.append(f"Line {i+1}: TAB 포함")
        if '\r' in line:
            id_problems.append(f"Line {i+1}: Windows 줄바꿈(\\r) 포함")
    
    if id_problems:
        for p in id_problems:
            print(f"  ❌ {p}")
    else:
        print(f"  ✅ 형식 OK (2컬럼, 공백 구분)")
    
    print(f"  첫 3줄: ")
    for line in id_lines[:3]:
        print(f"    {line}")

    # [진단 2] GRM 파일
    print(f"\n[2] GRM 파일 확인")
    print("-" * 40)
    
    expected_rows = n * (n + 1) // 2
    print(f"  기대되는 GRM 행 수: {expected_rows:,}")
    
    if not os.path.exists(grm_gz_file):
        print(f"  ❌ {grm_gz_file} 없음")
        # 다른 GRM 파일 찾기
        for f in os.listdir('.'):
            if '.grm' in f:
                print(f"    발견: {f}")
        return
    
    size_mb = os.path.getsize(grm_gz_file) / (1024**2)
    print(f"  파일: {grm_gz_file} ({size_mb:.2f} MB)")
    
    # GRM 구조 분석
    actual_rows = 0
    col_count = 0
    try:
        with gzip.open(grm_gz_file, 'rt') as f:
            for i, line in enumerate(f):
                if i == 0:
                    col_count = len(line.split())
                actual_rows += 1
                if i < 3:
                    print(f"  샘플 {i+1}행: {line.strip()}")
    except Exception as e:
        print(f"  ❌ 읽기 오류: {e}")
        return

    print(f"  열 개수: {col_count}")
    print(f"  행 개수: {actual_rows:,}")

    # GRM ID 파일
    if os.path.exists(grm_id_file):
        with open(grm_id_file, 'r') as f:
            grm_n = sum(1 for _ in f)
        print(f"  GRM ID 개체 수: {grm_n}")
        
        if grm_n != n:
            print(f"  ❌ ID 파일({n})과 GRM ID({grm_n}) 개체 수 불일치!")
        else:
            print(f"  ✅ ID 개체 수 일치")
    
    # [진단 3] GRM 형식 분석
    print(f"\n[3] GRM 형식 분석")
    print("-" * 40)
    
    if col_count == 4:
        print(f"  ✅ 4열 GCTA/PLINK 형식")
        print(f"  → MTG2 -zg {grm_gz_file} 로 직접 사용 가능 (변환 불필요)")
    elif col_count == 3:
        print(f"  ✅ 3열 MTG2 텍스트 형식")
        print(f"  → MTG2 -g 로 사용 가능")
    else:
        print(f"  ⚠️ 예상 외 열 수: {col_count}")
    
    if actual_rows != expected_rows:
        print(f"  ⚠️ 행 수 불일치: 실제 {actual_rows:,} vs 기대 {expected_rows:,}")
    else:
        print(f"  ✅ 행 수 정확")

    # [진단 4] 표현형 파일
    print(f"\n[4] 표현형 파일 확인")
    print("-" * 40)
    
    for trait_num in range(1, 5):
        # 여러 가능한 파일명 체크
        for fname in [f'trait{trait_num}_matched.pheno', f'trait{trait_num}.pheno']:
            if os.path.exists(fname):
                with open(fname, 'r') as f:
                    lines = f.readlines()
                
                pheno_n = len(lines)
                
                # 결측치 방식 확인
                na_count = sum(1 for l in lines if 'NA' in l.split()[-1] if l.strip())
                minus9_count = sum(1 for l in lines 
                                  if l.strip() and l.strip().split()[-1] in ['-9', '-9.0', '-9.000000'])
                
                status = "✅" if pheno_n == n else "❌"
                print(f"  {fname}: {pheno_n}행 {status}")
                
                if minus9_count > 0:
                    print(f"    ❌ -9 결측치 {minus9_count}개 발견 (MTG2는 NA 필요!)")
                if na_count > 0:
                    print(f"    ✅ NA 결측치 {na_count}개 (올바름)")
                
                break
        else:
            print(f"  trait{trait_num}: 파일 없음")

    # [진단 5] grm_list.txt (다중 GRM 사용 시)
    print(f"\n[5] grm_list.txt 확인")
    print("-" * 40)
    
    if os.path.exists('grm_list.txt'):
        with open('grm_list.txt', 'r') as f:
            content = f.read().strip()
        print(f"  내용: '{content}'")
        
        # -mzg 사용 시 .grm.gz 확장자 필요
        if not content.endswith('.grm.gz'):
            print(f"  ⚠️ -mzg 사용 시 .grm.gz 확장자가 필요합니다!")
            print(f"     수정: echo '{content}.grm.gz' > grm_list.txt")
            print(f"     또는 단일 GRM이면 -zg 직접 사용 권장:")
            print(f"     mtg2 ... -zg {grm_gz_file}")
        else:
            print(f"  ✅ -mzg 형식 OK")
    else:
        print(f"  없음 (단일 GRM → -zg 직접 사용 권장)")

    # [결론]
    print(f"\n{'='*60}")
    print("권장 MTG2 명령어")
    print("=" * 60)
    print(f"""
mtg2 -p {id_file} \\
     -d trait1_matched.pheno \\
     -zg {grm_gz_file} \\
     -out trait1_result \\
     -mod 1 \\
     -sv 0.3 0.7 \\
     -bv trait1_ebv
""")

if __name__ == "__main__":
    diagnose()
